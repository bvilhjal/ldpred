"""Unit tests for the Gibbs sampler kernel and its two backends.

The sampler has a compiled (Numba) sweep and an equivalent NumPy sweep. They
must produce the same numbers, and the packing helpers that feed them must
reproduce the window bookkeeping the original inline loop did. Nothing here
needs plink or HDF5, so it runs in under a second.
"""

import unittest

import numpy as np

from ldpred import LDpred_gibbs
from ldpred._numba import HAVE_NUMBA


def make_ld(m, ld_radius, seed=0, dtype=np.float32):
    rng = np.random.default_rng(seed)
    ld_dict = {}
    for i in range(m):
        start, stop = max(0, i - ld_radius), min(m, i + ld_radius + 1)
        v = (rng.standard_normal(stop - start) * 0.05).astype(dtype)
        v[min(ld_radius, i)] = 1.0          # unit diagonal
        ld_dict[i] = v
    return ld_dict


class PackingTests(unittest.TestCase):

    def test_flatten_matches_the_python_window_arithmetic(self):
        m, r = 200, 15
        ld_dict = make_ld(m, r)
        flat, offsets, lengths, starts = LDpred_gibbs.flatten_ld_dict(
            ld_dict, m, r, None)
        for i in range(m):
            self.assertEqual(starts[i], max(0, i - r))
            self.assertEqual(lengths[i], len(ld_dict[i]))
            np.testing.assert_allclose(
                flat[offsets[i]:offsets[i] + lengths[i]], ld_dict[i], rtol=0, atol=0)

    def test_flatten_honours_explicit_ld_boundaries(self):
        m, r = 60, 5
        ld_boundaries = [[max(0, i - 3), min(m, i + 4)] for i in range(m)]
        ld_dict = {i: np.ones(b[1] - b[0]) for i, b in enumerate(ld_boundaries)}
        _, _, lengths, starts = LDpred_gibbs.flatten_ld_dict(
            ld_dict, m, r, ld_boundaries)
        for i, (lo, hi) in enumerate(ld_boundaries):
            self.assertEqual(starts[i], lo)
            self.assertEqual(lengths[i], hi - lo)

    def test_constants_as_arrays_fixed_n(self):
        m = 50
        cd = LDpred_gibbs.prepare_constants(1000.0, None, m, 0.1, 0.3, 1.0)
        n, hdmpn, hh, c, d = LDpred_gibbs.constants_as_arrays(cd, m)
        for arr, key in ((n, 'n'), (hdmpn, 'hdmpn'), (hh, 'hdmp_hdmpn'),
                         (c, 'c_const'), (d, 'd_const')):
            self.assertEqual(arr.shape, (m,))
            np.testing.assert_allclose(arr, cd[key])

    def test_constants_as_arrays_varying_n(self):
        m = 50
        ns = np.linspace(800, 1200, m)
        cd = LDpred_gibbs.prepare_constants(None, ns, m, 0.1, 0.3, 1.0)
        n, hdmpn, hh, c, d = LDpred_gibbs.constants_as_arrays(cd, m)
        for i in range(m):
            entry = cd['snp_dict'][i]
            self.assertAlmostEqual(n[i], entry['n'])
            self.assertAlmostEqual(hdmpn[i], entry['hdmpn'])
            self.assertAlmostEqual(hh[i], entry['hdmp_hdmpn'])
            self.assertAlmostEqual(c[i], entry['c_const'])
            self.assertAlmostEqual(d[i], entry['d_const'])


class BackendEquivalenceTests(unittest.TestCase):
    """The compiled and NumPy sweeps must agree to floating point noise."""

    def _args(self, m=400, r=25, seed=1, skip=None):
        rng = np.random.default_rng(seed)
        ld_dict = make_ld(m, r, seed=seed)
        flat, offsets, lengths, starts = LDpred_gibbs.flatten_ld_dict(
            ld_dict, m, r, None)
        cd = LDpred_gibbs.prepare_constants(2e5, None, m, 0.05, 0.3, 1.0)
        consts = LDpred_gibbs.constants_as_arrays(cd, m)
        skip_mask = np.zeros(m, dtype=bool) if skip is None else skip
        return dict(
            beta_hats=np.ascontiguousarray(rng.standard_normal(m) * 2e-3),
            ld_flat=flat, ld_offsets=offsets, ld_lengths=lengths, ld_starts=starts,
            base_betas=rng.standard_normal(m) * 1e-3,
            rand_ps=rng.random(m), rand_norms=rng.standard_normal(m) * 1e-4,
            skip_mask=skip_mask, consts=consts, m=m)

    def _run(self, sweep, a):
        curr = a['base_betas'].copy()
        post = np.zeros(a['m'])
        sweep(a['beta_hats'], a['ld_flat'], a['ld_offsets'], a['ld_lengths'],
              a['ld_starts'], curr, post, a['rand_ps'], a['rand_norms'],
              a['skip_mask'], *a['consts'], 0.99, a['m'])
        return curr, post

    def test_backends_agree(self):
        a = self._args()
        c1, p1 = self._run(LDpred_gibbs._gibbs_sweep_jit, a)
        c2, p2 = self._run(LDpred_gibbs._gibbs_sweep_numpy, a)
        # Summation order differs (scalar loop vs BLAS), so allow last-bit noise.
        np.testing.assert_allclose(c1, c2, rtol=1e-11, atol=1e-18)
        np.testing.assert_allclose(p1, p2, rtol=1e-11, atol=1e-18)

    def test_backends_agree_with_skipped_variants(self):
        m = 400
        skip = np.zeros(m, dtype=bool)
        skip[100:150] = True
        a = self._args(m=m, skip=skip)
        c1, p1 = self._run(LDpred_gibbs._gibbs_sweep_jit, a)
        c2, p2 = self._run(LDpred_gibbs._gibbs_sweep_numpy, a)
        np.testing.assert_allclose(c1, c2, rtol=1e-11, atol=1e-18)
        np.testing.assert_allclose(p1, p2, rtol=1e-11, atol=1e-18)

    def test_skipped_variants_keep_their_starting_effect(self):
        """A skipped variant is frozen, not resampled.

        Its posterior mean stays zero, which is why it currently contributes a
        zero weight to the output rather than its LDpred-inf weight.
        """
        m = 300
        skip = np.zeros(m, dtype=bool)
        skip[50:80] = True
        a = self._args(m=m, skip=skip)
        curr, post = self._run(LDpred_gibbs._gibbs_sweep_numpy, a)
        np.testing.assert_array_equal(curr[skip], a['base_betas'][skip])
        np.testing.assert_array_equal(post[skip], np.zeros(int(skip.sum())))
        self.assertTrue(np.any(post[~skip] != 0.0))


class SamplerSmokeTests(unittest.TestCase):

    def test_ldpred_gibbs_runs_and_is_finite(self):
        m, r = 300, 20
        rng = np.random.default_rng(3)
        res = LDpred_gibbs.ldpred_gibbs(
            rng.standard_normal(m) * 2e-3,
            start_betas=rng.standard_normal(m) * 1e-3,
            h2=0.3, n=2e5, ns=None, ld_radius=r, num_iter=12, burn_in=4,
            p=0.05, ld_dict=make_ld(m, r), verbose=False, print_progress=False)
        self.assertEqual(res['betas'].shape, (m,))
        self.assertTrue(np.all(np.isfinite(res['betas'])))

    def test_skipped_variants_get_their_ldpred_inf_weight(self):
        """Long-range LD variants come back with their starting effect.

        They are held out of the sweep because the LD band is misspecified
        there, but they must still carry their LDpred-inf weight into the
        output, which is what --incl-long-range-ld's help text promises.  They
        used to come back as zeros, silently dropping ~5% of the genome.
        """
        m, r = 300, 20
        rng = np.random.default_rng(5)
        start_betas = rng.standard_normal(m) * 1e-3 + 0.05   # distinctive
        snp_lrld = np.zeros(m, dtype=bool)
        snp_lrld[100:130] = True
        res = LDpred_gibbs.ldpred_gibbs(
            rng.standard_normal(m) * 2e-3, start_betas=start_betas.copy(),
            h2=0.3, n=2e5, ns=None, ld_radius=r, num_iter=20, burn_in=5,
            p=0.05, ld_dict=make_ld(m, r), snp_lrld=snp_lrld,
            verbose=False, print_progress=False)
        np.testing.assert_allclose(res['betas'][snp_lrld],
                                   start_betas[snp_lrld], rtol=1e-12)
        self.assertTrue(np.all(res['betas'][snp_lrld] != 0.0))
        # Retained variants are still sampled, so they must not simply echo
        # their starting values.
        self.assertTrue(np.any(res['betas'][~snp_lrld] != start_betas[~snp_lrld]))

    def test_skipped_variants_are_unaffected_by_iteration_count(self):
        """Averaging a constant over the retained sweeps returns it unchanged."""
        m, r = 200, 15
        rng = np.random.default_rng(6)
        start_betas = rng.standard_normal(m) * 1e-3 + 0.02
        snp_lrld = np.zeros(m, dtype=bool)
        snp_lrld[50:70] = True
        kw = dict(h2=0.3, n=2e5, ns=None, ld_radius=r, p=0.05,
                  ld_dict=make_ld(m, r), snp_lrld=snp_lrld,
                  verbose=False, print_progress=False)
        beta_hats = rng.standard_normal(m) * 2e-3
        a = LDpred_gibbs.ldpred_gibbs(beta_hats.copy(), num_iter=15, burn_in=5,
                                      start_betas=start_betas.copy(), **kw)
        b = LDpred_gibbs.ldpred_gibbs(beta_hats.copy(), num_iter=40, burn_in=10,
                                      start_betas=start_betas.copy(), **kw)
        np.testing.assert_allclose(a['betas'][snp_lrld], b['betas'][snp_lrld],
                                   rtol=1e-12)

    def test_varying_sample_sizes_path_runs(self):
        m, r = 300, 20
        rng = np.random.default_rng(4)
        res = LDpred_gibbs.ldpred_gibbs(
            rng.standard_normal(m) * 2e-3,
            start_betas=rng.standard_normal(m) * 1e-3,
            h2=0.3, n=None, ns=rng.integers(150000, 250000, m).astype(float),
            ld_radius=r, num_iter=12, burn_in=4, p=0.05,
            ld_dict=make_ld(m, r), verbose=False, print_progress=False)
        self.assertTrue(np.all(np.isfinite(res['betas'])))

    def test_backend_reported_in_status(self):
        from ldpred._numba import numba_status
        status = numba_status()
        self.assertIsInstance(status, str)
        if HAVE_NUMBA:
            self.assertIn('enabled', status)
        else:
            self.assertIn('unavailable', status)


if __name__ == '__main__':
    unittest.main()
