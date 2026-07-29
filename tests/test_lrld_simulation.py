"""Regression tests for simulations/lrld_policy.py.

These keep the simulation runnable, they do not re-derive its conclusions. The
one scientific assertion is the qualitative effect that motivated it and that
holds with a large margin: with a strong long-range factor, sampling those
variants normally makes the Gibbs sampler diverge.
"""

import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, 'simulations'))

import lrld_policy  # noqa: E402


FAST = dict(
    m=800, ld_radius=20, lrld_start=300, lrld_width=200, rho=0.85,
    n_ref=1000, n_gwas=8000, n_test=1000, h2=0.5, n_causal=20,
    frac_causal_in_lrld=0.25, num_iter=20, burn_in=5,
)


class PanelTests(unittest.TestCase):

    def test_background_ld_decays_inside_the_band(self):
        rng = np.random.default_rng(0)
        X = lrld_policy.simulate_panel(rng, 800, 600, rho=0.85,
                                       lrld_slice=slice(0, 0), lrld_loading=0.0)
        r_near = abs(np.corrcoef(X[:, 100], X[:, 105])[0, 1])
        r_edge = abs(np.corrcoef(X[:, 100], X[:, 120])[0, 1])
        r_far = abs(np.corrcoef(X[:, 100], X[:, 300])[0, 1])
        self.assertGreater(r_near, r_edge)
        self.assertLess(r_far, 0.15)

    def test_long_range_factor_creates_distance_independent_ld(self):
        rng = np.random.default_rng(0)
        sl = slice(200, 500)
        X = lrld_policy.simulate_panel(rng, 1500, 600, rho=0.85,
                                       lrld_slice=sl, lrld_loading=1.0)
        # Two variants far apart inside the region stay strongly correlated.
        r_far_inside = abs(np.corrcoef(X[:, 220], X[:, 480])[0, 1])
        self.assertGreater(r_far_inside, 0.2)

    def test_band_share_falls_as_the_factor_strengthens(self):
        rng = np.random.default_rng(0)
        sl = slice(200, 500)
        shares = []
        for loading in (0.0, 1.5):
            X = lrld_policy.simulate_panel(rng, 1200, 600, rho=0.85,
                                           lrld_slice=sl, lrld_loading=loading)
            shares.append(lrld_policy.describe_ld(X, sl, 20)['band_share_of_ld_mass'])
        self.assertGreater(shares[0], shares[1])


class ReplicateTests(unittest.TestCase):

    def test_replicate_reports_every_policy(self):
        cfg = dict(FAST, loading=1.0)
        row = lrld_policy.one_replicate(0, cfg)
        for policy in lrld_policy.POLICIES:
            self.assertIn(policy + '_r2', row)
            self.assertIn(policy + '_sum_sq', row)
            self.assertTrue(np.isfinite(row[policy + '_r2']))
            self.assertGreaterEqual(row[policy + '_r2'], 0.0)

    def test_replicate_is_deterministic(self):
        cfg = dict(FAST, loading=1.0)
        a = lrld_policy.one_replicate(7, cfg)
        b = lrld_policy.one_replicate(7, cfg)
        for policy in lrld_policy.POLICIES:
            self.assertEqual(a[policy + '_r2'], b[policy + '_r2'])

    def test_sampling_long_range_variants_diverges(self):
        """The finding that motivates excluding them at all.

        Sum of squared effects should sit near the true h2 for a converged fit.
        Under 'include' it runs away by orders of magnitude, which is exactly
        what LDpred's own convergence warning reports.
        """
        cfg = dict(FAST, loading=1.0)
        row = lrld_policy.one_replicate(0, cfg)
        self.assertGreater(row['include_sum_sq'], 10 * cfg['h2'])
        for policy in ('zero', 'inf', 'drop'):
            self.assertLess(row[policy + '_sum_sq'], 5 * cfg['h2'])

    def test_excluding_beats_including_under_strong_long_range_ld(self):
        cfg = dict(FAST, loading=1.0)
        row = lrld_policy.one_replicate(0, cfg)
        for policy in ('zero', 'inf', 'drop'):
            self.assertGreater(row[policy + '_r2'], row['include_r2'])


if __name__ == '__main__':
    unittest.main()
