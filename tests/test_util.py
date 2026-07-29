"""Unit tests for ldpred.util.

The existing suite is made entirely of end-to-end golden-file comparisons on
simulated data whose SNP positions run from 0 to a few thousand base pairs.
That data can never enter a long-range LD region (the earliest starts at 25 Mb),
so the golden tests pass whether or not the long-range LD filter works at all.
These tests cover that logic directly.
"""

import unittest

import numpy as np

from ldpred import util


class LongRangeLDTests(unittest.TestCase):

    def setUp(self):
        self.lrld_dict = util.load_lrld_dict()

    def test_regions_are_loaded_into_reg_dict(self):
        # Regions used to be stored as siblings of 'reg_dict' rather than
        # inside it, which left every lookup reading an empty dict.
        total = sum(len(v['reg_dict']) for v in self.lrld_dict.values())
        self.assertGreater(total, 0)
        for chrom, entry in self.lrld_dict.items():
            self.assertEqual(set(entry.keys()), {'reg_dict'},
                             'chromosome %s has stray keys outside reg_dict' % chrom)

    def test_mhc_region_is_detected(self):
        # chr6:25.5-33.5 Mb is the MHC, the canonical long-range LD region.
        self.assertTrue(util.is_in_lrld(6, 30000000, self.lrld_dict))
        self.assertFalse(util.is_in_lrld(6, 10000000, self.lrld_dict))

    def test_some_position_is_flagged_genome_wide(self):
        flagged = sum(
            util.is_in_lrld(chrom, pos, self.lrld_dict)
            for chrom in util.ok_chromosomes
            for pos in range(0, 250000000, 1000000)
        )
        self.assertGreater(flagged, 0)

    def test_chromosome_may_be_a_group_name(self):
        # Callers iterate util.chromosomes_list, which yields 'chrom_<n>'
        # names; the counter they used to pass ran one past the last key.
        self.assertTrue(util.is_in_lrld('chrom_6', 30000000, self.lrld_dict))
        self.assertEqual(util.get_lrld_chrom_key('chrom_X'), 23)
        self.assertEqual(util.get_lrld_chrom_key('chrom_6'), 6)
        self.assertEqual(util.get_lrld_chrom_key(6), 6)

    def test_every_chromosome_in_the_iteration_order_resolves(self):
        for chrom_str in util.chromosomes_list:
            key = util.get_lrld_chrom_key(chrom_str)
            self.assertIn(key, self.lrld_dict)

    def test_status_mask_matches_scalar_lookup(self):
        rng = np.random.default_rng(0)
        positions = rng.integers(0, 250000000, 5000)
        for chrom in (1, 6, 11):
            mask = util.get_snp_lrld_status(chrom, positions, self.lrld_dict)
            self.assertEqual(mask.dtype, np.dtype(bool))
            expected = np.array([util.is_in_lrld(chrom, int(p), self.lrld_dict)
                                 for p in positions])
            np.testing.assert_array_equal(mask, expected)

    def test_status_mask_handles_empty_input(self):
        mask = util.get_snp_lrld_status(1, np.array([], dtype=int), self.lrld_dict)
        self.assertEqual(len(mask), 0)


class CalcAucTests(unittest.TestCase):

    @staticmethod
    def _brute_force_auc(y_true, y_hat):
        """P(case > control) + 0.5 P(tie), evaluated over every pair."""
        cases = y_hat[y_true == 1]
        controls = y_hat[y_true == 0]
        return np.mean([(c > k) + 0.5 * (c == k)
                        for c in cases for k in controls])

    def test_matches_brute_force_without_ties(self):
        rng = np.random.default_rng(0)
        y_true = (rng.random(200) < 0.3).astype(float)
        y_hat = y_true * 0.7 + rng.normal(size=200)
        self.assertAlmostEqual(util.calc_auc(y_true, y_hat),
                               self._brute_force_auc(y_true, y_hat), places=12)

    def test_matches_brute_force_with_ties(self):
        # The old threshold sweep double counted tied scores and returned
        # 0.889 here, where the correct answer is 0.722.
        y_true = np.array([0., 0., 1., 1., 0., 1.])
        y_hat = np.array([1., 2., 2., 3., 3., 3.])
        self.assertAlmostEqual(util.calc_auc(y_true, y_hat),
                               self._brute_force_auc(y_true, y_hat), places=12)

    def test_perfect_and_reversed_separation(self):
        y_true = np.array([0., 0., 1., 1.])
        self.assertAlmostEqual(util.calc_auc(y_true, np.array([1., 2., 3., 4.])), 1.0)
        self.assertAlmostEqual(util.calc_auc(y_true, np.array([4., 3., 2., 1.])), 0.0)

    def test_accepts_one_two_coded_phenotypes(self):
        y_hat = np.array([1., 2., 3., 4.])
        self.assertAlmostEqual(util.calc_auc(np.array([1., 1., 2., 2.]), y_hat),
                               util.calc_auc(np.array([0., 0., 1., 1.]), y_hat))


class ChromosomeListTests(unittest.TestCase):

    def test_chromosomes_list_has_no_duplicates(self):
        self.assertEqual(len(util.chromosomes_list),
                         len(set(util.chromosomes_list)))


if __name__ == '__main__':
    unittest.main()
