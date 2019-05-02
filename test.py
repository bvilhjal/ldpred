"""
A test file for LDpred.

Examples
--------
To run all tests:
$ python test.py

To run a specific test:
$ python -m unittest test.TestLDpred.test_ldpred_inf
"""

import LDpred
import cPickle
import filecmp
import glob
import gzip
import h5py
from ldpred import coord_genotypes
from ldpred import ld
from ldpred import sum_stats_parsers
import numpy as np
import os
import tempfile
import unittest

np.set_printoptions(linewidth=int(os.environ.get('COLUMNS', 100)))
TEST_DIR = os.path.dirname(os.path.abspath(__file__))

def run_test(mesg, cmd_str, error_mesg, *actual_and_golden_outputs):
  print(mesg)
  print(cmd_str + '\n')
  cmd_args = cmd_str.split()
  try:
    LDpred.main_with_args(cmd_args)
    for i in xrange(0, len(actual_and_golden_outputs), 2):
      actual_output = actual_and_golden_outputs[i]
      golden_output = os.path.join(TEST_DIR, actual_and_golden_outputs[i + 1])
      print('Diffing actual (%s) vs. golden (%s) outputs...' % (actual_output, golden_output))
      assert_files_equal(actual_output, golden_output)
      print('Diff passed!')
  except:
    print(error_mesg + '\n')
    raise


def h5_node_walker(h5_node, key_prefix=''):
  """Generator function that walks an hdf5 File or Group object.

  Args:
    h5_node: an h5py.File or h5py.Group object
    key_prefix: the '/' delimited string representing the name path of the
       node within the .hdf5 file.

  Yields:
    (child_key, child_value)
  """
  for k, v in h5_node.iteritems():
    v_type = type(v)
    v_path = key_prefix + '/' + k
    if v_type == h5py.Group:
      for nested_key, nested_value in h5_node_walker(v, v_path):
        yield nested_key, nested_value
    elif v_type == h5py.Dataset:
      yield v_path, v[:]
    else:
      assert False, 'Unexpected v_type: %s' % v_type


def h5_file_walker(h5_file):
  """Generator function that walks an hdf5 file.

  Args:
    h5_file: a string, the name of the .hdf5 file to walk.

  Yields:
    (child_key, child_value)
  """
  with h5py.File(h5_file, 'r') as h5_root_node:
    for k, v in h5_node_walker(h5_root_node):
      yield k, v


def pkl_node_walker(pkl_node, key_prefix=''):
  """Generator function that walks a Python pickle node (i.e. a dict).

  Args:
    pkl_node: A dict coming from a depickled object.
    key_prefix: the '/' delimited string representing the name path of the
       node within the pickle file.

  Yields:
    (child_key, child_value)
  """
  for k in sorted(pkl_node.keys()):
    v = pkl_node[k]
    v_type = type(v)
    v_path = key_prefix + '/' + str(k)
    if v_type == dict:
      for nested_key, nested_value in pkl_node_walker(v, v_path):
        yield nested_key, nested_value
    elif v_type == list:
      # Convert Python list to Numpy ndarray for assert_deep_equals.
      yield v_path, np.array(v)
    elif v_type in (float, np.float64, np.float32, int, str, np.ndarray):
      yield v_path, v
    else:
      assert False, 'Unexpected v_type: %s' % v_type


def pkl_file_walker(pkl_file):
  """Generator function that walks a Python pickle file.

  Args:
    pkl_file: a string, the name of the .pkl.gz file to walk.

  Yields:
    (child_key, child_value)
  """
  with gzip.open(pkl_file) as f:
    pkl_root_node = cPickle.load(f)
    for k, v in pkl_node_walker(pkl_root_node):
      yield k, v


def assert_deep_equals(walker1, walker2):
  """Test function that does a deep comparison of two structure walkers."""
  for (k1, v1), (k2, v2) in zip(walker1, walker2):
    assert k1 == k2, 'Key mismatch: %s vs. %s' % (k1, k2)
    assert type(v1) == type(v2), 'Type mismatch: %s vs. %s' % (type(v1), type(v2))
    if isinstance(v1, str) or isinstance(v1, int):
      assert v1 == v2, 'Value mismatch: %s vs. %s' % (v1, v2)
    elif isinstance(v1, float) or isinstance(v1, np.float32):
      assert np.isclose(v1, v2), 'Float mismatch: %s vs. %s' % (v1, v2)
    elif isinstance(v1, np.ndarray):
      assert v1.dtype == v2.dtype, 'dtype mismatch: %s vs. %s' % (v1.dtype, v2.dtype)
      if np.issubdtype(v1.dtype, np.number):
        assert np.allclose(v1, v2), 'ndarray number mismatch in key %s' % k1
      else:
        assert np.array_equal(v1, v2), 'ndarray non-number mismatch in key %s' % k1


def assert_files_equal(file1, file2):
  if file1.endswith('.hdf5'):
    assert_deep_equals(h5_file_walker(file1), h5_file_walker(file2))
  elif file1.endswith('.pkl.gz'):
    assert_deep_equals(pkl_file_walker(file1), pkl_file_walker(file2))
  else:
    assert filecmp.cmp(file1, file2), "Mismatch between: %s and %s" % (file1, file2)

def make_p_dict(*args):
  return vars(LDpred.parser.parse_args(args))

class TestLDpred(unittest.TestCase):
  @classmethod
  def setUpClass(cls):
    print('Testing LDpred.\n')
    print('Note that this test currently only tests the core functionality of LDpred.')
    print('Please report bugs on github (https://github.com/bvilhjal/ldpred) or to Bjarni J Vilhjalmsson (bjarni.vilhjalmsson@gmail.com).\n')

  def setUp(self):
    self.tf = tempfile.NamedTemporaryFile()
    self.tmp_file_prefix = next(tempfile._get_candidate_names())

  def tearDown(self):
    print('Cleaning up files.')
    cmd_str = 'rm -f %s*' % self.tmp_file_prefix
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems cleaning up test files!  Testing stopped'

  def test_parse_sum_stats(self):
    p_dict = {
        'ssf': os.path.join(TEST_DIR, 'test_data/coord_genotypes_ss.txt'),
        'ssf_format': 'STANDARD',
        'only_hm3': False,
        'N': 10000,
        'debug': True,
        'match_genomic_pos': False,
    }
    bimfile = os.path.join(TEST_DIR, 'test_data/LDpred_cc_data_p0.001_train_0.bim')
    summary_dict = {}
    out = '%s_parse_sum_stats.hdf5' % self.tmp_file_prefix
    with h5py.File(out, 'w') as h5f:
      sum_stats_parsers.parse_sum_stats(h5f, p_dict, bimfile, summary_dict)
      self.assertEqual(len(h5f['sum_stats']['chrom_1']['betas']), 10)

  def test_coord_genotypes(self):
    p_dict = make_p_dict(
        '--debug',
        'coord',
        '--gf=%s/test_data/LDpred_data_p0.001_train_0' % TEST_DIR,
        '--vgf=%s/test_data/LDpred_data_p0.001_test_0' % TEST_DIR,
        '--ssf=%s/test_data/coord_genotypes_ss.txt' % TEST_DIR,
        '--ssf-format=STANDARD',
        '--N=10000',
        '--out=%s_coord_genotypes.hdf5' % self.tmp_file_prefix,
    )
    summary_dict = coord_genotypes.main(p_dict)
    # summary_dict[11]['value'], if present, is the count of non-matching nts.
    # It should be 0.
    self.assertEqual(summary_dict.get(11, {}).get('value', 0), 0)
    with h5py.File(p_dict['out'], 'r') as h5f:
      self.assertEqual(len(h5f['sum_stats']['chrom_1']['betas']), 10)

  def test_ld_calculation(self):
    df = h5py.File('%s/test_data/goldens/golden.coord0.hdf5' % TEST_DIR, 'r')
    g = df['cord_data']['chrom_1']
    snps, n_raw_snps, n_snps = ld.extract_snps_from_cord_data_chrom(g)
    first_10_snps = snps[:10]
    self.assertEqual(len(first_10_snps), 10)
    ld_dict_and_scores = ld.get_LDpred_ld_tables(first_10_snps)
    ld_dict = ld_dict_and_scores['ld_dict']
    ld_mat = np.vstack([ld_dict[i] for i in xrange(10)])
    golden_ld_mat = np.load(os.path.join(TEST_DIR, 'test_data/ld_data.npz'))['ld']
    self.assertTrue(np.allclose(ld_mat, golden_ld_mat))

  def test_get_chromosome_herits(self):
    p_dict = make_p_dict(
        '--debug',
        'inf',
        '--cf=%s/test_data/goldens/golden.coord.hdf5' % TEST_DIR,
        '--ldr=100',
        '--ldf=' + self.tmp_file_prefix,
        '--N=10000',
        '--out=' + self.tmp_file_prefix,
    )
    summary_dict = {}
    ld_dict = ld.get_ld_dict_using_p_dict(p_dict, summary_dict)
    coord_file = os.path.join(TEST_DIR, 'test_data/goldens/golden.coord.hdf5')
    df = h5py.File(coord_file, 'r')
    herit_dict = ld.get_chromosome_herits(df['cord_data'], ld_dict['ld_scores_dict'], n=p_dict['N'], h2=None)
    print(herit_dict)
    self.assertAlmostEqual(herit_dict['chrom_1'], 0.0741219)
    self.assertAlmostEqual(herit_dict['gw_h2_ld_score_est'], 0.0741219)

  def test_ldpred_coord0(self):
    coord_file = self.tmp_file_prefix + '.coord0.hdf5'
    run_test(
        'Coordinating test data into file %s' % coord_file,
        'coord --gf=%s/test_data/LDpred_data_p0.001_train_0 --vgf=%s/test_data/LDpred_data_p0.001_test_0 --ssf=%s/test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --N=10000  --out=%s' % (TEST_DIR, TEST_DIR, TEST_DIR, coord_file),
        'Problems when coordinating data!',
        coord_file,
        'test_data/goldens/golden.coord0.hdf5'
    )

  def test_ldpred_coord(self):
    coord_file = self.tmp_file_prefix + '.coord.hdf5'
    run_test(
        'Coordinating test data into file %s' % coord_file,
        '--debug coord --gf=%s/test_data/LDpred_data_p0.001_train_0 --vbim=%s/test_data/LDpred_data_p0.001_test_0.bim --ssf=%s/test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --beta --N=10000  --out=%s' % (TEST_DIR, TEST_DIR, TEST_DIR, coord_file),
        'Problems when coordinating data!',
        coord_file,
        'test_data/goldens/golden.coord.hdf5')

  def test_ldpred_inf(self):
    run_test(
        'Running LDpred-inf with coordinated file prefix: %s ' % self.tmp_file_prefix,
        '--debug inf --cf=%s/test_data/goldens/golden.coord.hdf5  --ldr=100   --ldf=%s  --N=10000  --out=%s' % (TEST_DIR, self.tmp_file_prefix, self.tmp_file_prefix),
        'Problems when running LDpred_inf!',
        self.tmp_file_prefix + '.txt',
        'test_data/goldens/golden.txt',
        self.tmp_file_prefix + '_ldradius100.pkl.gz',
        'test_data/goldens/golden_ldradius100.pkl.gz')

  def test_ldpred_gibbs(self):
    np.random.seed(42)  # Set random seed to stabilize test results
    run_test(
        'Running LDpred with coordinated file prefix: %s ' % self.tmp_file_prefix,
        '--debug gibbs --cf=%s/test_data/goldens/golden.coord.hdf5  --ldr=100   --ldf=%s  --f=0.001 --N=10000  --out=%s' % (TEST_DIR, self.tmp_file_prefix, self.tmp_file_prefix),
        'Problems when running LDpred!',
        self.tmp_file_prefix + '_LDpred-inf.txt',
        'test_data/goldens/golden_LDpred-inf.txt',
        self.tmp_file_prefix + '_LDpred_p1.0000e-03.txt',
        'test_data/goldens/golden_LDpred_p1.0000e-03.txt')

  def test_ldpred_p_plus_t(self):
    run_test(
        'Running P+T with coordinated file prefix: %s ' % self.tmp_file_prefix,
        '--debug p+t --cf=%s/test_data/goldens/golden.coord.hdf5  --ldr=100  --p=0.001 --out=%s' % (TEST_DIR, self.tmp_file_prefix),
        'Problems when running P+T!',
        self.tmp_file_prefix + '_P+T_r0.20_p1.0000e-03.txt',
        'test_data/goldens/golden_P+T_r0.20_p1.0000e-03.txt')

  def test_ldpred_score_1(self):
    prs_file_prefix = self.tmp_file_prefix + 'prs'
    run_test(
        'Validating results with output file prefix: %s' % prs_file_prefix,
        '--debug score --gf=%s/test_data/LDpred_data_p0.001_test_0  --rf=%s/test_data/goldens/golden  --out=%s' % (TEST_DIR, TEST_DIR, prs_file_prefix),
        'Problems with the validation step!',
        prs_file_prefix + '_LDpred-inf.txt',
        'test_data/goldens/goldenprs_LDpred-inf.txt',
        prs_file_prefix + '_LDpred_p1.0000e-03.txt',
        'test_data/goldens/goldenprs_LDpred_p1.0000e-03.txt')

  def test_ldpred_score_2(self):
    prs_file_prefix2 = self.tmp_file_prefix + 'prs_2'
    run_test(
        'Validating results with output file prefix: %s' % self.tmp_file_prefix,
        'score --gf=%s/test_data/LDpred_data_p0.001_test_0  --rf=%s/test_data/goldens/golden  --out=%s' % (TEST_DIR, TEST_DIR, prs_file_prefix2),
        'Problems with the validation step!',
        prs_file_prefix2 + '_LDpred-inf.txt',
        'test_data/goldens/goldenprs_2_LDpred-inf.txt',
        prs_file_prefix2 + '_LDpred_p1.0000e-03.txt',
        'test_data/goldens/goldenprs_2_LDpred_p1.0000e-03.txt',)

  def test_ldpred_score_p_plus_t(self):
    prs_file_prefix2 = self.tmp_file_prefix + 'prs_2'
    run_test(
        'Validating results with output file prefix: %s' % self.tmp_file_prefix,
        'score --gf=%s/test_data/LDpred_data_p0.001_test_0  --rf=%s/test_data/goldens/golden  --rf-format=P+T --out=%s' % (TEST_DIR, TEST_DIR, prs_file_prefix2),
        'Problems with the P+T validation step!',
        prs_file_prefix2 + '_P+T_p1.0000e-03.txt',
        'test_data/goldens/goldenprs_2_P+T_p1.0000e-03.txt')

  def test_ldpred_score_3(self):
    prs_file_prefix3 = self.tmp_file_prefix + 'prs_3'
    run_test(
        'Validating results with output file prefix: %s' % self.tmp_file_prefix,
        'score --gf=%s/test_data/LDpred_data_p0.001_test_0  --only-score --rf=%s/test_data/goldens/golden  --rf-format=P+T --out=%s' % (TEST_DIR, TEST_DIR, prs_file_prefix3),
        'Problems with the P+T validation step!',
        prs_file_prefix3 + '_P+T_p1.0000e-03.txt',
        'test_data/goldens/goldenprs_3_P+T_p1.0000e-03.txt')


if __name__ == '__main__':
  unittest.main()
