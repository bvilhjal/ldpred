"""
A test file for LDpred.
"""

import LDpred
import cPickle
import filecmp
import glob
import gzip
import h5py
from ldpred import sum_stats_parsers
import numpy as np
import os
import tempfile
import unittest

def run_test(mesg, cmd_str, error_mesg, *actual_and_golden_outputs):
  print(mesg)
  print(cmd_str + '\n')
  cmd_args = cmd_str.split()
  try:
    np.random.seed(42)  # Set random seed to stabilize test results
    LDpred.main_with_args(cmd_args)
    for i in xrange(0, len(actual_and_golden_outputs), 2):
      actual_output = actual_and_golden_outputs[i]
      golden_output = actual_and_golden_outputs[i + 1]
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
  with h5py.File(h5_file) as h5_root_node:
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
    assert filecmp.cmp(file1, file2, "Mismatch between: %s and %s" % (file1, file2))


class TestLDPred(unittest.TestCase):
  def test_parse_sum_stats(self):
    with h5py.File('test_output.hdf5', 'w') as h5f:
      bimfile = 'test_data/LDpred_cc_data_p0.001_train_0.bim'
      p_dict = {
          'ssf': 'test_data/LDpred_cc_data_p0.001_ss_0.txt',
          'ssf_format': 'STANDARD',
          'only_hm3': False,
          'N': 10000,
          'debug': True,
          'match_genomic_pos': False,
         }
      sum_stats_parsers.parse_sum_stats(h5f, p_dict, bimfile, {})
      self.assertEqual(len(h5f['sum_stats']['chrom_1']['betas']), 9999)

  def test_ldpred(self):
    tf = tempfile.NamedTemporaryFile()
    tmp_file_prefix = next(tempfile._get_candidate_names())
    try:
        print('Testing LDpred.\n')
        print('Note that this test currently only tests the core functionality of LDpred.')
        print('Please report bugs on github (https://github.com/bvilhjal/ldpred) or to Bjarni J Vilhjalmsson (bjarni.vilhjalmsson@gmail.com).\n')

        coord_file = tmp_file_prefix + '.coord0.hdf5'
        run_test(
            'Coordinating test data into file %s' % coord_file,
            'coord --gf=./test_data/LDpred_data_p0.001_train_0 --vgf=./test_data/LDpred_data_p0.001_test_0 --ssf=./test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --N=10000  --out=%s' % coord_file,
            'Problems when coordinating data!',
            coord_file,
            'test_data/goldens/golden.coord0.hdf5'
        )

        coord_file = tmp_file_prefix + '.coord.hdf5'
        run_test(
            'Coordinating test data into file %s' % coord_file,
            '--debug coord --gf=./test_data/LDpred_data_p0.001_train_0 --vbim=./test_data/LDpred_data_p0.001_test_0.bim --ssf=./test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --beta --N=10000  --out=%s' % coord_file,
            'Problems when coordinating data!',
            coord_file,
            'test_data/goldens/golden.coord.hdf5')

        run_test(
            'Running LDpred-inf with coordinated file prefix: %s ' % tmp_file_prefix,
            '--debug inf --cf=%s  --ldr=100   --ldf=%s  --N=10000  --out=%s' % (coord_file, tmp_file_prefix, tmp_file_prefix),
            'Problems when running LDpred_inf!',
            tmp_file_prefix + '.txt',
            'test_data/goldens/golden.txt',
            tmp_file_prefix + '_ldradius100.pkl.gz',
            'test_data/goldens/golden_ldradius100.pkl.gz')

        run_test(
            'Running LDpred with coordinated file prefix: %s ' % tmp_file_prefix,
            '--debug gibbs --cf=%s  --ldr=100   --ldf=%s  --f=0.001 --N=10000  --out=%s' % (coord_file, tmp_file_prefix, tmp_file_prefix),
            'Problems when running LDpred!',
            tmp_file_prefix + '_LDpred-inf.txt',
            'test_data/goldens/golden_LDpred-inf.txt',
            tmp_file_prefix + '_LDpred_p1.0000e-03.txt',
            'test_data/goldens/golden_LDpred_p1.0000e-03.txt')

        run_test(
            'Running P+T with coordinated file prefix: %s ' % tmp_file_prefix,
            '--debug p+t --cf=%s  --ldr=100  --p=0.001 --out=%s' % (coord_file, tmp_file_prefix),
            'Problems when running P+T!',
            tmp_file_prefix + '_P+T_r0.20_p1.0000e-03.txt',
            'test_data/goldens/golden_P+T_r0.20_p1.0000e-03.txt')

        prs_file_prefix = tmp_file_prefix + 'prs'
        run_test(
            'Validating results with output file prefix: %s' % tmp_file_prefix,
            '--debug score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s' % (tmp_file_prefix, prs_file_prefix),
            'Problems with the validation step!',
            prs_file_prefix + '_LDpred-inf.txt',
            'test_data/goldens/goldenprs_LDpred-inf.txt',
            prs_file_prefix + '_LDpred_p1.0000e-03.txt',
            'test_data/goldens/goldenprs_LDpred_p1.0000e-03.txt')

        prs_file_prefix2 = tmp_file_prefix + 'prs_2'
        run_test(
            'Validating results with output file prefix: %s' % tmp_file_prefix,
            'score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s' % (tmp_file_prefix, prs_file_prefix2),
            'Problems with the validation step!',
            prs_file_prefix2 + '_LDpred-inf.txt',
            'test_data/goldens/goldenprs_2_LDpred-inf.txt',
            prs_file_prefix2 + '_LDpred_p1.0000e-03.txt',
            'test_data/goldens/goldenprs_2_LDpred_p1.0000e-03.txt',)

        run_test(
            'Validating results with output file prefix: %s' % tmp_file_prefix,
            'score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --rf-format=P+T --out=%s' % (tmp_file_prefix, prs_file_prefix),
            'Problems with the P+T validation step!',
            prs_file_prefix + '_P+T_p1.0000e-03.txt',
            'test_data/goldens/goldenprs_2_P+T_p1.0000e-03.txt')

        prs_file_prefix3 = tmp_file_prefix + 'prs_3'
        run_test(
            'Validating results with output file prefix: %s' % tmp_file_prefix,
            'score --gf=./test_data/LDpred_data_p0.001_test_0  --only-score --rf=%s  --rf-format=P+T --out=%s' % (tmp_file_prefix, prs_file_prefix3),
            'Problems with the P+T validation step!',
            prs_file_prefix3 + '_P+T_p1.0000e-03.txt',
            'test_data/goldens/goldenprs_3_P+T_p1.0000e-03.txt')
        print('Test finished successfully!')

    except Exception as e:
        print("Test failed: ",e)
        print('Cleaning up files.')
        cmd_str = 'rm %s*' % tmp_file_prefix
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems cleaning up test files!  Testing stopped'
        raise Exception('Test failed.')

    print('Cleaning up files.')
    cmd_str = 'rm %s*' % tmp_file_prefix
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems cleaning up test files!  Testing stopped'


if __name__ == '__main__':
  unittest.main()
