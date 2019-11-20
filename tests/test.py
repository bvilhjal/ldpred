"""
A test file for LDpred.

Examples
--------
To run all tests:
$ python -m tests.test

To run a specific test:
$ python -m unittest tests.test.SimpleTests.test_ldpred_inf
"""

import pickle
import filecmp
import gzip
import h5py
from ldpred import coord_genotypes
from ldpred import ld
from ldpred import sum_stats_parsers
from ldpred import run
import numpy as np
import os
import tempfile
import unittest
import sys

np.set_printoptions(linewidth=int(os.environ.get('COLUMNS', 100)))
TEST_DIR = os.path.dirname(os.path.abspath(__file__))

def run_test(mesg, cmd_str, error_mesg, *actual_and_golden_outputs):
    print(mesg)
    print(cmd_str + '\n')
    cmd_args = cmd_str.split()
    try:
        run.main_with_args(cmd_args)
        for i in range(0, len(actual_and_golden_outputs), 2):
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
    for k, v in h5_node.items():
        v_type = type(v)
        v_path = key_prefix + '/' + k
        if v_type == h5py.Group:
            for nested_key, nested_value in h5_node_walker(v, v_path):
                yield nested_key, nested_value
        elif v_type == h5py.Dataset:
            yield v_path, v[...]
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
    try:
        with gzip.open(pkl_file) as f:
            pkl_root_node = pickle.load(f)
    except UnicodeDecodeError as e:
        with gzip.open(pkl_file) as f:
            pkl_root_node = pickle.load(f,encoding='latin1')
    except Exception as e:
        print('Unable to load data ', pkl_file, ':', e)
        raise
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
                assert np.array_equal(v1, v2), 'ndarray non-number mismatch in key %s: v1=%s ; v2=%s' % (k1,str(v1),str(v2))


def assert_files_equal(file1, file2):
    if file1.endswith('.hdf5'):
        assert_deep_equals(h5_file_walker(file1), h5_file_walker(file2))
    elif file1.endswith('.pkl.gz'):
        assert_deep_equals(pkl_file_walker(file1), pkl_file_walker(file2))
    else:
        assert filecmp.cmp(file1, file2), "Mismatch between: %s and %s" % (file1, file2)

def make_p_dict(*args):
    return vars(run.parser.parse_args(args))

class SimpleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        print('Testing LDpred.\n')
        print('Note that this test currently only tests the core functionality of LDpred.')
        print('Please report bugs on github (https://github.com/bvilhjal/ldpred) or to Bjarni J Vilhjalmsson (bjarni.vilhjalmsson@gmail.com).\n')

    def setUp(self):
        self.tf = tempfile.NamedTemporaryFile()
        self.tmp_file_prefix = next(tempfile._get_candidate_names())

    def tearDown(self):
        print('Cleaning up files: %s* ' % self.tmp_file_prefix)
        cmd_str = 'rm -f %s*' % self.tmp_file_prefix
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems cleaning up test files!  Testing stopped'

    def test_parse_sum_stats(self):
        p_dict = {
            'ssf': os.path.join(TEST_DIR, 'test_data/sim1_0_ss.txt'),
            'ssf_format': 'LDPRED',
            'only_hm3': False,
            'N': 10000,
            'debug': True,
            'z_from_se':False,
            'match_genomic_pos': False,
            'eff_type':'LINREG'}
        bimfile = os.path.join(TEST_DIR, 'test_data/sim1_0_test.bim')
        summary_dict = {}
        out = '%s_parse_sum_stats.hdf5' % self.tmp_file_prefix
        with h5py.File(out, 'w') as h5f:
            sum_stats_parsers.parse_sum_stats(h5f, p_dict, bimfile, summary_dict)
            self.assertEqual(len(h5f['sum_stats']['chrom_1']['betas']), 2000)


        p_dict = {
            'ssf': os.path.join(TEST_DIR, 'test_data/sim4_0_ss.txt'),
            'ssf_format': 'LDPRED',
            'only_hm3': False,
            'N': None,
            'debug': True,
            'z_from_se':True,
            'match_genomic_pos': False,}
        bimfile = os.path.join(TEST_DIR, 'test_data/sim4_0_test.bim')
        summary_dict = {}
        out = '%s_parse_sum_stats.hdf5' % self.tmp_file_prefix
        with h5py.File(out, 'w') as h5f:
            sum_stats_parsers.parse_sum_stats(h5f, p_dict, bimfile, summary_dict)
            self.assertEqual(len(h5f['sum_stats']['chrom_1']['betas']), 2000)



    def test_coord_genotypes(self):
        p_dict = make_p_dict(
            '--debug',
            'coord',
            '--gf=%s/test_data/sim1_0_test' % TEST_DIR,
            '--vgf=%s/test_data/sim1_0_test' % TEST_DIR,
            '--ssf=%s/test_data/sim1_0_ss.txt' % TEST_DIR,
            '--ssf-format=LDPRED',
            '--out=%s_coord_genotypes.hdf5' % self.tmp_file_prefix,
        )
        summary_dict = coord_genotypes.main(p_dict)
        # summary_dict[11]['value'], if present, is the count of non-matching nts.
        # It should be 0.
        self.assertEqual(summary_dict.get(11, {}).get('value', 0), 0)
        with h5py.File(p_dict['out'], 'r') as h5f:
            self.assertEqual(len(h5f['sum_stats']['chrom_1']['betas']), 2000)



    def test_ld_calculation(self):
        df = h5py.File('%s/test_data/goldens/golden.coord0.hdf5' % TEST_DIR, 'r')
        g = df['cord_data']['chrom_1']
        snps, n_raw_snps, n_snps = ld.extract_snps_from_cord_data_chrom(g)
        first_10_snps = snps[:10]
        self.assertEqual(len(first_10_snps), 10)
        ld_dict_and_scores = ld.get_LDpred_ld_tables(first_10_snps)
        ld_dict = ld_dict_and_scores['ld_dict']
        ld_mat = np.vstack([ld_dict[i] for i in range(10)])
#         np.savez(os.path.join(TEST_DIR, 'test_data/goldens/ld_data'),ld=ld_mat)
        golden_ld_mat = np.load(os.path.join(TEST_DIR, 'test_data/goldens/ld_data.npz'))['ld']
        self.assertTrue(np.allclose(ld_mat, golden_ld_mat))

    def test_get_chromosome_herits(self):
        p_dict = make_p_dict(
            '--debug',
            'inf',
            '--cf=%s/test_data/goldens/golden.coord.hdf5' % TEST_DIR,
            '--ldr=100',
            '--ldf=' + self.tmp_file_prefix,
            '--N=4000',
            '--out=' + self.tmp_file_prefix,
        )
        summary_dict = {}
        ld_dict = ld.get_ld_dict_using_p_dict(p_dict, summary_dict)
        coord_file = os.path.join(TEST_DIR, 'test_data/goldens/golden.coord.hdf5')
        df = h5py.File(coord_file, 'r')
        herit_dict = ld.get_chromosome_herits(df['cord_data'], ld_dict['ld_scores_dict'], n=p_dict['N'])
        print(herit_dict)
        self.assertAlmostEqual(herit_dict['chrom_1']['h2'], 0.10640501626651437)
        self.assertAlmostEqual(herit_dict['gw_h2_ld_score_est'], 0.10640501626651437)

    def test_ldpred_coord0(self):
        coord_file = self.tmp_file_prefix + '.coord0.hdf5'
        run_test(
            'Coordinating test data into file %s' % coord_file,
            'coord --gf=%s/test_data/sim1_0_test --vgf=%s/test_data/sim1_0_test --ssf=%s/test_data/sim1_0_ss.txt --ssf-format=LDPRED --eff_type LINREG --out=%s' % (TEST_DIR, TEST_DIR, TEST_DIR, coord_file),
            'Problems when coordinating data!',
            coord_file,
            'test_data/goldens/golden.coord0.hdf5'
        )

    def test_ldpred_coord(self):
        coord_file = self.tmp_file_prefix + '.coord.hdf5'
        run_test(
            'Coordinating test data into file %s' % coord_file,
            '--debug coord --gf=%s/test_data/sim2_0_test --vbim=%s/test_data/sim2_0_test.bim --ssf=%s/test_data/sim2_0_ss.txt --ssf-format=LDPRED  --eff_type LINREG  --out=%s' % (TEST_DIR, TEST_DIR, TEST_DIR, coord_file),
            'Problems when coordinating data!',
            coord_file,
            'test_data/goldens/golden.coord.hdf5')

    def test_ldpred_inf(self):
        run_test(
            'Running LDpred-inf with output file prefix: %s ' % self.tmp_file_prefix,
            '--debug inf --cf=%s/test_data/goldens/golden.coord.hdf5  --ldr=100   --ldf=%s  --out=%s' % (TEST_DIR, self.tmp_file_prefix, self.tmp_file_prefix),
            'Problems when running LDpred_inf!',
            self.tmp_file_prefix + '_ldradius100.pkl.gz',
            'test_data/goldens/golden_inf_ldradius100.pkl.gz')

    def test_ldpred_fast(self):
        run_test(
            'Running LDpred-inf with output file prefix: %s ' % self.tmp_file_prefix,
            '--debug fast --cf=%s/test_data/goldens/golden.coord.hdf5 --f 0.3 0.1 0.03 0.01 --ldr=100   --ldf=%s  --out=%s' % (TEST_DIR, self.tmp_file_prefix, self.tmp_file_prefix),
            'Problems when running LDpred_fast!')

    def test_ldpred_gibbs(self):
        run_test(
            'Running LDpred with output file prefix: %s ' % self.tmp_file_prefix,
            '--debug gibbs --cf=%s/test_data/goldens/golden.coord.hdf5  --ldr=100   --ldf=%s  --f=0.001  --out=%s' % (TEST_DIR, self.tmp_file_prefix, self.tmp_file_prefix),
            'Problems when running LDpred!')

    def test_ldpred_p_plus_t(self):
        run_test(
            'Running P+T with coordinated file prefix: %s ' % self.tmp_file_prefix,
            '--debug p+t --cf=%s/test_data/goldens/golden.coord.hdf5  --ldr=100  --p=0.001 --out=%s' % (TEST_DIR, self.tmp_file_prefix),
            'Problems when running P+T!',
            self.tmp_file_prefix + '_P+T_r0.20_p1.0000e-03.txt',
            'test_data/goldens/golden_P+T_r0.20_p1.0000e-03.txt')

    def test_ldpred_score_1(self):
        prs_file_prefix = self.tmp_file_prefix 
        run_test(
            'Validating results with output file prefix: %s' % prs_file_prefix,
            '--debug score --gf=%s/test_data/sim2_0_test  --rf=%s/test_data/goldens/golden  --out=%s' % (TEST_DIR, TEST_DIR, prs_file_prefix),
            'Problems with the validation step!',
            prs_file_prefix + '_LDpred-inf.txt',
            'test_data/goldens/goldenprs_LDpred-inf.txt',
            prs_file_prefix + '_LDpred_p1.0000e-03.txt',
            'test_data/goldens/goldenprs_LDpred_p1.0000e-03.txt')

    def test_ldpred_score_2(self):
        prs_file_prefix = self.tmp_file_prefix 
        run_test(
            'Validating results with output file prefix: %s' % self.tmp_file_prefix,
            'score --gf=%s/test_data/sim2_0_test  --rf-format LDPRED --rf=%s/test_data/goldens/golden  --out=%s' % (TEST_DIR, TEST_DIR, prs_file_prefix),
            'Problems with the validation step!',
            prs_file_prefix + '_LDpred-inf.txt',
            'test_data/goldens/goldenprs_LDpred-inf.txt',
            prs_file_prefix + '_LDpred_p1.0000e-03.txt',
            'test_data/goldens/goldenprs_LDpred_p1.0000e-03.txt',)

    def test_ldpred_score_3(self):
        prs_file_prefix = self.tmp_file_prefix
        run_test(
            'Validating results with output file prefix: %s' % self.tmp_file_prefix,
            'score --gf=%s/test_data/sim2_0_test  --only-score --rf=%s/test_data/goldens/golden  --rf-format=P+T --out=%s' % (TEST_DIR, TEST_DIR, prs_file_prefix),
            'Problems with the P+T validation step!',
            prs_file_prefix + '_P+T_r0.20_p1.0000e-03.txt',
            'test_data/goldens/goldenprs_only_score_P+T_r0.20_p1.0000e-03.txt')

    def test_ldpred_score_4(self):
        prs_file_prefix = self.tmp_file_prefix
        run_test(
            'Validating results with output file prefix: %s' % self.tmp_file_prefix,
            'score --gf=%s/test_data/sim2_0_test  --rf=%s/test_data/goldens/golden  --rf-format=P+T --out=%s' % (TEST_DIR, TEST_DIR, prs_file_prefix),
            'Problems with the P+T validation step!',
            prs_file_prefix + '_P+T_r0.20_p1.0000e-03.txt',
            'test_data/goldens/goldenprs_P+T_r0.20_p1.0000e-03.txt')


class ComplexTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        print('Testing LDpred: Integration tests.\n')
        print('Note that this test currently only tests the core functionality of LDpred.')
        print('Please report bugs on github (https://github.com/bvilhjal/ldpred) or to Bjarni J Vilhjalmsson (bjarni.vilhjalmsson@gmail.com).\n')

    def setUp(self):
        self.tf = tempfile.NamedTemporaryFile()
        self.tmp_file_prefix = next(tempfile._get_candidate_names())

    def tearDown(self):
        print('Cleaning up files: %s* ' % self.tmp_file_prefix)
        cmd_str = 'rm -f %s*' % self.tmp_file_prefix
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems cleaning up test files!  Testing stopped'

    def test_mix1(self):
        t_i = 0
        label='mix1'
        for sim_i in range(1,6):
            td = '%s/test_data/sim%d'%(TEST_DIR, sim_i)
    
            file_prefix = '%s_%s_sim%d_%d'%(self.tmp_file_prefix,label,sim_i,t_i)
            df_prefix = '%s_%d'%(td,t_i)
            coord_file = file_prefix+'.hdf5'
             
            run_test(
                'Validating results with output file prefix: %s' % self.tmp_file_prefix,
                'coord --gf=%s_test --vbim=%s_test.bim --ssf=%s_ss.txt --ssf-format=LDPRED --out=%s' % (df_prefix,df_prefix,df_prefix,coord_file),
                'Problems when coordinating data!')
    
            ld_file = file_prefix+'.ld'
            weights_file = file_prefix+'.weights'
            
            run_test(
                'Running LDpred-fast with coordinated file prefix: %s ' % coord_file,
                '--debug fast --cf=%s --f 0.3 0.1 0.03 0.01 --ldr=100   --ldf=%s  --out=%s' % (coord_file, ld_file, weights_file),
                'Problems when running LDpred_fast!')

            run_test(
                'Running LDpred with coordinated file prefix: %s ' % coord_file,
                'gibbs --N 5500 --use-gw-h2  --n-burn-in 5 --n-iter 50 --cf=%s  --ldr=100   --ldf=%s  --f 1 0.3 0.1 --out=%s' % (coord_file, ld_file, weights_file),
                'Problems when running LDpred!')
    
            run_test(
                'Running P+T with coordinated file prefix: %s ' % coord_file,
                'p+t --cf=%s  --ldr=100  --p 1 0.3 0.1 --out=%s' % (coord_file, weights_file),
                'Problems when running P+T!')
    

            prs_file_prefix = file_prefix+'.prs'
            golden_prs_prefix = '%s/test_data/goldens/golden_%s_prs_%i_%i'%(TEST_DIR,label,sim_i,t_i)
            golden_summary_file = '%s.summary.txt'%golden_prs_prefix
            summary_file = file_prefix+'.summary.txt'
            run_test(
                'Validating results with output file prefix: %s' % prs_file_prefix,
                'score --gf=%s_test  --rf=%s  --out=%s --summary-file=%s' % (df_prefix, weights_file, prs_file_prefix, summary_file),
                'Problems with the validation step!',
                summary_file,golden_summary_file)
                
                
    def test_mix2(self):
        t_i = 0
        label='mix2'
        for sim_i in range(1,6):
            td = '%s/test_data/sim%d'%(TEST_DIR, sim_i)
    
            file_prefix = '%s_%s_sim%d_%d'%(self.tmp_file_prefix,label,sim_i,t_i)
            df_prefix = '%s_%d'%(td,t_i)
            coord_file = file_prefix+'.hdf5'
             
            run_test(
                'Validating results with output file prefix: %s' % self.tmp_file_prefix,
                'coord --gf=%s_test --vbim=%s_test.bim --z-from-se --ssf=%s_ss.txt --ssf-format=LDPRED --out=%s' % (df_prefix,df_prefix,df_prefix,coord_file),
                'Problems when coordinating data!')
    
            ld_file = file_prefix+'.ld'
            weights_file = file_prefix+'.weights'
            run_test(
                'Running LDpred-fast with coordinated file prefix: %s ' % coord_file,
                '--debug fast --cf=%s --f 0.3 0.1 0.03 0.01 0.001 --ldr=150   --ldf=%s  --out=%s' % (coord_file, ld_file, weights_file),
                'Problems when running LDpred_fast!')

            
            run_test(
                'Running LDpred with coordinated file prefix: %s ' % coord_file,
                'gibbs --n-burn-in 5 --n-iter 50 --cf=%s  --ldr=150   --ldf=%s  --f 1 0.1 0.01 0.001 --out=%s' % (coord_file, ld_file, weights_file),
                'Problems when running LDpred!')
    
         
            run_test(
                'Running P+T with coordinated file prefix: %s ' % coord_file,
                'p+t --cf=%s  --ldr=150  --r2 0.5 0.2 0.1 --p 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 0.00001 --out=%s' % (coord_file, weights_file),
                'Problems when running P+T!')
    
         
            prs_file_prefix = file_prefix+'.prs'
            golden_prs_prefix = '%s/test_data/goldens/golden_%s_prs_%i_%i'%(TEST_DIR,label,sim_i,t_i)
            golden_summary_file = '%s.summary.txt'%golden_prs_prefix
            summary_file = file_prefix+'.summary.txt'
            run_test(
                'Validating results with output file prefix: %s' % prs_file_prefix,
                'score --gf=%s_test  --r2 0.5 0.2 0.1 --p 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 0.00001 --rf=%s  --out=%s --summary-file=%s' % (df_prefix, weights_file, 
                                                                                              prs_file_prefix, summary_file),
                'Problems with the validation step!', 
                summary_file,golden_summary_file)
                

def update_golden_files_mix1():
    label = 'mix1'
    tf = tempfile.NamedTemporaryFile()
    tmp_file_prefix = next(tempfile._get_candidate_names())
    
    for sim_i in range(1,6):
        print('Updating golden results')
        coord_file = '%s_%i_coord.hdf5'%(tmp_file_prefix,sim_i)

        cmd_str = 'python -m ldpred --debug coord --gf %s/test_data/sim%i_0_test --vbim %s/test_data/sim%i_0_test.bim --ssf %s/test_data/sim%i_0_ss.txt --ssf-format LDPRED --out=%s' % (TEST_DIR,sim_i,TEST_DIR,sim_i,TEST_DIR,sim_i,coord_file)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'

        weights_prefix = '%s_%i_weights'%(tmp_file_prefix,sim_i)
        ld_prefix = '%s_%i'%(tmp_file_prefix,sim_i)

        cmd_str = 'python -m ldpred fast --cf %s  --ldr 100  --f 0.3 0.1 0.03 0.01 --ldf %s  --out %s' % (coord_file,ld_prefix,weights_prefix)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'

        cmd_str = 'python -m ldpred gibbs --N 5500 --use-gw-h2  --n-burn-in 5 --n-iter 50 --cf %s  --ldr 100   --ldf %s  --f 1 0.3 0.1 --out %s' % (coord_file,ld_prefix,weights_prefix)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'
   
        cmd_str = 'python -m ldpred p+t --cf %s  --ldr 100  --p 1 0.3 0.1 --out %s' % (coord_file,weights_prefix)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'

        prs_prefix = '%s_prs_%i_0'%(tmp_file_prefix,sim_i)
        golden_summary_file = '%s/test_data/goldens/golden_%s_prs_%i_0.summary.txt'%(TEST_DIR, label,sim_i) 
        cmd_str = 'python -m ldpred --debug score --gf %s/test_data/sim%i_0_test --rf %s  --out %s --summary-file %s' % (TEST_DIR, sim_i,weights_prefix,prs_prefix, golden_summary_file)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'

        print('Cleaning up files.')
        cmd_str = 'rm %s*' % tmp_file_prefix
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems cleaning up test files!  Testing stopped'

def update_golden_files_mix2():
    label = 'mix2'
    tf = tempfile.NamedTemporaryFile()
    tmp_file_prefix = next(tempfile._get_candidate_names())
    for sim_i in range(1,6):
        print('Updating golden results')
        coord_file = '%s_%i_coord.hdf5'%(tmp_file_prefix,sim_i)

        cmd_str = 'python -m ldpred coord --gf %s/test_data/sim%i_0_test --vbim %s/test_data/sim%i_0_test.bim --z-from-se --ssf %s/test_data/sim%i_0_ss.txt --ssf-format LDPRED --out=%s' % (TEST_DIR,sim_i,TEST_DIR,sim_i,TEST_DIR,sim_i,coord_file)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'

        weights_prefix = '%s_%i_weights'%(tmp_file_prefix,sim_i)
        ld_prefix = '%s_%i'%(tmp_file_prefix,sim_i)
        
        cmd_str = 'python -m ldpred fast --cf %s  --ldr 150  --f 0.3 0.1 0.03 0.01 0.001 --ldf %s  --out %s' % (coord_file,ld_prefix,weights_prefix)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'

        cmd_str = 'python -m ldpred gibbs --n-burn-in 5 --n-iter 50 --cf %s  --ldr 150   --ldf %s  --f 1 0.1 0.01 0.001 --out %s' % (coord_file,ld_prefix,weights_prefix)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'
   
        cmd_str = 'python -m ldpred p+t --cf %s  --ldr 150  --r2 0.5 0.2 0.1 --p 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 0.00001 --out %s' % (coord_file,weights_prefix)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'

        prs_prefix = '%s_prs_%i_0'%(tmp_file_prefix,sim_i)
        golden_summary_file = '%s/test_data/goldens/golden_%s_prs_%i_0.summary.txt'%(TEST_DIR, label,sim_i) 
        cmd_str = 'python -m ldpred score --gf %s/test_data/sim%i_0_test --r2 0.5 0.2 0.1 --p 1 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001 0.00001 --rf %s  --out %s --summary-file %s' % (TEST_DIR,sim_i,weights_prefix,prs_prefix, golden_summary_file)
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems when updating golden files'

        print('Cleaning up files.')
        cmd_str = 'rm %s*' % tmp_file_prefix
        print(cmd_str + '\n')
        assert os.system(cmd_str) == 0, 'Problems cleaning up test files!  Testing stopped'



def run_integration_tests():
    complex_suite = unittest.TestLoader().loadTestsFromTestCase(ComplexTests)
    unittest.TextTestRunner().run(complex_suite)
    
def run_unit_tests():
    simple_suite = unittest.TestLoader().loadTestsFromTestCase(SimpleTests)
    unittest.TextTestRunner().run(simple_suite)

if __name__ == '__main__':
    unittest.main()
    


