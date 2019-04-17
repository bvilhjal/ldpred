"""
A test file for LDpred.
"""

import LDpred
import tempfile
import os    

def run_test(mesg, cmd_str, error_mesg):
  print(mesg)
  print(cmd_str + '\n')
  cmd_args = cmd_str.split()
  try:
    LDpred.main_with_args(cmd_args)
  except:
    print(error_mesg + '\n')
    raise

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
        'Problems when coordinating data!')
    
    coord_file = tmp_file_prefix + '.coord.hdf5'
    run_test(
        'Coordinating test data into file %s' % coord_file,
        '--debug coord --gf=./test_data/LDpred_data_p0.001_train_0 --vbim=./test_data/LDpred_data_p0.001_test_0.bim --ssf=./test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --beta --N=10000  --out=%s' % coord_file,
        'Problems when coordinating data!')

    run_test(
        'Running LDpred-inf with coordinated file prefix: %s ' % tmp_file_prefix,
        '--debug inf --cf=%s  --ldr=100   --ldf=%s  --N=10000  --out=%s' % (coord_file, tmp_file_prefix, tmp_file_prefix),
        'Problems when running LDpred_inf!')
    
    run_test(
        'Running LDpred with coordinated file prefix: %s ' % tmp_file_prefix,
        '--debug gibbs --cf=%s  --ldr=100   --ldf=%s  --f=0.001 --N=10000  --out=%s' % (coord_file, tmp_file_prefix, tmp_file_prefix),
        'Problems when running LDpred!')
    
    run_test(
        'Running P+T with coordinated file prefix: %s ' % tmp_file_prefix,
        '--debug p+t --cf=%s  --ldr=100  --p=0.001 --out=%s' % (coord_file, tmp_file_prefix),
        'Problems when running P+T!')
    
    prs_file_prefix = tmp_file_prefix + 'prs'
    run_test(
        'Validating results with output file prefix: %s' % tmp_file_prefix,
        '--debug score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s' % (tmp_file_prefix, prs_file_prefix),
        'Problems with the validation step!')
    
    prs_file_prefix2 = tmp_file_prefix + 'prs_2'
    run_test(
        'Validating results with output file prefix: %s' % tmp_file_prefix,
        'score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s' % (tmp_file_prefix, prs_file_prefix2),
        'Problems with the validation step!')

    run_test(
        'Validating results with output file prefix: %s' % tmp_file_prefix,
        'score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --rf-format=P+T --out=%s' % (tmp_file_prefix, prs_file_prefix),
        'Problems with the P+T validation step!')

    prs_file_prefix3 = tmp_file_prefix + 'prs_3'
    run_test(
        'Validating results with output file prefix: %s' % tmp_file_prefix,
        'score --gf=./test_data/LDpred_data_p0.001_test_0  --only-score --rf=%s  --rf-format=P+T --out=%s' % (tmp_file_prefix, prs_file_prefix3),
        'Problems with the P+T validation step!')
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

