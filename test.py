"""
A test file for LDpred.
"""

import tempfile
import os    

tf = tempfile.NamedTemporaryFile()
tmp_file_prefix = next(tempfile._get_candidate_names())
try:
    print('Testing LDpred.\n')
    print('Note that this test currently only tests the core functionality of LDpred.')
    print('Please report bugs on github (https://github.com/bvilhjal/ldpred) or to Bjarni J Vilhjalmsson (bjarni.vilhjalmsson@gmail.com).\n')
    
    coord_file = tmp_file_prefix + '.coord0.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=./test_data/LDpred_data_p0.001_train_0 --vgf=./test_data/LDpred_data_p0.001_test_0 --ssf=./test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --N=10000  --out=%s' % coord_file
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'
    
    coord_file = tmp_file_prefix + '.coord.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=./test_data/LDpred_data_p0.001_train_0 --vbim=./test_data/LDpred_data_p0.001_test_0.bim --ssf=./test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --beta --N=10000  --out=%s' % coord_file
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'

    print('Running LDpred-inf with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug inf --cf=%s  --ldr=100   --ldf=%s  --N=10000  --out=%s' % (coord_file, tmp_file_prefix, tmp_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running LDpred_inf!'
    
    print('Running LDpred with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug gibbs --cf=%s  --ldr=100   --ldf=%s  --f=0.001 --N=10000  --out=%s' % (coord_file, tmp_file_prefix, tmp_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running LDpred!'
    
    print('Running P+T with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug p+t --cf=%s  --ldr=100  --p=0.001 --out=%s' % (coord_file, tmp_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running P+T!'
    
    prs_file_prefix = tmp_file_prefix + 'prs'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s' % (tmp_file_prefix, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the validation step!'
    
    prs_file_prefix2 = tmp_file_prefix + 'prs_2'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s' % (tmp_file_prefix, prs_file_prefix2)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the validation step!'
    print('Test finished successfully!')

    prs_file_prefix = tmp_file_prefix + 'prs_P+T'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --rf-format=P+T --out=%s' % (tmp_file_prefix, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the P+T validation step!'
    print('Test finished successfully!')


except Exception as e:
    print("Test failed: ",e)

print('Cleaning up files.')
cmd_str = 'rm %s*' % tmp_file_prefix
print(cmd_str + '\n')
assert os.system(cmd_str) == 0, 'Problems cleaning up test files!  Testing stopped'

