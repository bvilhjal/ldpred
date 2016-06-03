"""
A test file for LDpred.
"""

import sys
import tempfile
import os

tf = tempfile.NamedTemporaryFile()
tmp_file_prefix = next(tempfile._get_candidate_names())

print '\033[1mTesting LDpred.'
print 'Note that this test currently only tests the core functionality of LDpred.'
print 'Please report bugs in bitbucket (https://bitbucket.org/bjarni_vilhjalmsson/ldpred) or to Bjarni J Vilhjalmsson (bjarni.vilhjalmsson@gmail.com).\n'

coord_file = tmp_file_prefix+'.coord.hdf5'
print '\033[93mCoordinating test data into file %s'%coord_file
cmd_str = 'python coord_genotypes.py --gf=test_data/LDpred_data_p0.001_train_0 --vgf=test_data/LDpred_data_p0.001_test_0 --ssf=test_data/LDpred_data_p0.001_ss_0.txt --N=10000  --out=%s'%coord_file
print cmd_str+'\033[0m'
assert os.system(cmd_str)==0, '\033[91mProblems when coordinating data!  Testing stopped'

out_file = tmp_file_prefix+'.res'
print '\033[93mCoordinating test data with LD file and results file prefix: %s '%tmp_file_prefix
cmd_str = 'python LDpred.py --coord=%s  --ld_radius=100   --local_ld_file_prefix=%s  --PS=0.001 --N=10000  --out=%s'%(coord_file,tmp_file_prefix,tmp_file_prefix)
print cmd_str+'\033[0m'
assert os.system(cmd_str)==0, '\033[91mProblems when running LDpred!  Testing stopped'

out_file = tmp_file_prefix+'.res'
print '\033[93mValidating results with output file prefix: %s'%tmp_file_prefix
cmd_str = 'python validate.py --vgf=test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s'%(tmp_file_prefix,tmp_file_prefix)
print cmd_str+'\033[0m'
assert os.system(cmd_str)==0, '\033[91mProblems with the validation step!  Testing stopped'


print '\033[93mCleaning up files.'
cmd_str = 'rm %s*'%tmp_file_prefix
print cmd_str+'\033[0m'
assert os.system(cmd_str)==0, '\033[91mProblems cleaning up test files!  Testing stopped'

print '\033[1mTest finished successfully!\033[0m'