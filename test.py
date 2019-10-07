"""
A test file for LDpred.
"""

import tempfile
import os    

tf = tempfile.NamedTemporaryFile()
tmp_file_prefix = next(tempfile._get_candidate_names())

print('Testing LDpred.\n')
print('Note that this test currently only tests the core functionality of LDpred.')
print('Please report bugs on github (https://github.com/bvilhjal/ldpred) or to Bjarni J Vilhjalmsson (bjarni.vilhjalmsson@gmail.com).\n')

def test_coord(label='coord'):
    print("Testing P+T Workflow")
    
    coord_file = tmp_file_prefix+label + '.coord0.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=./test_data/LDpred_cc_data_p0.001_train_0 --vgf=./test_data/LDpred_cc_data_p0.001_test_0 --ssf=./test_data/LDpred_cc_data_p0.001_ss_0.txt --ssf-format=STANDARD  --N=8000  --out=%s' % coord_file
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'

    coord_file = tmp_file_prefix+label + '.coord1.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --z-from-se --gf=./test_data/LDpred_cc_data_p0.001_train_0 --vgf=./test_data/LDpred_cc_data_p0.001_test_0 --ssf=./test_data/LDpred_cc_data_p0.001_ss_0.txt --ssf-format=STANDARD  --N=8000  --out=%s' % coord_file
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'


def test_simple_pt(label='pt'):
    print("Testing P+T Workflow")
    
    coord_file = tmp_file_prefix+label + '.coord0.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=./test_data/LDpred_data_p0.001_train_0 --vgf=./test_data/LDpred_data_p0.001_test_0 --ssf=./test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --N=8000  --out=%s' % coord_file
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'
    
    weights_file = tmp_file_prefix+label+'.weights'
    print('Running P+T with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug p+t --cf=%s  --ldr=100  --p=0.001 --out=%s' % (coord_file, weights_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running P+T!'

    prs_file_prefix = tmp_file_prefix + label+'.prs'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --rf-format=P+T --out=%s' % (weights_file, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the P+T validation step!'
    print('Test finished successfully!')

    # Expected output:
    # =============================== Scoring Summary ================================
    # Validation genotype file (prefix):                
    #                                            ./test_data/LDpred_data_p0.001_test_0
    # Input weight file(s) (prefix):                                            OjrZ35
    # Output scores file(s) (prefix):                                        OjrZ35prs
    # ---------------------------------- Phenotypes ----------------------------------
    # Phenotype file (plink format):                    
    #                                        ./test_data/LDpred_data_p0.001_test_0.fam
    # Individuals with phenotype information:                                     2000
    # Running time for parsing phenotypes:                         0 min and 0.00 secs
    # ----------------------------------- Scoring ------------------------------------
    # Best P+T (r2=0.20, p=1.00e-03) (unadjusted) R2:                           0.0675
    # Running time for calculating scores:                         0 min and 0.09 secs
    # --------------------------- Optimal polygenic score ----------------------------
    # Method with highest (unadjusted) Pearson R2:                     P+T_p1.0000e-03
    # Best (unadjusted) Pearson R2:                                             0.0675
    # ================================================================================




def test_LD_pred_inf(label='inf'):
    print("Testing LDpred-inf Workflow")
    
    coord_file = tmp_file_prefix + label+'.coord.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=./test_data/LDpred_data_p0.001_train_0 --vbim=./test_data/LDpred_data_p0.001_test_0.bim --ssf=./test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --beta --N=8000  --out=%s' % coord_file
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'

    ld_file = tmp_file_prefix+label+'.ld'
    weights_file = tmp_file_prefix+label+'.weights'
    print('Running LDpred-inf with coordinated file prefix: %s ' % coord_file)
    cmd_str = 'python LDpred.py --debug inf --cf=%s  --ldr=100   --ldf=%s  --N=10000  --out=%s' % (coord_file, ld_file, weights_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running LDpred_inf!'

    prs_file_prefix = tmp_file_prefix + label+'.prs'
    print('Validating results with output file prefix: %s' % weights_file)
    cmd_str = 'python LDpred.py --debug score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s' % (weights_file, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the validation step!'


def test_gibbs(label='gibbs'):
    print("Testing LDpred-gibbs Workflow")

    coord_file = tmp_file_prefix + label+'.coord.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=./test_data/LDpred_data_p0.001_train_0 --vbim=./test_data/LDpred_data_p0.001_test_0.bim --ssf=./test_data/LDpred_data_p0.001_ss_0.txt --ssf-format=STANDARD  --beta --N=8000  --out=%s' % coord_file
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'

    ld_file = tmp_file_prefix+label+'.ld'
    weights_file = tmp_file_prefix+label+'.weights'
   
    print('Running LDpred with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug gibbs --cf=%s  --ldr=100   --ldf=%s  --f=0.001 --N=10000  --out=%s' % (coord_file, ld_file, weights_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running LDpred!'
    
    # Expected output (approximately):
    # =========================== Summary of LDpred Gibbs ============================
    # Coordinated data filename                                      OjrZ35.coord.hdf5
    # SNP weights output file (prefix)                                          OjrZ35
    # LD data filename (prefix)                                                 OjrZ35
    # LD radius used                                                               100
    # -------------------------------- LD information --------------------------------
    # Genome-wide (LDscore) estimated heritability:                             0.0741
    # Chi-square lambda (inflation statistic).                                  1.7022
    # Running time for calculating LD information:                 0 min and 1.17 secs
    # ----------------------------- LDpred Gibbs sampler -----------------------------
    # Gibbs sampler fractions used                                             [0.001]
    # Convergence issues (for each fraction)                                    ['No']
    # Running time for Gibbs sampler(s):                           0 min and 6.62 secs
    # ================================================================================

    prs_file_prefix = tmp_file_prefix + label+'.prs'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py score --gf=./test_data/LDpred_data_p0.001_test_0  --rf=%s  --out=%s' % (weights_file, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the validation step!'
    print('Test finished successfully!')

    # Expected output (approximately):
    # =============================== Scoring Summary ================================
    # Validation genotype file (prefix):                
    #                                            ./test_data/LDpred_data_p0.001_test_0
    # Input weight file(s) (prefix):                                            OjrZ35
    # Output scores file(s) (prefix):                                      OjrZ35prs_2
    # ---------------------------------- Phenotypes ----------------------------------
    # Phenotype file (plink format):                    
    #                                        ./test_data/LDpred_data_p0.001_test_0.fam
    # Individuals with phenotype information:                                     2000
    # Running time for parsing phenotypes:                         0 min and 0.00 secs
    # ----------------------------------- Scoring ------------------------------------
    # LDpred_inf (unadjusted) Pearson R2:                                       0.0379
    # Best LDpred (f=1.00e-03) (unadjusted) R2:                                 0.0864
    # Best P+T (r2=0.20, p=1.00e-03) (unadjusted) R2:                           0.0675
    # Running time for calculating scores:                         0 min and 3.12 secs
    # --------------------------- Optimal polygenic score ----------------------------
    # Method with highest (unadjusted) Pearson R2:                  LDpred_p1.0000e-03
    # Best (unadjusted) Pearson R2:                                             0.0864
    # ================================================================================


def test_mix(label='mix', td='./test_data/LDpred_data_p0.10'):
    print("Testing mixed LDpred Workflow")

    coord_file = tmp_file_prefix + label+'.coord.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=%s_train_0 --vbim=%s_test_0.bim --ssf=%s_ss_0.txt --ssf-format=LDPRED --out=%s' % (td,td,td,coord_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'

    ld_file = tmp_file_prefix+label+'.ld'
    weights_file = tmp_file_prefix+label+'.weights'
   
    print('Running LDpred with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug gibbs --cf=%s  --ldr=100   --ldf=%s  --f 1 0.3 0.1 --N=8000  --out=%s' % (coord_file, ld_file, weights_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running LDpred!'

    print('Running P+T with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug p+t --cf=%s  --ldr=100  --p 1 0.3 0.1 --out=%s' % (coord_file, weights_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running P+T!'

    prs_file_prefix = tmp_file_prefix + label+'.prs'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug score --gf=%s_test_0  --rf=%s  --out=%s' % (td, weights_file, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the validation step!'
    print('Test finished successfully!')

    prs_file_prefix2 = tmp_file_prefix + label+'.prs2'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug score --gf=%s_test_0  --only-score --rf=%s  --rf-format=P+T --out=%s' % (td, weights_file, prs_file_prefix2)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the P+T validation step!'
    print('Test finished successfully!')

    
try:
    #test_coord()
    #test_simple_pt()
    #test_LD_pred_inf()    
    #test_gibbs()
    test_mix()


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


