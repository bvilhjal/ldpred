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

def test_coord(label='coord', td_prefix='./test_data/sim5_0'):
    print("Testing P+T Workflow")
    
    coord_file = tmp_file_prefix+label + '.coord0.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=%s_train --vgf=%s_test --ssf=%s_ss.txt --ssf-format=LDPRED  --N=8000  --out=%s' %(td_prefix,td_prefix,td_prefix, coord_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'

    coord_file = tmp_file_prefix+label + '.coord1.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --z-from-se --gf=%s_train --vgf=%s_test --ssf=%s_ss.txt --ssf-format=LDPRED  --N=8000  --out=%s' %(td_prefix,td_prefix,td_prefix, coord_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'


def test_simple_pt(label='pt', td_prefix='./test_data/sim2_0'):
    print("Testing P+T Workflow")
    
    coord_file = tmp_file_prefix+label + '.coord0.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=%s_train --vgf=%s_test --ssf=%s_ss.txt --ssf-format=LDPRED  --N=8000  --out=%s'%(td_prefix,td_prefix,td_prefix, coord_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'
    
    weights_file = tmp_file_prefix+label+'.weights'
    print('Running P+T with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug p+t --cf=%s  --ldr=100  --p=0.001 --out=%s' % (coord_file, weights_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running P+T!'

    prs_file_prefix = tmp_file_prefix + label+'.prs'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py score --gf=%s_test  --rf=%s  --rf-format=P+T --out=%s' % (td_prefix, weights_file, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the P+T validation step!'
    print('Test finished successfully!')




def test_LD_pred_inf(label='inf', td_prefix='./test_data/sim1_0'):
    print("Testing LDpred-inf Workflow")
    
    coord_file = tmp_file_prefix + label+'.coord.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=%s_train --vbim=%s_test.bim --ssf=%s_ss.txt --ssf-format=LDPRED  --beta --N=8000  --out=%s' %(td_prefix,td_prefix,td_prefix, coord_file)
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
    cmd_str = 'python LDpred.py --debug score --gf=%s_test  --rf=%s  --out=%s' % (td_prefix, weights_file, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the validation step!'


def test_gibbs(label='gibbs',td_prefix='./test_data/sim2_0'):
    print("Testing LDpred-gibbs Workflow")

    coord_file = tmp_file_prefix + label+'.coord.hdf5'
    print('Coordinating test data into file %s' % coord_file)
    cmd_str = 'python LDpred.py --debug coord --gf=%s_train --vbim=%s_test.bim --ssf=%s_ss.txt --ssf-format=LDPRED  --beta --N=8000  --out=%s' %(td_prefix,td_prefix,td_prefix, coord_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when coordinating data!'

    ld_file = tmp_file_prefix+label+'.ld'
    weights_file = tmp_file_prefix+label+'.weights'
   
    print('Running LDpred with coordinated file prefix: %s ' % tmp_file_prefix)
    cmd_str = 'python LDpred.py --debug gibbs --cf=%s  --ldr=100   --ldf=%s  --f=0.001 --N=10000  --out=%s' % (coord_file, ld_file, weights_file)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems when running LDpred!'
    

    prs_file_prefix = tmp_file_prefix + label+'.prs'
    print('Validating results with output file prefix: %s' % tmp_file_prefix)
    cmd_str = 'python LDpred.py score --gf=%s_test --rf=%s  --out=%s' % (td_prefix, weights_file, prs_file_prefix)
    print(cmd_str + '\n')
    assert os.system(cmd_str) == 0, 'Problems with the validation step!'
    print('Test finished successfully!')


def test_mix(label='mix', td_prefix='./test_data/',num_traits=1):
    print("Testing mixed LDpred Workflow")
    
    for sim_i in range(1,6):
        td = './test_data/sim%d'%sim_i
        for t_i in range(num_traits):
            file_prefix = '%s_%s_sim%d_%d'%(tmp_file_prefix,label,sim_i,t_i)
            df_prefix = '%s_%d'%(td,t_i)
            
            coord_file = file_prefix+'.hdf5'
            print('Coordinating test data into file %s' % coord_file)
            cmd_str = 'python LDpred.py --debug coord --gf=%s_train --vbim=%s_test.bim --ssf=%s_ss.txt --ssf-format=LDPRED --out=%s' % (df_prefix,df_prefix,df_prefix,coord_file)
            print(cmd_str + '\n')
            assert os.system(cmd_str) == 0, 'Problems when coordinating data!'
        
            ld_file = file_prefix+'.ld'
            weights_file = file_prefix+'.weights'
            print('Running LDpred with coordinated file prefix: %s ' % coord_file)
            cmd_str = 'python LDpred.py --debug gibbs --n-burn-in 3 --n-iter 30 --cf=%s  --ldr=100   --ldf=%s  --f 1 0.3 0.1 --out=%s' % (coord_file, ld_file, weights_file)
            print(cmd_str + '\n')
            assert os.system(cmd_str) == 0, 'Problems when running LDpred!'
        
            print('Running P+T with coordinated file prefix: %s ' % coord_file)
            cmd_str = 'python LDpred.py --debug p+t --cf=%s  --ldr=100  --p 1 0.3 0.1 --out=%s' % (coord_file, weights_file)
            print(cmd_str + '\n')
            assert os.system(cmd_str) == 0, 'Problems when running P+T!'
        
            prs_file_prefix = file_prefix+'.prs'
            print('Validating results with output file prefix: %s' % prs_file_prefix)
            cmd_str = 'python LDpred.py --debug score --gf=%s_test  --rf=%s  --out=%s' % (df_prefix, weights_file, prs_file_prefix)
            print(cmd_str + '\n')
            assert os.system(cmd_str) == 0, 'Problems with the validation step!'
            print('Test finished successfully!')
        
            prs_file_prefix2 = file_prefix+'.prs2'
            print('Validating results with output file prefix: %s' % prs_file_prefix2)
            cmd_str = 'python LDpred.py --debug score --gf=%s_test  --only-score --rf=%s  --rf-format=P+T --out=%s' % (df_prefix, weights_file, prs_file_prefix2)
            print(cmd_str + '\n')
            assert os.system(cmd_str) == 0, 'Problems with the P+T validation step!'
            print('Test finished successfully!')

    
try:
    test_coord()
    test_simple_pt()
    test_LD_pred_inf()    
    test_gibbs()
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


