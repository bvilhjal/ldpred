"""
Implementation of Joel Mefford's BLUP-based shrink (Mefford, PhD thesis 2018).

"""


import sys
import time

import scipy as sp
from scipy import stats
import h5py

from ldpred import LDpred_inf
from ldpred import util
from ldpred import ld
from ldpred import reporting
from ldpred import coord_genotypes




def get_LDpred_sample_size(n,ns,verbose):
    if n is None:
        #If coefficient of variation is very small, then use one N nevertheless.
        n_cv = sp.std(ns)/sp.mean(ns)
        if n_cv<0.01:
            ldpred_n = sp.mean(ns)
            if verbose:
                print ("Sample size does not vary much (CV=%0.4f).  Using a fixed sample size of %0.2f"%(n_cv,ldpred_n))
        else:
            if verbose:
                print ("Using varying sample sizes")
                print ("Sample size ranges between %d and %d"%(min(ns),max(ns)))
                print ("Average sample size is %0.2f "%(sp.mean(ns)))
        ldpred_inf_n = sp.mean(ns)
        ldpred_n = None
    else:
        ldpred_n = float(n)
        if verbose:
            print ("Using the given fixed sample size of %d"%(n))
        ldpred_inf_n = float(n)
    return ldpred_n,ldpred_inf_n
        
        
def ldpred_fast(betas_blup, h2=None, Ns=None, p_c=None):
    """
    LDpred-fast
    
    Joel Mefford's sparsified BLUP shrink. 
    """
    M = len(betas_blup)
    h2mn = h2+M/Ns
    h2mpn = h2+(M*p_c)/Ns
    C1 = h2mn/h2mpn
    C2 = h2/h2mn
    V_nc = (C2/Ns)*(1-C2*h2)
    V_c = C2*(h2/(M*p_c))+V_nc
    P_beta_ncaus = stats.norm.pdf(betas_blup,scale=sp.sqrt(V_nc))
    P_beta_caus = stats.norm.pdf(betas_blup,scale=sp.sqrt(V_c))
    P_caus_beta =  P_beta_caus *p_c/(P_beta_caus *p_c + P_beta_ncaus *(1-p_c))
    beta_jms = betas_blup * C1 * P_caus_beta 
    return beta_jms


def get_herit_dict(df, eff_type, n, h2, use_gw_h2= False, ld_dict=None, summary_dict=None, debug=False):
    if eff_type != 'BLUP':
        ld_scores_dict = ld_dict['ld_scores_dict']
        cord_data_g = df['cord_data']
    
        #Calculating genome-wide heritability using LD score regression, and partition heritability by chromsomes
        return ld.get_chromosome_herits(cord_data_g, ld_scores_dict, n, h2=h2, use_gw_h2=use_gw_h2, summary_dict=summary_dict, debug=debug)
    else:
        raise NotImplementedError
    
def init_res_dict(out_file_prefix,cord_data_g):
    results_cdict = {}
    
    for chrom_str in util.chromosomes_list:
        if chrom_str in cord_data_g:
            g = cord_data_g[chrom_str]
            positions = g['positions'][...]
            
            rd = {'raw_betas': g['log_odds'][...], 
                  'out_sids': (g['sids'][...]).astype(util.sids_u_dtype),
                  'chromosomes': [chrom_str] * len(positions), 
                  'out_positions': positions, 
                  'out_nts': (g['nts'][...]).astype(util.nts_u_dtype), 
                  'ldpred_fast_betas_dict':{}}
            results_cdict[chrom_str] = rd
    return results_cdict
    
    
def write_res_dict_2_file(out_file_prefix, results_cdict, ps):
    for p_c in ps:
        p_str = '%0.4f' % p_c
        weights_out_file = '%s_LDpred_fast_p%0.4e.txt' % (out_file_prefix, p_c)
        with open(weights_out_file, 'w') as f:
            f.write('chrom    pos    sid    nt1    nt2    raw_beta     ldpred_beta\n')
            for chrom_str in util.chromosomes_list:
                if chrom_str in results_cdict:
                    rd = results_cdict[chrom_str]
                    for chrom, pos, sid, nt, raw_beta, ldpred_beta in zip(rd['chromosomes'], 
                                                                          rd['out_positions'], 
                                                                          rd['out_sids'], 
                                                                          rd['out_nts'], 
                                                                          rd['raw_betas'], 
                                                                          rd['ldpred_fast_betas_dict'][p_str]):
                        f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n' % (chrom, pos, sid, nt[0], nt[1], raw_beta, ldpred_beta))

        
def ldpred_fast_genomewide(data_file=None, ld_radius=None, ld_dict=None, out_file_prefix=None, 
                           ps=None,n=None, h2=None, use_gw_h2=False, 
                           eff_type=None, summary_dict=None, debug=False):
    """
    Calculate LDpred for a genome
    """    
    print('Applying LDpred-fast')

    df = h5py.File(data_file, 'r')
    cord_data_g = df['cord_data']
    mean_n = coord_genotypes.get_mean_sample_size(n, cord_data_g)
    herit_dict = get_herit_dict(df, eff_type, mean_n, h2, ld_dict=ld_dict, summary_dict=summary_dict, debug=debug)
    results_cdict = init_res_dict(out_file_prefix,cord_data_g)
    blup_betas_chrom_dict = {}

    if eff_type != 'BLUP':
        chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
        print('Calculating LDpred-inf weights')
        for chrom_str in util.chromosomes_list:
            if chrom_str in cord_data_g:
                if debug:    
                    print('Calculating LDpred-inf weights for Chromosome %s' % ((chrom_str.split('_'))[1]))           
                g = cord_data_g[chrom_str]
    
                # Filter monomorphic SNPs
                snp_stds = g['snp_stds_ref'][...]
                snp_stds = snp_stds.flatten()
                ok_snps_filter = snp_stds > 0
                pval_derived_betas = g['betas'][...]
                pval_derived_betas = pval_derived_betas[ok_snps_filter]      
                                  
                h2_chrom = herit_dict[chrom_str]['h2']
                
                blup_betas = LDpred_inf.ldpred_inf(pval_derived_betas, genotypes=None, reference_ld_mats=chrom_ref_ld_mats[chrom_str],
                                                    h2=h2_chrom, n=mean_n, ld_window_size=2 * ld_radius, verbose=debug)

                snp_stds = g['snp_stds_ref'][...]
                snp_stds = snp_stds.flatten()
                ldpred_inf_betas = blup_betas / snp_stds
                results_cdict[chrom_str]['ldpred_inf_betas'] = ldpred_inf_betas
                blup_betas_chrom_dict[chrom_str]=blup_betas
                
    else:
        raise NotImplementedError
        
    
    chrom_i = 0
    num_chrom = len(cord_data_g.keys())
    for chrom_str in util.chromosomes_list:
        chrom_i+=1
        if chrom_str in cord_data_g:
            g = cord_data_g[chrom_str]
            h2_chrom = herit_dict[chrom_str]['h2']
            ns = g['ns'][...]
            blup_betas = blup_betas_chrom_dict[chrom_str]

            if debug:
                print('Calculating LDpred-fast weights for Chromosome %s' % ((chrom_str.split('_'))[1]))
            
            for p_c in ps:
                p_str = '%0.4f' % p_c
                ldpred_fast_betas = ldpred_fast(blup_betas,h2=h2_chrom, Ns=ns, p_c=p_c)
                results_cdict[chrom_str]['ldpred_fast_betas_dict'][p_str] = ldpred_fast_betas / snp_stds

            if not debug:
                sys.stdout.write('\r%0.2f%%' % (100.0 * (min(1, float(chrom_i) / num_chrom))))
                sys.stdout.flush()
    
    write_res_dict_2_file(out_file_prefix, results_cdict, ps)

    if not debug:
        sys.stdout.write('\r%0.2f%%\n' % (100.0))
        sys.stdout.flush()
        
    summary_dict[2.0]={'name':'LDpred-fast fractions used','value':str(ps)}

def get_eff_type(cf):
    df = h5py.File(cf, 'r')
    eff_type = df['parameters']['eff_type'][...]
    df.close()
    return eff_type



def main(p_dict):
    #Check parameters
    summary_dict = {}
    summary_dict[0]={'name':'Coordinated data filename','value':p_dict['cf']}
    summary_dict[0.1]={'name':'SNP weights output file (prefix)', 'value':p_dict['out']}
    
    eff_type = get_eff_type(p_dict['cf'])
    #If already BLUP betas, then skip LD calculation
    if eff_type!='BLUP':
        summary_dict[0.2]={'name':'LD data filename (prefix)', 'value':p_dict['ldf']}
        summary_dict[1.01]={'name':'LD radius used','value':str(p_dict['ldr'])}
        summary_dict[1]={'name':'dash', 'value':'LD information'}
        t0 = time.time()
        ld_dict = ld.get_ld_dict_using_p_dict(p_dict, summary_dict)
        t1 = time.time()
        t = (t1 - t0)
        summary_dict[1.2]={'name':'Running time for calculating LD information:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
        t0 = time.time()
    
    summary_dict[1.9]={'name':'dash', 'value':'LDpred-fast'}
    ldpred_fast_genomewide(data_file=p_dict['cf'], out_file_prefix=p_dict['out'], ps=p_dict['f'], ld_radius=p_dict['ldr'],
                      ld_dict=ld_dict, n=p_dict['N'], h2=p_dict['h2'], use_gw_h2=p_dict['use_gw_h2'], eff_type = eff_type,
                      summary_dict=summary_dict, debug=p_dict['debug'],)
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[3]={'name':'Running time for LDpred-fast:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    reporting.print_summary(summary_dict, 'Summary of LDpred-fast')

        

        
