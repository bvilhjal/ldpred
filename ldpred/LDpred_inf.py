import sys
import time
import h5py
import numpy as np
from scipy import linalg
from ldpred import util
from ldpred import ld
from ldpred import reporting
from ldpred import coord_genotypes
 


def ldpred_inf(beta_hats, h2=0.1, n=1000, inf_shrink_matrices=None, 
               reference_ld_mats=None, genotypes=None, ld_window_size=100, verbose=False):
    """
    Apply the infinitesimal shrink w LD (which requires LD information).
    
    If reference_ld_mats are supplied, it uses those, otherwise it uses the LD in the genotype data.
    
    If genotypes are supplied, then it assumes that beta_hats and the genotypes are synchronized.

    """
    n = float(n)
    if verbose:
        print('Doing LD correction')
    t0 = time.time()
    m = len(beta_hats)
    updated_betas = np.empty(m)

    for i, wi in enumerate(range(0, m, ld_window_size)):
        start_i = wi
        stop_i = min(m, wi + ld_window_size)
        curr_window_size = stop_i - start_i
        if inf_shrink_matrices is not None:
            A_inv = inf_shrink_matrices[i]
            updated_betas[start_i: stop_i] = np.dot(A_inv * n , beta_hats[start_i: stop_i])  # Adjust the beta_hats
        else:
            if reference_ld_mats is not None:
                D = reference_ld_mats[i]
            else:
                if genotypes is not None:
                    X = genotypes[start_i: stop_i]
                    num_indivs = X.shape[1]
                    D = np.dot(X, X.T) / num_indivs
                else:
                    raise NotImplementedError
            A = ((m / h2) * np.eye(curr_window_size) + (n / (1.0)) * D)
            # Solve A x = n * beta_hat rather than inverting A.  A is
            # (m/h2) I + n D, so it is symmetric and normally positive
            # definite, and a Cholesky solve is several times cheaper than the
            # SVD behind pinv.  Shrinking D can cost positive definiteness, so
            # fall back to the general symmetric solver and then to pinv.
            b = n * beta_hats[start_i: stop_i]
            try:
                updated_betas[start_i: stop_i] = linalg.solve(A, b, assume_a='pos')
            except linalg.LinAlgError:
                try:
                    updated_betas[start_i: stop_i] = linalg.solve(A, b, assume_a='sym')
                except linalg.LinAlgError:
                    updated_betas[start_i: stop_i] = np.dot(linalg.pinv(A), b)

        if verbose:
            sys.stdout.write('\r%0.2f%%' % (100.0 * (min(1, float(wi + ld_window_size) / m))))
            sys.stdout.flush()

    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('\nIt took %d minutes and %0.2f seconds to perform the Infinitesimal LD shrink' % (t / 60, t % 60))
    return updated_betas


def ldpred_inf_genomewide(data_file=None, ld_radius = None, ld_dict=None, out_file_prefix=None,
                          n=None, h2=None, use_gw_h2=False, verbose=False, summary_dict=None):
    """
    Calculate LDpred for a genome
    """    
    
    df = h5py.File(data_file,'r')
    has_phenotypes = 'y' in df
    if has_phenotypes:
        y = df['y'][...]  # Phenotype
        num_individs = len(y)
        risk_scores_pval_derived = np.zeros(num_individs)

    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
        
    print('Applying LDpred-inf with LD radius: %d' % ld_radius)
    results_dict = {}
    cord_data_g = df['cord_data']

    mean_n = coord_genotypes.get_mean_sample_size(n, cord_data_g)

    #Calculating genome-wide heritability using LD score regression, and partition heritability by chromsomes
    herit_dict = ld.get_chromosome_herits(cord_data_g, ld_scores_dict, mean_n, h2=h2, use_gw_h2=use_gw_h2, 
                                          debug=verbose,summary_dict=summary_dict)

    if out_file_prefix:
        #Preparing output files
        raw_effect_sizes = []
        ldpred_effect_sizes = []
        sids = []
        chromosomes = []
        positions = []
        nts = []
        
    for chrom_str in util.chromosomes_list:
        if chrom_str in cord_data_g:
            g = cord_data_g[chrom_str]
            if has_phenotypes:
                if 'raw_snps_val' in g:
                    raw_snps = g['raw_snps_val'][...]
                else:
                    raw_snps = g['raw_snps_ref'][...]
            
            pval_derived_betas = g['betas'][...]
            if out_file_prefix:
                chromosomes.extend([chrom_str]*len(pval_derived_betas))
                positions.extend(g['positions'][...])
                sids_arr = (g['sids'][...]).astype(util.sids_u_dtype)
                sids.extend(sids_arr)
                raw_effect_sizes.extend(g['log_odds'][...])
                nts_arr = (g['nts'][...]).astype(util.nts_u_dtype)
                nts.extend(nts_arr)
        
            h2_chrom = herit_dict[chrom_str]['h2'] 
            updated_betas = ldpred_inf(pval_derived_betas, genotypes=None, reference_ld_mats=chrom_ref_ld_mats[chrom_str], 
                                                h2=h2_chrom, n=mean_n, ld_window_size=2*ld_radius, verbose=False)
                    
            if verbose:
                print('Calculating scores for Chromosome %s'%((chrom_str.split('_'))[1]))
            snp_stds = g['snp_stds_ref'][...]
            updated_betas = updated_betas / (snp_stds.flatten())
            ldpred_effect_sizes.extend(updated_betas)
            if has_phenotypes:
                prs = np.dot(updated_betas, raw_snps)
                risk_scores_pval_derived += prs
                corr = np.corrcoef(y, prs)[0, 1]
                r2 = corr ** 2
                if verbose:
                    print('The R2 prediction accuracy of PRS using %s was: %0.4f' %(chrom_str, r2))

                
    if has_phenotypes:
        num_indivs = len(y)
        results_dict['y']=y
        results_dict['risk_scores_pd']=risk_scores_pval_derived
        print('Prediction accuracy was assessed using %d individuals.'%(num_indivs))

        corr = np.corrcoef(y, risk_scores_pval_derived)[0, 1]
        r2 = corr ** 2
        results_dict['r2_pd']=r2
        print('The  R2 prediction accuracy (observed scale) for the whole genome was: %0.4f (%0.6f)' % (r2, ((1-r2)**2)/num_indivs))

        if corr<0:
            risk_scores_pval_derived = -1* risk_scores_pval_derived
        auc = util.calc_auc(y,risk_scores_pval_derived)
        print('AUC for the whole genome was: %0.4f'%auc)

        #Now calibration                               
        denominator = np.dot(risk_scores_pval_derived.T, risk_scores_pval_derived)
        y_norm = (y-np.mean(y))/np.std(y)
        numerator = np.dot(risk_scores_pval_derived.T, y_norm)
        regression_slope = (numerator / denominator)
        print('The slope for predictions with P-value derived  effects is: %0.4f'%regression_slope)
        results_dict['slope_pd']=regression_slope
    
    weights_out_file = '%s_LDpred-inf.txt' % (out_file_prefix)
    with open(weights_out_file,'w') as f:
        f.write('chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta\n')
        for chrom, pos, sid, nt, raw_beta, ldpred_beta in zip(chromosomes, positions, sids, nts, raw_effect_sizes, ldpred_effect_sizes):
            nt1,nt2 = nt[0],nt[1]
            f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n'%(chrom, pos, sid, nt1, nt2, raw_beta, ldpred_beta))



def main(p_dict):
    
    summary_dict = {}
    summary_dict[0]={'name':'Coordinated data filename','value':p_dict['cf']}
    summary_dict[0.1]={'name':'SNP weights output file (prefix)', 'value':p_dict['out']}
    summary_dict[0.2]={'name':'LD data filename (prefix)', 'value':p_dict['ldf']}
    summary_dict[1]={'name':'LD radius used','value':str(p_dict['ldr'])}
    t0 = time.time()
    summary_dict[1.09]={'name':'dash', 'value':'LD information'}
    ld_dict = ld.get_ld_dict_using_p_dict(p_dict, summary_dict={})
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[1.2]={'name':'Running time for calculating LD information:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    t0 = time.time()
    summary_dict[1.9]={'name':'dash', 'value':'LDpred infinitesimal model'}

    ldpred_inf_genomewide(data_file=p_dict['cf'], out_file_prefix=p_dict['out'], ld_radius=p_dict['ldr'], 
                          ld_dict = ld_dict, n=p_dict['N'], h2=p_dict['h2'], use_gw_h2=p_dict['use_gw_h2'],
                          verbose=p_dict['debug'], summary_dict=summary_dict)
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[2.2]={'name':'Running time for LDpred-inf:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    reporting.print_summary(summary_dict, 'Summary of LDpred-inf')
    
        
            
                    

