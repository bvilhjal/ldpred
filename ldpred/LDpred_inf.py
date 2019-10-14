import sys
import time
import h5py
import scipy as sp
from scipy import linalg
from ldpred import util
from ldpred import ld

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
    updated_betas = sp.empty(m)

    for i, wi in enumerate(range(0, m, ld_window_size)):
        start_i = wi
        stop_i = min(m, wi + ld_window_size)
        curr_window_size = stop_i - start_i
        if inf_shrink_matrices!=None:
            A_inv = inf_shrink_matrices[i]
        else:
            if reference_ld_mats != None:
                D = reference_ld_mats[i]
            else:
                if genotypes != None:
                    X = genotypes[start_i: stop_i]
                    num_indivs = X.shape[1]
                    D = sp.dot(X, X.T) / num_indivs
                else:
                    raise NotImplementedError
            A = ((m / h2) * sp.eye(curr_window_size) + (n / (1.0)) * D)
            A_inv = linalg.pinv(A)
        updated_betas[start_i: stop_i] = sp.dot(A_inv * n , beta_hats[start_i: stop_i])  # Adjust the beta_hats

        if verbose:
            sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(wi + 1.0) / m))))
            sys.stdout.flush()

    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('\nIt took %d minutes and %0.2f seconds to perform the Infinitesimal LD shrink' % (t / 60, t % 60))
    return updated_betas


def ldpred_inf_genomewide(data_file=None, ld_radius = None, ld_dict=None, out_file_prefix=None,
                          n=None, h2=None, verbose=False):
    """
    Calculate LDpred for a genome
    """    
    
    df = h5py.File(data_file,'r')
    has_phenotypes=False
    if 'y' in df:
        'Validation phenotypes found.'
        y = df['y'][...]  # Phenotype
        num_individs = len(y)
        risk_scores_pval_derived = sp.zeros(num_individs)
        has_phenotypes=True

    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
        
    print('Applying LDpred-inf with LD radius: %d' % ld_radius)
    results_dict = {}
    cord_data_g = df['cord_data']

    #Calculating genome-wide heritability using LD score regression, and partition heritability by chromsomes
    herit_dict = ld.get_chromosome_herits(cord_data_g, ld_scores_dict, n, h2=h2)

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
            
            snp_stds = g['snp_stds_ref'][...]
            pval_derived_betas = g['betas'][...]
            if out_file_prefix:
                chromosomes.extend([chrom_str]*len(pval_derived_betas))
                positions.extend(g['positions'][...])
                sids_arr = (g['sids'][...]).astype(util.sids_u_dtype)
                sids.extend(sids_arr)
                raw_effect_sizes.extend(g['log_odds'][...])
                nts_arr = (g['nts'][...]).astype(util.nts_u_dtype)
                nts.extend(nts_arr)
        
            h2_chrom = herit_dict[chrom_str] 
            updated_betas = ldpred_inf(pval_derived_betas, genotypes=None, reference_ld_mats=chrom_ref_ld_mats[chrom_str], 
                                                h2=h2_chrom, n=n, ld_window_size=2*ld_radius, verbose=False)
                    
            print('Calculating scores for Chromosome %s'%((chrom_str.split('_'))[1]))
            updated_betas = updated_betas / (snp_stds.flatten())
            ldpred_effect_sizes.extend(updated_betas)
            if has_phenotypes:
                prs = sp.dot(updated_betas, raw_snps)
                risk_scores_pval_derived += prs
                corr = sp.corrcoef(y, prs)[0, 1]
                r2 = corr ** 2
                print('The R2 prediction accuracy of PRS using %s was: %0.4f' %(chrom_str, r2))

                
    if has_phenotypes:
        num_indivs = len(y)
        results_dict['y']=y
        results_dict['risk_scores_pd']=risk_scores_pval_derived
        print('Prediction accuracy was assessed using %d individuals.'%(num_indivs))

        corr = sp.corrcoef(y, risk_scores_pval_derived)[0, 1]
        r2 = corr ** 2
        results_dict['r2_pd']=r2
        print('The  R2 prediction accuracy (observed scale) for the whole genome was: %0.4f (%0.6f)' % (r2, ((1-r2)**2)/num_indivs))

        if corr<0:
            risk_scores_pval_derived = -1* risk_scores_pval_derived
        auc = util.calc_auc(y,risk_scores_pval_derived)
        print('AUC for the whole genome was: %0.4f'%auc)

        #Now calibration                               
        denominator = sp.dot(risk_scores_pval_derived.T, risk_scores_pval_derived)
        y_norm = (y-sp.mean(y))/sp.std(y)
        numerator = sp.dot(risk_scores_pval_derived.T, y_norm)
        regression_slope = (numerator / denominator)
        print('The slope for predictions with P-value derived  effects is: %0.4f'%regression_slope)
        results_dict['slope_pd']=regression_slope
    
    weights_out_file = '%s.txt'%(out_file_prefix)
    with open(weights_out_file,'w') as f:
        f.write('chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta\n')
        for chrom, pos, sid, nt, raw_beta, ldpred_beta in zip(chromosomes, positions, sids, nts, raw_effect_sizes, ldpred_effect_sizes):
            nt1,nt2 = nt[0],nt[1]
            f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n'%(chrom, pos, sid, nt1, nt2, raw_beta, ldpred_beta))



def main(p_dict):
    ld_dict = ld.get_ld_dict_using_p_dict(p_dict, summary_dict={})
    
    ldpred_inf_genomewide(data_file=p_dict['cf'], out_file_prefix=p_dict['out'], ld_radius=p_dict['ldr'], 
                          ld_dict = ld_dict, n=p_dict['N'], h2=p_dict['h2'], verbose=p_dict['debug'])
            
                    

