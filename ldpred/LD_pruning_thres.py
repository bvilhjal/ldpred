import h5py
import scipy as sp
import time
from ldpred import ld
from ldpred import util

def smart_ld_pruning(beta_hats, ld_table, pvalues=None, max_ld=0.2, verbose=False):
    """
    Smart LD pruning.
    """    
    if verbose:
        print('Doing smart LD pruning')
    t0 = time.time()
    if pvalues is not None:
        pruning_vector = ld.smart_ld_pruning(pvalues, ld_table, max_ld=max_ld, verbose=verbose)    
    else:
        pruning_vector = ld.smart_ld_pruning(beta_hats ** 2, ld_table, max_ld=max_ld, verbose=verbose, reverse=True)    

    if verbose:
        if sp.sum(pruning_vector) == 0:
            print('No SNPs, skipping chromosome')
    shrunk_betas = beta_hats * pruning_vector
    
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('\nIt took %d minutes and %0.2f seconds to perform the LD shrink' % (t / 60, t % 60))
    return shrunk_betas, pruning_vector



def ld_pruning(data_file=None, ld_radius = None, out_file_prefix=None, p_thres=None, 
              verbose=False, max_r2s=[1,0.2]):
    """
    LD pruning + P-value thresholding 
    """    
    df = h5py.File(data_file,'r')
    has_phenotypes=False
    if 'y' in list(df.keys()):
        print('Validation phenotypes found.')
        y = df['y'][...]  # Phenotype
        num_individs = len(y)
        has_phenotypes=True
        
    for max_r2 in max_r2s:
        if has_phenotypes:
            risk_scores = sp.zeros(num_individs)    
            
        print('')
        if max_r2<1:
            print('Applying LD-pruning + P-value thresholding with p-value threshold of %0.2e, a LD radius of %d SNPs, and a max r2 of %0.2f' %(p_thres, ld_radius, max_r2))
        else:
            if p_thres<1:
                print('Applying P-value thresholding with p-value threshold of %0.2e' %(p_thres))
            else:
                print('Calculating polygenic risk score using all SNPs')        
        results_dict = {}
        num_snps = 0
        cord_data_g = df['cord_data']
    
        chromsomes = []
        for chrom_str in list(cord_data_g.keys()):
            g = cord_data_g[chrom_str]
            betas = g['betas'][...]
            n_snps = len(betas)
            num_snps += n_snps
            chromsomes.append(int((chrom_str.split('_'))[1]))
            
        chromsomes.sort()
        p_str = '%0.4f'%p_thres
        results_dict[p_str]={}
    
        if out_file_prefix:
            #Preparing output files
            raw_effect_sizes = []
            raw_pval_effect_sizes = []
            updated_effect_sizes = []
            updated_pval_effect_sizes = []
            sids = []
            chromosomes = []
            positions = []
            nts = []
            
        tot_num_snps = 0
        num_snps_used = 0
        for chrom in chromsomes:
            chrom_str = 'chrom_%d'%chrom
            g = cord_data_g[chrom_str]
            pvalues = g['ps'][...]
            snp_filter = pvalues < p_thres
            num_snps = sp.sum(snp_filter)
            if num_snps == 0:
                continue
            tot_num_snps += num_snps
    
            pvalues = pvalues[snp_filter]
            if 'raw_snps_val' in list(g.keys()):
                raw_snps = g['raw_snps_val'][...][snp_filter]
    
            else:
                raw_snps = g['raw_snps_ref'][...][snp_filter]
                    
            snp_means = g['snp_means_ref'][...][snp_filter]
            snp_stds = g['snp_stds_ref'][...][snp_filter]
            raw_betas = g['log_odds'][...][snp_filter]
            pval_derived_betas = g['betas'][...][snp_filter]
            if out_file_prefix:
                chromosomes.extend([chrom_str]*len(pval_derived_betas))
                positions.extend(g['positions'][...][snp_filter])
                sids_arr = (g['sids'][...]).astype(util.sids_u_dtype)
                sids.extend(sids_arr[snp_filter])
                raw_effect_sizes.extend(raw_betas)
                raw_pval_effect_sizes.extend(pval_derived_betas)
                nts_arr = (g['nts'][...]).astype(util.nts_u_dtype)
                nts.extend(nts_arr[snp_filter])
    
            if max_r2<1:
                snp_means.shape = (len(snp_means),1)   
                snp_stds.shape = (len(snp_means),1)   
                #Normalize SNPs..
                norm_ref_snps = sp.array((raw_snps - snp_means)/snp_stds,dtype='float32') 
                ld_table = ld.calc_ld_table(norm_ref_snps, max_ld_dist=ld_radius, min_r2=max_r2, verbose=verbose)
                
                updated_raw_betas, pruning_vector = smart_ld_pruning(raw_betas, ld_table, pvalues=pvalues, max_ld=max_r2, verbose=verbose)
                updated_pval_derived_betas = pval_derived_betas * pruning_vector
                num_snps_used += sp.sum(pruning_vector)
            else:
                updated_raw_betas = sp.copy(raw_betas)
                updated_pval_derived_betas = sp.copy(pval_derived_betas) 
                updated_pval_derived_betas = updated_pval_derived_betas / (snp_stds.flatten())
                pruning_vector = sp.ones(len(pval_derived_betas))
                num_snps_used += sp.sum(pruning_vector)
    
            if out_file_prefix:
                updated_effect_sizes.extend(updated_raw_betas)
                updated_pval_effect_sizes.extend(updated_pval_derived_betas)
    
                    
            if has_phenotypes:
                print('Calculating scores for Chromosome %s'%chrom_str) 
                prs = sp.dot(updated_raw_betas, raw_snps)
                risk_scores += prs
                corr = sp.corrcoef(y, prs)[0, 1]
                r2 = corr ** 2
                print('The R2 prediction accuracy of PRS using %s was: %0.4f' %(chrom_str, r2))
    
                    
        print('There were %d (SNP) effects after p-value thresholding' % tot_num_snps)
        print('After LD-pruning %d SNPs had non-zero effects'%num_snps_used)
        if has_phenotypes:
            results_dict[p_str]['y']=y
            results_dict[p_str]['risk_scores']=risk_scores
            print('Prediction accuracy was assessed using %d individuals.'%(num_individs))
    
            corr = sp.corrcoef(y, risk_scores)[0, 1]
            r2 = corr ** 2
            results_dict[p_str]['r2_pd']=r2
            print('The  R2 prediction accuracy (observed scale) for the whole genome was: %0.4f (%0.6f)' % (r2, ((1-r2)**2)/num_individs))
    
            if corr<0:
                risk_scores = -1* risk_scores
    
            #Now calibration                               
            denominator = sp.dot(risk_scores.T, risk_scores)
            y_norm = (y-sp.mean(y))/sp.std(y)
            numerator = sp.dot(risk_scores.T, y_norm)
            regression_slope = (numerator / denominator)
            print('The slope for predictions with P-value derived  effects is: %0.4f' %regression_slope)
            results_dict[p_str]['slope_pd']=regression_slope
        
        weights_out_file = '%s_P+T_r%0.2f_p%0.4e.txt'%(out_file_prefix, max_r2, p_thres)
        
        with open(weights_out_file,'w') as f:
            f.write('chrom    pos    sid    nt1    nt2    raw_beta    raw_pval_beta    updated_beta    updated_pval_beta \n')
            for chrom, pos, sid, nt, raw_beta, raw_pval_beta, upd_beta, upd_pval_beta in zip(chromosomes, positions, sids, nts, raw_effect_sizes, raw_pval_effect_sizes, updated_effect_sizes, updated_pval_effect_sizes):
                nt1,nt2 = nt[0],nt[1]
                f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e    %0.4e    %0.4e\n'%(chrom, pos, sid, nt1, nt2, raw_beta, raw_pval_beta, upd_beta, upd_pval_beta))
    df.close()

def main(p_dict):
    for p_thres in reversed(p_dict['p']):
        ld_pruning(data_file=p_dict['cf'], out_file_prefix=p_dict['out'], p_thres=p_thres, ld_radius=p_dict['ldr'],
                   max_r2s=p_dict['r2'])

    