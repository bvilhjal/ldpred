import h5py
import scipy as sp
import time
import sys
from ldpred import ld
from ldpred import reporting
from ldpred import util

def ld_clumping(beta_hats, ld_table, pruning_stat =None, pvalues=None, max_ld=0.2, verbose=False):
    """
    Smart LD pruning.  (Also known as LD clumping)
    """    
    if verbose:
        print('LD clumping')
    t0 = time.time()
    if pruning_stat is not None: 
        pruning_vector = ld.smart_ld_pruning(pruning_stat, ld_table, max_ld=max_ld, verbose=verbose, reverse=True)    
    elif pvalues is not None:
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
                
                updated_raw_betas, pruning_vector = ld_clumping(raw_betas, ld_table, pruning_stat=pval_derived_betas**2, 
                                                                max_ld=max_r2, verbose=verbose)
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
    
        if verbose:            
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

def run_pt(p_dict,summary_dict={}):
    pi=0
    for p_thres in reversed(p_dict['p']):
        max_r2s=p_dict['r2']
        if p_dict['debug']:
            print('')
            if len(max_r2s)>1 or max_r2s[0]<1:
                r2s_string = ', '.join(['%0.2f'%r2 for r2 in max_r2s])
                print('Applying LD-pruning + P-value thresholding with p-value threshold of %0.2e, a LD radius of %d SNPs, and a max r2(s) of %s' %(p_thres, p_dict['ldr'], r2s_string))
            elif p_thres<1:
                print('Applying P-value thresholding with p-value threshold of %0.2e' %(p_thres))
            else:
                print('Calculating polygenic risk score using all SNPs')        
        else:
            sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (float(pi) / len(p_dict['p']))))
            sys.stdout.flush()

        ld_pruning(data_file=p_dict['cf'], out_file_prefix=p_dict['out'], p_thres=p_thres, ld_radius=p_dict['ldr'],
                   max_r2s=max_r2s)
        pi+=1
    if not p_dict['debug']:
        sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%\n' % (100.0))
        sys.stdout.flush()
    summary_dict[1.2]={'name':'P-value thresholds used','value':str(p_dict['p'])}
    summary_dict[1.3]={'name':'Maximum r2 LD thresholds used','value':str(max_r2s)}
   

def main(p_dict):
    summary_dict = {}
    summary_dict[0]={'name':'Coordinated data filename','value':p_dict['cf']}
    summary_dict[0.1]={'name':'SNP weights output file (prefix)', 'value':p_dict['out']}
    summary_dict[1]={'name':'LD radius used','value':str(p_dict['ldr'])}
    t0 = time.time()
    summary_dict[1.1]={'name':'dash', 'value':'LD-pruning + Thresholding'}
    run_pt(p_dict,summary_dict)
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[2]={'name':'Running time for calculating P+T:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    
    reporting.print_summary(summary_dict, 'Summary of LD-pruning + Tresholding')

        
            


    