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
        
        
def prepare_constants(ldpred_n,ns,m,p,h2,sampl_var_shrink_factor):
    Mp = m * p
    hdmp = (h2 / Mp)
    const_dict = {'Mp':Mp, 'hdmp':hdmp}
    rv_scalars = sp.zeros(m)
    if ldpred_n is not None:
        hdmpn = hdmp + 1.0 / ldpred_n
        hdmp_hdmpn = (hdmp / hdmpn)
        c_const = (p / sp.sqrt(hdmpn))
        d_const = (1.0 - p) / (sp.sqrt(1.0 / ldpred_n))
        const_dict['n']=ldpred_n
        const_dict['hdmpn']=hdmpn
        const_dict['hdmp_hdmpn']=hdmp_hdmpn
        const_dict['c_const']=c_const
        const_dict['d_const']=d_const
        rv_scalars[:]=sampl_var_shrink_factor* sp.sqrt((hdmp_hdmpn) * (1.0 / ldpred_n))
    else:
        snp_dict = {}
        for i in range(m):
            ni = ns[i]
            hdmpn_i = hdmp + 1.0 / ni
            hdmp_hdmpn_i = (hdmp / hdmpn_i)
            c_const_i = (p / sp.sqrt(hdmpn_i))
            d_const_i = (1.0 - p) / (sp.sqrt(1.0 / ni))
            snp_dict[i]={'n':ni, 'hdmpn':hdmpn_i, 
                           'hdmp_hdmpn':hdmp_hdmpn_i, 
                           'c_const':c_const_i, 
                           'd_const':d_const_i}
            rv_scalars[i]=sampl_var_shrink_factor* sp.sqrt((hdmp_hdmpn_i) * (1.0 / ni))
        const_dict['snp_dict']=snp_dict
    const_dict['rv_scalars']=rv_scalars
    return const_dict


def get_constants(snp_i,const_dict):
    if 'snp_dict' in const_dict:
        return const_dict['snp_dict'][snp_i]
    else:
        return const_dict
        
            
def ldpred_gibbs(beta_hats, genotypes=None, start_betas=None, h2=None, n=None, ns= None, ld_radius=100,
                 num_iter=60, burn_in=10, p=None, zero_jump_prob=0.01, sampl_var_shrink_factor=0.9, 
                 tight_sampling=False,ld_dict=None, reference_ld_mats=None, ld_boundaries=None, 
                 snp_lrld=None, verbose=False, print_progress=True):
    """
    LDpred (Gibbs Sampler) 
    """
    # Set random seed to stabilize results
    sp.random.seed(42) 

    t0 = time.time()
    m = len(beta_hats)


    ldpred_n, ldpred_inf_n = get_LDpred_sample_size(n,ns,verbose)
    
    # If no starting values for effects were given, then use the infinitesimal model starting values.
    if start_betas is None and verbose:
        print('Initializing LDpred effects with posterior mean LDpred-inf effects.')
        print('Calculating LDpred-inf effects.')
        start_betas = LDpred_inf.ldpred_inf(beta_hats, genotypes=genotypes, reference_ld_mats=reference_ld_mats,
                                            h2=h2, n=ldpred_inf_n, ld_window_size=2 * ld_radius, verbose=False)
    curr_betas = sp.copy(start_betas)
    assert len(curr_betas)==m,'Betas returned by LDpred_inf do not have the same length as expected.'
    curr_post_means = sp.zeros(m)
    avg_betas = sp.zeros(m)

    # Iterating over effect estimates in sequential order
    iter_order = sp.arange(m)
    
    # Setting up the marginal Bayes shrink
    const_dict = prepare_constants(ldpred_n,ns,m,p,h2,sampl_var_shrink_factor)


    for k in range(num_iter):  # Big iteration
        h2_est = max(0.00001, sp.sum(curr_betas ** 2))
        if tight_sampling:
            # Force an alpha shrink if estimates are way off compared to heritability estimates.  
            #(May improve MCMC convergence.)
            alpha = min(1.0 - zero_jump_prob, 1.0 / h2_est, (h2 + 1.0 / sp.sqrt(ldpred_n)) / h2_est)
        else:
            alpha = 1.0 - zero_jump_prob

        rand_ps = sp.random.random(m)
        
        rand_norms = stats.norm.rvs(0.0, 1, size=m)*const_dict['rv_scalars']

        for i, snp_i in enumerate(iter_order):
            if ld_boundaries is None:
                start_i = max(0, snp_i - ld_radius)
                focal_i = min(ld_radius, snp_i)
                stop_i = min(m, snp_i + ld_radius + 1)
            else:
                start_i = ld_boundaries[snp_i][0]
                stop_i = ld_boundaries[snp_i][1]
                focal_i = snp_i - start_i
                
            if snp_lrld is not None:
                if snp_lrld[snp_i]:
                    continue
                
            #Figure out what sample size and constants to use
            cd = get_constants(snp_i,const_dict)

            # Local LD matrix
            D_i = ld_dict[snp_i]
            
            # Local (most recently updated) effect estimates
            local_betas = curr_betas[start_i: stop_i]
            
            # Calculate the local posterior mean, used when sampling.
            local_betas[focal_i] = 0.0
            res_beta_hat_i = beta_hats[snp_i] - sp.dot(D_i , local_betas)
            b2 = res_beta_hat_i ** 2


            d_const_b2_exp = cd['d_const'] * sp.exp(-b2 * cd['n'] / 2.0)
            if sp.isreal(d_const_b2_exp):
                numerator = cd['c_const'] * sp.exp(-b2 / (2.0 * cd['hdmpn']))
                if sp.isreal(numerator):
                    if numerator == 0.0:
                        postp = 0.0
                    else:
                        postp = numerator / (numerator + d_const_b2_exp)
                        assert sp.isreal(postp), 'The posterior mean is not a real number?  Possibly due to problems with summary stats, LD estimates, or parameter settings.' 
                else:
                    postp = 0.0
            else:
                postp = 1.0
            curr_post_means[snp_i] = cd['hdmp_hdmpn'] * postp * res_beta_hat_i

            if rand_ps[i] < postp * alpha:
                # Sample from the posterior Gaussian dist.
                proposed_beta = rand_norms[snp_i] + cd['hdmp_hdmpn'] * res_beta_hat_i

            else:
                # Sample 0
                proposed_beta = 0.0

            curr_betas[snp_i] = proposed_beta  # UPDATE BETA
             
        if verbose and print_progress:
            sys.stdout.write('\r%0.2f%%' % (100.0 * (min(1, float(k + 1) / num_iter))))
            sys.stdout.flush()

        if k >= burn_in:
            avg_betas += curr_post_means  # Averaging over the posterior means instead of samples.
    if verbose and print_progress:
        sys.stdout.write('\r%0.2f%%\n' % (100.0))
        sys.stdout.flush()

    avg_betas = avg_betas / float(num_iter - burn_in)
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('Took %d minutes and %0.2f seconds' % (t / 60, t % 60))
    return {'betas':avg_betas, 'inf_betas':start_betas}


def ldpred_genomewide(data_file=None, ld_radius=None, ld_dict=None, out_file_prefix=None, 
                      summary_dict=None, ps=None,n=None, h2=None, use_gw_h2=False, 
                      sampl_var_shrink_factor=1, incl_long_range_ld=False,
                      num_iter=None, verbose=False, zero_jump_prob=0.01, 
                      burn_in=5):
    """
    Calculate LDpred for a genome
    """    
    print('Applying LDpred with LD radius: %d' % ld_radius)

    df = h5py.File(data_file, 'r')
    has_phenotypes = False
    if 'y' in df:
        y = df['y'][...]  # Phenotype
        num_individs = len(y)
        risk_scores_pval_derived = sp.zeros(num_individs)
        has_phenotypes = True


    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ld_dict = ld_dict['chrom_ld_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
        
    cord_data_g = df['cord_data']
    mean_n = coord_genotypes.get_mean_sample_size(n, cord_data_g)

    #Calculating genome-wide heritability using LD score regression, and partition heritability by chromsomes
    herit_dict = ld.get_chromosome_herits(cord_data_g, ld_scores_dict, mean_n, h2=h2, use_gw_h2=use_gw_h2, 
                                          debug=verbose, summary_dict=summary_dict)


    if herit_dict['gw_h2_ld_score_est']>ld_radius/10.0:
        print ('\033[93m Warning: LD radius seems small in comparison to the average LD score. '
               'Please consider a larger one, or a smaller number of SNPs used in the analysis. \033[0m')

    LDpred_inf_chrom_dict = {}
    print('Calculating LDpred-inf weights')
    for chrom_str in util.chromosomes_list:
        if chrom_str in cord_data_g:
            if verbose:    
                print('Calculating SNP weights for Chromosome %s' % ((chrom_str.split('_'))[1]))           
            g = cord_data_g[chrom_str]

            # Filter monomorphic SNPs
            snp_stds = g['snp_stds_ref'][...]
            snp_stds = snp_stds.flatten()
            ok_snps_filter = snp_stds > 0
            pval_derived_betas = g['betas'][...]
            pval_derived_betas = pval_derived_betas[ok_snps_filter]      
                              
            h2_chrom = herit_dict[chrom_str]['h2']
            start_betas = LDpred_inf.ldpred_inf(pval_derived_betas, genotypes=None, reference_ld_mats=chrom_ref_ld_mats[chrom_str],
                                                h2=h2_chrom, n=mean_n, ld_window_size=2 * ld_radius, verbose=verbose)
            LDpred_inf_chrom_dict[chrom_str] = start_betas

    if not incl_long_range_ld:
        lrld_dict = util.load_lrld_dict()
        num_snps_in_lrld = 0
    
    
    results_dict = {}
    convergence_report = {}
    for p in ps:
        convergence_report[p] = False
        print('Starting LDpred gibbs with f=%0.4f' % p)
        p_str = '%0.4f' % p
        results_dict[p_str] = {}
        
        if out_file_prefix:
            # Preparing output files
            raw_effect_sizes = []
            ldpred_effect_sizes = []
            ldpred_inf_effect_sizes = []
            out_sids = []
            chromosomes = []
            out_positions = []
            out_nts = []
            
        chrom_i = 0
        num_chrom = len(cord_data_g.keys())
        for chrom_str in util.chromosomes_list:
            chrom_i+=1
            if chrom_str in cord_data_g:
                g = cord_data_g[chrom_str]
                
                if out_file_prefix:
                    positions = g['positions'][...]
                    sids = (g['sids'][...]).astype(util.sids_u_dtype)
                    log_odds = g['log_odds'][...]
                    nts = (g['nts'][...]).astype(util.nts_u_dtype)
                    
                    chromosomes.extend([chrom_str] * len(positions))
                    out_positions.extend(positions)
                    out_sids.extend(sids)
                    raw_effect_sizes.extend(log_odds)
                    out_nts.extend(nts)
        
                pval_derived_betas = g['betas'][...]
                ns = g['ns'][...]
                h2_chrom = herit_dict[chrom_str]['h2']
                
                snp_lrld = None
                if not incl_long_range_ld:
                    snp_lrld = util.get_snp_lrld_status(chrom_i, positions, lrld_dict)      
                    num_snps_in_lrld +=sp.sum(snp_lrld)
                    
                ld_boundaries = None
                if 'chrom_ld_boundaries' in ld_dict:
                    ld_boundaries = ld_dict['chrom_ld_boundaries'][chrom_str]
                if verbose:
                    print('Calculating SNP weights for Chromosome %s' % ((chrom_str.split('_'))[1]))
                res_dict = ldpred_gibbs(pval_derived_betas,h2=h2_chrom, n=n, ns=ns, p=p, ld_radius=ld_radius,
                                        verbose=verbose, num_iter=num_iter, burn_in=burn_in, ld_dict=chrom_ld_dict[chrom_str],
                                        start_betas=LDpred_inf_chrom_dict[chrom_str], ld_boundaries=ld_boundaries,
                                        zero_jump_prob=zero_jump_prob,sampl_var_shrink_factor=sampl_var_shrink_factor,
                                        snp_lrld=snp_lrld, print_progress=False)
                
                updated_betas = res_dict['betas']
                updated_inf_betas = res_dict['inf_betas']
                sum_sqr_effects = sp.sum(updated_betas ** 2)
                if sum_sqr_effects > herit_dict['gw_h2_ld_score_est']:
                    print('Sum of squared updated effects estimates seems too large: %0.4f'% sum_sqr_effects)
                    print('This suggests that the Gibbs sampler did not convergence.')
                    convergence_report[p] = True
                
                snp_stds = g['snp_stds_ref'][...]
                snp_stds = snp_stds.flatten()
                updated_betas = updated_betas / snp_stds
                updated_inf_betas = updated_inf_betas / snp_stds
                ldpred_effect_sizes.extend(updated_betas)
                ldpred_inf_effect_sizes.extend(updated_inf_betas)
                
                if not verbose:
                    sys.stdout.write('\r%0.2f%%' % (100.0 * (min(1, float(chrom_i) / num_chrom))))
                    sys.stdout.flush()

                else:
                    if has_phenotypes:
                        if 'raw_snps_val' in g:
                            raw_snps = g['raw_snps_val'][...]
                        else:
                            raw_snps = g['raw_snps_ref'][...]
                        prs = sp.dot(updated_betas, raw_snps)
                        risk_scores_pval_derived += prs
                        corr = sp.corrcoef(y, prs)[0, 1]
                        r2 = corr ** 2
                        print('The R2 prediction accuracy of PRS using %s was: %0.4f' % (chrom_str, r2))

        if not incl_long_range_ld:
            summary_dict[1.3]={'name':'SNPs in long-range LD regions','value':'%d'%num_snps_in_lrld}
        
        if not verbose:
            sys.stdout.write('\r%0.2f%%\n' % (100.0))
            sys.stdout.flush()
        if verbose and has_phenotypes:
            num_indivs = len(y)
            results_dict[p_str]['y'] = y
            results_dict[p_str]['risk_scores_pd'] = risk_scores_pval_derived
            print('Prediction accuracy was assessed using %d individuals.' % (num_indivs))
    
            corr = sp.corrcoef(y, risk_scores_pval_derived)[0, 1]
            r2 = corr ** 2
            results_dict[p_str]['r2_pd'] = r2
            print('The  R2 prediction accuracy (observed scale) for the whole genome was: %0.4f (%0.6f)' % (r2, ((1 - r2) ** 2) / num_indivs))
    
            if corr < 0:
                risk_scores_pval_derived = -1 * risk_scores_pval_derived
            auc = util.calc_auc(y, risk_scores_pval_derived)
            print('AUC for the whole genome was: %0.4f' % auc)
    
            # Now calibration                               
            denominator = sp.dot(risk_scores_pval_derived.T, risk_scores_pval_derived)
            y_norm = (y - sp.mean(y)) / sp.std(y)
            numerator = sp.dot(risk_scores_pval_derived.T, y_norm)
            regression_slope = (numerator / denominator)  # [0][0]
            print('The slope for predictions with P-value derived  effects is: %0.4f' % regression_slope)
            results_dict[p_str]['slope_pd'] = regression_slope
            
        
        weights_out_file = '%s_LDpred_p%0.4e.txt' % (out_file_prefix, p)
        with open(weights_out_file, 'w') as f:
            f.write('chrom    pos    sid    nt1    nt2    raw_beta     ldpred_beta\n')
            for chrom, pos, sid, nt, raw_beta, ldpred_beta in zip(chromosomes, out_positions, out_sids, out_nts, raw_effect_sizes, ldpred_effect_sizes):
                nt1, nt2 = nt[0], nt[1]
                f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n' % (chrom, pos, sid, nt1, nt2, raw_beta, ldpred_beta))

    weights_out_file = '%s_LDpred-inf.txt' % (out_file_prefix)
    with open(weights_out_file, 'w') as f:
        f.write('chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta \n')
        for chrom, pos, sid, nt, raw_beta, ldpred_inf_beta in zip(chromosomes, out_positions, out_sids, out_nts, raw_effect_sizes, ldpred_inf_effect_sizes):
            nt1, nt2 = nt[0], nt[1]
            f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n' % (chrom, pos, sid, nt1, nt2, raw_beta, ldpred_inf_beta))

    summary_dict[2.0]={'name':'Gibbs sampler fractions used','value':str(ps)}
    ['Yes' if convergence_report[p] else 'No' for p in ps]
    summary_dict[2.1]={'name':'Number of burn-iterations used','value':'%i'%burn_in}
    summary_dict[2.2]={'name':'Number of iterations used','value':'%i'%num_iter}
    summary_dict[2.3]={'name':'Convergence issues (for each fraction)','value':str(['Yes' if convergence_report[p] else 'No' for p in ps])}


def main(p_dict):
    #Check parameters
    summary_dict = {}
    summary_dict[0]={'name':'Coordinated data filename','value':p_dict['cf']}
    summary_dict[0.1]={'name':'SNP weights output file (prefix)', 'value':p_dict['out']}
    summary_dict[0.2]={'name':'LD data filename (prefix)', 'value':p_dict['ldf']}
    summary_dict[1.01]={'name':'LD radius used','value':str(p_dict['ldr'])}
    t0 = time.time()
    summary_dict[1]={'name':'dash', 'value':'LD information'}
    ld_dict = ld.get_ld_dict_using_p_dict(p_dict, summary_dict)
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[1.2]={'name':'Running time for calculating LD information:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    t0 = time.time()
    summary_dict[1.9]={'name':'dash', 'value':'LDpred Gibbs sampler'}
    ldpred_genomewide(data_file=p_dict['cf'], out_file_prefix=p_dict['out'], ps=p_dict['f'], ld_radius=p_dict['ldr'],
                      ld_dict=ld_dict, n=p_dict['N'], num_iter=p_dict['n_iter'], burn_in=p_dict['n_burn_in'], 
                      h2=p_dict['h2'], use_gw_h2=p_dict['use_gw_h2'], incl_long_range_ld=p_dict['incl_long_range_ld'], 
                      sampl_var_shrink_factor=1, verbose=p_dict['debug'], summary_dict=summary_dict)
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[3]={'name':'Running time for Gibbs sampler(s):','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    reporting.print_summary(summary_dict, 'Summary of LDpred Gibbs')

        
