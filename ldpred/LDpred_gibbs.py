import sys
import time

import scipy as sp
from scipy import stats
import h5py

from ldpred import LDpred_inf
from ldpred import util
from ldpred import ld
from ldpred import reporting



        
def ldpred_gibbs(beta_hats, genotypes=None, start_betas=None, h2=None, n=1000, ld_radius=100,
                 num_iter=60, burn_in=10, p=None, zero_jump_prob=0.05, tight_sampling=False,
                 ld_dict=None, reference_ld_mats=None, ld_boundaries=None, verbose=False,
                 print_progress=True):
    """
    LDpred (Gibbs Sampler) 
    """
    t0 = time.time()
    m = len(beta_hats)
    n = float(n)
    
    # If no starting values for effects were given, then use the infinitesimal model starting values.
    if start_betas is None and verbose:
        print('Initializing LDpred effects with posterior mean LDpred-inf effects.')
        print('Calculating LDpred-inf effects.')
        start_betas = LDpred_inf.ldpred_inf(beta_hats, genotypes=genotypes, reference_ld_mats=reference_ld_mats,
                                            h2=h2, n=n, ld_window_size=2 * ld_radius, verbose=False)
    curr_betas = sp.copy(start_betas)
    assert len(curr_betas)==m,'Betas returned by LDpred_inf do not have the same length as expected.'
    curr_post_means = sp.zeros(m)
    avg_betas = sp.zeros(m)

    # Iterating over effect estimates in sequential order
    iter_order = sp.arange(m)
    
    # Setting up the marginal Bayes shrink
    Mp = m * p
    hdmp = (h2 / Mp)
    hdmpn = hdmp + 1.0 / n
    hdmp_hdmpn = (hdmp / hdmpn)
    c_const = (p / sp.sqrt(hdmpn))
    d_const = (1.0 - p) / (sp.sqrt(1.0 / n))

    for k in range(num_iter):  # Big iteration

        # Force an alpha shrink if estimates are way off compared to heritability estimates.  (Improves MCMC convergence.)
        h2_est = max(0.00001, sp.sum(curr_betas ** 2))
        if tight_sampling:
            alpha = min(1.0 - zero_jump_prob, 1.0 / h2_est, (h2 + 1.0 / sp.sqrt(n)) / h2_est)
        else:
            alpha = 1.0 - zero_jump_prob

        rand_ps = sp.random.random(m)
        rand_norms = stats.norm.rvs(0.0, (hdmp_hdmpn) * (1.0 / n), size=m)

        if ld_boundaries is None:
            for i, snp_i in enumerate(iter_order):
                start_i = max(0, snp_i - ld_radius)
                focal_i = min(ld_radius, snp_i)
                stop_i = min(m, snp_i + ld_radius + 1)
                
                # Local LD matrix
                D_i = ld_dict[snp_i]
                
                # Local (most recently updated) effect estimates
                local_betas = curr_betas[start_i: stop_i]
                
                # Calculate the local posterior mean, used when sampling.
                local_betas[focal_i] = 0.0
                res_beta_hat_i = beta_hats[snp_i] - sp.dot(D_i , local_betas)
                b2 = res_beta_hat_i ** 2
    
                d_const_b2_exp = d_const * sp.exp(-b2 * n / 2.0)
                if sp.isreal(d_const_b2_exp):
                    numerator = c_const * sp.exp(-b2 / (2.0 * hdmpn))
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
                curr_post_means[snp_i] = hdmp_hdmpn * postp * res_beta_hat_i
    
                if rand_ps[i] < postp * alpha:
                    # Sample from the posterior Gaussian dist.
                    proposed_beta = rand_norms[i] + hdmp_hdmpn * res_beta_hat_i
    
                else:
                    # Sample 0
                    proposed_beta = 0.0
    
                curr_betas[snp_i] = proposed_beta  # UPDATE BETA
        else:
            for i, snp_i in enumerate(iter_order):
                start_i = ld_boundaries[snp_i][0]
                stop_i = ld_boundaries[snp_i][1]
                focal_i = snp_i - start_i
                
                # Local LD matrix
                D_i = ld_dict[snp_i]
                
                # Local (most recently updated) effect imates
                local_betas = curr_betas[start_i: stop_i]
                
                # Calculate the local posterior mean, used when sampling.
                local_betas[focal_i] = 0.0
                res_beta_hat_i = beta_hats[snp_i] - sp.dot(D_i , local_betas)
                b2 = res_beta_hat_i ** 2
    
                d_const_b2_exp = d_const * sp.exp(-b2 * n / 2.0)
                if sp.isreal(d_const_b2_exp):
                    numerator = c_const * sp.exp(-b2 / (2.0 * hdmpn))
                    if sp.isreal(numerator):
                        if numerator == 0.0:
                            postp = 0.0
                        else:
                            postp = numerator / (numerator + d_const_b2_exp)
                            assert sp.isreal(postp), 'Posterior mean is not a real number? Possibly due to problems with summary stats, LD estimates, or parameter settings.' 
                    else:
                        postp = 0.0
                else:
                    postp = 1.0
                curr_post_means[snp_i] = hdmp_hdmpn * postp * res_beta_hat_i
    
                if rand_ps[i] < postp * alpha:
                    # Sample from the posterior Gaussian dist.
                    proposed_beta = rand_norms[i] + hdmp_hdmpn * res_beta_hat_i
    
                else:
                    # Sample 0
                    proposed_beta = 0.0
    
                curr_betas[snp_i] = proposed_beta  # UPDATE BETA                
        if verbose and print_progress:
            sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(k + 1) / num_iter))))
            sys.stdout.flush()

        if k >= burn_in:
            avg_betas += curr_post_means  # Averaging over the posterior means instead of samples.

    avg_betas = avg_betas / float(num_iter - burn_in)
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('Took %d minutes and %0.2f seconds' % (t / 60, t % 60))
    return {'betas':avg_betas, 'inf_betas':start_betas}


def ldpred_genomewide(data_file=None, ld_radius=None, ld_dict=None, out_file_prefix=None, 
                      summary_dict=None, ps=None,
                      n=None, h2=None, num_iter=None, 
                      verbose=False, zero_jump_prob=0.05, burn_in=5):
    """
    Calculate LDpred for a genome
    """    
    df = h5py.File(data_file, 'r')
    has_phenotypes = False
    if 'y' in df:
        'Validation phenotypes found.'
        y = df['y'][...]  # Phenotype
        num_individs = len(y)
        risk_scores_pval_derived = sp.zeros(num_individs)
        has_phenotypes = True

    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ld_dict = ld_dict['chrom_ld_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
        
    print('Applying LDpred with LD radius: %d' % ld_radius)
    results_dict = {}
    cord_data_g = df['cord_data']

    #Calculating genome-wide heritability using LD score regression, and partition heritability by chromsomes
    herit_dict = ld.get_chromosome_herits(cord_data_g, ld_scores_dict, n, h2=h2, debug=verbose,summary_dict=summary_dict)

    LDpred_inf_chrom_dict = {}
    print('Calculating LDpred-inf weights')
    for chrom_str in util.chromosomes_list:
        if chrom_str in cord_data_g:
            print('Calculating SNP weights for Chromosome %s' % ((chrom_str.split('_'))[1]))           
            g = cord_data_g[chrom_str]

            # Filter monomorphic SNPs
            snp_stds = g['snp_stds_ref'][...]
            snp_stds = snp_stds.flatten()
            ok_snps_filter = snp_stds > 0
            pval_derived_betas = g['betas'][...]
            pval_derived_betas = pval_derived_betas[ok_snps_filter]            
            h2_chrom = herit_dict[chrom_str]
            start_betas = LDpred_inf.ldpred_inf(pval_derived_betas, genotypes=None, reference_ld_mats=chrom_ref_ld_mats[chrom_str],
                                                h2=h2_chrom, n=n, ld_window_size=2 * ld_radius, verbose=False)
            LDpred_inf_chrom_dict[chrom_str] = start_betas

    
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
        num_chrom = len(util.chromosomes_list)
        for chrom_str in util.chromosomes_list:
            chrom_i+=1
            if chrom_str in cord_data_g:
                g = cord_data_g[chrom_str]
                if verbose and has_phenotypes:
                    if 'raw_snps_val' in g:
                        raw_snps = g['raw_snps_val'][...]
                    else:
                        raw_snps = g['raw_snps_ref'][...]
                
                # Filter monomorphic SNPs
                snp_stds = g['snp_stds_ref'][...]
                snp_stds = snp_stds.flatten()
                pval_derived_betas = g['betas'][...]
                positions = g['positions'][...]
                sids = (g['sids'][...]).astype(util.sids_u_dtype)
                log_odds = g['log_odds'][...]
                nts = (g['nts'][...]).astype(util.nts_u_dtype)
                ok_snps_filter = snp_stds > 0
                if not sp.all(ok_snps_filter):
                    snp_stds = snp_stds[ok_snps_filter]
                    pval_derived_betas = pval_derived_betas[ok_snps_filter]
                    positions = positions[ok_snps_filter]
                    sids = sids[ok_snps_filter]
                    log_odds = log_odds[ok_snps_filter]
                    nts = nts[ok_snps_filter]
                    if verbose and has_phenotypes:
                        raw_snps = raw_snps[ok_snps_filter]


                if out_file_prefix:
                    chromosomes.extend([chrom_str] * len(pval_derived_betas))
                    out_positions.extend(positions)
                    out_sids.extend(sids)
                    raw_effect_sizes.extend(log_odds)
                    out_nts.extend(nts)
        
                
                h2_chrom = herit_dict[chrom_str]
                if 'chrom_ld_boundaries' in ld_dict:
                    ld_boundaries = ld_dict['chrom_ld_boundaries'][chrom_str]
                    res_dict = ldpred_gibbs(pval_derived_betas, h2=h2_chrom, n=n, p=p, ld_radius=ld_radius,
                                            verbose=verbose, num_iter=num_iter, burn_in=burn_in, ld_dict=chrom_ld_dict[chrom_str],
                                            start_betas=LDpred_inf_chrom_dict[chrom_str], ld_boundaries=ld_boundaries,
                                            zero_jump_prob=zero_jump_prob,
                                            print_progress=False)
                else:
                    res_dict = ldpred_gibbs(pval_derived_betas, h2=h2_chrom, n=n, p=p, ld_radius=ld_radius,
                                            verbose=verbose, num_iter=num_iter, burn_in=burn_in, ld_dict=chrom_ld_dict[chrom_str],
                                            start_betas=LDpred_inf_chrom_dict[chrom_str], zero_jump_prob=zero_jump_prob,
                                            print_progress=False)
                
                updated_betas = res_dict['betas']
                updated_inf_betas = res_dict['inf_betas']
                sum_sqr_effects = sp.sum(updated_betas ** 2)
                if sum_sqr_effects > herit_dict['gw_h2_ld_score_est']:
                    print('Sum of squared updated effects estimates seems too large: %0.4f'% sum_sqr_effects)
                    print('This suggests that the Gibbs sampler did not convergence.')
                    convergence_report[p] = True
                
                if verbose:
                    print('Calculating SNP weights for Chromosome %s' % ((chrom_str.split('_'))[1]))
                elif verbose:
                    print('Pprogress: %0.2f%%' % (100.0 * (min(1, float(chrom_i + 1) / num_chrom))))
                    sys.stdout.flush()
                else:
                    sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(chrom_i + 1) / num_chrom))))
                    sys.stdout.flush()

                updated_betas = updated_betas / (snp_stds.flatten())
                updated_inf_betas = updated_inf_betas / (snp_stds.flatten())
                ldpred_effect_sizes.extend(updated_betas)
                ldpred_inf_effect_sizes.extend(updated_inf_betas)
                if verbose and has_phenotypes:
                    prs = sp.dot(updated_betas, raw_snps)
                    risk_scores_pval_derived += prs
                    corr = sp.corrcoef(y, prs)[0, 1]
                    r2 = corr ** 2
                    print('The R2 prediction accuracy of PRS using %s was: %0.4f' % (chrom_str, r2))
        
                    
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
    summary_dict[2.1]={'name':'Convergence issues (for each fraction)','value':str(['Yes' if convergence_report[p] else 'No' for p in ps])}


def main(p_dict):
    summary_dict = {}
    summary_dict[0]={'name':'Coordinated data filename','value':p_dict['cf']}
    summary_dict[0.1]={'name':'SNP weights output file (prefix)', 'value':p_dict['out']}
    summary_dict[0.2]={'name':'LD data filename (prefix)', 'value':p_dict['ldf']}
    summary_dict[1]={'name':'LD radius used','value':str(p_dict['ldr'])}
    t0 = time.time()
    ld_dict = ld.get_ld_dict(p_dict['cf'], p_dict['ldf'], p_dict['ldr'], verbose=p_dict['debug'],
                              compressed=not p_dict['no_ld_compression'], use_hickle=p_dict['hickle_ld'], summary_dict=summary_dict)    
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[1.2]={'name':'Running time for calculating LD information:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    t0 = time.time()
    ldpred_genomewide(data_file=p_dict['cf'], out_file_prefix=p_dict['out'], ps=p_dict['f'], ld_radius=p_dict['ldr'],
                      ld_dict=ld_dict, n=p_dict['N'], num_iter=p_dict['n_iter'], burn_in=p_dict['n_burn_in'], 
                      h2=p_dict['h2'], verbose=p_dict['debug'], summary_dict=summary_dict)
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[2.2]={'name':'Running time for Gibbs sampler(s):','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    reporting.print_summary(summary_dict, 'Summary of LDpred Gibbs')

        

        
