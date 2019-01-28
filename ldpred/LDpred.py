#!/usr/bin/env python
"""
Implements LDpred, an approximate Gibbs sampler that calculate posterior means of effects, conditional on LD information.  
The method requires the user to have generated a coordinated dataset using coord_genotypes.py


Usage: 
ldpred --coord=COORD_DATA_FILE  --ld_radius=LD_RADIUS   --local_ld_file_prefix=LD_FILE_NAME  --PS=FRACTIONS_CAUSAL    
                          --N=SAMPLE_SIZE  --out=OUTPUT_FILE_PREFIX  [ --num_iter=NUM_ITER  --H2=HERTIABILITY  --gm_ld_radius=GEN_MAP_RADIUS]
    
 - COORD_DATA_FILE: The HDF4 file obtained by running the coord_genotypes.py
 
 - LD_RADIUS: An integer number which denotes the number of SNPs on each side of the focal SNP for which LD should be adjusted.  
              A value corresponding M/3000, where M is the number of SNPs in the genome is recommended.
 
 - LD_FILE_NAME: A path and filename prefix for the LD file.  If it doesn't exist, it will be generated.  This can take up to several hours, 
                 depending on LD radius number of SNPs, etc.  If it does exits, that file will be used.
    
 - FRACTIONS_CAUSAL: A list of comma separated (without space) values between 1 and 0, excluding 0.  1 corresponds to the infinitesimal model and will yield results 
                     similar to LDpred-inf.  Default is --PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001 
 
 - N: This is the sample size which LDpred assumes was used to calculate the GWAS summary statistics.

 - OUTPUT_FILE_PREFIX:  The prefix of output file.  

 - NUM_ITER (optional): The number of iterations used by the Gibbs sampler. The default is 60, and burn-in is fixed to 5.
 
 - HERTIABILITY (optional): The heritability assumed by LDpred.  By default it estimates the heritability from the GWAS summary statistics.
 
 - GEN_MAP_RADIUS (optional):  If this option is set, then a genetic map will be used to calculate LD-radius.  A value around 1 is arguably reasonable.
 
  
 2015 (c) Bjarni J Vilhjalmsson: bjarni.vilhjalmsson@gmail.com
 
 """

import getopt
import sys
import time

from scipy import stats
import h5py

import LDpred_inf
import itertools as it
import ld
import scipy as sp

import util

__version__ = '0.9.1'

def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """

    long_options_list = ['coord=', 'ld_radius=', 'local_ld_file_prefix=', 'PS=', 'out=', 'N=',
                         'num_iter=', 'H2=', 'gm_ld_radius=', 'h', 'help']

    p_dict = {'coord':None, 'ld_radius':None, 'local_ld_file_prefix':None, 'PS':[1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001], 'out':None,
              'N':None, 'num_iter': 60, 'H2':None, 'gm':None, 'gm_ld_radius':None}

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)
    
        except:
            print("Some problems with parameters.  Please read the usage documentation carefully.")
            print("Use the -h option for usage information.")
            sys.exit(2)
    
        for opt, arg in opts:
            if opt == "-h" or opt == "--h" or opt == '--help':
                print(__doc__)
                sys.exit(0)
            elif opt == "--coord": p_dict['coord'] = arg
            elif opt == "--ld_radius": p_dict['ld_radius'] = int(arg)
            elif opt == "--local_ld_file_prefix": p_dict['local_ld_file_prefix'] = arg
            elif opt == "--PS": p_dict['PS'] = list(map(float, arg.split(',')))
            elif opt == "--out": p_dict['out'] = arg
            elif opt == "--N": p_dict['N'] = int(arg)
            elif opt == "--num_iter": p_dict['num_iter'] = int(arg)
            elif opt == "--H2": p_dict['H2'] = float(arg)
            elif opt == "--gm_ld_radius": p_dict['gm_ld_radius'] = float(arg)
            else:
                print("Unkown option:", opt)
                print("Use -h option for usage information.")
                sys.exit(2)
    else:
        print(__doc__)
        sys.exit(0)

    if _validate_parameters_(p_dict):
        print ('Use -h flag to print documentation.')
        sys.exit(0)
   
    return p_dict

def _validate_parameters_(p_dict):
    any_missing = False
    if p_dict['coord'] is None:
        print('--coord flag missing:  Please provide a coordinated data file (generated using "python coord_genotypes.py").')
        any_missing = True
    if p_dict['ld_radius'] is None:
        print('--ld_radius flag missing: Please provide a LD radius size (in number of SNPs).')
        any_missing = True
    if p_dict['local_ld_file_prefix'] is None:
        print('--local_ld_file_prefix flag missing: Please provide a LD file (prefix).')
        any_missing = True
    if p_dict['N'] is None:
        print('--N flag missing: Please provide estimated sample size used to calculate the GWAS summary statistics.')
        any_missing = True
    if p_dict['out'] is None:
        print('--out flag missing: Please provide an output filename (prefix).')
        any_missing = True
    return any_missing



def ldpred_genomewide(data_file=None, ld_radius=None, ld_dict=None, out_file_prefix=None, ps=None,
               n=None, h2=None, num_iter=None, verbose=False, zero_jump_prob=0.05, burn_in=5):
    """
    Calculate LDpred for a genome
    """    
    
    df = h5py.File(data_file, 'r')
    has_phenotypes = False
    if 'y' in list(df.keys()):
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
    herit_dict = ld.get_chromosome_herits(cord_data_g, ld_scores_dict, n, h2)

    LDpred_inf_chrom_dict = {}
    print('Calculating LDpred-inf weights')
    for chrom_str in util.chromosomes_list:
        if chrom_str in list(cord_data_g.keys()):
            print('Calculating scores for Chromosome %s' % ((chrom_str.split('_'))[1]))           
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
    
    
    for p in ps:
        print('Starting LDpred with p=%0.4f' % p)
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
            
        for chrom_str in util.chromosomes_list:
            if chrom_str in list(cord_data_g.keys()):
                g = cord_data_g[chrom_str]
                if has_phenotypes:
                    if 'raw_snps_val' in list(g.keys()):
                        raw_snps = g['raw_snps_val'][...]
                    else:
                        raw_snps = g['raw_snps_ref'][...]
                
                # Filter monomorphic SNPs
                snp_stds = g['snp_stds_ref'][...]
                snp_stds = snp_stds.flatten()
                ok_snps_filter = snp_stds > 0
                snp_stds = snp_stds[ok_snps_filter]
                pval_derived_betas = g['betas'][...]
                pval_derived_betas = pval_derived_betas[ok_snps_filter]
                positions = g['positions'][...]
                positions = positions[ok_snps_filter]
                sids = g['sids'][...]
                sids = sids[ok_snps_filter]
                log_odds = g['log_odds'][...]
                log_odds = log_odds[ok_snps_filter]
                nts = g['nts'][...]
                nts = nts[ok_snps_filter]


                if out_file_prefix:
                    chromosomes.extend([chrom_str] * len(pval_derived_betas))
                    out_positions.extend(positions)
                    out_sids.extend(sids)
                    raw_effect_sizes.extend(log_odds)
                    out_nts.extend(nts)
        
                
                h2_chrom = herit_dict[chrom_str]
                if 'chrom_ld_boundaries' in list(ld_dict.keys()):
                    ld_boundaries = ld_dict['chrom_ld_boundaries'][chrom_str]
                    res_dict = ldpred_gibbs(pval_derived_betas, h2=h2_chrom, n=n, p=p, ld_radius=ld_radius,
                                            verbose=verbose, num_iter=num_iter, burn_in=burn_in, ld_dict=chrom_ld_dict[chrom_str],
                                            start_betas=LDpred_inf_chrom_dict[chrom_str], ld_boundaries=ld_boundaries,
                                            zero_jump_prob=zero_jump_prob)
                else:
                    res_dict = ldpred_gibbs(pval_derived_betas, h2=h2_chrom, n=n, p=p, ld_radius=ld_radius,
                                            verbose=verbose, num_iter=num_iter, burn_in=burn_in, ld_dict=chrom_ld_dict[chrom_str],
                                            start_betas=LDpred_inf_chrom_dict[chrom_str], zero_jump_prob=zero_jump_prob)
                
                updated_betas = res_dict['betas']
                updated_inf_betas = res_dict['inf_betas']
                sum_sqr_effects = sp.sum(updated_betas ** 2)
                if sum_sqr_effects > herit_dict['gw_h2_ld_score_est']:
                    print('Sum of squared updated effects estimates seems too large: %0.4f'% sum_sqr_effects)
                    print('This suggests that the Gibbs sampler did not convergence.')
                
                print('Calculating scores for Chromosome %s' % ((chrom_str.split('_'))[1]))
                updated_betas = updated_betas / (snp_stds.flatten())
                updated_inf_betas = updated_inf_betas / (snp_stds.flatten())
                ldpred_effect_sizes.extend(updated_betas)
                ldpred_inf_effect_sizes.extend(updated_inf_betas)
                if has_phenotypes:
                    prs = sp.dot(updated_betas, raw_snps)
                    risk_scores_pval_derived += prs
                    corr = sp.corrcoef(y, prs)[0, 1]
                    r2 = corr ** 2
                    print('The R2 prediction accuracy of PRS using %s was: %0.4f' % (chrom_str, r2))
        
                    
        if has_phenotypes:
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
            for chrom, pos, sid, nt, raw_beta, ldpred_beta in it.izip(chromosomes, out_positions, out_sids, out_nts, raw_effect_sizes, ldpred_effect_sizes):
                nt1, nt2 = nt[0], nt[1]
                f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n' % (chrom, pos, sid, nt1, nt2, raw_beta, ldpred_beta))

    weights_out_file = '%s_LDpred-inf.txt' % (out_file_prefix)
    with open(weights_out_file, 'w') as f:
        f.write('chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta \n')
        for chrom, pos, sid, nt, raw_beta, ldpred_inf_beta in it.izip(chromosomes, out_positions, out_sids, out_nts, raw_effect_sizes, ldpred_inf_effect_sizes):
            nt1, nt2 = nt[0], nt[1]
            f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n' % (chrom, pos, sid, nt1, nt2, raw_beta, ldpred_inf_beta))


        
def ldpred_gibbs(beta_hats, genotypes=None, start_betas=None, h2=None, n=1000, ld_radius=100,
                 num_iter=60, burn_in=10, p=None, zero_jump_prob=0.05,
                 ld_dict=None, reference_ld_mats=None, ld_boundaries=None, verbose=False):
    """
    LDpred (Gibbs Sampler) 
    """
    t0 = time.time()
    m = len(beta_hats)
    
    # If no starting values for effects were given, then use the infinitesimal model starting values.
    if start_betas is None:
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
        alpha = min(1 - zero_jump_prob, 1.0 / h2_est, (h2 + 1.0 / sp.sqrt(n)) / h2_est)

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
        if verbose:
            sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(k + 1) / num_iter))))
            sys.stdout.flush()

        if k >= burn_in:
            avg_betas += curr_post_means  # Averaging over the posterior means instead of samples.

    avg_betas = avg_betas / float(num_iter - burn_in)
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('\nTook %d minutes and %0.2f seconds' % (t / 60, t % 60))
    return {'betas':avg_betas, 'inf_betas':start_betas}


def main():
    p_dict = parse_parameters()
    
    print("""
Note: For maximal accuracy all SNPs with LDpred weights should be included in the validation data set.
If they are a subset of the validation data set, then we suggest recalculate LDpred for the overlapping SNPs. 
""")
    ld_dict = ld.get_ld_dict(p_dict['coord'], p_dict['local_ld_file_prefix'], p_dict['ld_radius'], p_dict['gm_ld_radius'])
        
    ldpred_genomewide(data_file=p_dict['coord'], out_file_prefix=p_dict['out'], ps=p_dict['PS'], ld_radius=p_dict['ld_radius'],
                      ld_dict=ld_dict, n=float(p_dict['N']), num_iter=p_dict['num_iter'], h2=float(p_dict['H2']), verbose=False)
            
        

if __name__ == '__main__':
    main()

        
