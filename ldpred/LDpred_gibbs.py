import sys
import time

import numpy as np
from scipy import stats
import h5py

from ldpred import LDpred_inf
from ldpred import util
from ldpred import ld
from ldpred import reporting
from ldpred import coord_genotypes
from ldpred._numba import HAVE_NUMBA, _jit, numba_status


# What to do with variants inside the long-range LD regions (Price et al. AJHG
# 2008; 34 regions, ~5% of the genome).  LDpred models LD as a band of
# +/- ld_radius, which cannot capture LD in those regions, so sampling them
# normally can stop the Gibbs sampler from converging.
#
#   include   sample them like any other variant
#   inf       hold them out of the sampler, report their LDpred-inf weight
#   exclude   hold them out with a zero effect, so they leave the model
#             entirely, and report a zero weight
#
# simulations/lrld_policy.py compares these.  'inf' is the default because it
# is the only one that is safe in both regimes: it never diverges, and unlike
# 'exclude' it costs nothing when a listed region turns out to behave normally.
# 'exclude' is the better choice when a region really does carry strong
# long-range LD.
LRLD_INCLUDE = 'include'
LRLD_INF = 'inf'
LRLD_EXCLUDE = 'exclude'
LRLD_POLICIES = (LRLD_INCLUDE, LRLD_INF, LRLD_EXCLUDE)
LRLD_DEFAULT = LRLD_INF


def resolve_lrld_policy(p_dict):
    """Resolve the long-range LD policy from parsed command line arguments.

    --incl-long-range-ld is the historical boolean spelling of
    ``--long-range-ld include``.  It still works, but warns.
    """
    policy = p_dict.get('long_range_ld') or LRLD_DEFAULT
    if p_dict.get('incl_long_range_ld'):
        if policy != LRLD_DEFAULT:
            raise ValueError(
                '--incl-long-range-ld cannot be combined with --long-range-ld '
                '%s; use --long-range-ld include instead.' % policy)
        print('\033[93m Warning: --incl-long-range-ld is deprecated. '
              'Use --long-range-ld include instead. \033[0m')
        policy = LRLD_INCLUDE
    if policy not in LRLD_POLICIES:
        raise ValueError('--long-range-ld must be one of %s, got %r'
                         % (', '.join(LRLD_POLICIES), policy))
    return policy


def prepare_constants(ldpred_n,ns,m,p,h2,sampl_var_shrink_factor):
    Mp = m * p
    hdmp = (h2 / Mp)
    const_dict = {'Mp':Mp, 'hdmp':hdmp}
    rv_scalars = np.zeros(m)
    if ldpred_n is not None:
        hdmpn = hdmp + 1.0 / ldpred_n
        hdmp_hdmpn = (hdmp / hdmpn)
        c_const = (p / np.sqrt(hdmpn))
        d_const = (1.0 - p) / (np.sqrt(1.0 / ldpred_n))
        const_dict['n']=ldpred_n
        const_dict['hdmpn']=hdmpn
        const_dict['hdmp_hdmpn']=hdmp_hdmpn
        const_dict['c_const']=c_const
        const_dict['d_const']=d_const
        rv_scalars[:]=sampl_var_shrink_factor* np.sqrt((hdmp_hdmpn) * (1.0 / ldpred_n))
    else:
        snp_dict = {}
        for i in range(m):
            ni = ns[i]
            hdmpn_i = hdmp + 1.0 / ni
            hdmp_hdmpn_i = (hdmp / hdmpn_i)
            c_const_i = (p / np.sqrt(hdmpn_i))
            d_const_i = (1.0 - p) / (np.sqrt(1.0 / ni))
            snp_dict[i]={'n':ni, 'hdmpn':hdmpn_i, 
                           'hdmp_hdmpn':hdmp_hdmpn_i, 
                           'c_const':c_const_i, 
                           'd_const':d_const_i}
            rv_scalars[i]=sampl_var_shrink_factor* np.sqrt((hdmp_hdmpn_i) * (1.0 / ni))
        const_dict['snp_dict']=snp_dict
    const_dict['rv_scalars']=rv_scalars
    return const_dict


def get_constants(snp_i,const_dict):
    if 'snp_dict' in const_dict:
        return const_dict['snp_dict'][snp_i]
    else:
        return const_dict


def constants_as_arrays(const_dict, m):
    """Flatten the per-SNP shrink constants into arrays for the sweep kernel.

    prepare_constants returns either one set of scalars (fixed N) or a dict of
    per-SNP dicts (varying N).  The kernel wants plain arrays either way, so the
    fixed-N case is broadcast.  Values are unchanged.
    """
    if 'snp_dict' in const_dict:
        snp_dict = const_dict['snp_dict']
        n_arr = np.empty(m)
        hdmpn_arr = np.empty(m)
        hdmp_hdmpn_arr = np.empty(m)
        c_const_arr = np.empty(m)
        d_const_arr = np.empty(m)
        for i in range(m):
            cd = snp_dict[i]
            n_arr[i] = cd['n']
            hdmpn_arr[i] = cd['hdmpn']
            hdmp_hdmpn_arr[i] = cd['hdmp_hdmpn']
            c_const_arr[i] = cd['c_const']
            d_const_arr[i] = cd['d_const']
    else:
        n_arr = np.full(m, const_dict['n'], dtype=np.float64)
        hdmpn_arr = np.full(m, const_dict['hdmpn'], dtype=np.float64)
        hdmp_hdmpn_arr = np.full(m, const_dict['hdmp_hdmpn'], dtype=np.float64)
        c_const_arr = np.full(m, const_dict['c_const'], dtype=np.float64)
        d_const_arr = np.full(m, const_dict['d_const'], dtype=np.float64)
    return n_arr, hdmpn_arr, hdmp_hdmpn_arr, c_const_arr, d_const_arr


def flatten_ld_dict(ld_dict, m, ld_radius, ld_boundaries):
    """Pack the per-SNP LD vectors into one contiguous array.

    ld_dict maps a SNP index to its LD vector over the surrounding window.  The
    kernel cannot index a Python dict, so the vectors are concatenated and
    addressed by (offset, length) with the window's left edge recorded
    separately.  The left edge is derived exactly as the Python loop derives it,
    from ld_boundaries when a genetic map was used and from ld_radius otherwise.
    """
    lengths = np.empty(m, dtype=np.int64)
    starts = np.empty(m, dtype=np.int64)
    for snp_i in range(m):
        lengths[snp_i] = len(ld_dict[snp_i])
        if ld_boundaries is None:
            starts[snp_i] = max(0, snp_i - ld_radius)
        else:
            starts[snp_i] = ld_boundaries[snp_i][0]
    offsets = np.zeros(m, dtype=np.int64)
    if m > 1:
        offsets[1:] = np.cumsum(lengths)[:-1]
    flat = np.empty(int(lengths.sum()), dtype=np.float64)
    for snp_i in range(m):
        o = offsets[snp_i]
        flat[o:o + lengths[snp_i]] = ld_dict[snp_i]
    return flat, offsets, lengths, starts


def _gibbs_sweep(beta_hats, ld_flat, ld_offsets, ld_lengths, ld_starts,
                 curr_betas, curr_post_means, rand_ps, rand_norms, skip_mask,
                 n_arr, hdmpn_arr, hdmp_hdmpn_arr, c_const_arr, d_const_arr,
                 alpha, m):
    """One Gibbs sweep over all SNPs; mutates curr_betas / curr_post_means.

    This is a transcription of the loop that used to live inline in
    ldpred_gibbs, arranged so Numba can compile it.  Zeroing the focal element
    before the dot product reproduces the original code exactly: it took a view
    of curr_betas and set the focal entry to 0, and the entry is overwritten
    with the new draw at the end of the iteration either way.

    The random draws are passed in rather than generated here so the random
    stream is byte-for-byte what the pure Python version consumed.
    """
    for snp_i in range(m):
        if skip_mask[snp_i]:
            continue

        offset = ld_offsets[snp_i]
        length = ld_lengths[snp_i]
        start_i = ld_starts[snp_i]

        curr_betas[snp_i] = 0.0
        res_beta_hat_i = beta_hats[snp_i]
        for t in range(length):
            res_beta_hat_i -= ld_flat[offset + t] * curr_betas[start_i + t]
        b2 = res_beta_hat_i * res_beta_hat_i

        d_const_b2_exp = d_const_arr[snp_i] * np.exp(-b2 * n_arr[snp_i] / 2.0)
        numerator = c_const_arr[snp_i] * np.exp(-b2 / (2.0 * hdmpn_arr[snp_i]))
        if numerator == 0.0:
            postp = 0.0
        else:
            postp = numerator / (numerator + d_const_b2_exp)
        curr_post_means[snp_i] = hdmp_hdmpn_arr[snp_i] * postp * res_beta_hat_i

        if rand_ps[snp_i] < postp * alpha:
            curr_betas[snp_i] = (rand_norms[snp_i]
                                 + hdmp_hdmpn_arr[snp_i] * res_beta_hat_i)
        else:
            curr_betas[snp_i] = 0.0


_gibbs_sweep_jit = _jit(_gibbs_sweep)


def _gibbs_sweep_numpy(beta_hats, ld_flat, ld_offsets, ld_lengths, ld_starts,
                       curr_betas, curr_post_means, rand_ps, rand_norms, skip_mask,
                       n_arr, hdmpn_arr, hdmp_hdmpn_arr, c_const_arr, d_const_arr,
                       alpha, m):
    """Fallback sweep for when Numba is unavailable.

    Identical arithmetic to :func:`_gibbs_sweep`, but the inner product is a
    single np.dot rather than an element-by-element loop.  Running the kernel's
    scalar loop in the interpreter would be far slower than the code this
    replaced, so the two paths deliberately differ in how the dot is spelled and
    in nothing else.
    """
    for snp_i in range(m):
        if skip_mask[snp_i]:
            continue

        offset = ld_offsets[snp_i]
        length = ld_lengths[snp_i]
        start_i = ld_starts[snp_i]

        curr_betas[snp_i] = 0.0
        res_beta_hat_i = beta_hats[snp_i] - np.dot(
            ld_flat[offset:offset + length],
            curr_betas[start_i:start_i + length])
        b2 = res_beta_hat_i * res_beta_hat_i

        d_const_b2_exp = d_const_arr[snp_i] * np.exp(-b2 * n_arr[snp_i] / 2.0)
        numerator = c_const_arr[snp_i] * np.exp(-b2 / (2.0 * hdmpn_arr[snp_i]))
        if numerator == 0.0:
            postp = 0.0
        else:
            postp = numerator / (numerator + d_const_b2_exp)
        curr_post_means[snp_i] = hdmp_hdmpn_arr[snp_i] * postp * res_beta_hat_i

        if rand_ps[snp_i] < postp * alpha:
            curr_betas[snp_i] = (rand_norms[snp_i]
                                 + hdmp_hdmpn_arr[snp_i] * res_beta_hat_i)
        else:
            curr_betas[snp_i] = 0.0


def ldpred_gibbs(beta_hats, genotypes=None, start_betas=None, h2=None, n=None, ns= None, ld_radius=100,
                 num_iter=60, burn_in=10, p=None, zero_jump_prob=0.01, sampl_var_shrink_factor=0.9, 
                 tight_sampling=False,ld_dict=None, reference_ld_mats=None, ld_boundaries=None, 
                 snp_lrld=None, lrld_policy=LRLD_DEFAULT, verbose=False, print_progress=True):
    """
    LDpred (Gibbs Sampler)

    snp_lrld marks variants in long-range LD regions; lrld_policy says what to
    do with them, and is ignored when snp_lrld is None.  See LRLD_POLICIES.
    """
    if lrld_policy not in LRLD_POLICIES:
        raise ValueError('lrld_policy must be one of %s, got %r'
                         % (', '.join(LRLD_POLICIES), lrld_policy))
    # Set random seed to stabilize results
    np.random.seed(42) 

    t0 = time.time()
    m = len(beta_hats)


    ldpred_n, ldpred_inf_n = util.get_LDpred_sample_size(n,ns,verbose)
    
    # If no starting values for effects were given, then use the infinitesimal model starting values.
    if start_betas is None and verbose:
        print('Initializing LDpred effects with posterior mean LDpred-inf effects.')
        print('Calculating LDpred-inf effects.')
        start_betas = LDpred_inf.ldpred_inf(beta_hats, genotypes=genotypes, reference_ld_mats=reference_ld_mats,
                                            h2=h2, n=ldpred_inf_n, ld_window_size=2 * ld_radius, verbose=False)
    # float64 throughout. The Python loop got this for free because np.dot
    # upcast a float32 LD vector against float64 effects; the kernel is
    # explicitly typed, so pin the dtype rather than let it depend on whatever
    # the caller passed in.
    beta_hats = np.ascontiguousarray(beta_hats, dtype=np.float64)
    curr_betas = np.ascontiguousarray(start_betas, dtype=np.float64).copy()
    assert len(curr_betas)==m,'Betas returned by LDpred_inf do not have the same length as expected.'
    curr_post_means = np.zeros(m)
    avg_betas = np.zeros(m)

    # Setting up the marginal Bayes shrink
    const_dict = prepare_constants(ldpred_n,ns,m,p,h2,sampl_var_shrink_factor)

    # Pack everything the sweep touches into arrays once, outside the iteration
    # loop, so the per-SNP work is pure arithmetic on contiguous memory.
    n_arr, hdmpn_arr, hdmp_hdmpn_arr, c_const_arr, d_const_arr = constants_as_arrays(const_dict, m)
    ld_flat, ld_offsets, ld_lengths, ld_starts = flatten_ld_dict(
        ld_dict, m, ld_radius, ld_boundaries)
    if snp_lrld is None:
        skip_mask = np.zeros(m, dtype=bool)
    else:
        skip_mask = np.asarray(snp_lrld, dtype=bool)

    # Variants in long-range LD regions are held out of the sweep: the LD band
    # is badly misspecified there, and sampling them makes the chain diverge.
    # Skipped variants stay frozen at whatever value they start with, so the
    # starting value decides both what they contribute to their neighbours'
    # residuals and, via curr_post_means, what weight they are reported with.
    # See simulations/lrld_policy.py for the comparison behind the default.
    if skip_mask.any():
        if lrld_policy == 'exclude':
            # Out of the model entirely: a zero effect contributes nothing to
            # any neighbour, and the posterior mean stays at its initial zero,
            # so the reported weight is zero too.
            curr_betas[skip_mask] = 0.0
        else:                                     # 'inf'
            # Held out of sampling but kept in the model at the LDpred-inf
            # effect, and reported with it. The sweep never writes these
            # entries, so averaging the seeded value over the retained
            # iterations returns it unchanged.
            curr_post_means[skip_mask] = curr_betas[skip_mask]

    sweep = _gibbs_sweep_jit if HAVE_NUMBA else _gibbs_sweep_numpy

    for k in range(num_iter):  # Big iteration
        h2_est = max(0.00001, np.sum(curr_betas ** 2))
        if tight_sampling:
            # Force an alpha shrink if estimates are way off compared to heritability estimates.
            #(May improve MCMC convergence.)
            alpha = min(1.0 - zero_jump_prob, 1.0 / h2_est, (h2 + 1.0 / np.sqrt(ldpred_n)) / h2_est)
        else:
            alpha = 1.0 - zero_jump_prob

        rand_ps = np.random.random(m)

        rand_norms = stats.norm.rvs(0.0, 1, size=m)*const_dict['rv_scalars']

        sweep(beta_hats, ld_flat, ld_offsets, ld_lengths, ld_starts,
              curr_betas, curr_post_means, rand_ps, rand_norms, skip_mask,
              n_arr, hdmpn_arr, hdmp_hdmpn_arr, c_const_arr, d_const_arr,
              alpha, m)

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
                      sampl_var_shrink_factor=1, lrld_policy=LRLD_DEFAULT,
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
        risk_scores_pval_derived = np.zeros(num_individs)
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

    if lrld_policy != 'include':
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
                positions = g['positions'][...]

                if out_file_prefix:
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
                if lrld_policy != 'include':
                    # Key off the group name, not the loop counter: chromosomes_list
                    # runs to 24 entries while the long-range LD table is keyed 1-23.
                    snp_lrld = util.get_snp_lrld_status(chrom_str, positions, lrld_dict)
                    num_snps_in_lrld +=np.sum(snp_lrld)
                    
                ld_boundaries = None
                if 'chrom_ld_boundaries' in ld_dict:
                    ld_boundaries = ld_dict['chrom_ld_boundaries'][chrom_str]
                if verbose:
                    print('Calculating SNP weights for Chromosome %s' % ((chrom_str.split('_'))[1]))
                res_dict = ldpred_gibbs(pval_derived_betas,h2=h2_chrom, n=n, ns=ns, p=p, ld_radius=ld_radius,
                                        verbose=verbose, num_iter=num_iter, burn_in=burn_in, ld_dict=chrom_ld_dict[chrom_str],
                                        start_betas=LDpred_inf_chrom_dict[chrom_str], ld_boundaries=ld_boundaries,
                                        zero_jump_prob=zero_jump_prob,sampl_var_shrink_factor=sampl_var_shrink_factor,
                                        snp_lrld=snp_lrld, lrld_policy=lrld_policy, print_progress=False)
                
                updated_betas = res_dict['betas']
                updated_inf_betas = res_dict['inf_betas']
                sum_sqr_effects = np.sum(updated_betas ** 2)
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
                        prs = np.dot(updated_betas, raw_snps)
                        risk_scores_pval_derived += prs
                        corr = np.corrcoef(y, prs)[0, 1]
                        r2 = corr ** 2
                        print('The R2 prediction accuracy of PRS using %s was: %0.4f' % (chrom_str, r2))

        if lrld_policy != LRLD_INCLUDE:
            summary_dict[1.3]={'name':'SNPs in long-range LD regions','value':'%d'%num_snps_in_lrld}
            summary_dict[1.31]={
                'name':'Long-range LD handling',
                'value':('held out, LDpred-inf weights' if lrld_policy == LRLD_INF
                         else 'held out, zero weights')}
        else:
            summary_dict[1.31]={'name':'Long-range LD handling',
                                'value':'sampled like any other SNP'}
        
        if not verbose:
            sys.stdout.write('\r%0.2f%%\n' % (100.0))
            sys.stdout.flush()
        if verbose and has_phenotypes:
            num_indivs = len(y)
            results_dict[p_str]['y'] = y
            results_dict[p_str]['risk_scores_pd'] = risk_scores_pval_derived
            print('Prediction accuracy was assessed using %d individuals.' % (num_indivs))
    
            corr = np.corrcoef(y, risk_scores_pval_derived)[0, 1]
            r2 = corr ** 2
            results_dict[p_str]['r2_pd'] = r2
            print('The  R2 prediction accuracy (observed scale) for the whole genome was: %0.4f (%0.6f)' % (r2, ((1 - r2) ** 2) / num_indivs))
    
            if corr < 0:
                risk_scores_pval_derived = -1 * risk_scores_pval_derived
            auc = util.calc_auc(y, risk_scores_pval_derived)
            print('AUC for the whole genome was: %0.4f' % auc)
    
            # Now calibration                               
            denominator = np.dot(risk_scores_pval_derived.T, risk_scores_pval_derived)
            y_norm = (y - np.mean(y)) / np.std(y)
            numerator = np.dot(risk_scores_pval_derived.T, y_norm)
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
    summary_dict[2.15]={'name':'Sampler backend','value':numba_status()}
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
                      h2=p_dict['h2'], use_gw_h2=p_dict['use_gw_h2'], lrld_policy=resolve_lrld_policy(p_dict), 
                      sampl_var_shrink_factor=1, verbose=p_dict['debug'], summary_dict=summary_dict)
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[3]={'name':'Running time for Gibbs sampler(s):','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    reporting.print_summary(summary_dict, 'Summary of LDpred Gibbs')

        
