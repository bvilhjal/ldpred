"""
Various useful LD functions.


"""
import scipy as sp    
import sys, os, gzip
import time
import h5py
import pickle
from ldpred import util

hickle_available = False
try:
    import hickle
    hickle_available = True
except:
    # hickle not found, using pickle instead
    pass

from scipy import linalg 

        
def shrink_r2_mat(D_i,n):
    D_i = sp.clip(D_i, -1, 1)
    D_i[(1.0/(n-1))>sp.absolute(D_i)]=0
    return D_i

def shrink_r2_mat_2(D_i,n):
    D_i = sp.clip(D_i, 0, 1)
    for j in range(len(D_i)):
        D_i[j]=shrink_r2(D_i[j],n)
    return D_i

def shrink_r2(r,n):
    abs_r = sp.absolute(r)
    if abs_r<0.0001:
        return 0
    else:
        return r * sp.sqrt(max(0,(1.0+1.0/float(n-2))*(r**2-1.0/(1.0+float(n-2)))))/abs_r 



def _calculate_D(snps, snp_i, snp, start_i, stop_i, ld_dict, ld_scores):
    m, n = snps.shape
    X = snps[start_i: stop_i]
    D_i = sp.dot(snp, X.T) / n
    D_i_shrunk = shrink_r2_mat(D_i,n)
    r2s = D_i ** 2
    ld_dict[snp_i] = D_i_shrunk
    lds_i = sp.sum(r2s - (1 - r2s) / (n - 2), dtype='float32')
    ld_scores[snp_i] = lds_i


def get_LDpred_ld_tables(snps, ld_radius=100, ld_window_size=0, h2=None, n_training=None, gm=None, gm_ld_radius=None):
    """
    Calculates LD tables, and the LD score in one go...
    """
    
    ld_dict = {}
    m, n = snps.shape
    ld_scores = sp.ones(m)
    ret_dict = {}
    if gm_ld_radius is None:
        for snp_i, snp in enumerate(snps):
            start_i = max(0, snp_i - ld_radius)
            stop_i = min(m, snp_i + ld_radius + 1)
            _calculate_D(snps, snp_i, snp, start_i, stop_i, ld_dict, ld_scores)
    else:
        assert gm is not None, 'Genetic map is missing.'
        window_sizes = []
        ld_boundaries = []
        for snp_i, snp in enumerate(snps):
            curr_cm = gm[snp_i] 
            
            # Now find lower boundary
            start_i = snp_i
            min_cm = gm[snp_i]
            while start_i > 0 and min_cm > curr_cm - gm_ld_radius:
                start_i = start_i - 1
                min_cm = gm[start_i]
            
            # Now find the upper boundary
            stop_i = snp_i
            max_cm = gm[snp_i]
            while stop_i > 0 and max_cm < curr_cm + gm_ld_radius:
                stop_i = stop_i + 1
                max_cm = gm[stop_i]
            
            ld_boundaries.append([start_i, stop_i])    
            curr_ws = stop_i - start_i
            window_sizes.append(curr_ws)
            assert curr_ws > 0, 'Some issues with the genetic map'
            _calculate_D(snps, snp_i, snp, start_i, stop_i, ld_dict, ld_scores)
        
        avg_window_size = sp.mean(window_sizes)
        print('Average # of SNPs in LD window was %0.2f' % avg_window_size)
        if ld_window_size == 0:
            ld_window_size = avg_window_size * 2
        ret_dict['ld_boundaries'] = ld_boundaries
    ret_dict['ld_dict'] = ld_dict
    ret_dict['ld_scores'] = ld_scores
    
    if ld_window_size > 0:
        ref_ld_matrices = []
        inf_shrink_matrices = []
        for wi in range(0, m, ld_window_size):
            start_i = wi
            stop_i = min(m, wi + ld_window_size)
            curr_window_size = stop_i - start_i
            X = snps[start_i: stop_i]
            D = sp.dot(X, X.T) / n
            D = shrink_r2_mat(D,n)
            ref_ld_matrices.append(D)
            if h2 != None and n_training != None:
                A = ((m / h2) * sp.eye(curr_window_size) + (n_training / (1)) * D)
                A_inv = linalg.pinv(A)
                inf_shrink_matrices.append(A_inv)
        ret_dict['ref_ld_matrices'] = ref_ld_matrices
        if h2 != None and n_training != None:
            ret_dict['inf_shrink_matrices'] = inf_shrink_matrices
    return ret_dict


def calc_ld_table(snps, max_ld_dist=2000, min_r2=0.2, verbose=True, normalize=False):
    """
    Calculate LD between all SNPs using a sliding LD square
    
    This function only retains r^2 values above the given threshold
    """
    # Normalize SNPs (perhaps not necessary, but cheap)
    if normalize:
        snps = snps.T
        snps = (snps - sp.mean(snps, 0)) / sp.std(snps, 0)
        snps = snps.T

    
    if verbose:
        print('Calculating LD table')
    t0 = time.time()
    num_snps, num_indivs = snps.shape    
    ld_table = {}
    for i in range(num_snps):
        ld_table[i] = {}

    a = min(max_ld_dist, num_snps)
    num_pairs = (a * (num_snps - 1)) - a * (a + 1) * 0.5
    if verbose:
        print('Correlation between %d pairs will be tested' % num_pairs)
    num_stored = 0
    for i in range(0, num_snps - 1):
        start_i = i + 1
        end_i = min(start_i + max_ld_dist, num_snps)
        ld_vec = sp.dot(snps[i], sp.transpose(snps[start_i:end_i])) / float(num_indivs)
        ld_vec = sp.array(ld_vec).flatten()
        for k in range(start_i, end_i):
            ld_vec_i = k - start_i
            if ld_vec[ld_vec_i] ** 2 > min_r2:
                ld_table[i][k] = ld_vec[ld_vec_i]
                ld_table[k][i] = ld_vec[ld_vec_i]
                num_stored += 1
        if verbose:
            if i % 1000 == 0:
                sys.stdout.write('.')
                sys.stdout.write('\r%0.2f%%' % (100.0 * (min(1, float(i + 1) / (num_snps - 1)))))
                sys.stdout.flush()
    if verbose:
        sys.stdout.write('Done.\n')
        if num_pairs > 0:
            print('Stored %d (%0.4f%%) correlations that made the cut (r^2>%0.3f).' % (num_stored, 100 * (num_stored / float(num_pairs)), min_r2))
        else:
            print('-')
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('\nIt took %d minutes and %0.2f seconds to calculate the LD table' % (t / 60, t % 60))
    del snps
    return ld_table


def extract_snps_from_cord_data_chrom(g):
    raw_snps = g['raw_snps_ref'][...]
    snp_stds = g['snp_stds_ref'][...]
    snp_means = g['snp_means_ref'][...]
    n_raw_snps = len(raw_snps)

    # Filter monomorphic SNPs
    ok_snps_filter = snp_stds > 0
    ok_snps_filter = ok_snps_filter.flatten()
    raw_snps = raw_snps[ok_snps_filter]
    snp_means = snp_means[ok_snps_filter]
    snp_stds = snp_stds[ok_snps_filter]

    n_snps = len(raw_snps)
    snp_means.shape = (n_snps, 1)
    snp_stds.shape = (n_snps, 1)

    # Normalize SNPs..
    snps = sp.array((raw_snps - snp_means) / snp_stds, dtype='float32')
    assert snps.shape == raw_snps.shape, 'Problems normalizing SNPs (array shape mismatch).'
    return snps, n_raw_snps, n_snps


def smart_ld_pruning(scores, ld_table, max_ld=0.5, verbose=False, reverse=False):
    """
    Prunes SNPs in LD, but with smaller scores (p-values or betas)
    
    If using betas, set reversed to True.
    """
    if verbose:
        print('Performing smart LD pruning')
    t0 = time.time()
    if type(scores) == type([]):
        l = list(zip(scores, list(range(len(scores)))))
    else:
        l = list(zip(scores.tolist(), list(range(len(scores)))))
    l.sort(reverse=reverse)
    l = list(map(list, list(zip(*l))))
    rank_order = l[1]
    indices_to_keep = []
    remaining_indices = set(rank_order)
    for i in rank_order:
        if len(remaining_indices) == 0:
            break
        elif not (i in remaining_indices):
            continue
        else:
            indices_to_keep.append(i)
            for j in ld_table[i]:
                if ld_table[i][j] > max_ld and j in remaining_indices:
                    remaining_indices.remove(j)
    indices_to_keep.sort()
    pruning_vector = sp.zeros(len(scores), dtype='bool8')
    pruning_vector[indices_to_keep] = 1
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('\nIt took %d minutes and %0.2f seconds to LD-prune' % (t / 60, t % 60))
    return pruning_vector             


def get_ld_dict_using_p_dict(p_dict, summary_dict):
    return get_ld_dict(p_dict['cf'], p_dict['ldf'], p_dict['ldr'],
                       verbose=p_dict['debug'],
                       compressed=not p_dict['no_ld_compression'],
                       use_hickle=p_dict['hickle_ld'],
                       summary_dict=summary_dict)


def get_ld_dict(cord_data_file, local_ld_file_prefix, ld_radius, 
                gm_ld_radius=None, verbose=False, compressed=True, 
                use_hickle=False, summary_dict=None):
    """
    Returns the LD dictionary.  Creates a new LD file, if the file doesn't already exist.
    """
    if use_hickle:
        if not hickle_available:
            print ('Unable to find hickle on your system, using pickle instead.')
            print ('See http://telegraphic.github.io/hickle/ for how to install.')
            use_hickle = False
    
    local_ld_dict_file = _get_ld_filename_(local_ld_file_prefix, ld_radius, compressed=compressed)
    
    if not os.path.isfile(local_ld_dict_file):             
        t0 = time.time()
        chrom_ld_scores_dict = {}
        chrom_ld_dict = {}
        chrom_ref_ld_mats = {}
        if gm_ld_radius is not None:
            chrom_ld_boundaries = {}
        ld_score_sum = 0
        num_snps = 0
        num_raw_snps=0
        print('Calculating LD information w. radius %d' % ld_radius)

        df = h5py.File(cord_data_file, 'r')
        cord_data_g = df['cord_data']

        for chrom_str in cord_data_g:
            if verbose:
                print('Calculating LD for chromosome %s' % chrom_str)
            g = cord_data_g[chrom_str]
            snps, n_raw_snps, n_snps = extract_snps_from_cord_data_chrom(g)
            num_raw_snps += n_raw_snps
            num_snps += n_snps

            if gm_ld_radius is not None:
                assert 'genetic_map' in g, 'Genetic map is missing.'
                gm = g['genetic_map'][...]
                ret_dict = get_LDpred_ld_tables(snps, gm=gm, gm_ld_radius=gm_ld_radius)
                chrom_ld_boundaries[chrom_str] = ret_dict['ld_boundaries']
            else:
                ret_dict = get_LDpred_ld_tables(snps, ld_radius=ld_radius, ld_window_size=2 * ld_radius)
            chrom_ld_dict[chrom_str] = ret_dict['ld_dict']
            chrom_ref_ld_mats[chrom_str] = ret_dict['ref_ld_matrices']
            ld_scores = ret_dict['ld_scores']
            chrom_ld_scores_dict[chrom_str] = {'ld_scores':ld_scores, 'avg_ld_score':sp.mean(ld_scores)}
            ld_score_sum += sp.sum(ld_scores)
        avg_gw_ld_score = ld_score_sum / float(num_snps)
        ld_scores_dict = {'avg_gw_ld_score': avg_gw_ld_score, 'chrom_dict':chrom_ld_scores_dict, 
                          'num_snps':num_snps, 'num_raw_snps':num_raw_snps}    
        summary_dict[1.1]={'name':'Average LD score:','value':avg_gw_ld_score}
        
        t1 = time.time()
        t = (t1 - t0)
        if verbose:
            print('\nIt took %d minutes and %0.2f seconds to calculate LD information' % (t / 60, t % 60))
            print('Done calculating the LD table and LD score, writing to file: %s'%local_ld_dict_file)
        ld_dict = {'ld_scores_dict':ld_scores_dict, 'chrom_ld_dict':chrom_ld_dict, 'chrom_ref_ld_mats':chrom_ref_ld_mats}
        if gm_ld_radius is not None:
            ld_dict['chrom_ld_boundaries'] = chrom_ld_boundaries 
        _serialize_ld_info_(local_ld_dict_file, ld_dict, verbose=verbose, compressed=compressed)
        if verbose:
            print('LD information has now been serialized (written to disk).')
    else:
        if verbose:
            print('Loading LD information from file: %s' % local_ld_dict_file)
        ld_dict = _load_ld_info_(local_ld_dict_file, verbose=verbose, compressed=compressed)
        
        num_raw_snps=0
        #Verify LD data
        df = h5py.File(cord_data_file, 'r')
        cord_data_g = df['cord_data']
        for chrom_str in cord_data_g:
            g = cord_data_g[chrom_str]
            num_raw_snps += len(g['raw_snps_ref'])
        df.close()
        assert num_raw_snps == ld_dict['ld_scores_dict']['num_raw_snps'], 'LD reference stored in the provided file, does not seem to match the coordinated SNP data.  Perhaps you want to delete the file and try again.'
        
    return ld_dict


def get_chromosome_herits(cord_data_g, ld_scores_dict, n, h2=None, max_h2=1, use_gw_h2=False, debug=False, summary_dict={}):
    """
    Calculating genome-wide heritability using LD score regression, and partition heritability by chromosome
    """
    tot_num_snps = 0
    tot_sum_beta2s = 0
    herit_dict = {}
    for chrom_str in util.chromosomes_list:
        if chrom_str in cord_data_g:
            g = cord_data_g[chrom_str]
            betas = g['betas'][...]
            n_snps = len(betas)
            tot_num_snps += n_snps
            sum_beta2s = sp.sum(betas ** 2)
            tot_sum_beta2s += sum_beta2s
            chr_chi_square_lambda = sp.mean(n * sum_beta2s / float(n_snps))
            chr_avg_ld_score = ld_scores_dict['chrom_dict'][chrom_str]['avg_ld_score']
            chr_h2_ld_score_est = max(0.0001, (max(1, chr_chi_square_lambda) - 1) / (n * (chr_avg_ld_score / n_snps)))
            herit_dict[chrom_str]= {'n_snps':n_snps,'h2':chr_h2_ld_score_est}
            if debug:
                print('%s heritability LD score estimate: %0.4f'% (chrom_str,chr_h2_ld_score_est))
            
    
    print('%d SNP effects were found' % tot_num_snps)

    L = ld_scores_dict['avg_gw_ld_score']
    chi_square_lambda = sp.mean(n * tot_sum_beta2s / float(tot_num_snps))
    gw_h2_ld_score_est = max(0.0001, (max(1, chi_square_lambda) - 1) / (n * (L / tot_num_snps)))
    if debug:
        print('Genome-wide lambda inflation: %0.4f'% chi_square_lambda)
        print('Genome-wide mean LD score: %0.4f'% L)
        print('LD-score estimated genome-wide heritability: %0.4f'% gw_h2_ld_score_est)
    summary_dict[1.11]={'name':'Genome-wide (LDscore) estimated heritability:','value':'%0.4f'%gw_h2_ld_score_est}
    summary_dict[1.12]={'name':'Chi-square lambda (inflation statistic).','value':'%0.4f'%chi_square_lambda}
    assert chi_square_lambda>1, 'Something is wrong with the GWAS summary statistics, parsing of them, or the given GWAS sample size (N). Lambda (the mean Chi-square statistic) is too small.  '

    #Only use LD score heritability if it is not given as argument. 
    if h2==None:
        if gw_h2_ld_score_est>1:
            print ('LD-score estimated heritability is suspiciously large, suggesting that the given sample size is wrong, or that SNPs are enriched for heritability (e.g. using p-value thresholding). If the SNPs are enriched for heritability we suggest using the --h2 flag to provide a more reasonable heritability estimate.')
        h2 = min(gw_h2_ld_score_est,max_h2)
    if debug:
        print('Heritability used for inference: %0.4f'%h2)
        
    #Distributing heritabilities among chromosomes.
    if use_gw_h2:
        for chrom_str in herit_dict:
            herit_dict[chrom_str]['h2'] = h2 * herit_dict[chrom_str]['n_snps']/float(tot_num_snps)
        
    herit_dict['gw_h2_ld_score_est'] = gw_h2_ld_score_est

    return herit_dict

def _get_ld_filename_(local_ld_file_prefix, ld_radius, compressed=True, use_hickle=False):
    if use_hickle:
        if compressed:
            return '%s_ldradius%d_gzip.hkl'%(local_ld_file_prefix, ld_radius)
        else:
            return '%s_ldradius%d.hkl'%(local_ld_file_prefix, ld_radius)

    else:
        if compressed:
            return '%s_ldradius%d.pkl.gz'%(local_ld_file_prefix, ld_radius)
        else:
            return '%s_ldradius%d.pkl'%(local_ld_file_prefix, ld_radius)


def _serialize_ld_info_(local_ld_dict_file, ld_dict, verbose=False, compressed=True, use_hickle=False):
    t0 = time.time()
    if use_hickle:
        f = h5py.File(local_ld_dict_file, 'w')
        if compressed:
            print('Storing compressed LD information to hdf5 file')
            hickle.dump(ld_dict, f, compression='gzip')
        else:
            hickle.dump(ld_dict, f)
        f.close()
    else:
        if compressed:
            print('Storing LD information to compressed pickle file')
            f = gzip.open(local_ld_dict_file, 'wb')
        else:
            f = open(local_ld_dict_file, 'wb')            
        pickle.dump(ld_dict, f, protocol=-1)
        f.close()
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('\nIt took %d minutes and %0.2f seconds to write LD information to disk.' % (t / 60, t % 60))
        print('LD information file size on disk: %0.4f Mb' % float(os.path.getsize(local_ld_dict_file)/1000000.0))


def _load_ld_info_(local_ld_dict_file, verbose=True, compressed=True, use_hickle=False):
    t0 = time.time()
    if use_hickle:
        f = h5py.File(local_ld_dict_file, 'r')
        ld_dict = hickle.load(f)
        f.close()
    else:
        if compressed:
            f = gzip.open(local_ld_dict_file, 'r')
        else:
            f = open(local_ld_dict_file, 'r')            
        ld_dict = pickle.load(f)
        f.close()
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print('\nIt took %d minutes and %0.2f seconds to load LD information from disk.' % (t / 60, t % 60))
    return ld_dict
