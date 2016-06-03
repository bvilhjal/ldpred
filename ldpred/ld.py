"""
Various useful LD functions.


"""
try: 
    import scipy as sp
except Exception:
    print 'Using Numpy instead of Scipy.'
    import numpy as sp
    
import sys
import time
from numpy import linalg 
import pdb
import os
import cPickle
import plinkio
from plinkio import plinkfile
import itertools as it
import random



def get_LDpred_ld_tables(snps, ld_radius=100, ld_window_size=0, h2=None, n_training=None, gm=None, gm_ld_radius=None):
    """
    Calculates LD tables, and the LD score in one go...
    """
    
    ld_dict = {}
    m,n = snps.shape
    print m,n
    ld_scores = sp.ones(m)
    ret_dict = {}
    if gm_ld_radius is None:
        for snp_i, snp in enumerate(snps):
            # Calculate D
            start_i = max(0, snp_i - ld_radius)
            stop_i = min(m, snp_i + ld_radius + 1)
            X = snps[start_i: stop_i]
            D_i = sp.dot(snp, X.T) / n
            r2s = D_i ** 2
            ld_dict[snp_i] = D_i
            lds_i = sp.sum(r2s - (1-r2s) / (n-2),dtype='float32')
            #lds_i = sp.sum(r2s - (1-r2s)*empirical_null_r2)
            ld_scores[snp_i] =lds_i
    else:
        assert gm is not None, 'Genetic map is missing.'
        window_sizes = []
        ld_boundaries =[]
        for snp_i, snp in enumerate(snps):
            curr_cm = gm[snp_i] 
            
            #Now find lower boundary
            start_i = snp_i
            min_cm = gm[snp_i]
            while start_i>0 and min_cm>curr_cm-gm_ld_radius:
                start_i = start_i -1
                min_cm = gm[start_i]
            
            #Now find the upper boundary
            stop_i = snp_i
            max_cm = gm[snp_i]
            while stop_i>0 and max_cm<curr_cm+gm_ld_radius:
                stop_i = stop_i +1
                max_cm = gm[stop_i]
            
            ld_boundaries.append([start_i,stop_i])    
            curr_ws = stop_i-start_i
            window_sizes.append(curr_ws)
            assert curr_ws>0, 'Some issues with the genetic map'

            X = snps[start_i: stop_i]
            D_i = sp.dot(snp, X.T) / n
            r2s = D_i ** 2
            ld_dict[snp_i] = D_i
            lds_i = sp.sum(r2s - (1-r2s) / (n-2),dtype='float32')
            #lds_i = sp.sum(r2s - (1-r2s)*empirical_null_r2)
            ld_scores[snp_i] =lds_i
        
        avg_window_size=sp.mean(window_sizes)
        print 'Average # of SNPs in LD window was %0.2f'%avg_window_size
        if ld_window_size==0:
            ld_window_size = avg_window_size*2
        ret_dict['ld_boundaries'] = ld_boundaries
    ret_dict['ld_dict']=ld_dict
    ret_dict['ld_scores']=ld_scores
    
    if ld_window_size>0:
        ref_ld_matrices = []
        inf_shrink_matrices = []
        for i, wi in enumerate(range(0, m, ld_window_size)):
            start_i = wi
            stop_i = min(m, wi + ld_window_size)
            curr_window_size = stop_i - start_i
            X = snps[start_i: stop_i]
            D = sp.dot(X, X.T) / n
            ref_ld_matrices.append(D)
            if h2!=None and n_training!=None:
                A = ((m / h2) * sp.eye(curr_window_size) + (n_training / (1)) * D)
                A_inv = linalg.pinv(A)
                inf_shrink_matrices.append(A_inv)
        ret_dict['ref_ld_matrices']=ref_ld_matrices
        if h2!=None and n_training!=None:
            ret_dict['inf_shrink_matrices']=inf_shrink_matrices
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
        print 'Calculating LD table'
    t0 = time.time()
    num_snps, num_indivs = snps.shape    
    ld_table = {}
    for i in range(num_snps):
        ld_table[i] = {}

    a = min(max_ld_dist, num_snps)
    num_pairs = (a * (num_snps - 1)) - a * (a + 1) * 0.5
    if verbose:
        print 'Correlation between %d pairs will be tested' % num_pairs
    num_stored = 0
    for i in range(0, num_snps - 1):
        start_i = i + 1
        end_i = min(start_i + max_ld_dist, num_snps)
        ld_vec = sp.dot(snps[i], sp.transpose(snps[start_i:end_i])) / float(num_indivs)
        ld_vec = sp.array(ld_vec).flatten()
        for k in range(start_i, end_i):
            ld_vec_i = k - start_i
            if ld_vec[ld_vec_i] > min_r2:
                ld_table[i][k] = ld_vec[ld_vec_i]
                ld_table[k][i] = ld_vec[ld_vec_i]
                num_stored += 1
        if verbose:
            if i % 1000 == 0:
                sys.stdout.write('.')
#                 sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(i + 1) / (num_snps - 1)))))
                sys.stdout.flush()
    if verbose:
        sys.stdout.write('Done.\n')
        if num_pairs>0:
            print 'Stored %d (%0.4f%%) correlations that made the cut (r^2>%0.3f).' % (num_stored, 100 * (num_stored / float(num_pairs)), min_r2)
        else:
            print '-'
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to calculate the LD table' % (t / 60, t % 60)
    del snps
    return ld_table



def calc_full_ld_table(ld_mat, min_ld=0.2, verbose=True):
    """
    Calculate LD between all SNPs using an efficient sliding LD square
    """
    if verbose:
        print 'Calculating LD table'
    t0 = time.time()
    m = ld_mat.shape[0]
    ld_table = {}
    for i in range(m):
        ld_table[i] = {}

    D2 = sp.array(ld_mat) ** 2
    # print D2
    for i in range(m - 1):
        for j in range(i + 1, m):
            if D2[i, j] > min_ld:
                ld_table[j][i] = D2[i, j]
                ld_table[i][j] = D2[i, j]
        if verbose:
            sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(i + 1) / (m)))))
            sys.stdout.flush()
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to calculate the LD table' % (t / 60, t % 60)
    return ld_table




def ml_LD_shrink(beta_hats, genotypes=None, reference_ld_mats=None, window_method='sliding',
                              ld_window_size=100, ld_radius=20, verbose=False):
    """
    Calculate the joint least square estimates using LD information.
    
    Two methods are implemented, i.e. 'sliding', and 'tiling'. 
    
    If reference_ld_mats are supplied, it uses those, otherwise it uses the LD in the genotype data.
    """    
    if verbose:
        print 'Doing LD correction'
    t0 = time.time()
    num_betas = len(beta_hats)
    updated_betas = sp.empty(num_betas)
    m = len(beta_hats)
    
    if window_method == 'tiling':
        for i, wi in enumerate(range(0, num_betas, ld_window_size)):
            start_i = wi 
            stop_i = min(num_betas, wi + ld_window_size)
            curr_window_size = stop_i - start_i
            if reference_ld_mats != None:
                D = reference_ld_mats[i]
            else:
                if genotypes != None:
                    X = genotypes[start_i: stop_i]
                    num_indivs = X.shape[1]
                    D = sp.dot(X, X.T) / num_indivs
                else:
                    raise NotImplementedError
            D_inv = linalg.pinv(D)
            updated_betas[start_i: stop_i] = sp.dot(D_inv , beta_hats[start_i: stop_i])  # Adjust the beta_hats    
            if verbose:
                sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(wi + 1) / num_betas))))
                sys.stdout.flush()
    
    elif window_method == 'sliding':
        if genotypes == None:
            raise NotImplementedError
        X = genotypes[0: ld_radius]
        num_indivs = X.shape[1]
        D_prev = sp.dot(X, X.T) / num_indivs
        D_inv = linalg.pinv(D_prev)                        
        updated_betas[0: ld_radius] = sp.dot(D_inv , beta_hats[0: ld_radius])
        for i in range(1, num_betas - ld_radius):
            D = sp.empty((ld_radius, ld_radius))
            D[0:ld_radius - 1, 0:ld_radius - 1] = D_prev[1:ld_radius, 1:ld_radius]
            X = genotypes[i:i + ld_radius]
            X_add = genotypes[i + ld_radius]
            last_row = sp.dot(X, X_add.T) / num_indivs
            for j in range(ld_radius):
                D[ld_radius - 1, j] = last_row[j]
                D[j, ld_radius - 1] = last_row[j]
            D_inv = linalg.pinv(D)                        
            i_focal = i + ld_radius / 2
            updated_betas[i_focal] = sp.dot(D_inv[ld_radius / 2] , beta_hats[i:i + ld_radius])
            if verbose:
                sys.stdout.write('\b\b\b\b\b\b\b\b\b%0.1f%%' % (100.0 * (min(1, float(i + 1) / (num_betas - ld_radius)))))
                sys.stdout.flush()     
        updated_betas[-ld_radius:] = sp.dot(D_inv , beta_hats[-ld_radius:])
        

    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to perform the Infinitesimal LD shrink' % (t / 60, t % 60)
    return updated_betas


 
def ml_iter(beta_hats, genotypes, ld_radius=20,
            verbose=False, iter_percentile=0.05,
            ld_dict_file_prefix=None, ld_dict=None):
    """
    Yang et al. iterative scheme.
    
    Idea:
    - While # selected betas is < threshold
    -     Sort betas
    -     Pick the largest beta that hasn't been selected. 
    -     For each marker in window around selected marker:
    -         Update LD matrix
    -         Invert LD matrix
    -         Update beta

    """

    if verbose:
        print 'Performing iterative approach'
    t0 = time.time()
    
    #Ordering the beta_hats, and storing the order
    m = len(beta_hats)
    indices = range(m)
    beta_hats = beta_hats.tolist()

    ld_table = {}
    for i in range(m):
        ld_table[i] = {'ld_partners':[i], 'beta_hat_list':[beta_hats[i]], 'D' : sp.array([1.0]), 'D_inv':sp.array([1.0])}
    genotypes = sp.array(genotypes)
    n_test = len(genotypes[0])
    genotypes = sp.array(genotypes)
    assert len(genotypes) == m, 'WTF?'
    if verbose:
        print '# SNPs: %d' % m

    max_num_selected = int(m * iter_percentile)
    selected_indices = set()
    updated_beta_hats = beta_hats[:]
    while len(selected_indices) < max_num_selected:
        #Sort and select beta
        l = zip((sp.array(updated_beta_hats) ** 2).tolist(), range(m))
        l.sort(reverse=True)
        for beta_hat, beta_i in l:
            if not beta_i in selected_indices:
                selected_indices.add(beta_i)
                break
        #Iterating over the window around the selected beta
        start_i = max(0, beta_i - ld_radius)
        end_i = min(beta_i + ld_radius, m)
        for i in range(start_i, end_i):
            if i == beta_i:
                continue
            else:
                #Updaate the LD matrix
                d = ld_table[i]
                ld_partners = d['ld_partners']
                beta_hat_list = d['beta_hat_list']
                ld_partners.append(beta_i)
                beta_hat_list.append(beta_hats[beta_i])
                bs = sp.array(beta_hat_list)
                X = genotypes[ld_partners]
                D = sp.dot(X, X.T) / float(n_test)
                ld_table['D'] = D
                D_inv = linalg.pinv(D)
                ld_table['D_inv'] = D_inv
                updated_beta = sp.dot(D_inv[0], bs)
                #assert updated_beta != beta_hats[i]
                updated_beta_hats[i] = updated_beta
    for i in range(m):
        if not i in selected_indices:
            updated_beta_hats[i] = 0
    updated_beta_hats = sp.array(updated_beta_hats)
    return updated_beta_hats
    


def smart_ld_pruning(scores, ld_table, max_ld=0.5, verbose=False, reverse=False):
    """
    Prunes SNPs in LD, but with smaller scores (p-values or betas)
    
    If using betas, set reversed to True.
    """
    if verbose:
        print 'Performing smart LD pruning'
    t0 = time.time()
    if type(scores) == type([]):
        l = zip(scores, range(len(scores)))
    else:
        l = zip(scores.tolist(), range(len(scores)))
    l.sort(reverse=reverse)
    l = map(list, zip(*l))
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
        print '\nIt took %d minutes and %0.2f seconds to LD-prune' % (t / 60, t % 60)
    return pruning_vector             



def ld_pruning(ld_table, max_ld=0.5, verbose=False):
    """
    Prunes SNPs in LD, in random order. 
    """
    import random
    if verbose:
        print 'Calculating LD table'
    t0 = time.time()
    indices_to_keep = []
    num_snps = len(ld_table)
    indices = sp.random.permutation(num_snps)
    remaining_indices = set(indices)
    for i in indices:
        if len(remaining_indices) == 0:
            break
        elif not (i in remaining_indices):
            continue
        else:
            indices_to_keep.append(i)
            for j in ld_table[i]:
                if ld_table[i][j] > max_ld and j in remaining_indices:
                    remaining_indices.remove(j)
    filter_vector = sp.zeros(num_snps, dtype='bool')
    filter_vector[indices_to_keep] = 1
    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to LD-prune' % (t / 60, t % 60)
    return filter_vector


