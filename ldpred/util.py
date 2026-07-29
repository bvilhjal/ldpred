"""
Various general utility functions.

"""
import numpy as np
from scipy import stats
import pickle
import gzip
import os
from itertools import takewhile
from itertools import repeat
import sys
import re

# LDpred currently ignores the Y and MT chromosomes.
ok_chromosomes = set(range(1, 24))

chrom_name_map = {'X':23,'chr_X':23,'chrom_X':23}
for chrom in ok_chromosomes:
    chrom_name_map['%d' % (chrom)]=chrom
    chrom_name_map['chrom_%d' % (chrom)]=chrom
    chrom_name_map['chr_%d' % (chrom)]=chrom


def get_chrom_num(chrom):
    return chrom_name_map.get(re.sub("chr", "", chrom),0)

#Various auxiliary variables
ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
opp_strand_dict = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}

valid_nts = set(['A', 'T', 'C', 'G'])


# LDpred currently ignores the Y and MT chromosomes.  Note that X is stored as
# chromosome 23 in coordinated data (see chrom_name_map), so 'chrom_X' never
# matches a group; it is kept here only so the ordering is unchanged.
valid_chromosomes = ['%d' % (x) for x in range(1, 24)]
valid_chromosomes.append('X')

chromosomes_list = ['chrom_%s' % (chrom) for chrom in valid_chromosomes]


#Conversion sizes for strings (necessary for using h5py and python 3)
fids_dtype = '|S64'
iids_dtype = '|S64'
sids_dtype = "|S30" 
nts_dtype = "|S1"

sids_u_dtype = '<U30'
nts_u_dtype = '<U1'


my_path = os.path.dirname(os.path.abspath(__file__))
hm3_file = os.path.join(my_path, 'reference','hm3_sids.txt.gz')
lrld_file = os.path.join(my_path, 'reference','long-range-ld-price-2008hg38.txt')

def check_chromosomes(missing_chromosomes):
        if len(missing_chromosomes) > 0:
            print('Ignored chromosomes:', ','.join(list(missing_chromosomes)))
            print('Please note that only data on chromosomes 1-23, and X is parsed.')


def calc_auc(y_true, y_hat, show_plot=False):
    """
    Calculate the Area Under the Curve (AUC) for a predicted and observed case-control phenotype.
    """
    y_true = np.copy(y_true)
    if len(np.unique(y_true)) == 2:
        y_min = y_true.min()
        y_max = y_true.max()
        if y_min != 0 or y_max != 1:
            print('Transforming back to a dichotomous trait')
            y_true[y_true == y_min] = 0
            y_true[y_true == y_max] = 1
        
    else:
        print('Warning: Calculating AUC for a quantitative phenotype.')
        y_mean = np.mean(y_true)
        zero_filter = y_true <= y_mean
        one_filter = y_true > y_mean
        y_true[zero_filter] = 0
        y_true[one_filter] = 1

    num_cases = np.sum(y_true == 1)
    num_controls = np.sum(y_true == 0)
    assert num_cases + num_controls == len(y_true), 'The phenotype is not defined as expected. It is not binary (0 1 case-control status).'
    print('%d cases, %d controls' % (num_cases, num_controls)) 
    
    y_hat = np.asarray(y_hat)

    # The AUC is the probability that a randomly drawn case scores above a
    # randomly drawn control, which is the cases' rank sum rescaled (the
    # Mann-Whitney U statistic).  rankdata averages tied ranks, so ties
    # contribute a half each, and the whole thing is one O(n log n) sort.  The
    # previous version swept every observed score as a threshold and rescanned
    # the full prediction vector at each one, which is O(n^2).
    ranks = stats.rankdata(y_hat)
    auc = ((ranks[y_true == 1].sum() - num_cases * (num_cases + 1) / 2.0)
           / float(num_cases * num_controls))
    print('AUC: %0.4f' % auc)
    if show_plot:
        # Sweep the scores from high to low and accumulate the case/control
        # counts to trace the ROC in one pass.
        order = np.argsort(-y_hat, kind='mergesort')
        y_sorted = y_true[order]
        roc_y = np.concatenate([[0.0], np.cumsum(y_sorted == 1) / float(num_cases)])
        roc_x = np.concatenate([[0.0], np.cumsum(y_sorted == 0) / float(num_controls)])
        import pylab
        pylab.plot(roc_x, roc_y)
        pylab.show()
    return auc



def obs_h2_to_liab(R2_osb,K=0.01,P=0.5):
    """
    Transformation from observed to liability scale.
    
    Lee et al. AJHG 2011 conversion? 
    
    For heritability only
    """
    t = stats.norm.ppf(1-K)
    z = stats.norm.pdf(t)
    c = P*(1-P)*z**2/(K**2*(1-K)**2)
    R2_liab = R2_osb/c
    return R2_liab


def obs_r2_to_liab(R2_osb,K=0.01,P=0.5):
    """
    Lee et al., Gen Epi 2012 conversion
    
    For R2 only

    """
    t = stats.norm.ppf(K)
    z = stats.norm.pdf(t)
    m = z/K
    C = (K*(1-K))**2/((z**2)*(P*(1-P)))
    d =  m*((P-K)/(1-K))
    theta =d**2 - d*t
    R2_liab_cc =  (R2_osb*C)/(1+(R2_osb*C*theta))
    return R2_liab_cc
    
def get_LDpred_sample_size(n,ns,verbose):
    if n is None:
        #If coefficient of variation is very small, then use one N nevertheless.
        n_cv = np.std(ns)/np.mean(ns)
        if n_cv<0.01:
            ldpred_n = np.mean(ns)
            if verbose:
                print ("Sample size does not vary much (CV=%0.4f).  Using a fixed sample size of %0.2f"%(n_cv,ldpred_n))
        else:
            if verbose:
                print ("Using varying sample sizes")
                print ("Sample size ranges between %d and %d"%(min(ns),max(ns)))
                print ("Average sample size is %0.2f "%(np.mean(ns)))
        ldpred_inf_n = np.mean(ns)
        ldpred_n = None
    else:
        ldpred_n = float(n)
        if verbose:
            print ("Using the given fixed sample size of %d"%(n))
        ldpred_inf_n = float(n)
    return ldpred_n,ldpred_inf_n


def load_hapmap_SNPs():
    f = gzip.open(hm3_file, 'r')
    hm3_sids = pickle.load(f)
    f.close()
    return hm3_sids


def load_lrld_dict():
    #Load Price et al. AJHG 2008 long range LD table.
    d = {}
    for chrom in ok_chromosomes:
        d[chrom] = {'reg_dict':{}}
    with open(lrld_file, 'r') as f:
        for line in f:
            l = line.split()
            d[chrom_name_map[l[0]]]['reg_dict'][l[3]] = {'start_pos':int(l[1]), 'end_pos':int(l[2])}
    return d


def get_lrld_chrom_key(chromosome):
    """Resolve a chromosome label to the key used by :func:`load_lrld_dict`.

    Accepts the 'chrom_<n>' group names used in coordinated data, a bare
    number, or an int.  X maps to 23, matching how coordinated data stores it.
    """
    if isinstance(chromosome, str):
        return chrom_name_map[chromosome]
    return int(chromosome)


def is_in_lrld(chrom, pos, lrld_dict):
    for lrld_reg in lrld_dict[get_lrld_chrom_key(chrom)]['reg_dict'].values():
        if lrld_reg['start_pos'] < pos < lrld_reg['end_pos']:
            return True
    return False


def get_snp_lrld_status(chromosome, positions, lrld_dict):
    """Boolean mask marking which positions fall in a long-range LD region."""
    positions = np.asarray(positions)
    snp_lrld = np.zeros(len(positions), dtype=bool)
    for lrld_reg in lrld_dict[get_lrld_chrom_key(chromosome)]['reg_dict'].values():
        snp_lrld |= ((positions > lrld_reg['start_pos'])
                     & (positions < lrld_reg['end_pos']))
    return snp_lrld
    
def is_gz(name):
    return name.lower().endswith(('.gz', '.gzip'))



def count_lines(filename):
    if sys.version_info >= (3,0):
        return count_lines_fast(filename)
    else:
        return count_lines_slow(filename)


def count_lines_fast(filename):
    opener = open
    if is_gz(filename):
        opener = gzip.open
    try:
        with opener(filename, 'rb') as f:
            bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024) for _ in repeat(None)))
            num_lines =sum( buf.count(b'\n') for buf in bufgen )
    except Exception:
        num_lines = -1
    return num_lines 
            

def count_lines_slow(filename):
    opener = open
    if is_gz(filename):
        opener = gzip.open    
    try:
        with opener(filename, 'rb') as f:
            num_lines = sum(1 for line in f)
    except Exception:
        num_lines=-1
    return num_lines
