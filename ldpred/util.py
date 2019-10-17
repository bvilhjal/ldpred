"""
Various general utility functions.

"""
import scipy as sp
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
chromosomes_list = ['chrom_%d' % (chrom) for chrom in ok_chromosomes]
chrom_name_map = {'X':'23','chr_X':'23','chrom_X':'23'}
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


# LDpred currently ignores the Y and MT chromosomes.
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

def check_chromosomes(missing_chromosomes):
        if len(missing_chromosomes) > 0:
            print('Ignored chromosomes:', ','.join(list(missing_chromosomes)))
            print('Please note that only data on chromosomes 1-23, and X is parsed.')


def calc_auc(y_true, y_hat, show_plot=False):
    """
    Calculate the Area Under the Curve (AUC) for a predicted and observed case-control phenotype.
    """
    y_true = sp.copy(y_true)
    if len(sp.unique(y_true)) == 2:
        y_min = y_true.min()
        y_max = y_true.max()
        if y_min != 0 or y_max != 1:
            print('Transforming back to a dichotomous trait')
            y_true[y_true == y_min] = 0
            y_true[y_true == y_max] = 1
        
    else:
        print('Warning: Calculating AUC for a quantitative phenotype.')
        y_mean = sp.mean(y_true)
        zero_filter = y_true <= y_mean
        one_filter = y_true > y_mean
        y_true[zero_filter] = 0
        y_true[one_filter] = 1

    num_cases = sp.sum(y_true == 1)
    num_controls = sp.sum(y_true == 0)
    assert num_cases + num_controls == len(y_true), 'The phenotype is not defined as expected. It is not binary (0 1 case-control status).'
    print('%d cases, %d controls' % (num_cases, num_controls)) 
    
    num_indivs = float(len(y_true))
    tot_num_pos = float(sp.sum(y_true))
    tot_num_neg = float(num_indivs - tot_num_pos)
        
    l = y_hat.tolist()
    l.sort(reverse=True)
    roc_x = []
    roc_y = []
    auc = 0.0
    prev_fpr = 0.0
    for thres in l:
        thres_filter = y_hat >= thres
        y_t = y_true[thres_filter]
        n = len(y_t)
        tp = sp.sum(y_t)
        fp = n - tp
        
        fpr = fp / tot_num_neg
        tpr = tp / tot_num_pos
        roc_x.append(fpr)
        roc_y.append(tpr)
        delta_fpr = fpr - prev_fpr
        auc += tpr * delta_fpr
        prev_fpr = fpr
    print('AUC: %0.4f' % auc)
    if show_plot:
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
    
def load_hapmap_SNPs():
    f = gzip.open(hm3_file, 'r')
    hm3_sids = pickle.load(f)
    f.close()
    return hm3_sids

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
