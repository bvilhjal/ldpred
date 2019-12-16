"""
Code for handling plink files.

Uses plinkio.

Some of the code is borrowed from plink-pandas

"""
import scipy as sp
from plinkio import plinkfile

import pandas_plink as pdp
import pandas as pd
from pandas.api.types import CategoricalDtype
from collections import OrderedDict as odict
from numpy import int64, float64


def read_bim(fn):
    fn = fn+'.bim'

    header = odict(
        [
            ("chrom", bytes),
            ("snp", bytes),
            ("cm", float64),
            ("pos", int64),
            ("a0", bytes),
            ("a1", bytes),
        ]
    )
    df = pd.read_csv(fn,
                  delim_whitespace=True,
                  header=None,
                  names=header.keys(),
                  dtype=header,
                  compression=None,
                  engine="c",
                  )

    df["chrom"] = df["chrom"].astype(CategoricalDtype(ordered=True))
    df["a0"] = df["a0"].astype("category")
    df["a1"] = df["a1"].astype("category")
    df["i"] = range(df.shape[0])
    df["pos"] = df["pos"].astype(int64)
    df["cm"] = df["cm"].astype(float64)
    return df

    
def read_fam(fn):
    fn = fn+'.fam'

    header = odict(
        [
            ("fid", str),
            ("iid", str),
            ("father", str),
            ("mother", str),
            ("gender", bytes),
            ("trait", float64),
        ]
    )
    
    df = pd.read_csv(fn,
                  delim_whitespace=True,
                  header=None,
                  names=header.keys(),
                  dtype=header,
                  compression=None,
                  engine="c",
                  )


    df["gender"] = df["gender"].astype("category")
    df["i"] = range(df.shape[0])
    df["trait"] = df["trait"].astype(float64)
    return df


def get_chrom_dict(bim_df, chromosomes):
    chr_dict = {}
    for chrom in chromosomes:
        chr_str = 'chrom_%s' % chrom
        chr_dict[chr_str] = bim_df.query('chrom=="%s"'%chrom)
    return chr_dict



def parse_snps(pf,snp_mask, imp_type='mode'):
    bedf = "%s.bed"%(pf)
    bimf = "%s.bim"%(pf)
    famf = "%s.fam"%(pf)
    G = pdp.read_plink1_bin(bedf, bimf, famf, verbose=False)
    #Excluding SNPs not in mask (which is assumed to be ordered)
    G = G[:,snp_mask]
    num_indivs,num_snps = G.shape
    G_T = G.T
    snps=sp.zeros(G_T.shape,dtype='int8')
    imp_frac = 0
    for snp_i, snp in enumerate(G_T.values):
        snp_is_nan = sp.isnan(snp)
        if sp.any(snp_is_nan):
            num_nans = sp.sum(snp_is_nan)
            imp_frac += num_nans
            if imp_type=='mode':
                snp[snp_is_nan]=3
                snp = sp.array(snp,dtype='int8')
                bin_counts = sp.bincount(snp,minlength=4)
                mode_v = sp.argmax(bin_counts[:2])
                snp[snp == 3] = mode_v
            elif imp_type=='random':
                snp[snp_is_nan] = sp.random.choice(snp[~snp_is_nan],size=num_nans)
                snp = sp.array(snp,dtype='int8')
        else:
            snp = sp.array(snp,dtype='int8')
        snps[snp_i]=snp
        
    imp_frac = imp_frac/float(num_indivs*num_snps)
    freqs = sp.sum(snps, 1, dtype='float32') / (2 * float(num_indivs))
    return snps, freqs, imp_frac

# def parse_plink_snps(genotype_file, snp_indices):
#     plinkf = plinkfile.PlinkFile(genotype_file)
#     samples = plinkf.get_samples()
#     num_individs = len(samples)
#     num_snps = len(snp_indices)
#     raw_snps = sp.empty((num_snps, num_individs), dtype='int8')
# 
#     # If these indices are not in order then we place them in the right place while parsing SNPs.
#     snp_order = sp.argsort(snp_indices)
#     ordered_snp_indices = list(snp_indices[snp_order])
#     ordered_snp_indices.reverse()
#     # Iterating over file to load SNPs
#     snp_i = 0
#     next_i = ordered_snp_indices.pop()
#     line_i = 0
#     max_i = ordered_snp_indices[0]
#     imp_frac = 0
#     while line_i <= max_i: 
#         if line_i < next_i:
#             next(plinkf)
#         elif line_i == next_i:
#             line = next(plinkf)
#             snp = sp.array(line, dtype='int8')
#             bin_counts = line.allele_counts()
#             if bin_counts[-1] > 0:
#                 imp_frac += bin_counts[-1]
#                 mode_v = sp.argmax(bin_counts[:2])
#                 snp[snp == 3] = mode_v
#             s_i = snp_order[snp_i]
#             raw_snps[s_i] = snp
#             if line_i < max_i:
#                 next_i = ordered_snp_indices.pop()
#             snp_i += 1
#         line_i += 1
#     plinkf.close()
#     imp_frac = imp_frac/float(num_individs*num_snps)
#     assert snp_i == len(raw_snps), 'Parsing SNPs from plink file failed.'
#     num_indivs = len(raw_snps[0])
#     freqs = sp.sum(raw_snps, 1, dtype='float32') / (2 * float(num_indivs))
#     return raw_snps, freqs, imp_frac

def get_num_indivs(genotype_file):
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    plinkf.close()
    return len(samples)

