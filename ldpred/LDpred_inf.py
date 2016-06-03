#!/usr/bin/env python
"""
Implements LDpred-inf.  The method requires the user to have generated a coordinated dataset using coord_genotypes.py


Usage: 
ldpred-inf --coord=COORD_DATA_FILE  --ld_radius=LD_RADIUS   --local_ld_file_prefix=LD_FILE_NAME
                          --N=SAMPLE_SIZE  --out=OUTPUT_FILE_PREFIX  [--H2=HERTIABILITY ]
    
 - COORD_DATA_FILE: The HDF4 file obtained by running the coord_genotypes.py
 
 - LD_RADIUS: An integer number which denotes the number of SNPs on each side of the focal SNP for which LD should be adjusted.  
              A value corresponding M/3000, where M is the number of SNPs in the genome is recommended.
 
 - LD_FILE_NAME: A path and filename prefix for the LD file.  If it doesn't exist, it will be generated.  This can take up to several hours, 
                 depending on LD radius number of SNPs, etc.  If it does exits, that file will be used.
     
 - N: This is the sample size which LDpred assumes was used to calculate the GWAS summary statistics.

 - OUTPUT_FILE_PREFIX:  The prefix of output file.  
 
 - HERTIABILITY (optional): The heritability assumed by LDpred.  By default it estimates the heritability from the GWAS summary statistics
                             using LDscore regression.
 
 
 
 (c) Bjarni J Vilhjalmsson: bjarni.vilhjalmsson@gmail.com
 
 """

import getopt
import sys
import traceback
import h5py
import scipy as sp
from scipy import linalg
import ld
import cPickle
import time
import os
import gzip
import itertools as it
import LDpred 

def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
#    if len(sys.argv) == 1:
#        print __doc__
#        sys.exit(2)

    long_options_list = ['coord=', 'ld_radius=', 'ld_prefix=', 'out=', 'N=', 'H2=',]

    p_dict = {'coord':None, 'ld_radius':None, 'ld_prefix':None, 'out':None, 'N':None, 'H2':None}

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)
    
        except:
            print "Some problems with usage.  Please read the usage documentation carefully."
            traceback.print_exc()
            print __doc__
            sys.exit(2)
    
        for opt, arg in opts:
            if opt == "-h" or opt=="--h" or opt=='--help':
                print __doc__
                sys.exit(0)
            elif opt in ("--coord"): p_dict['coord'] = arg
            elif opt in ("--ld_radius"): p_dict['ld_radius'] = int(arg)
            elif opt in ("--ld_prefix"): p_dict['ld_prefix'] = arg
            elif opt in ("--out"): p_dict['out'] = arg
            elif opt in ("--N"): p_dict['N'] = int(arg)
            elif opt in ("--H2"): p_dict['H2'] = float(arg)
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict




def ldpred_inf(beta_hats, h2=0.1, n=1000, inf_shrink_matrices=None, 
               reference_ld_mats=None, genotypes=None, ld_window_size=100, verbose=False):
    """
    Apply the infinitesimal shrink w LD (which requires LD information).
    
    If reference_ld_mats are supplied, it uses those, otherwise it uses the LD in the genotype data.
    
    If genotypes are supplied, then it assumes that beta_hats and the genotypes are synchronized.

    """
    if verbose:
        print 'Doing LD correction'
    t0 = time.time()
    num_betas = len(beta_hats)
    updated_betas = sp.empty(num_betas)
    m = len(beta_hats)

    for i, wi in enumerate(range(0, num_betas, ld_window_size)):
        start_i = wi
        stop_i = min(num_betas, wi + ld_window_size)
        curr_window_size = stop_i - start_i
        if inf_shrink_matrices!=None:
            A_inv = inf_shrink_matrices[i]
        else:
            if reference_ld_mats != None:
                D = reference_ld_mats[i]
            else:
                if genotypes != None:
                    X = genotypes[start_i: stop_i]
                    num_indivs = X.shape[1]
                    D = sp.dot(X, X.T) / num_indivs
                else:
                    raise NotImplementedError
            A = ((m / h2) * sp.eye(curr_window_size) + (n / (1)) * D)
            A_inv = linalg.pinv(A)
        updated_betas[start_i: stop_i] = sp.dot(A_inv * n , beta_hats[start_i: stop_i])  # Adjust the beta_hats

        if verbose:
            sys.stdout.write('\b\b\b\b\b\b\b%0.2f%%' % (100.0 * (min(1, float(wi + 1) / num_betas))))
            sys.stdout.flush()

    t1 = time.time()
    t = (t1 - t0)
    if verbose:
        print '\nIt took %d minutes and %0.2f seconds to perform the Infinitesimal LD shrink' % (t / 60, t % 60)
    return updated_betas


def ldpred_inf_genomewide(data_file=None, ld_radius = None, ld_dict=None, out_file_prefix=None,
                          n=None, h2=None, verbose=False):
    """
    Calculate LDpred for a genome
    """    
    
    df = h5py.File(data_file,'r')
    has_phenotypes=False
    if 'y' in df.keys():
        'Validation phenotypes found.'
        y = df['y'][...]  # Phenotype
        num_individs = len(y)
        risk_scores_pval_derived = sp.zeros(num_individs)
        has_phenotypes=True

    ld_scores_dict = ld_dict['ld_scores_dict']
    chrom_ref_ld_mats = ld_dict['chrom_ref_ld_mats']
        
    print 'Applying LDpred-inf with LD radius: %d' % ld_radius
    results_dict = {}
    num_snps = 0
    sum_beta2s = 0
    cord_data_g = df['cord_data']

    for chrom_str in cord_data_g.keys():
        g = cord_data_g[chrom_str]
        betas = g['betas'][...]
        n_snps = len(betas)
        num_snps += n_snps
        sum_beta2s += sp.sum(betas ** 2)
        
    L = ld_scores_dict['avg_gw_ld_score']
    chi_square_lambda = sp.mean(n * sum_beta2s / float(num_snps))
    print 'Genome-wide lambda inflation:', chi_square_lambda,
    print 'Genome-wide mean LD score:', L
    gw_h2_ld_score_est = max(0.0001, (max(1, chi_square_lambda) - 1) / (n * (L / num_snps)))
    print 'Estimated genome-wide heritability:', gw_h2_ld_score_est
    
    assert chi_square_lambda>1, 'Something is wrong with the GWAS summary statistics, parsing of them, or the given GWAS sample size (N). Lambda (the mean Chi-square statistic) is too small.  '
    

    if out_file_prefix:
        #Preparing output files
        raw_effect_sizes = []
        ldpred_effect_sizes = []
        sids = []
        chromosomes = []
        positions = []
        nts = []
        
    for chrom_str in cord_data_g.keys():
        g = cord_data_g[chrom_str]
        if has_phenotypes:
            if 'raw_snps_val' in g.keys():
                raw_snps = g['raw_snps_val'][...]
            else:
                raw_snps = g['raw_snps_ref'][...]
        
        snp_stds = g['snp_stds_ref'][...]
        pval_derived_betas = g['betas'][...]
        if out_file_prefix:
            chromosomes.extend([chrom_str]*len(pval_derived_betas))
            positions.extend(g['positions'][...])
            sids.extend(g['sids'][...])
            raw_effect_sizes.extend(g['log_odds'][...])
            nts.extend(g['nts'][...])

        n_snps = len(pval_derived_betas)

        h2_chrom = gw_h2_ld_score_est * (n_snps / float(num_snps))
        #print 'Prior parameters: p=%0.3f, n=%d, m=%d, h2_chrom=%0.4f' % (p, n, n_snps, h2_chrom)
        updated_betas = ldpred_inf(pval_derived_betas, genotypes=None, reference_ld_mats=chrom_ref_ld_mats[chrom_str], 
                                            h2=h2_chrom, n=n, ld_window_size=2*ld_radius, verbose=False)

                
        print 'Calculating scores for Chromosome %s'%((chrom_str.split('_'))[1])
        updated_betas = updated_betas / (snp_stds.flatten())
        ldpred_effect_sizes.extend(updated_betas)
        if has_phenotypes:
            prs = sp.dot(updated_betas, raw_snps)
            risk_scores_pval_derived += prs
            corr = sp.corrcoef(y, prs)[0, 1]
            r2 = corr ** 2
            print 'The R2 prediction accuracy of PRS using %s was: %0.4f' %(chrom_str, r2)

                
    print 'There were %d (SNP) effects' % num_snps
    if has_phenotypes:
        num_indivs = len(y)
        results_dict['y']=y
        results_dict['risk_scores_pd']=risk_scores_pval_derived
        print 'Prediction accuracy was assessed using %d individuals.'%(num_indivs)

        corr = sp.corrcoef(y, risk_scores_pval_derived)[0, 1]
        r2 = corr ** 2
        results_dict['r2_pd']=r2
        print 'The  R2 prediction accuracy (observed scale) for the whole genome was: %0.4f (%0.6f)' % (r2, ((1-r2)**2)/num_indivs)

        if corr<0:
            risk_scores_pval_derived = -1* risk_scores_pval_derived
        auc = LDpred.calc_auc(y,risk_scores_pval_derived)
        print 'AUC for the whole genome was: %0.4f'%auc

        #Now calibration                               
        denominator = sp.dot(risk_scores_pval_derived.T, risk_scores_pval_derived)
        y_norm = (y-sp.mean(y))/sp.std(y)
        numerator = sp.dot(risk_scores_pval_derived.T, y_norm)
        regression_slope = (numerator / denominator)#[0][0]
        print 'The slope for predictions with P-value derived  effects is:',regression_slope
        results_dict['slope_pd']=regression_slope
    
    weights_out_file = '%s.txt'%(out_file_prefix)
    with open(weights_out_file,'w') as f:
        f.write('chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta\n')
        for chrom, pos, sid, nt, raw_beta, ldpred_beta in it.izip(chromosomes, positions, sids, nts, raw_effect_sizes, ldpred_effect_sizes):
            nt1,nt2 = nt[0],nt[1]
            f.write('%s    %d    %s    %s    %s    %0.4e    %0.4e\n'%(chrom, pos, sid, nt1, nt2, raw_beta, ldpred_beta))



def main():
    p_dict = parse_parameters()

    #Use the same LD file as LDpred
    local_ld_dict_file = '%s_ldradius%d.pickled.gz'%(p_dict['ld_prefix'], p_dict['ld_radius'])
    
    print """
Note: For maximal accuracy all SNPs with LDpred weights should be included in the validation data set.
If they are a subset of the validation data set, then we suggest recalculate LDpred for the overlapping SNPs. 
"""
    if not os.path.isfile(local_ld_dict_file):
        df = h5py.File(p_dict['coord'])
                 
        chrom_ld_scores_dict = {}
        chrom_ld_dict = {}
        chrom_ref_ld_mats = {}
        ld_score_sum = 0
        num_snps = 0
        print 'Calculating LD information w. radius %d'% p_dict['ld_radius']

        cord_data_g = df['cord_data']

        for chrom_str in cord_data_g.keys():
            print 'Working on %s'%chrom_str
            g = cord_data_g[chrom_str]
            if 'raw_snps_ref' in g.keys():
                raw_snps = g['raw_snps_ref'][...]
                snp_stds = g['snp_stds_ref'][...]
                snp_means = g['snp_means_ref'][...]
            
            n_snps = len(raw_snps)
            snp_means.shape = (n_snps,1)   
            snp_stds.shape = (n_snps,1)   
            
            # Normalize SNPs..
            snps = sp.array((raw_snps - snp_means)/snp_stds,dtype='float32')
            ret_dict = ld.get_LDpred_ld_tables(snps, ld_radius=p_dict['ld_radius'], ld_window_size=2*p_dict['ld_radius'])
            chrom_ld_dict[chrom_str] = ret_dict['ld_dict']
            chrom_ref_ld_mats[chrom_str] = ret_dict['ref_ld_matrices']
            ld_scores = ret_dict['ld_scores']
            chrom_ld_scores_dict[chrom_str] = {'ld_scores':ld_scores, 'avg_ld_score':sp.mean(ld_scores)}
            ld_score_sum += sp.sum(ld_scores)
            num_snps += n_snps
        avg_gw_ld_score = ld_score_sum / float(num_snps)
        ld_scores_dict = {'avg_gw_ld_score': avg_gw_ld_score, 'chrom_dict':chrom_ld_scores_dict}    
        
        print 'Done calculating the LD table and LD score, writing to file:', local_ld_dict_file
        print 'Genome-wide average LD score was:', ld_scores_dict['avg_gw_ld_score']
        ld_dict = {'ld_scores_dict':ld_scores_dict, 'chrom_ld_dict':chrom_ld_dict, 'chrom_ref_ld_mats':chrom_ref_ld_mats}
        with gzip.open(local_ld_dict_file, 'wb') as f:
            cPickle.dump(ld_dict, f, protocol=2)
        print 'LD information is now pickled.'
    else:
        print 'Loading LD information from file: %s'%local_ld_dict_file
        with gzip.open(local_ld_dict_file, 'r') as f:
            ld_dict = cPickle.load(f)
    
    ldpred_inf_genomewide(data_file=p_dict['coord'], out_file_prefix=p_dict['out'], ld_radius=p_dict['ld_radius'], 
                          ld_dict = ld_dict, n=p_dict['N'], h2=p_dict['H2'], verbose=False)
            
        

if __name__ == '__main__':
    main()
            

