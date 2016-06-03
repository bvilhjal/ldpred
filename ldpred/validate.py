#!/usr/bin/env python
"""
    Takes LDpred.py (or LD_pruning_thres.py) effect estimates, and (validation) genotypes in PLINK bed format as input.  
    The script then works out overlap and outputs predictions or risk scores as well as some prediction 
    accuracy statistics.
    
    Note that for maximal accuracy all SNPs with LDpred weights should be included in the validation dataset.
    If they are a subset of the validation dataset, then we suggest recalculate LDpred for the overlapping SNPs.
    
    
    

Usage: 
validate --vgf=PLINK_VAL_GENOTYPE_FILE  --rf=RESULT_FILE_PREFIX  --out=OUTPUT_FILE_PREFIX  [--res_format=LDPRED 
                    --split_by_chrom --pf=PHEN_FILE --pf_format=STANDARD --cov_file=COVARIATE_FILE --pcs_file=PCS_FILE 
                    --PS=FRACTIONS_CAUSAL  --TS=PVAL_THRESHOLDS --adjust_for_sex  --adjust_for_covariates  --adjust_for_pcs]
    
 
 - PLINK_VAL_GENOTYPE_FILE: PLINK formatted genotypes for which we want to calculate risk scores.
 
 - RESULT_FILE_PREFIX: SNP weights file, e.g. LDpred SNP weights.

 - OUTPUT_FILE_PREFIX:  The prefix of output file.  

 - RESULT_FILE_FORMAT: The format to expect the results to be in.  The default format is LDPRED, which refers to the format which
   running LDpred output. LDPRED-INF and P+T (LD-pruning + p-value thresholding) are also implemented.
 
 - PHEN_FILE: Is a file with individual IDs and phenotypes
 
 - PVAL_THRESHOLDS: This option is only valid if a P+T result file prefix is supplied.  It's a list of p-value thresholds, 
                    separated by a comma (without space), to be used for LDpred. Default values are 
                    --TS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,3E-5,1E-5,1E-6,1E-7,1E-8

 - FRACTIONS_CAUSAL: This option is only valid if a LDPRED result file prefix is supplied.  A list of comma separated 
                     (without space) values between 1 and 0, excluding 0.  1 corresponds to the infinitesimal model 
                     and will yield results similar to LDpred-inf.  Default values are 
                     --PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001 
 
 2015 (c) Bjarni J Vilhjalmsson: bjarni.vilhjalmsson@gmail.com
 
 """

import getopt
import sys
import os
import traceback
import scipy as sp
from scipy import linalg
from plinkio import plinkfile
import itertools as it
import time 
import h5py

ok_nts = ['A','T','C','G']
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}



def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
    long_options_list = ['vgf=', 'rf=', 'res_format=', 'out=', 'indiv_filter=', 'split_by_chrom', 'pf=', 'pf_format=', 'cov_file=', 
                         'pcs_file=', 'PS=', 'TS=', 'adjust_for_sex','adjust_for_covariates', 'adjust_for_pcs', 'h','help']
    
    p_dict = {'vgf':None, 'rf':None, 'out':None, 'res_format':'LDPRED', 'indiv_filter':None, 'split_by_chrom':False, 
              'pf':None, 'pf_format':'STANDARD', 'cov_file':None, 'pcs_file':None, 'PS':[1,0.3,0.1,0.03,0.01,0.003,0.001],
              'TS':[1,0.3,0.1,0.03,0.01,0.003,0.001,3*1E-4,1E-4,3*1E-5,1E-5,1E-6,1E-7,1E-8],
              'adjust_for_sex':False, 'adjust_for_covariates':False, 'adjust_for_pcs':False}

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)
    
        except:
            print "Some problems with parameters.  Please read the usage documentation carefully."
            print "Use the -h option for usage information."
#             traceback.print_exc()
#             print __doc__
            sys.exit(2)
    
        for opt, arg in opts:
            if opt == "-h" or opt=="--h" or opt=='--help':
                print __doc__
                sys.exit(0)
            elif opt in ("--vgf"): p_dict['vgf'] = arg
            elif opt in ("--rf"): p_dict['rf'] = arg
            elif opt in ("--res_format"): p_dict['res_format'] = arg
            elif opt in ("--indiv_filter"): p_dict['indiv_filter'] = arg
            elif opt in ("--out"): p_dict['out'] = arg
            elif opt in ("--split_by_chrom"): p_dict['split_by_chrom'] = True
            elif opt in ("--PS"): p_dict['PS'] = map(float,arg.split(','))
            elif opt in ("--TS"): p_dict['TS'] = map(float,arg.split(','))
            elif opt in ("--pf"): p_dict['pf'] = arg
            elif opt in ("--pf_format"): p_dict['pf_format'] = arg
            elif opt in ("--cov_file"): p_dict['cov_file'] = arg
            elif opt in ("--pcs_file"): p_dict['pcs_file'] = arg
            elif opt in ("--adjust_for_sex"): p_dict['adjust_for_sex'] = True
            elif opt in ("--adjust_for_covariates"): p_dict['adjust_for_covariates'] = True
            elif opt in ("--adjust_for_pcs"): p_dict['adjust_for_pcs'] = True
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict


def get_prs(genotype_file, rs_id_map, phen_map=None):
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    
    #1. Figure out indiv filter and get true phenotypes
    indiv_filter=sp.zeros(len(samples),dtype='bool8')
    true_phens = []
    iids = []
    if phen_map is not None:
        pcs = []
        sex = []
        covariates = []
        phen_iids = set(phen_map.keys())
        for samp_i, sample in enumerate(samples):
            if sample.iid in phen_iids:
                indiv_filter[samp_i] = True
                true_phens.append(phen_map[sample.iid]['phen'])
                iids.append(sample.iid)
                if 'pcs' in phen_map[sample.iid].keys():
                    pcs.append(phen_map[sample.iid]['pcs'])
                if 'sex' in phen_map[sample.iid].keys():
                    sex.append(phen_map[sample.iid]['sex'])
                if 'covariates' in phen_map[sample.iid].keys():
                    #Temp hack...
#                     if phen_map[sample.iid]['sex']==1:
#                         covariates.append([phen_map[sample.iid]['covariates'][0],0])
#                     else:
#                         covariates.append([0,phen_map[sample.iid]['covariates'][0]])
                    covariates.append(phen_map[sample.iid]['covariates'])
        if len(pcs)>0:
            assert len(pcs)==len(true_phens), 'PC information missing for some individuals with phenotypes'
        if len(sex)>0:
            assert len(sex)==len(true_phens), 'Sex information missing for some individuals with phenotypes'
        if len(covariates)>0:
            assert len(covariates)==len(true_phens), 'Covariates missing for some individuals with phenotypes'
    else:
        for samp_i, sample in enumerate(samples):
            if sample.affection!=2:
                indiv_filter[samp_i] = True
                true_phens.append(sample.affection)
                iids.append(sample.iid)
                
    num_individs = sp.sum(indiv_filter)
    assert num_individs>0, 'Issues in parsing the phenotypes and/or PCs?'

    assert not sp.any(sp.isnan(true_phens)),'Phenotypes appear to have some NaNs, or parsing failed.'
                
    print '%d individuals have phenotype and genotype information.'%num_individs
    
    num_non_matching_nts = 0
    num_flipped_nts = 0
    
    raw_effects_prs = sp.zeros(num_individs)
    pval_derived_effects_prs = sp.zeros(num_individs)
    #If these indices are not in order then we place them in the right place while parsing SNPs.
    print 'Iterating over BED file to calculate risk scores.'
    locus_list = plinkf.get_loci()
    snp_i = 0
    
    
    for locus, row in it.izip( locus_list, plinkf):
        upd_pval_beta = 0
        try:
            #Check rs-ID
#             sid = '%d_%d'%(locus.chromosome,locus.bp_position)
            sid = locus.name
            rs_info = rs_id_map[sid]
        except Exception: #Move on if rsID not found.
            continue
        
        if rs_info['upd_pval_beta']==0:
            continue
        
        #Check whether the nucleotides are OK, and potentially flip it.
        ss_nt = rs_info['nts']
        g_nt =  [locus.allele1,locus.allele2]
        flip_nts = False
        os_g_nt = sp.array([opp_strand_dict[g_nt[0]], opp_strand_dict[g_nt[1]]])
        if not (sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)):
            # Opposite strand nucleotides
            flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
            if flip_nts:
                raw_beta = -rs_info['raw_beta']
                upd_pval_beta = -rs_info['upd_pval_beta']
                num_flipped_nts+=1
            else:
                #print "Nucleotides don't match after all?: sid=%s, g_nt=%s, ss_nt=%s" % (locus.name, str(g_nt), str(ss_nt))
                num_non_matching_nts += 1
                continue
        else:
            raw_beta = rs_info['raw_beta']
            upd_pval_beta = rs_info['upd_pval_beta']
        
        #Parse SNP, and fill in the blanks if necessary.
        snp = sp.array(row, dtype='int8')[indiv_filter]
        bin_counts = row.allele_counts()
        if bin_counts[-1]>0:
            mode_v = sp.argmax(bin_counts[:2])
            snp[snp==3] = mode_v
        
        #Normalize SNP 
#         n_snp = (snp - sp.mean(snp))/sp.std(snp)

        
        #Update scores and move on.
        raw_effects_prs += snp*raw_beta
        assert not sp.any(sp.isnan(raw_effects_prs)),'Raw effects PRS is corrupted'
        pval_derived_effects_prs += snp*upd_pval_beta
        assert not sp.any(sp.isnan(pval_derived_effects_prs)),'Weighted effects PRS is corrupted'
        
        if snp_i>0 and snp_i%100000==0:
            print snp_i 
            print 'Number of non-matching NTs: %d'%num_non_matching_nts
            raw_eff_r2 = (sp.corrcoef(raw_effects_prs, true_phens)[0,1])**2
            pval_eff_r2  = (sp.corrcoef(pval_derived_effects_prs, true_phens)[0,1])**2
            print 'Raw effects PRS r2: %0.4f'%raw_eff_r2
            print 'Weigted effects PRS r2: %0.4f'%pval_eff_r2
        
        snp_i +=1

    plinkf.close()

    
    print "DONE!"
    print 'Number of non-matching NTs: %d'%num_non_matching_nts
    print 'Number of flipped NTs: %d'%num_flipped_nts
    raw_eff_corr = sp.corrcoef(raw_effects_prs, true_phens)[0,1]
    raw_eff_r2 = raw_eff_corr**2
    pval_eff_corr = sp.corrcoef(pval_derived_effects_prs, true_phens)[0,1]
    pval_eff_r2  = pval_eff_corr**2
    
    print 'Raw effects PRS correlation: %0.4f'%raw_eff_corr
    print 'Raw effects PRS r2: %0.4f'%raw_eff_r2
    print 'Weigted effects PRS correlation: %0.4f'%pval_eff_corr
    print 'Weigted effects PRS r2: %0.4f'%pval_eff_r2

    ret_dict = {'raw_effects_prs':raw_effects_prs.copy(), 'pval_derived_effects_prs':pval_derived_effects_prs.copy(), 
                'true_phens':true_phens[:], 'iids':iids}

    if len(pcs)>0:
        ret_dict['pcs'] = pcs
    if len(sex)>0:
        ret_dict['sex'] = sex
    if len(covariates)>0:
        ret_dict['covariates'] = covariates

    return ret_dict




def parse_phen_file(pf, pf_format):
    print pf
    phen_map ={}
    if pf!=None:
        if pf_format=='FAM':
            """
            Individual's family ID ('FID')
            Individual's within-family ID ('IID'; cannot be '0')
            Within-family ID of father ('0' if father isn't in dataset)
            Within-family ID of mother ('0' if mother isn't in dataset)
            Sex code ('1' = male, '2' = female, '0' = unknown)
            Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
            """
            print 'Parsing phenotypes'
            with open(pf,'r') as f:
                for line in f:
                    l = line.split()                     
                    iid = l[1] 
                    #iid = iid[:-4]
                    sex = int(l[4])
                    phen = float(l[5])  
                    if sex!=0 and phen !=-9:
                        phen_map[iid] = {'phen':phen, 'sex':sex}
                    
            iids = set(phen_map.keys())

        if pf_format=='STANDARD':
            """
            IID   PHE
            """
            print 'Parsing phenotypes'
            with open(pf,'r') as f:
                for line in f:
                    l = line.split()                     
                    iid = l[0]
                    phen = float(l[1])                     
                    phen_map[iid] = {'phen':phen}
                    
            iids = set(phen_map.keys())
#         elif pf_format=='Other':
#             print 'Parse phenotype file: %s'%p_dict['pf']
#             phen_map ={}
#             """
#             FID     IID     SEX     PHE     SCZ     BIP     BOTH
#             00C04395        00C04395        2       SCZ     2       -9      2
#             00C04941        00C04941        2       SCZ     2       -9      2
#             01C05110        01C05110        2       SCZ     2       -9      2
#             01C05278        01C05278        2       SCZ     2       -9      2
#             01C05402        01C05402        1       SCZ     2       -9      2
#             01C05566        01C05566        1       BIP     -9      2       2
#      
#             """
#             with open(p_dict['pf'],'r') as f:
#                 print f.next()
#                 for line in f:
#                     l = line.split()
#                     if p_dict['phen']=='SCZ':
#                         phen = float(l[4])
#                     elif p_dict['phen']=='BIP':
#                         phen = float(l[5])
#                     elif p_dict['phen']=='BOTH':
#                         phen = float(l[6])
#                      
#                     if phen ==-9:
#                         continue
#                      
#                     fid = l[0]
#                     iid = l[1]
#                     sex = int(l[2])
#                     phen_map[iid] = {'fid':fid, 'sex':sex, 'phen':phen}
        
        elif pf_format=='S2':
            """
            IID Age Sex Height_Inches
            """
            with open(pf,'r') as f:
                print f.next()
                for line in f:
                    l = line.split()                     
                    iid = l[0]
                    age = float(l[1])                     
                    if l[2]=='Male':
                        sex = 1 
                    elif l[2]=='Female':
                        sex = 2
                    else:
                        raise Exception('Sex missing')
                    phen = float(l[3])                     
                    phen_map[iid] = {'phen':phen, 'age':age, 'sex':sex}
                    
    return phen_map


def parse_ldpred_res(file_name):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta    ldpred_beta
    1, 798959, rs11240777, C, T, -1.1901e-02, 3.2443e-03, 2.6821e-04
    """
    with open(file_name,'r') as f:
        f.next()
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1,nt2]
            raw_beta = float(l[5])
            upd_pval_beta = float(l[6])
            rs_id_map[rs_id] = {'chrom':chrom, 'pos':pos, 'nts':nts, 'raw_beta':raw_beta, 
                                'upd_pval_beta':upd_pval_beta}
    return rs_id_map

def parse_pt_res(file_name):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    raw_pval_beta    upd_beta    upd_pval_beta
    1    798959    rs11240777    C    T    -1.1901e-02    -1.1901e-02    2.6821e-04    2.6821e-04
    """
    with open(file_name,'r') as f:
        f.next()
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1,nt2]
            raw_beta = float(l[5])
            upd_pval_beta = float(l[8])
            if raw_beta!=0:
                rs_id_map[rs_id] = {'chrom':chrom, 'pos':pos, 'nts':nts, 'raw_beta':raw_beta, 
                                    'upd_pval_beta':upd_pval_beta}
                non_zero_chromosomes.add(chrom)
                
    return rs_id_map

def calc_risk_scores(bed_file, rs_id_map, phen_map, out_file=None, split_by_chrom=False, adjust_for_sex=False,
                     adjust_for_covariates=False, adjust_for_pcs=False):
    print 'Parsing PLINK bed file: %s'%bed_file
    num_individs = len(phen_map)
    assert num_individs>0, 'No individuals found.  Problems parsing the phenotype file?' 
         
    if split_by_chrom:
        raw_effects_prs = sp.zeros(num_individs)
        pval_derived_effects_prs = sp.zeros(num_individs)
        
        for i in range(1,23):
            if i in non_zero_chromosomes:
                genotype_file=bed_file+'_%i_keep'%i
                if os.path.isfile(genotype_file+'.bed'):
                    print 'Working on chromosome %d'%i
                    prs_dict = get_prs(genotype_file, rs_id_map, phen_map)
                    
                    raw_effects_prs += prs_dict['raw_effects_prs']
                    pval_derived_effects_prs += prs_dict['pval_derived_effects_prs']
    #                 raw_eff_r2 = (sp.corrcoef(raw_effects_prs, prs_dict['true_phens'])[0,1])**2
    #                 pval_eff_r2  = (sp.corrcoef(pval_derived_effects_prs, prs_dict['true_phens'])[0,1])**2
    #                 print 'Overall raw effects PRS r2: %0.4f'%raw_eff_r2
    #                 print 'Overall weigted effects PRS r2: %0.4f'%pval_eff_r2
            else:
                print 'Skipping chromosome'

                    
                
    else:
        prs_dict = get_prs(bed_file, rs_id_map, phen_map)
        raw_effects_prs = prs_dict['raw_effects_prs']
        pval_derived_effects_prs = prs_dict['pval_derived_effects_prs']
        true_phens = prs_dict['true_phens']
     
    #Report prediction accuracy
    raw_eff_corr = sp.corrcoef(raw_effects_prs, prs_dict['true_phens'])[0,1]
    raw_eff_r2 = raw_eff_corr**2
    pval_eff_corr = sp.corrcoef(pval_derived_effects_prs, prs_dict['true_phens'])[0,1]
    pval_eff_r2  = pval_eff_corr**2
    
    print 'Final raw effects PRS correlation: %0.4f'%raw_eff_corr
    print 'Final raw effects PRS r2: %0.4f'%raw_eff_r2
    print 'Final weighted effects PRS correlation: %0.4f'%pval_eff_corr
    print 'Final weighted effects PRS r2: %0.4f'%pval_eff_r2

    res_dict = {'pred_r2':pval_eff_r2}

    raw_effects_prs.shape = (len(raw_effects_prs), 1)
    pval_derived_effects_prs.shape = (len(pval_derived_effects_prs), 1)
    true_phens = sp.array(true_phens)
    true_phens.shape = (len(true_phens), 1)
    
    #Store covariate weights, slope, etc.
    weights_dict = {}
    
    #Store Adjusted predictions
    adj_pred_dict = {}

    #Direct effect
    Xs = sp.hstack([pval_derived_effects_prs, sp.ones((len(true_phens), 1))])
    (betas, rss00, r, s) = linalg.lstsq(sp.ones((len(true_phens), 1)), true_phens)
    (betas, rss, r, s) = linalg.lstsq(Xs, true_phens)
    pred_r2 = 1 - rss / rss00
#     print 'Fitted effects (betas) for PRS, and intercept on true phenotype:',betas
    weights_dict['unadjusted']={'Intercept':betas[1][0],'ldpred_prs_effect':betas[0][0]}
#     print pred_r2
    
    #Adjust for sex
    if adjust_for_sex and 'sex' in prs_dict and len(prs_dict['sex'])>0:
        sex = sp.array(prs_dict['sex'])
        sex.shape = (len(sex),1)
        (betas, rss0, r, s) = linalg.lstsq(sp.hstack([sex, sp.ones((len(true_phens), 1))]), true_phens)
        (betas, rss, r, s) = linalg.lstsq(sp.hstack([raw_effects_prs, sex, sp.ones((len(true_phens), 1))]), true_phens)
        Xs = sp.hstack([pval_derived_effects_prs, sex, sp.ones((len(true_phens), 1))])
        (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
        weights_dict['sex_adj']={'Intercept':betas[2][0],'ldpred_prs_effect':betas[0][0], 'sex':betas[1][0]}
        print 'Fitted effects (betas) for PRS, sex, and intercept on true phenotype:',betas
        adj_pred_dict['sex_adj'] = sp.dot(Xs,betas)
        pred_r2 = 1 - rss / rss0
        print 'Sex adjusted prediction accuracy (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        pred_r2 = 1 - rss / rss00
        print 'Sex adjusted prediction + Sex (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        pred_r2 = 1 - rss_pd / rss0
        print 'Sex adjusted prediction accuracy (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        res_dict['PC_adj_pred_r2'] = pred_r2
        pred_r2 = 1 - rss_pd / rss00
        print 'Sex adjusted prediction + Sex (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        res_dict['PC_adj_pred_r2+PC'] = pred_r2
    
    #Adjust for PCs
    if adjust_for_pcs and 'pcs' in prs_dict and len(prs_dict['pcs'])>0:
        pcs = prs_dict['pcs']
        (betas, rss0, r, s) = linalg.lstsq(sp.hstack([pcs, sp.ones((len(true_phens), 1))]), true_phens)
        (betas, rss, r, s) = linalg.lstsq(sp.hstack([raw_effects_prs, pcs, sp.ones((len(true_phens), 1))]), true_phens)
        Xs = sp.hstack([pval_derived_effects_prs, sp.ones((len(true_phens), 1)),pcs])
        (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
        weights_dict['pc_adj']={'Intercept':betas[1][0],'ldpred_prs_effect':betas[0][0], 'pcs':betas[2][0]}
        adj_pred_dict['pc_adj'] = sp.dot(Xs,betas)
        pred_r2 = 1 - rss / rss0
        print 'PC adjusted prediction accuracy (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        pred_r2 = 1 - rss / rss00
        print 'PC adjusted prediction + PCs (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        pred_r2 = 1 - rss_pd / rss0
        print 'PC adjusted prediction accuracy (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        res_dict['PC_adj_pred_r2'] = pred_r2
        pred_r2 = 1 - rss_pd / rss00
        print 'PC adjusted prediction + PCs (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        res_dict['PC_adj_pred_r2+PC'] = pred_r2
        
        #Adjust for both PCs and Sex
        if adjust_for_sex and 'sex' in prs_dict and len(prs_dict['sex'])>0:
            sex = sp.array(prs_dict['sex'])
            sex.shape = (len(sex),1)
            (betas, rss0, r, s) = linalg.lstsq(sp.hstack([sex, pcs, sp.ones((len(true_phens), 1))]), true_phens)
            (betas, rss, r, s) = linalg.lstsq(sp.hstack([raw_effects_prs, sex, pcs, sp.ones((len(true_phens), 1))]), true_phens)
            Xs = sp.hstack([pval_derived_effects_prs, sex, sp.ones((len(true_phens), 1)), pcs])
            (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
            weights_dict['sex_pc_adj']={'Intercept':betas[2][0],'ldpred_prs_effect':betas[0][0], 'sex':betas[1][0], 'pcs':betas[3][0]}
            adj_pred_dict['sex_pc_adj'] = sp.dot(Xs,betas)
            pred_r2 = 1 - rss / rss0
            print 'PCs+Sex adjusted prediction accuracy (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
            pred_r2 = 1 - rss / rss00
            print 'PCs+Sex adjusted prediction and PCs+Sex (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
            pred_r2 = 1 - rss_pd / rss0
            print 'PCs+Sex adjusted prediction accuracy (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
            res_dict['PC_Sex_adj_pred_r2'] = pred_r2
            pred_r2 = 1 - rss_pd / rss00
            print 'PCs+Sex adjusted prediction and PCs+Sex (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
            res_dict['PC_Sex_adj_pred_r2+PC_Sex'] = pred_r2
            

    #Adjust for covariates
    if adjust_for_covariates and 'covariates' in prs_dict and len(prs_dict['covariates'])>0:
        covariates = prs_dict['covariates']
        (betas, rss0, r, s) = linalg.lstsq(sp.hstack([covariates, sp.ones((len(true_phens), 1))]), true_phens)
        (betas, rss, r, s) = linalg.lstsq(sp.hstack([raw_effects_prs, covariates, sp.ones((len(true_phens), 1))]), true_phens)
        Xs = sp.hstack([pval_derived_effects_prs, covariates, sp.ones((len(true_phens), 1))])
        (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
        adj_pred_dict['cov_adj'] = sp.dot(Xs,betas)
        pred_r2 = 1 - rss / rss0
        print 'Cov adjusted prediction accuracy (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        pred_r2 = 1 - rss / rss00
        print 'Cov adjusted prediction + Cov (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        pred_r2 = 1 - rss_pd / rss0
        print 'Cov adjusted prediction accuracy (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        res_dict['Cov_adj_pred_r2'] = pred_r2
        pred_r2 = 1 - rss_pd / rss00
        print 'Cov adjusted prediction + Cov (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
        res_dict['Cov_adj_pred_r2+Cov'] = pred_r2
        
        if  adjust_for_pcs and 'pcs' in prs_dict and len(prs_dict['pcs']) and 'sex' in prs_dict and len(prs_dict['sex'])>0:
            pcs = prs_dict['pcs']
            sex = sp.array(prs_dict['sex'])
            sex.shape = (len(sex),1)
            (betas, rss0, r, s) = linalg.lstsq(sp.hstack([covariates, sex, pcs, sp.ones((len(true_phens), 1))]), true_phens)
            (betas, rss, r, s) = linalg.lstsq(sp.hstack([raw_effects_prs, covariates, sex, pcs, sp.ones((len(true_phens), 1))]), true_phens)
            Xs = sp.hstack([pval_derived_effects_prs, covariates, sex, pcs, sp.ones((len(true_phens), 1))])
            (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
            adj_pred_dict['cov_sex_pc_adj'] = sp.dot(Xs,betas)
            pred_r2 = 1 - rss / rss0
            print 'Cov+PCs+Sex adjusted prediction accuracy (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
            pred_r2 = 1 - rss / rss00
            print 'Cov+PCs+Sex adjusted prediction and PCs+Sex (R^2) for the whole genome PRS with raw effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
            pred_r2 = 1 - rss_pd / rss0
            print 'Cov+PCs+Sex adjusted prediction accuracy (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
            res_dict['Cov_PC_Sex_adj_pred_r2'] = pred_r2
            pred_r2 = 1 - rss_pd / rss00
            print 'Cov+PCs+Sex adjusted prediction and PCs+Sex (R^2) for the whole genome PRS with weighted effects was: %0.4f (%0.6f)' % (pred_r2, (1-pred_r2)/sp.sqrt(num_individs))
            res_dict['Cov_PC_Sex_adj_pred_r2+Cov_PC_Sex'] = pred_r2

#     print sp.corrcoef(true_phens.T,adj_pred_dict['cov_sex_pc_adj'].T)**2

    #Now calibration
    y_norm = (true_phens-sp.mean(true_phens))/sp.std(true_phens)
    denominator = sp.dot(raw_effects_prs.T, raw_effects_prs)
    numerator = sp.dot(raw_effects_prs.T, y_norm)
    regression_slope = (numerator / denominator)[0][0]
    print 'The slope for predictions with raw effects is:',regression_slope
                           
    denominator = sp.dot(pval_derived_effects_prs.T, pval_derived_effects_prs)
    numerator = sp.dot(pval_derived_effects_prs.T, y_norm)
    regression_slope = (numerator / denominator)[0][0]
    print 'The slope for predictions with weighted effects is:',regression_slope
 
     
#     print sp.corrcoef(prs_dict['raw_effects_prs'], prs_dict['true_phens'])[0,1]
#     print sp.corrcoef(prs_dict['pval_derived_effects_prs'], prs_dict['true_phens'])[0,1]
    num_individs = len(prs_dict['pval_derived_effects_prs'])
     
    #Write PRS out to file.
    if out_file!=None:
        with open(out_file,'w') as f:
            out_str = 'IID, true_phens, raw_effects_prs, pval_derived_effects_prs'
            if 'sex' in prs_dict:
                out_str = out_str+', sex'                
            if 'pcs' in prs_dict:            
                pcs_str = ', '.join(['PC%d'%(1+pc_i) for pc_i in range(len(prs_dict['pcs'][0]))])
                out_str = out_str+', '+pcs_str
            out_str += '\n'
            f.write(out_str)
            for i in range(num_individs):
                out_str = '%s, %0.6e, %0.6e, %0.6e, '%(prs_dict['iids'][i], prs_dict['true_phens'][i],raw_effects_prs[i],
                                                       pval_derived_effects_prs[i])
                if 'sex' in prs_dict:
                    out_str = out_str + '%d, '%prs_dict['sex'][i] 
                if 'pcs' in prs_dict:
                    pcs_str = ', '.join(map(str,prs_dict['pcs'][i]))
                    out_str = out_str + pcs_str 
                out_str += '\n'
                f.write(out_str)
 
        if len(adj_pred_dict.keys())>0:
            with open(out_file+'.adj','w') as f:
                adj_prs_labels = adj_pred_dict.keys()
                out_str = 'IID, true_phens, raw_effects_prs, pval_derived_effects_prs, '+', '.join(adj_prs_labels)
                out_str += '\n'
                f.write(out_str)
                for i in range(num_individs):
                    out_str = '%s, %0.6e, %0.6e, %0.6e'%(prs_dict['iids'][i], prs_dict['true_phens'][i],raw_effects_prs[i],
                                                           pval_derived_effects_prs[i])
                    for adj_prs in adj_prs_labels:
                        out_str += ', %0.4f'%adj_pred_dict[adj_prs][i]
                    out_str += '\n'
                    f.write(out_str)
        if weights_dict!=None:
            oh5f = h5py.File(out_file+'.weights.hdf5','w')
            for k1 in weights_dict.keys():
                kg = oh5f.create_group(k1)
                for k2 in weights_dict[k1]:
                    kg.create_dataset(k2,data=sp.array(weights_dict[k1][k2]))
            oh5f.close()
    return res_dict


def main():
    p_dict = parse_parameters()
    non_zero_chromosomes =set()

    #Parse phenotypes
    if p_dict['pf'] is None:
        if p_dict['vgf'] is not None:
            phen_map = parse_phen_file(p_dict['vgf']+'.fam', 'FAM')
        else:
            raise Exception('Validation phenotypes were not found.')
    else:    
        phen_map = parse_phen_file(p_dict['pf'], p_dict['pf_format'])
    iids = set(phen_map.keys())            
    
                        

    if p_dict['cov_file']!=None:
        print 'Parsing additional covariates'

        with open(p_dict['cov_file'],'r') as f:
            num_missing = 0
            for line in f:
                l = line.split()
                iid = l[0]
                if iid in phen_map:
                    covariates = map(float,l[1:])
                    phen_map[iid]['covariates']=covariates
                else:
                    num_missing +=1
            if num_missing>0:
                print 'Unable to find %d iids in phen file!'%num_missing
                    
    if p_dict['pcs_file']:
        print 'Parsing PCs'

        with open(p_dict['pcs_file'],'r') as f:
            num_missing = 0
            for line in f:
                l = line.split()
                iid = l[1]
                if iid in phen_map:
                    pcs = map(float,l[2:])
                    phen_map[iid]['pcs']=pcs
                else:
                    num_missing +=1
            if num_missing>0:
                print 'Unable to find %d iids in phen file!'%num_missing    
    
    
    num_individs = len(phen_map)
    assert num_individs>0, 'No phenotypes were found!' 
    
    res_dict = {}    
    if p_dict['res_format']=='LDPRED':
        weights_file = '%s_LDpred-inf.txt'%(p_dict['rf'])
        if os.path.isfile(weights_file):
            print ''
            print 'Calculating LDpred-inf risk scores'
            rs_id_map = parse_ldpred_res(weights_file)       
            out_file = '%s_LDpred-inf.txt'%(p_dict['out'])
            calc_risk_scores(p_dict['vgf'], rs_id_map, phen_map, out_file=out_file, split_by_chrom=p_dict['split_by_chrom'], 
                             adjust_for_sex=p_dict['adjust_for_sex'], adjust_for_covariates=p_dict['adjust_for_covariates'], 
                             adjust_for_pcs=p_dict['adjust_for_pcs'])        
        
        for p in p_dict['PS']:
            weights_file = '%s_LDpred_p%0.4e.txt'%(p_dict['rf'], p)
            if os.path.isfile(weights_file):
                print ''
                print 'Calculating LDpred risk scores using p=%0.3e'%p
                rs_id_map = parse_ldpred_res(weights_file)       
                out_file = '%s_LDpred_p%0.4e.txt'%(p_dict['out'], p)
                method_str = 'LDpred_p%0.4e'%(p)
                res_dict[method_str]=calc_risk_scores(p_dict['vgf'], rs_id_map, phen_map, out_file=out_file, 
                                                      split_by_chrom=p_dict['split_by_chrom'], adjust_for_sex=p_dict['adjust_for_sex'], 
                                                      adjust_for_covariates=p_dict['adjust_for_covariates'], 
                                                      adjust_for_pcs=p_dict['adjust_for_pcs'])        

        #Plot results?

    elif p_dict['res_format']=='P+T':
        weights_file = '%s_all_snps.txt'%(p_dict['rf'])
        if os.path.isfile(weights_file):
            print ''
            print 'Calculating risk scores using all SNPs'
            rs_id_map = parse_ldpred_res(weights_file)       
            out_file = '%s_all_snps.txt'%(p_dict['out'])
            res_dict['all_snps'] = calc_risk_scores(p_dict['vgf'], rs_id_map, phen_map, out_file=out_file, 
                                                    split_by_chrom=p_dict['split_by_chrom'], adjust_for_sex=p_dict['adjust_for_sex'], 
                                                    adjust_for_covariates=p_dict['adjust_for_covariates'], 
                                                    adjust_for_pcs=p_dict['adjust_for_pcs'])        
        
        for p_thres in p_dict['TS']:
            weights_file = '%s_P+T_p%0.4e.txt'%(p_dict['rf'], p_thres)
            print weights_file
            if os.path.isfile(weights_file):
                print ''
                print 'Calculating P+T risk scores using p-value threshold of %0.3e'%p_thres
                rs_id_map = parse_pt_res(weights_file)       
                out_file = '%s_P+T_p%0.4e.txt'%(p_dict['out'], p_thres)
                method_str = 'P+T_p%0.4e'%(p_thres)
                res_dict[method_str] = calc_risk_scores(p_dict['vgf'], rs_id_map, phen_map, out_file=out_file, 
                                                        split_by_chrom=p_dict['split_by_chrom'], adjust_for_sex=p_dict['adjust_for_sex'], 
                                                        adjust_for_covariates=p_dict['adjust_for_covariates'], 
                                                        adjust_for_pcs=p_dict['adjust_for_pcs'])        
        
        #Plot results?
    else:
        raise NotImplementedError('Results file format missing or unknown: %s'%p_dict['res_format'])
    
    




if __name__ == '__main__':
    main()

