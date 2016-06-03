#!/usr/bin/env python
"""
Coordinate genotypes and summary statistics datasets for calculating polygenic risk scores.  
Only SNPs that overlap between the two (or three) different datasets are retained.  

Usage: 
coord --gf=PLINK_LD_REF_GENOTYPE_FILE --ssf=SUM_STATS_FILE --N=SS_SAMPLE_SIZE  --out=OUT_COORD_FILE 
                            [--vgf=PLINK_VAL_GENOTYPE_FILE --vbim=VAL_PLINK_BIM_FILE  --ssf_format=SSF_FORMAT --gmdir=GENETIC_MAP_DIR
                             --gf_format=GENOTYPE_FILE_FORMAT --indiv_list=INDIV_LIST_FILE  --maf=MAF_THRES  --skip_coordination  --check_mafs] 

 - PLINK_LD_REF_GENOTYPE_FILE (and PLINK_VAL_GENOTYPE_FILE) should be a (full path) filename prefix to a standard PLINK bed file 
   (without .bed) Make sure that the fam and bim files with same names are in the same directory.  PLINK_LD_REF_GENOTYPE_FILE refers 
   LD reference genotypes, and PLINK_VAL_GENOTYPE_FILE refers to validation genotypes. It is not necessary to have LD validation 
   genotypes at this stage.

 - SUM_STATS_FILE should be a (full path) filename prefix for a text file with the GWAS summary statistics.  The STANDARD format 
   (see below) of this file should be:
     chr     pos     ref     alt     reffrq  info    rs       pval    effalt
    chr1    1020428 C       T       0.85083 0.98732 rs6687776    0.0587  -0.0100048507289348
    chr1    1020496 G       A       0.85073 0.98751 rs6678318    0.1287  -0.00826075392985992
    ..
    ..
    
 - SS_SAMPLE_SIZE: This is the approximate number of individuals used for calculating the GWAS summary statistics.

 - OUT_COORD_FILE: The output file.  This file will follow a HDF5 format and contain both LD-reference genotypes and summary 
   statistics.
 
 - VAL_PLINK_BIM_FILE (optional): This is a PLINK BIM file which can be used to filter the set of SNPs down to the set of validation SNPs.  To 
   maximize accuracy, it's best to calculate the LDpred weights for the SNPs that are used to calculate the risk scores.
 
 - SSF_FORMAT (optional): This is the format type of the summary statistics file.  Currently there are two implementations, "STANDARD", "BASIC", 
   "GIANT", "PGC", and "PGC_large".  The standard format is described above.
   
   "BASIC" format, which contains of the basic required information, is as follows:
    hg19chrc    snpid    a1    a2    bp    or    p       
    chr1    rs4951859    C    G    729679    0.97853    0.2083  
    chr1    rs142557973    T    C    731718    1.01949    0.3298  
    ..
    ..

 - GENOTYPE_FILE_FORMAT (optional): The expected genotype format.  The standard format is PLINK.  Other formats implemented is 
   DECODE format.  If the DECODE format is used, then the program assumes that the data directory is supplied in the --gf flag.
   
 - INDIV_LIST_FILE (optional): List of individuals to include in the analysis.  Currently required for the DECODE format.
 
 - MAF_THRES (optional): Minor allele frequency threshold, i.e. SNPs with MAF smaller than threshold will be excluded from analysis.

 - --skip_cordination flag assumes that the alleles have already been coordinated between LD reference, validation samples,
   and the summary statistics files.

 - GENETIC_MAP_DIR (optional): The directory of genetic a genetic map. 


 2015 (c) Bjarni J Vilhjalmsson: bjarni.vilhjalmsson@gmail.com
"""


import getopt
import sys
import os
import traceback
import h5py
import scipy as sp
from scipy import stats
from plinkio import plinkfile
import itertools as it
import gzip
import random

ambig_nts = set([('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')])
# recode_dict = {'1':'A', '2':'C', '3':'G', '4':'T'}
opp_strand_dict = {'A':'T', 'G':'C', 'T':'A', 'C':'G'}

valid_nts = set(['A','T','C','G'])




def parse_parameters():
    """
    Parse the parameters into a dict, etc.
    """
#    if len(sys.argv) == 1:
#        print __doc__
#        sys.exit(2)

    long_options_list = ['gf=', 'vgf=', 'ssf=', 'out=', 'vbim=', 'N=', 'ssf_format=','gmdir=', 'indiv_list=',  
                         'gf_format=', 'maf=','skip_coordination','check_mafs','h','help', 'debug']

    p_dict = {'gf':None, 'vgf':None, 'ssf':None, 'out':None, 'vbim':None, 'N':None, 'ssf_format':'STANDARD', 'gmdir':None, 
              'indiv_list':None, 'gf_format':'PLINK', 'maf':0.01, 'skip_coordination':False, 'debug':False, 'check_mafs':False}

    if len(sys.argv) == 1:
        print __doc__
    elif len(sys.argv) > 1:
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
            elif opt in ("--gf"): p_dict['gf'] = arg
            elif opt in ("--vgf"): p_dict['vgf'] = arg
            elif opt in ("--ssf"): p_dict['ssf'] = arg
            elif opt in ("--out"): p_dict['out'] = arg
            elif opt in ("--vbim"): p_dict['vbim'] = arg
            elif opt in ("--gmdir"): p_dict['gmdir'] = arg
            elif opt in ("--indiv_list"): p_dict['indiv_list'] = arg
            elif opt in ("--gf_format"): p_dict['gf_format'] = arg
            elif opt in ("--skip_coordination"): p_dict['skip_coordination'] = True
            elif opt in ("--check_mafs"): p_dict['check_mafs'] = True
            elif opt in ("--maf"): p_dict['maf'] = float(arg)
            elif opt in ("--debug"): p_dict['debug'] = True
            elif opt in ("--ssf_format"): p_dict['ssf_format'] = arg
            elif opt in ("--N"): p_dict['N'] = int(arg)
            else:
                print "Unkown option:", opt
                print "Use -h option for usage information."
                sys.exit(2)
    else:
        print __doc__
        sys.exit(0)
    return p_dict
    

def _get_chrom_dict_(loci, chromosomes):
    chr_dict = {}
    for chrom in chromosomes:
        chr_str = 'chrom_%d'%chrom
        chr_dict[chr_str] = {'sids':[],'snp_indices':[],'positions':[], 'nts':[]}
     
    for i, l in enumerate(loci):
        chrom = l.chromosome
        pos = l.bp_position
        chr_str = 'chrom_%d'%chrom
        chr_dict[chr_str]['sids'].append(l.name)
#         chr_dict[chr_str]['sids'].append('%d_%d'%(chrom,pos))
        chr_dict[chr_str]['snp_indices'].append(i)
        chr_dict[chr_str]['positions'].append(pos)
        chr_dict[chr_str]['nts'].append([l.allele1,l.allele2])
     
    print 'Genotype dictionary filled'
    return chr_dict


# def _get_chrom_dict_bim_(bim_file, chromosomes):
#     chr_dict = {}
#     for chrom in chromosomes:
#         chr_str = 'chrom_%d'%chrom
#         chr_dict[chr_str] = {'sids':[],'snp_indices':[],'positions':[], 'nts':[]}
#      
#     with open(bim_file) as f:
#         for i, line in enumerate(f):
#             l = line.split()
#             chrom = int(l[0])
#             chr_str = 'chrom_%d'%chrom
#             sid_str = l[1]
#             sid_list = sid_str.split(':')
#             sid = sid_list[0]
# #             if sid[:2]=='rs':
#             chr_dict[chr_str]['sids'].append(sid)
#             chr_dict[chr_str]['snp_indices'].append(i)
#             chr_dict[chr_str]['positions'].append(int(l[3]))
#             chr_dict[chr_str]['nts'].append([l[4],l[5]])
#      
#     print 'Genotype dictionary filled'
#     return chr_dict

def _parse_plink_snps_(genotype_file, snp_indices):
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    num_individs = len(samples)
    num_snps = len(snp_indices)
    raw_snps = sp.empty((num_snps,num_individs),dtype='int8')
    #If these indices are not in order then we place them in the right place while parsing SNPs.
    snp_order = sp.argsort(snp_indices)
    ordered_snp_indices = list(snp_indices[snp_order])
    ordered_snp_indices.reverse()
    print 'Iterating over file to load SNPs'
    snp_i = 0
    next_i = ordered_snp_indices.pop()
    line_i = 0
    max_i = ordered_snp_indices[0]
    while line_i <= max_i:
        if line_i < next_i:
            plinkf.next()
        elif line_i==next_i:
            line = plinkf.next()
            snp = sp.array(line, dtype='int8')
            bin_counts = line.allele_counts()
            if bin_counts[-1]>0:
                mode_v = sp.argmax(bin_counts[:2])
                snp[snp==3] = mode_v
            s_i = snp_order[snp_i]
            raw_snps[s_i]=snp
            if line_i < max_i:
                next_i = ordered_snp_indices.pop()
            snp_i+=1
        line_i +=1
    plinkf.close()
    assert snp_i==len(raw_snps), 'Failed to parse SNPs?'
    num_indivs = len(raw_snps[0])
    freqs = sp.sum(raw_snps,1, dtype='float32')/(2*float(num_indivs))
    return raw_snps, freqs


def _parse_decode_genotypes_(decode_file, sids, pns, ocg):
    ih5f = h5py.File(decode_file,'r')

    #Determine individual filter
    pns1 = ih5f['PNs'][...]
    assert len(sp.unique(pns1))==len(pns1), 'WTF?'
    pn_sort_indices = sp.argsort(pns1)
    pn_overlap_filter = sp.in1d(pns1,pns)
    pn_sort_indices = pn_sort_indices[pn_overlap_filter]
    assert sp.all(pns==pns1[pn_sort_indices]), 'Re-ordering of individuals failed?'

    ocg.create_dataset('indivs',data=pns)

    #Determine marker filter
    mns1 = ih5f['Marker-names'][...]
    assert len(sp.unique(mns1))==len(mns1), 'WTF?'
    mn_filter = sp.in1d(sids, mns1)
    mns = sids[mn_filter]

    indices = range(len(mns1))
    mn_indices_dict = dict((key, value) for (key, value) in it.izip(mns1,indices))
    
    mn_indices = []
    for mn in mns:
        mn_indices.append(mn_indices_dict[mn])
#         mn_indices = sp.array(mn_indices)
    mn_indices = sp.array(mn_indices)
    assert sp.all(mns==mns1[mn_indices]), 'Re-ordering failed?'
    

    #Pull off position from marker name
    positions = sp.array([int(mn.split(':')[1]) for mn in mns])
    
    #Sort by position
    order = sp.argsort(positions)
    positions = positions[order]
    mn_indices = mn_indices[order]
    mns = mns[order]

    #Get nucleotide
    alleles= ih5f["Alleles"][...]
    a = sp.arange(len(alleles))
    even_map = a % 2 == 0
    odd_map = a % 2 == 1
    nts = (sp.vstack([alleles[even_map],alleles[odd_map]])).T
    nts = nts[mn_indices]
    
    

    n_snps = len(mns)
    n_indivs = len(pns)
    print 'Parsing SNPs (%d x %d matrix)'%(n_snps,n_indivs)
    snps = ocg.create_dataset('raw_snps_ref',shape=(n_snps,n_indivs),dtype='single',compression='lzf')
    freqs = sp.zeros(len(mn_indices))
    snp_means = sp.zeros(len(mn_indices))
    for i, m_i in enumerate(mn_indices):
        if i%1000==0:
            print "Reached %d'th SNP"%i
        probs = ih5f["Probabilities2"][m_i,pn_sort_indices]
        pat_snp = sp.array(map(lambda x: x[0], probs),'float32')
        mat_snp = sp.array(map(lambda x: x[1], probs),'float32')
        snp = pat_snp+mat_snp
        ok_filter = (pat_snp>=0)*(mat_snp>=0)
        if not sp.all(ok_filter):
            #impute missing
            mean_gt = sp.mean(snp[ok_filter])
            snp[~ok_filter]=mean_gt

        snps[i] = snp        
        freq = mean_gt/2.0
        snp_means[i] = mean_gt
        freqs[i] = freq
    
        
    #Calculate stds
    snp_stds = sp.sqrt(2*freqs*(1-freqs)) #sp.std(raw_snps, 1)

    ocg.create_dataset('snp_stds_ref', data=snp_stds)
    ocg.create_dataset('snp_means_ref', data=snp_means)
    ocg.create_dataset('freqs_ref', data=freqs)
    ocg.create_dataset('positions', data=positions)
    ocg.create_dataset('nts', data=nts)
    ocg.create_dataset('sids',data=mns)

    return {'ss_filter':mn_filter, 'ss_order':order}





lc_2_cap_map = {'a':'A', 'c':'C', 'g':'G', 't':'T'}


def parse_sum_stats_standard(filename=None,
                    bimfile =None,
                    hdf5_file=None, 
                    n=None):
    """
    Input format:
    
     chr     pos     ref     alt     reffrq  info    rs       pval    effalt
    chr1    1020428 C       T       0.85083 0.98732 rs6687776    0.0587  -0.0100048507289348
    chr1    1020496 G       A       0.85073 0.98751 rs6678318    0.1287  -0.00826075392985992

    """

    if bimfile is not None:
        print 'Parsing SNP list'
        valid_sids = set()
        print 'Parsing bim file: %s'%bimfile
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])
        print len(valid_sids)

    chrom_dict = {}
#     for chrom_i in range(1, 23):
#         chrom = 'chr%d'%chrom_i
#         chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'infos':[], 'freqs':[],
#                              'betas':[], 'nts': [], 'sids': [], 'positions': []}

    

    print 'Parsing the file: %s' % filename
    with open(filename) as f:
        print f.next()
        for line in f:
            l = (line.strip()).split()
            chrom = int(l[0][3:])
            pos = int(l[1])
            sid = l[6]
            if sid in valid_sids:
                if not chrom in chrom_dict.keys():
                    chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'infos':[], 'freqs':[],
                             'betas':[], 'nts': [], 'sids': [], 'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                chrom_dict[chrom]['freqs'].append(float(l[4]))
                info = l[5]
                chrom_dict[chrom]['infos'].append(info)
                pval = float(l[7])
                chrom_dict[chrom]['ps'].append(pval)
                nt = [l[2], l[3]]
                chrom_dict[chrom]['nts'].append(nt)                
                raw_beta = float(l[8])
                chrom_dict[chrom]['log_odds'].append(raw_beta)
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                
                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
            

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not 'sum_stats' in hdf5_file.keys(), 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict.keys():
        print '%d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['infos'],
                 chrom_dict[chrom]['freqs'], chrom_dict[chrom]['ps'])
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        infos = []
        freqs=[]
        prev_pos = -1
        for pos, sid, nt, beta, lo, info, frq, p in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            infos.append(info)
            freqs.append(frq)
        print 'Still %d SNPs on chromosome %d'%(len(ps),chrom)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('freqs', data=freqs)
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps +=len(log_odds)
        g.create_dataset('infos', data=infos)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print '%d SNPs parsed from summary statistics file.'% num_snps
    


def parse_sum_stats_giant(filename=None,
                          bimfile =None,
                          hdf5_file=None,
                          debug=False):
    """
    Input format:
    
    MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  b       SE      p       N
    MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
    rs10 a c 0.0333 0.8826 78380
    rs1000000 a g 0.3667 0.1858 133822

    """
    lc_CAPs_dict = {'a':'A', 'c':'C', 'g':'G', 't':'T'}

    snps_pos_map = {}
    assert bimfile is not None, 'BIM file is required'
    print 'Parsing SNP list'
    valid_sids = set()
    print 'Parsing bim file: %s'%bimfile
    with open(bimfile) as f:
        for line in f:
            l = line.split()
#             sid_str = l[1]
#             sid_list = sid_str.split(':')
#             sid = sid_list[0]
#             if sid[:2]=='rs':
#             valid_sids.add(sid[0])
            valid_sids.add(l[1])
            #valid_sids.add("%d_%d"%(l[0],l[3]))
            snps_pos_map[l[1]]={'pos':int(l[3]), 'chrom':int(l[0])}
    print len(valid_sids)

    chrom_dict = {}

    print 'Parsing the file: %s' % filename
    with open(filename) as f:
        print f.next()
        for line in f:
            l = (line.strip()).split()
            sid = l[0]
            if debug and random.random()>0.1:
                continue
            if sid in valid_sids:                        
                pos = snps_pos_map[sid]['pos']
                chrom = snps_pos_map[sid]['chrom']
                if not chrom in chrom_dict.keys():
                    chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'freqs':[],
                             'betas':[], 'nts': [], 'sids': [], 'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
#                 chrom_dict[chrom]['sids'].append('%d_%d'%(chrom,pos))
                chrom_dict[chrom]['positions'].append(pos)
                if l[3]=='.' or l[3]=='NA':
                    freq = -1
                else:
                    freq = float(l[3])
                chrom_dict[chrom]['freqs'].append(freq)
                raw_beta = float(l[4])
                pval = float(l[6])
                n = float(l[7])
                chrom_dict[chrom]['ps'].append(pval)
                #nt = [lc_CAPs_dict[l[1]], lc_CAPs_dict[l[2]]]
                nt = [l[1], l[2]]
                chrom_dict[chrom]['nts'].append(nt)                
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
                chrom_dict[chrom]['log_odds'].append(beta / sp.sqrt(n))

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not 'sum_stats' in hdf5_file.keys(), 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict.keys():
        print '%d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'],
                 chrom_dict[chrom]['freqs'], chrom_dict[chrom]['ps'])
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        freqs=[]
        prev_pos = -1
        for pos, sid, nt, beta, lo, frq, p in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            freqs.append(frq)
        print 'Still %d SNPs on chromosome %d'%(len(ps),chrom)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('freqs', data=freqs)
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps +=len(log_odds)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print '%d SNPs parsed from summary statistics file.'% num_snps
    
    
    
def parse_sum_stats_decode(filename=None,
                          hdf5_file=None,
                          debug_filter=1):
    """
    Input format:
    
    MarkerName      rs      Allele1 Allele2 Freq1   Weight  Zscore  pval    SE      Beta    chr     pos     ref     alt     mn      iref    a1      a2      frq1    wt1

    """
    lc_CAPs_dict = {'a':'A', 'c':'C', 'g':'G', 't':'T'}

    chrom_dict = {}

    print 'Parsing the file: %s' % filename
    with open(filename) as f:
        print f.next()
        for line in f:
#             if random.random()>debug_filter:
#                 continue
            l = (line.strip()).split()
            sid = l[14]
            pos = l[11]
            chrom = int(l[10][3:])
            if not chrom in chrom_dict.keys():
                chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'freqs':[],
                         'betas':[], 'nts': [], 'sids': [], 'positions': []}
            chrom_dict[chrom]['sids'].append(sid)
            chrom_dict[chrom]['positions'].append(pos)
            freq = float(l[18])
            chrom_dict[chrom]['freqs'].append(freq)
            pval = float(l[7])
            raw_beta = float(l[19])
            n = int(float(l[5]))
            chrom_dict[chrom]['ps'].append(pval)
            nt = [l[16], l[17]]
            chrom_dict[chrom]['nts'].append(nt)                
            chrom_dict[chrom]['log_odds'].append(raw_beta)
            beta = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
            chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not 'sum_stats' in hdf5_file.keys(), 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict.keys():
        print '%d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'],
                 chrom_dict[chrom]['freqs'], chrom_dict[chrom]['ps'])
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        freqs=[]
        prev_pos = -1
        for pos, sid, nt, beta, lo, frq, p in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            freqs.append(frq)
        print 'Still %d SNPs on chromosome %d'%(len(ps),chrom)
        g = ssg.create_group('chrom_%d' % chrom)
        ps = sp.array(ps)
        assert not sp.any(sp.isnan(ps)),'Some p-values are nan'
        g.create_dataset('betas', data=betas)
        g.create_dataset('ps', data=ps)
        g.create_dataset('freqs', data=freqs)
        betas = sp.array(betas)
        assert not sp.any(sp.isnan(betas)),'Some betas are nan'
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps +=len(log_odds)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print '%d SNPs parsed from summary statistics file.'% num_snps
       
    
    
    
def parse_sum_stats_pgc(filename=None,
                        bimfile =None,
                        hdf5_file=None, 
                        n=None):
    """
    Input format:

    CHR    SNP    BP    A1    A2    FRQ_A_30232    FRQ_U_40578    INFO    OR    SE    P    ngt    Direction    HetISqt    HetChiSq    HetDf    HetPVa
    ...
    
    """
    
    if bimfile is not None:
        print 'Parsing SNP list'
        valid_sids = set()
        print 'Parsing bim file: %s'%bimfile
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])
        print len(valid_sids)

    chrom_dict = {}

    

    print 'Parsing the file: %s' % filename
    denom = float(30542+40629)
    a_scalar = 30542/denom
    u_scalar = 40629/denom
    with open(filename) as f:
        print f.next()
        for line in f:
            l = (line.strip()).split()
            chrom = int(l[0])
            pos = int(l[2])
            sid = l[1]
            if sid in valid_sids:
                if not chrom in chrom_dict.keys():
                    chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'infos':[], 'freqs':[],
                             'betas':[], 'nts': [], 'sids': [], 'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                freq_a = float(l[5]) 
                freq_u = float(l[6])
                freq =  freq_a*a_scalar+freq_u*u_scalar
                chrom_dict[chrom]['freqs'].append(freq)
                info = l[7]
                chrom_dict[chrom]['infos'].append(info)
                pval = float(l[10])
                chrom_dict[chrom]['ps'].append(pval)
                nt = [l[3], l[4]]
                chrom_dict[chrom]['nts'].append(nt)                
                raw_beta = sp.log(float(l[8]))
                chrom_dict[chrom]['log_odds'].append(raw_beta)
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                
                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
            

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not 'sum_stats' in hdf5_file.keys(), 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict.keys():
        print 'Parsed summary stats for %d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['infos'],
                 chrom_dict[chrom]['freqs'], chrom_dict[chrom]['ps'])
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        infos = []
        freqs=[]
        prev_pos = -1
        for pos, sid, nt, beta, lo, info, frq, p in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            infos.append(info)
            freqs.append(frq)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('freqs', data=freqs)
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps +=len(log_odds)
        g.create_dataset('infos', data=infos)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print 'In all, %d SNPs parsed from summary statistics file.'%num_snps
    
    
    
def parse_sum_stats_pgc_small(filename=None,
                              bimfile =None,
                              hdf5_file=None, 
                              n=None):
    """
    Input format:

    hg19chrc        snpid   a1      a2      bp      info    or      se      p       ngt
    chr1    rs4951859       C       G       729679  0.631   0.97853 0.0173  0.2083  0
    chr1    rs142557973     T       C       731718  0.665   1.01949 0.0198  0.3298  0
    ...
    
    """
    
    if bimfile is not None:
        print 'Parsing SNP list'
        valid_sids = set()
        print 'Parsing bim file: %s'%bimfile
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])
        print len(valid_sids)

    chrom_dict = {}

    

    print 'Parsing the file: %s' % filename
    with open(filename) as f:
        print f.next()
        for line in f:
            l = (line.strip()).split()
            chrom_str = l[0]
            chrom = chrom_str[3:]
            if chrom.isdigit():
                chrom = int(chrom)
                pos = int(l[4])
                sid = l[1]
                if sid in valid_sids:
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'infos':[],
                                             'betas':[], 'nts': [], 'sids': [], 
                                             'positions': []}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    info = l[5]
                    chrom_dict[chrom]['infos'].append(info)
                    pval = float(l[8])
                    chrom_dict[chrom]['ps'].append(pval)
                    nt = [l[2], l[3]]
                    chrom_dict[chrom]['nts'].append(nt)                
                    raw_beta = sp.log(float(l[6]))
                    chrom_dict[chrom]['log_odds'].append(raw_beta)
                    beta = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    
                    chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
            

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not 'sum_stats' in hdf5_file.keys(), 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict.keys():
        print 'Parsed summary stats for %d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['infos'],
                 chrom_dict[chrom]['ps'])
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        infos = []
        prev_pos = -1
        for pos, sid, nt, beta, lo, info, p in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            infos.append(info)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps +=len(log_odds)
        g.create_dataset('infos', data=infos)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print 'In all, %d SNPs parsed from summary statistics file.'%num_snps
        


def parse_sum_stats_basic(filename=None,
                              bimfile =None,
                              hdf5_file=None, 
                              n=None):
    """
    Input format:

    hg19chrc    snpid    a1    a2    bp    or    p       
    chr1    rs4951859    C    G    729679    0.97853    0.2083  
    chr1    rs142557973    T    C    731718    1.01949    0.3298  
    ...
    
    """
    
    if bimfile is not None:
        print 'Parsing SNP list'
        valid_sids = set()
        print 'Parsing bim file: %s'%bimfile
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])
        print len(valid_sids)
    chrom_dict = {}


    print 'Parsing the file: %s' % filename
    with open(filename) as f:
        print f.next()
        for line in f:
            l = (line.strip()).split()
            chrom_str = l[0]
            chrom = chrom_str[3:]
            if chrom.isdigit():
                chrom = int(chrom)
                pos = int(l[4])
                sid = l[1]
                if sid in valid_sids:
                    if not chrom in chrom_dict.keys():
                        chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'infos':[],
                                             'betas':[], 'nts': [], 'sids': [], 
                                             'positions': []}
                    chrom_dict[chrom]['sids'].append(sid)
                    chrom_dict[chrom]['positions'].append(pos)
                    pval = float(l[6])
                    chrom_dict[chrom]['ps'].append(pval)
                    nt = [l[2], l[3]]
                    chrom_dict[chrom]['nts'].append(nt)                
                    raw_beta = sp.log(float(l[5]))
                    chrom_dict[chrom]['log_odds'].append(raw_beta)
                    beta = sp.sign(raw_beta) * stats.norm.ppf(pval/2.0)
                    chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
     
            

    print 'SS file loaded, now sorting and storing in HDF5 file.'
    assert not 'sum_stats' in hdf5_file.keys(), 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict.keys():
        print 'Parsed summary stats for %d SNPs on chromosome %d'%(len(chrom_dict[chrom]['positions']),chrom)
        sl = zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['ps'])
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        prev_pos = -1
        for pos, sid, nt, beta, lo, p in sl:
            if pos == prev_pos:
                print 'duplicated position %d' % pos
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
        g = ssg.create_group('chrom_%d' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps +=len(log_odds)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print 'In all, %d SNPs parsed from summary statistics file.'%num_snps
            


def coordinate_decode_genot_ss(genotype_file=None,
                                hdf5_file=None,
                                genetic_map_dir=None, indiv_file=None):
    """
    Assumes deCODE genotype files.  Imputes missing genotypes.
    """
    pns = []
    with open(indiv_file) as f:
        for line in f:
            pns.append(line.strip())
    print 'Parsed IDs for %d individuals.'%len(pns)
    pns = sp.array(pns)
    
    #Figure out overlap in individuals, and order them
#     for chrom in range(1,23):
#         fn = "%s/chr%d.hdf5"%(genotype_file,chrom)
#         ih5f = h5py.File(fn,"r")
#         pns1 = ih5f['PNs'][...]
#         pns = sp.intersect1d(pns, pns1)
#         ih5f.close()
#     print 'Found %d overlapping pns'%len(pns)
#     pns = sp.sort(pns)

    hdf5_file.create_dataset('fids', data=pns)
    hdf5_file.create_dataset('iids', data=pns)
    ssf = hdf5_file['sum_stats']
    cord_data_g = hdf5_file.create_group('cord_data')

    #Figure out chromosomes and positions by looking at SNPs.  
    chromosomes = ssf.keys()
    num_common_snps = 0
    for chr_str in chromosomes:
        chrom = int(chr_str.split('_')[1])
        print 'Working on chromsome: %s'%chr_str
        try:
            ssg = ssf['chrom_%d' % chrom]
        except Exception, err_str:
            print err_str
            print 'Did not find chromsome in SS dataset.'
            print 'Continuing.'
            continue
        ss_sids = ssg['sids'][...]
        ss_sid_set = set(ss_sids)
        assert len(ss_sid_set) == len(ss_sids), 'The summary statistics contain some duplicates?'
        
        ofg = cord_data_g.create_group('chrom_%d' % chrom)
        decode_file='%s/chr%d.hdf5'%(genotype_file,chrom)
        ret_d = _parse_decode_genotypes_(decode_file, ss_sids, pns, ofg)        
        
        ss_filter = ret_d['ss_filter']
        ss_order = ret_d['ss_order']
        betas = ssg['betas'][...]
        betas = (betas[ss_filter])[ss_order]
        log_odds = ssg['log_odds'][...]
        log_odds = (log_odds[ss_filter])[ss_order]
        ps = ssg['ps'][...]
        ps = (ps[ss_filter])[ss_order]
        
        assert not sp.any(sp.isnan(betas)), 'Some of the effect estimates are missing (parsed as NaN).'
        assert not sp.any(sp.isinf(betas)), 'Some of the effect estimates are missing (parsed as Inf).'
#        ofg.create_dataset('genetic_map', data=genetic_map)
        ofg.create_dataset('ps', data=ps)
        ofg.create_dataset('betas', data=betas)
        ofg.create_dataset('log_odds', data=log_odds)
        
        num_common_snps += len(betas)
    print 'There were %d SNPs in common' % num_common_snps
    print 'Done coordinating genotypes and summary statistics datasets.'


def coordinate_genot_ss(genotype_file=None,
                        hdf5_file=None,
                        genetic_map_dir=None,
                        check_mafs=False,
                        min_maf =0.01):
    """
    Assumes plink BED files.  Imputes missing genotypes.
    """
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    num_individs = len(samples)
#        num_individs = len(gf['chrom_1']['snps'][:, 0])
#     Y = sp.array(gf['indivs']['phenotype'][...] == 'Case', dtype='int8')
    Y = [s.phenotype for s in samples]
    fids = [s.fid for s in samples]
    iids = [s.iid for s in samples]
    unique_phens = sp.unique(Y)
    if len(unique_phens)==1:
        print 'Unable to find phenotype values.'
        has_phenotype=False
    elif len(unique_phens)==2:
        cc_bins = sp.bincount(Y)
        assert len(cc_bins)==2, 'Problems with loading phenotype'
        print 'Loaded %d controls and %d cases'%(cc_bins[0], cc_bins[1])
        has_phenotype=True
    else:
        print 'Found quantitative phenotype values'
        has_phenotype=True
    risk_scores = sp.zeros(num_individs)
    rb_risk_scores = sp.zeros(num_individs)
    num_common_snps = 0
    corr_list = []
    rb_corr_list = []

    if has_phenotype:
        hdf5_file.create_dataset('y', data=Y)
    
    hdf5_file.create_dataset('fids', data=fids)
    hdf5_file.create_dataset('iids', data=iids)
    ssf = hdf5_file['sum_stats']
    cord_data_g = hdf5_file.create_group('cord_data')

    #Figure out chromosomes and positions by looking at SNPs.  
    loci = plinkf.get_loci()
    plinkf.close()
    gf_chromosomes = [l.chromosome for l in loci] 

    chromosomes = sp.unique(gf_chromosomes)
    chromosomes.sort()
    chr_dict = _get_chrom_dict_(loci, chromosomes)
    
    tot_num_non_matching_nts = 0
    for chrom in chromosomes:
        chr_str = 'chrom_%d'%chrom
        print 'Working on chromsome: %s'%chr_str
        
        chrom_d = chr_dict[chr_str]
        try:
            ssg = ssf['chrom_%d' % chrom]
        except Exception, err_str:
            print err_str
            print 'Did not find chromsome in SS dataset.'
            print 'Continuing.'
            continue

        g_sids = chrom_d['sids']
        g_sid_set = set(g_sids)
        assert len(g_sid_set) == len(g_sids), 'Some duplicates?'
        ss_sids = ssg['sids'][...]
        ss_sid_set = set(ss_sids)
        assert len(ss_sid_set) == len(ss_sids), 'Some duplicates?'

        #Figure out filters:
        g_filter = sp.in1d(g_sids,ss_sids)
        ss_filter = sp.in1d(ss_sids,g_sids)

        #Order by SNP IDs
        g_order = sp.argsort(g_sids)
        ss_order = sp.argsort(ss_sids)

        g_indices = []
        for g_i in g_order:
            if g_filter[g_i]:
                g_indices.append(g_i)

        ss_indices = []
        for ss_i in ss_order:
            if ss_filter[ss_i]:
                ss_indices.append(ss_i)

        g_nts = chrom_d['nts']
        snp_indices = chrom_d['snp_indices']
        ss_nts = ssg['nts'][...]
        betas = ssg['betas'][...]
        log_odds = ssg['log_odds'][...]
        assert not sp.any(sp.isnan(betas)), 'WTF?'
        assert not sp.any(sp.isinf(betas)), 'WTF?'

        num_non_matching_nts = 0
        num_ambig_nts = 0
        ok_nts = []
        print 'Found %d SNPs present in both datasets'%(len(g_indices))

        if 'freqs' in ssg.keys():
            ss_freqs = ssg['freqs'][...]
            ss_freqs_list=[]
        
        ok_indices = {'g':[], 'ss':[]}
        for g_i, ss_i in it.izip(g_indices, ss_indices):
            
            #Is the nucleotide ambiguous?
            #g_nt = [recode_dict[g_nts[g_i][0]],recode_dict[g_nts[g_i][1]]
            g_nt = [g_nts[g_i][0],g_nts[g_i][1]]
            if tuple(g_nt) in ambig_nts:
                num_ambig_nts +=1
                tot_num_non_matching_nts += 1
                continue
            
            #First check if nucleotide is sane?
            if (not g_nt[0] in valid_nts) or (not g_nt[1] in valid_nts):
                num_non_matching_nts += 1
                tot_num_non_matching_nts += 1                
                continue

            ss_nt = ss_nts[ss_i]
            #Are the nucleotides the same?
            flip_nts = False
            os_g_nt = sp.array([opp_strand_dict[g_nt[0]], opp_strand_dict[g_nt[1]]])
            if not (sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)):
                # Opposite strand nucleotides
                flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
                if flip_nts:
                    betas[ss_i] = -betas[ss_i]
                    log_odds[ss_i] = -log_odds[ss_i]
                    if 'freqs' in ssg.keys():
                        ss_freqs[ss_i] = 1-ss_freqs[ss_i]
                else:
#                     print "Nucleotides don't match after all?: g_sid=%s, ss_sid=%s, g_i=%d, ss_i=%d, g_nt=%s, ss_nt=%s" % \
#                         (g_sids[g_i], ss_sids[ss_i], g_i, ss_i, str(g_nt), str(ss_nt))
                    num_non_matching_nts += 1
                    tot_num_non_matching_nts += 1
                        
                    continue

            # everything seems ok.
            ok_indices['g'].append(g_i)
            ok_indices['ss'].append(ss_i)
            ok_nts.append(g_nt)

        print '%d SNPs were excluded due to ambiguous nucleotides.' % num_ambig_nts
        print '%d SNPs were excluded due to non-matching nucleotides.' % num_non_matching_nts

        #Resorting by position
        positions = sp.array(chrom_d['positions'])[ok_indices['g']]
        order = sp.argsort(positions)
        ok_indices['g'] = list(sp.array(ok_indices['g'])[order])
        ok_indices['ss'] = list(sp.array(ok_indices['ss'])[order])
        positions = positions[order]
        
        #Parse SNPs
        snp_indices = sp.array(chrom_d['snp_indices'])
        snp_indices = snp_indices[ok_indices['g']] #Pinpoint where the SNPs are in the file.
        raw_snps, freqs = _parse_plink_snps_(genotype_file, snp_indices)
        print 'raw_snps.shape=', raw_snps.shape

        snp_stds = sp.sqrt(2*freqs*(1-freqs)) #sp.std(raw_snps, 1) 
        snp_means = freqs*2 #sp.mean(raw_snps, 1)

        betas = betas[ok_indices['ss']]
        log_odds = log_odds[ok_indices['ss']]
        ps = ssg['ps'][...][ok_indices['ss']]
        nts = sp.array(ok_nts)[order]
        sids = ssg['sids'][...][ok_indices['ss']]

        #Check SNP frequencies..
        if check_mafs and 'freqs' in ssg.keys():
            ss_freqs = ss_freqs[ok_indices['ss']]
            freq_discrepancy_snp = sp.absolute(ss_freqs-(1-freqs))>0.15
            if sp.any(freq_discrepancy_snp):
                print 'Warning: %d SNPs appear to have high frequency discrepancy between summary statistics and validation sample'%sp.sum(freq_discrepancy_snp)
                print freqs[freq_discrepancy_snp]
                print ss_freqs[freq_discrepancy_snp]
                
                #Filter freq_discrepancy_snps
                ok_freq_snps = sp.negative(freq_discrepancy_snp)
                raw_snps = raw_snps[ok_freq_snps]
                snp_stds = snp_stds[ok_freq_snps]
                snp_means = snp_means[ok_freq_snps]
                freqs = freqs[ok_freq_snps]
                ps = ps[ok_freq_snps]
                positions = positions[ok_freq_snps]
                nts = nts[ok_freq_snps]
                sids = sids[ok_freq_snps]
                betas = betas[ok_freq_snps]
                log_odds = log_odds[ok_freq_snps]
                 
        
        #Filter minor allele frequency SNPs.
        maf_filter = (freqs>min_maf)*(freqs<(1-min_maf))
        maf_filter_sum = sp.sum(maf_filter)
        n_snps = len(maf_filter)
        assert maf_filter_sum<=n_snps, "WTF?"
        if sp.sum(maf_filter)<n_snps:
            raw_snps = raw_snps[maf_filter]
            snp_stds = snp_stds[maf_filter]
            snp_means = snp_means[maf_filter]
            freqs = freqs[maf_filter]
            ps = ps[maf_filter]
            positions = positions[maf_filter]
            nts = nts[maf_filter]
            sids = sids[maf_filter]
            betas = betas[maf_filter]
            log_odds = log_odds[maf_filter]
            
            
            print '%d SNPs with MAF < %0.3f were filtered'%(n_snps-maf_filter_sum,min_maf)

        print '%d SNPs were retained on chromosome %d.' % (maf_filter_sum, chrom)
        
        rb_prs = sp.dot(sp.transpose(raw_snps), log_odds)
        if has_phenotype:
            print 'Normalizing SNPs'
            snp_means.shape = (len(raw_snps),1)
            snp_stds.shape = (len(raw_snps),1)
            snps = (raw_snps - snp_means) / snp_stds
            assert snps.shape==raw_snps.shape, 'Aha!'
            snp_stds = snp_stds.flatten()
            snp_means = snp_means.flatten()
            prs = sp.dot(sp.transpose(snps), betas)
            corr = sp.corrcoef(Y, prs)[0, 1]
            corr_list.append(corr)
            print 'PRS correlation for chromosome %d was %0.4f' % (chrom, corr)
            rb_corr = sp.corrcoef(Y, rb_prs)[0, 1]
            rb_corr_list.append(rb_corr)
            print 'Raw effect sizes PRS correlation for chromosome %d was %0.4f' % (chrom, rb_corr)
        
        sid_set = set(sids)
        if genetic_map_dir is not None:
            genetic_map = [] 
            with gzip.open(genetic_map_dir+'chr%d.interpolated_genetic_map.gz'%chrom) as f:
                for line in f:
                    l = line.split()
                    if l[0] in sid_set:
                        genetic_map.append(l[0])
        
        print 'Now storing coordinated data to HDF5 file.'
        ofg = cord_data_g.create_group('chrom_%d' % chrom)
        ofg.create_dataset('raw_snps_ref', data=raw_snps, compression='lzf')
        ofg.create_dataset('snp_stds_ref', data=snp_stds)
        ofg.create_dataset('snp_means_ref', data=snp_means)
        ofg.create_dataset('freqs_ref', data=freqs)
        ofg.create_dataset('ps', data=ps)
        ofg.create_dataset('positions', data=positions)
        ofg.create_dataset('nts', data=nts)
        ofg.create_dataset('sids', data=sids)
        if genetic_map_dir is not None:
            ofg.create_dataset('genetic_map', data=genetic_map)
#         print 'Sum of squared effect sizes:', sp.sum(betas ** 2)
#         print 'Sum of squared log odds:', sp.sum(log_odds ** 2)
        ofg.create_dataset('betas', data=betas)
        ofg.create_dataset('log_odds', data=log_odds)
        ofg.create_dataset('log_odds_prs', data=rb_prs)
        if has_phenotype:
            risk_scores += prs
        rb_risk_scores += rb_prs
        num_common_snps += len(betas)

    if has_phenotype:
        # Now calculate the prediction r^2
        corr = sp.corrcoef(Y, risk_scores)[0, 1]
        rb_corr = sp.corrcoef(Y, rb_risk_scores)[0, 1]
        print 'PRS R2 prediction accuracy for the whole genome was %0.4f (corr=%0.4f)' % (corr ** 2,corr)
        print 'Log-odds (effects) PRS R2 prediction accuracy for the whole genome was %0.4f (corr=%0.4f)' % (rb_corr ** 2, rb_corr)
    print 'There were %d SNPs in common' % num_common_snps
    print 'In all, %d SNPs were excluded due to nucleotide issues.' % tot_num_non_matching_nts
    print 'Done coordinating genotypes and summary statistics datasets.'





def coordinate_genotypes_ss_w_ld_ref(genotype_file = None,
                                    reference_genotype_file = None,
                                    hdf5_file = None,
                                    genetic_map_dir=None,
                                    check_mafs=False,
                                    min_maf=0.01):
#   recode_dict = {1:'A', 2:'T', 3:'C', 4:'G'} #1K genomes recoding..
    print 'Coordinating things w genotype file: %s \nref. genot. file: %s'%(genotype_file, reference_genotype_file) 
    plinkf = plinkfile.PlinkFile(genotype_file)
    
    #Loads only the individuals... (I think?)
    samples = plinkf.get_samples()
    num_individs = len(samples)
    Y = [s.phenotype for s in samples]
    fids = [s.fid for s in samples]
    iids = [s.iid for s in samples]
    
    unique_phens = sp.unique(Y)
    if len(unique_phens)==1:
        print 'Unable to find phenotype values.'
        has_phenotype=False
    elif len(unique_phens)==2:
        cc_bins = sp.bincount(Y)
        assert len(cc_bins)==2, 'Problems with loading phenotype'
        print 'Loaded %d controls and %d cases'%(cc_bins[0], cc_bins[1])
        has_phenotype=True
    else:
        print 'Found quantitative phenotype values'
        has_phenotype=True

    #Figure out chromosomes and positions.  
    print 'Parsing validation genotype bim file'
    loci = plinkf.get_loci()
    plinkf.close()
    gf_chromosomes = [l.chromosome for l in loci] 

    chromosomes = sp.unique(gf_chromosomes)
    chromosomes.sort()
    
    chr_dict = _get_chrom_dict_(loci, chromosomes)

    print 'Parsing LD reference genotype bim file'
    plinkf_ref = plinkfile.PlinkFile(reference_genotype_file)
    loci_ref = plinkf_ref.get_loci()
    plinkf_ref.close()
    
    chr_dict_ref = _get_chrom_dict_(loci_ref, chromosomes)
#     chr_dict_ref = _get_chrom_dict_bim_(reference_genotype_file+'.bim', chromosomes)
    
    #Open HDF5 file and prepare out data
    assert not 'iids' in hdf5_file.keys(), 'Something is wrong with the HDF5 file?'
    if has_phenotype:
        hdf5_file.create_dataset('y', data=Y)
    
    hdf5_file.create_dataset('fids', data=fids)
    hdf5_file.create_dataset('iids', data=iids)
    ssf = hdf5_file['sum_stats']
    cord_data_g = hdf5_file.create_group('cord_data')

    maf_adj_risk_scores = sp.zeros(num_individs)
    num_common_snps = 0
    #corr_list = []
    
    tot_g_ss_nt_concord_count = 0
    tot_rg_ss_nt_concord_count = 0
    tot_g_rg_nt_concord_count = 0    
    tot_num_non_matching_nts = 0
   
    #Now iterate over chromosomes
    for chrom in chromosomes:
        ok_indices = {'g':[], 'rg':[], 'ss':[]}
        
        chr_str = 'chrom_%d'%chrom
        print 'Working on chromsome: %s'%chr_str
        
        chrom_d = chr_dict[chr_str]
        chrom_d_ref = chr_dict_ref[chr_str]
        try:
            ssg = ssf['chrom_%d' % chrom]
        except Exception, err_str:
            print err_str
            print 'Did not find chromsome in SS dataset.'
            print 'Continuing.'
            continue

        ssg = ssf['chrom_%d' % chrom]
        g_sids = chrom_d['sids']
        rg_sids = chrom_d_ref['sids']
        ss_sids = ssg['sids'][...]
        print 'Found %d SNPs in validation data, %d SNPs in LD reference data, and %d SNPs in summary statistics.'%(len(g_sids), len(rg_sids), len(ss_sids))
        common_sids = sp.intersect1d(ss_sids, g_sids)
        common_sids = sp.intersect1d(common_sids, rg_sids)
        print 'Found %d SNPs on chrom %d that were common across all datasets'%(len(common_sids), chrom)

        ss_snp_map = []
        g_snp_map = []
        rg_snp_map = []
        
        ss_sid_dict = {}
        for i, sid in enumerate(ss_sids):
            ss_sid_dict[sid]=i

        g_sid_dict = {}
        for i, sid in enumerate(g_sids):
            g_sid_dict[sid]=i

        rg_sid_dict = {}
        for i, sid in enumerate(rg_sids):
            rg_sid_dict[sid]=i
            
        for sid in common_sids:
            g_snp_map.append(g_sid_dict[sid])
        
        #order by positions
        g_positions = sp.array(chrom_d['positions'])[g_snp_map]
        order = sp.argsort(g_positions)
        #order = order.tolist()
        g_snp_map = sp.array(g_snp_map)[order]
        g_snp_map = g_snp_map.tolist()
        common_sids = sp.array(common_sids)[order]

        #Get the other two maps
        for sid in common_sids:
            rg_snp_map.append(rg_sid_dict[sid])
        
        for sid in common_sids:
            ss_snp_map.append(ss_sid_dict[sid])
            
        
        g_nts = sp.array(chrom_d['nts'])
        rg_nts = sp.array(chrom_d_ref['nts'])
        rg_nts_ok = sp.array(rg_nts)[rg_snp_map]
#         rg_nts_l = []
#         for nt in rg_nts_ok:
#             rg_nts_l.append([recode_dict[nt[0]],recode_dict[nt[1]]])
#         rg_nts_ok = sp.array(rg_nts_l)
        ss_nts = ssg['nts'][...]
        betas = ssg['betas'][...]
        log_odds = ssg['log_odds'][...]

        if 'freqs' in ssg.keys():
            ss_freqs = ssg['freqs'][...]

        g_ss_nt_concord_count = sp.sum(g_nts[g_snp_map] == ss_nts[ss_snp_map])/2.0
        rg_ss_nt_concord_count = sp.sum(rg_nts_ok == ss_nts[ss_snp_map])/2.0
        g_rg_nt_concord_count = sp.sum(g_nts[g_snp_map] == rg_nts_ok)/2.0
        print 'Nucleotide concordance counts out of %d genotypes: vg-g: %d, vg-ss: %d, g-ss: %d'%(len(g_snp_map),g_rg_nt_concord_count, g_ss_nt_concord_count, rg_ss_nt_concord_count)
        tot_g_ss_nt_concord_count += g_ss_nt_concord_count
        tot_rg_ss_nt_concord_count += rg_ss_nt_concord_count
        tot_g_rg_nt_concord_count += g_rg_nt_concord_count


        num_non_matching_nts = 0
        num_ambig_nts = 0


        #Identifying which SNPs have nucleotides that are ok..
        ok_nts = []
        for g_i, rg_i, ss_i in it.izip(g_snp_map, rg_snp_map, ss_snp_map):
            
            #To make sure, is the SNP id the same?
            assert g_sids[g_i]==rg_sids[rg_i]==ss_sids[ss_i], 'Some issues with coordinating the genotypes.'
            
            g_nt = g_nts[g_i]
            rg_nt = rg_nts[rg_i]
#             rg_nt = [recode_dict[rg_nts[rg_i][0]],recode_dict[rg_nts[rg_i][1]]]
            ss_nt = ss_nts[ss_i]

            #Is the nucleotide ambiguous.
            g_nt = [g_nts[g_i][0],g_nts[g_i][1]]
            if tuple(g_nt) in ambig_nts:
                num_ambig_nts +=1
                tot_num_non_matching_nts += 1                
                continue
            
            #First check if nucleotide is sane?
            if (not g_nt[0] in valid_nts) or (not g_nt[1] in valid_nts):
                num_non_matching_nts += 1
                tot_num_non_matching_nts += 1                
                continue
            
            os_g_nt = sp.array([opp_strand_dict[g_nt[0]], opp_strand_dict[g_nt[1]]])

            flip_nts = False
            if not ((sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)) and (sp.all(g_nt == rg_nt) or sp.all(os_g_nt == rg_nt))):
                if sp.all(g_nt == rg_nt) or sp.all(os_g_nt == rg_nt):
                    flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
                    #Try flipping the SS nt
                    if flip_nts:
                        betas[ss_i] = -betas[ss_i]                        
                        log_odds[ss_i] = -log_odds[ss_i]    
                        if 'freqs' in ssg.keys():
                            ss_freqs[ss_i] = 1-ss_freqs[ss_i]
                    else:
                        print "Nucleotides don't match after all?: g_sid=%s, ss_sid=%s, g_i=%d, ss_i=%d, g_nt=%s, ss_nt=%s" % \
                            (g_sids[g_i], ss_sids[ss_i], g_i, ss_i, str(g_nt), str(ss_nt))
                        num_non_matching_nts += 1
                        tot_num_non_matching_nts += 1
                        continue

                    
                else:
                    num_non_matching_nts += 1
                    tot_num_non_matching_nts += 1
                    continue
                    # Opposite strand nucleotides
            
           
            # everything seems ok.
            ok_indices['g'].append(g_i)
            ok_indices['rg'].append(rg_i)
            ok_indices['ss'].append(ss_i)

            ok_nts.append(g_nt)
#             if flip_nts:
#                 ok_nts.append([ss_nt[1],ss_nt[0]])
#             else:
#                 ok_nts.append(ss_nt)                

                        
        #print '%d SNPs in LD references to be flipped.'%((len(ref_snp_directions)-sp.sum(ref_snp_directions))/2.0)
        print '%d SNPs had ambiguous nucleotides.' % num_ambig_nts 
        print '%d SNPs were excluded due to nucleotide issues.' % num_non_matching_nts 
        print '%d SNPs were retained on chromosome %d.' % (len(ok_indices['g']), chrom)

        #Resorting by position
        positions = sp.array(chrom_d['positions'])[ok_indices['g']]
#         order = sp.argsort(positions)
#         sorted_positions = positions[order]
#         assert sp.all(sorted_positions==positions), 'Perhaps something is wrong here?'
#         ok_indices['g'] = list(sp.array(ok_indices['g'])[order])
#         ok_indices['ss'] = list(sp.array(ok_indices['ss'])[order])

        
        #Now parse SNPs ..
        snp_indices = sp.array(chrom_d['snp_indices'])
        snp_indices = snp_indices[ok_indices['g']] #Pinpoint where the SNPs are in the file.
        raw_snps,freqs = _parse_plink_snps_(genotype_file, snp_indices)
        
        snp_indices_ref = sp.array(chrom_d_ref['snp_indices'])
        snp_indices_ref = snp_indices_ref[ok_indices['rg']] #Pinpoint where the SNPs are in the file.
        raw_ref_snps, freqs_ref = _parse_plink_snps_(reference_genotype_file, snp_indices_ref)
        
        
        snp_stds_ref = sp.sqrt(2*freqs_ref*(1-freqs_ref)) 
        snp_means_ref = freqs_ref*2

        snp_stds = sp.sqrt(2*freqs*(1-freqs)) 
        snp_means = freqs*2
        
        betas = betas[ok_indices['ss']]  # * sp.sqrt(freqs * (1 - freqs))
        log_odds = log_odds[ok_indices['ss']]  # * sp.sqrt(freqs * (1 - freqs))

        ps = ssg['ps'][...][ok_indices['ss']]
        nts = sp.array(ok_nts)#[order]
        sids = ssg['sids'][...][ok_indices['ss']]

        #For debugging...
#         g_sids = sp.array(chrom_d['sids'])[ok_indices['g']]
#         rg_sids = sp.array(chrom_d_ref['sids'])[ok_indices['rg']]
#         ss_sids = ssg['sids'][...][ok_indices['ss']]
#         assert sp.all(g_sids==rg_sids) and sp.all(rg_sids==ss_sids), 'WTF?'
        
        #Check SNP frequencies..
        if check_mafs and 'freqs' in ssg.keys():
            ss_freqs = ss_freqs[ok_indices['ss']]
            freq_discrepancy_snp = sp.absolute(ss_freqs-(1-freqs))>0.15
            if sp.any(freq_discrepancy_snp):
                print 'Warning: %d SNPs were filtered due to high allele frequency discrepancy between summary statistics and validation sample'%sp.sum(freq_discrepancy_snp)
#                 print freqs[freq_discrepancy_snp]
#                 print ss_freqs[freq_discrepancy_snp]
                 
                #Filter freq_discrepancy_snps
                ok_freq_snps = sp.negative(freq_discrepancy_snp)
                raw_snps = raw_snps[ok_freq_snps]
                snp_stds = snp_stds[ok_freq_snps]
                snp_means = snp_means[ok_freq_snps]
                raw_ref_snps = raw_ref_snps[ok_freq_snps]
                snp_stds_ref = snp_stds_ref[ok_freq_snps]
                snp_means_ref = snp_means_ref[ok_freq_snps]
                freqs = freqs[ok_freq_snps]
                freqs_ref = freqs_ref[ok_freq_snps]
                ps = ps[ok_freq_snps]
                positions = positions[ok_freq_snps]
                nts = nts[ok_freq_snps]
                sids = sids[ok_freq_snps]
                betas = betas[ok_freq_snps]
                log_odds = log_odds[ok_freq_snps]
                #For debugging...
#         if sp.any(freq_discrepancy_snp):
#             g_sids = g_sids[ok_freq_snps]
#             rg_sids = rg_sids[ok_freq_snps]
#             ss_sids = ss_sids[ok_freq_snps]
#         assert sp.all(g_sids==rg_sids) and sp.all(rg_sids==ss_sids), 'WTF?'

        
        
        #Filter minor allele frequency SNPs.
        maf_filter = (freqs>min_maf)*(freqs<(1-min_maf))
        maf_filter_sum = sp.sum(maf_filter)
        n_snps = len(maf_filter)
        assert maf_filter_sum<=n_snps, "WTF?"
        if sp.sum(maf_filter)<n_snps:
            raw_snps = raw_snps[maf_filter]
            snp_stds = snp_stds[maf_filter]
            snp_means = snp_means[maf_filter]
            raw_ref_snps = raw_ref_snps[maf_filter]
            snp_stds_ref = snp_stds_ref[maf_filter]
            snp_means_ref = snp_means_ref[maf_filter]
            freqs = freqs[maf_filter]
            freqs_ref = freqs_ref[maf_filter]
            ps = ps[maf_filter]
            positions = positions[maf_filter]
            nts = nts[maf_filter]
            sids = sids[maf_filter]
            betas = betas[maf_filter]
            log_odds = log_odds[maf_filter]
#         if sp.sum(maf_filter)<n_snps:
#             g_sids = g_sids[maf_filter]
#             rg_sids = rg_sids[maf_filter]
#             ss_sids = ss_sids[maf_filter]
#         assert sp.all(g_sids==rg_sids) and sp.all(rg_sids==ss_sids), 'WTF?'
        
        
        
        maf_adj_prs = sp.dot(log_odds, raw_snps)
        if has_phenotype:
            maf_adj_corr = sp.corrcoef(Y, maf_adj_prs)[0, 1]
            print 'Log odds, per genotype PRS correlation w phenotypes for chromosome %d was %0.4f' % (chrom, maf_adj_corr)

        genetic_map = [] 
        if genetic_map_dir is not None:
            with gzip.open(genetic_map_dir+'chr%d.interpolated_genetic_map.gz'%chrom) as f:
                for line in f:
                    l = line.split()
                    if l[0] in sid_set:
                        genetic_map.append(l[0])
        
        
        print 'Now storing coordinated data to HDF5 file.'
        ofg = cord_data_g.create_group('chrom_%d' % chrom)
        ofg.create_dataset('raw_snps_val', data=raw_snps, compression='lzf')
        ofg.create_dataset('snp_stds_val', data=snp_stds)
        ofg.create_dataset('snp_means_val', data=snp_means)
        ofg.create_dataset('freqs_val', data=freqs)
        ofg.create_dataset('raw_snps_ref', data=raw_ref_snps, compression='lzf')
        ofg.create_dataset('snp_stds_ref', data=snp_stds_ref)
        ofg.create_dataset('snp_means_ref', data=snp_means_ref)
        ofg.create_dataset('freqs_ref', data=freqs_ref)
        ofg.create_dataset('nts', data=nts)
        ofg.create_dataset('ps', data=ps)
        ofg.create_dataset('positions', data=positions)
        ofg.create_dataset('sids', data=sids)
        if genetic_map_dir is not None:
            ofg.create_dataset('genetic_map', data=genetic_map)
        ofg.create_dataset('betas', data=betas)
        ofg.create_dataset('log_odds', data=log_odds)
        ofg.create_dataset('log_odds_prs', data=maf_adj_prs)
#         print 'Sum betas', sp.sum(betas ** 2)
        #ofg.create_dataset('prs', data=prs)
        
        
        #risk_scores += prs
        maf_adj_risk_scores += maf_adj_prs
        num_common_snps += len(betas)
        
    # Now calculate the prediction r^2
    if has_phenotype:
        maf_adj_corr = sp.corrcoef(Y, maf_adj_risk_scores)[0, 1]
        #print 'PRS correlation for the whole genome was %0.4f (r^2=%0.4f)' % (corr, corr ** 2)
        print 'Log odds, per PRS correlation for the whole genome was %0.4f (r^2=%0.4f)' % (maf_adj_corr, maf_adj_corr ** 2)
    print 'Overall nucleotide concordance counts: g_rg: %d, g_ss: %d, rg_ss: %d'%(tot_g_rg_nt_concord_count, tot_g_ss_nt_concord_count, tot_rg_ss_nt_concord_count)
    print 'There were %d SNPs in common' % num_common_snps
    print 'In all, %d SNPs were excluded due to nucleotide issues.' % tot_num_non_matching_nts 
    print 'Done!'

def main():
    p_dict = parse_parameters()
    print """
    Note: For maximal accuracy all SNPs with LDpred weights should be included in the validation data set.
    If they are a subset of the validation data set, then we suggest recalculate LDpred for the overlapping SNPs. 
    You can coordinate across the three data sets by either using the same LD reference and the validation data, or using 
    the --vbim argument, and supply the validation data set PLINK formatted bim file. 
    """
    if p_dict['N'] is None:
        print 'Please specify an integer value for the sample size used to calculate the GWAS summary statistics.'
    print  'Preparing to parse summary statistics'
    if p_dict['vbim'] is not None:
        bimfile = p_dict['vbim']
    elif p_dict['vgf'] is not None:
        bimfile = p_dict['vgf']+'.bim'
    elif p_dict['gf'] is not None:
        bimfile = p_dict['gf']+'.bim'
    else:
        print 'Set of validation SNPs is missing!  Please specify either a validation PLINK genotype file, or a PLINK BIM file with the SNPs of interest.'
    if os.path.isfile(p_dict['out']):
        print 'Output file (%s) already exists!  Delete, rename it, or use a different output file.'%(p_dict['out'])
        raise Exception('Output file already exists!')
        
    h5f = h5py.File(p_dict['out'],'w')
    if p_dict['ssf_format']=='STANDARD':
        parse_sum_stats_standard(filename=p_dict['ssf'], bimfile = bimfile, hdf5_file=h5f, n=p_dict['N'])
    elif p_dict['ssf_format']=='PGC':
        parse_sum_stats_pgc_small(filename=p_dict['ssf'], bimfile = bimfile, hdf5_file=h5f, n=p_dict['N'])
    elif p_dict['ssf_format']=='PGC_large':
        parse_sum_stats_pgc(filename=p_dict['ssf'], bimfile = bimfile, hdf5_file=h5f, n=p_dict['N'])
    elif p_dict['ssf_format']=='BASIC':
        parse_sum_stats_basic(filename=p_dict['ssf'], bimfile = bimfile, hdf5_file=h5f, n=p_dict['N'])
    elif p_dict['ssf_format']=='GIANT':
        parse_sum_stats_giant(filename=p_dict['ssf'], bimfile = bimfile, hdf5_file=h5f, debug=p_dict['debug'])
    elif p_dict['ssf_format']=='DECODE':
        parse_sum_stats_decode(filename=p_dict['ssf'], hdf5_file=h5f)
    if not p_dict['vgf'] == None:
        assert p_dict['gf_format']=='PLINK', 'The validation genotype option currently only works with the PLINK format'
        coordinate_genotypes_ss_w_ld_ref(genotype_file=p_dict['vgf'], reference_genotype_file=p_dict['gf'], 
                                         genetic_map_dir=p_dict['gmdir'], check_mafs = p_dict['check_mafs'], 
                                         hdf5_file=h5f, min_maf=p_dict['maf'], skip_coordination=p_dict['skip_coordination'])
    else:
        if p_dict['gf_format']=='PLINK':
            coordinate_genot_ss(genotype_file=p_dict['gf'],  genetic_map_dir=p_dict['gmdir'], check_mafs = p_dict['check_mafs'], 
                                hdf5_file=h5f, min_maf=p_dict['maf'], skip_coordination=p_dict['skip_coordination'])
        elif p_dict['gf_format']=='DECODE':
            if not p_dict['skip_coordination']:
                raise Exception('This option requires you to skip coordination of nucleotides and some QC.  Please confirm with --skip_coordination flag.')
            coordinate_decode_genot_ss(genotype_file=p_dict['gf'],  genetic_map_dir=p_dict['gmdir'], indiv_file=p_dict['indiv_list'],  
                                       hdf5_file=h5f)
        else:
            raise Exception('Unknown genotype file format: %s'%p_dict['gf_format'])
    
    h5f.close()

if __name__ == '__main__':
    main()
        