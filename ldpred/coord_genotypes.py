#!/usr/bin/env python
"""
Coordinate genotypes and summary statistics datasets for calculating polygenic risk scores.
Only SNPs that overlap between the two (or three) different datasets are retained.

Usage:
coord --gf=PLINK_LD_REF_GENOTYPE_FILE --ssf=SUM_STATS_FILE --N=SS_SAMPLE_SIZE  --out=OUT_COORD_FILE
                            [--vgf=PLINK_VAL_GENOTYPE_FILE --vbim=VAL_PLINK_BIM_FILE  --ssf_format=SSF_FORMAT --gmdir=GENETIC_MAP_DIR
                             --gf_format=GENOTYPE_FILE_FORMAT --indiv_list=INDIV_LIST_FILE  --maf=MAF_THRES  --skip_coordination  --check_mafs
                             --chr=ChR_COLUMN_HEADER --pos=POSITION_COLUMN_HEADER --ref=REF_ALLELE_COLUMN_HEADER --alt=ALT_ALLELE_COLUMN_HEADER
                             --reffreq=REF_FRE_COLUMN_HEADER --info=INFO_COLUMN_HEADER --rs=RSID_COLUMN_HEADER --pval=PVAL_COLUMN_HEADER
                             --eff=EFFECT_SIZE_COLUMN_HEADER --ncol=SAMPLE_SIZE --beta]

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


 2019 (c) Bjarni J Vilhjalmsson: bjarni.vilhjalmsson@gmail.com
"""

import argparse
import os
import scipy as sp
import itertools as it
import gzip
import plinkfiles
import h5py
import util
import sum_stats_parsers as ssp





parser = argparse.ArgumentParser()
parser.add_argument('--gf', type=str, required=True,
                    help='LD Reference Genotype File. '
                         'Should be a (full path) filename prefix to a standard PLINK bed file (without .bed). '
                         'Make sure that the fam and bim files with same names are in the same directory. ')
parser.add_argument('--ssf', type=str, required=True,
                    help='Summary Statistic File. '
                         'Filename prefix for a text file with the GWAS summary statistics')
parser.add_argument('--N', type=int, default=None,
                    help='Number of Individuals in Summary Statistic File.  Required for most summary '
                         'statistics formats.')
parser.add_argument('--out', type=str, required=True,
                    help='Output Prefix')
parser.add_argument('--vbim', type=str, default=None,
                    help='Validation SNP file. '
                         'This is a PLINK BIM file which can be used to filter the set of SNPs down '
                         'to the set of validation SNPs. To maximize accuracy, we recommend calculating LDpred '
                         'weights for the subset of SNPs that are used to calculate the risk scores in the '
                         'target (validation) sample.')
parser.add_argument('--vgf', type=str, default=None,
                    help='Validation genotype file. '
                         'This is a PLINK BIM file which can be used to filter the set of SNPs down to the '
                         'set of validation SNPs.  To maximize accuracy, we recommend calculating LDpred '
                         'weights for the subset of SNPs that are used to calculate the risk scores in the '
                         'target (validation) sample.')
parser.add_argument('--gf-format', type=str, default="PLINK",
                    help='The expected genotype format. The standard format is PLINK. '
                         'Other formats implemented is DECODE format. '
                         'If the DECODE format is used, then the program assumes that '
                         'the data directory is supplied in the --gf flag.')
parser.add_argument('--indiv-list', type=str,
                    help='List of individuals to include in the analysis. '
                         'Currently required for the DECODE format.', default=None)
parser.add_argument('--gmdir', type=str,
                    help='The directory of genetic map.', default=None)
parser.add_argument('--skip-coordination', default=False, action='store_true',
                    help="Assumes that the alleles have already been coordinated between LD reference, "
                         "validation samples, and the summary statistics files")
parser.add_argument('--beta', default=False, action='store_true',
                    help="Assumes the summary statistics are BETA instead of OR")
parser.add_argument('--debug', default=False, action='store_true',
                    help="Activate debugging")
parser.add_argument('--check-maf', default=False, action='store_true',
                    help="Perform MAF checking")
parser.add_argument('--maf', type=float, default=0.01,
                    help='MAF filtering threshold')
parser.add_argument('--ssf-format', type=str, default="CUSTOM", 
                    help='This is the format type of the summary statistics file. '
                    'Currently there are two implementations, "STANDARD", "BASIC", "GIANT", '
                    '"PGC", and "PGC_large".  The standard format is described above.')
parser.add_argument('--rs', type=str, default="SNP",
                    help="Column header of SNP ID")
parser.add_argument('--A1', type=str, default="A1",
                    help="Column header containing the effective allele. "
                         "There isn't any standardized label for the effective allele, "
                         "therefore extra care must be taken to ensure the correct label is provided, "
                         "otherwise, the effect will be flipped.")
parser.add_argument('--A2', type=str, default="A2",
                    help="Column header containing non-effective allele.")
parser.add_argument('--pos', type=str, default="BP",
                    help="Column header containing the coordinate of SNPs.")
parser.add_argument('--info', type=str, default="INFO",
                    help="Column header containing the INFO score.")
parser.add_argument('--chr', type=str, default="CHR",
                    help="Column header containing the chromosome information.")
parser.add_argument('--reffreq', type=str, default="MAF",
                    help="Column header containing the reference MAF")
parser.add_argument('--pval', type=str, default="P",
                    help="Column header containing the P-value information")
parser.add_argument('--eff', type=str, default="OR",
                    help="Column header containing effect size information")
parser.add_argument('--ncol', type=str, default="N",
                    help="Column header containing sample size information")



def get_chrom_dict_bim(bim_file, chromosomes):
    chr_dict = {}
    for chrom in chromosomes:
        chr_str = 'chrom_%d' % chrom
        chr_dict[chr_str] = {'sids': [],
                             'snp_indices': [], 'positions': [], 'nts': []}

    with open(bim_file) as f:
        for i, line in enumerate(f):
            l = line.split()
            chrom = int(l[0])
            chr_str = 'chrom_%d' % chrom
            sid_str = l[1]
            sid_list = sid_str.split(':')
            sid = sid_list[0]
            chr_dict[chr_str]['sids'].append(sid)
            chr_dict[chr_str]['snp_indices'].append(i)
            chr_dict[chr_str]['positions'].append(int(l[3]))
            chr_dict[chr_str]['nts'].append([l[4], l[5]])

    print 'Genotype dictionary filled'
    return chr_dict




def coordinate_genot_ss(genotype_file=None,
                        hdf5_file=None,
                        genetic_map_dir=None,
                        check_mafs=False,
                        min_maf=0.01,
                        skip_coordination=False):
    """
    Assumes plink BED files.  Imputes missing genotypes.
    """
    from plinkio import plinkfile
    plinkf = plinkfile.PlinkFile(genotype_file)
    plinkf_dict = plinkfiles.get_phenotypes(plinkf)
    num_individs = plinkf_dict['num_individs']
    risk_scores = sp.zeros(num_individs)
    rb_risk_scores = sp.zeros(num_individs)
    num_common_snps = 0
    corr_list = []
    rb_corr_list = []

    if plinkf_dict['has_phenotype']:
        hdf5_file.create_dataset('y', data=plinkf_dict['phenotypes'])

    hdf5_file.create_dataset('fids', data=plinkf_dict['fids'])
    hdf5_file.create_dataset('iids', data=plinkf_dict['iids'])
    ssf = hdf5_file['sum_stats']
    cord_data_g = hdf5_file.create_group('cord_data')

    # Figure out chromosomes and positions by looking at SNPs.
    loci = plinkf.get_loci()
    plinkf.close()
    gf_chromosomes = [l.chromosome for l in loci]

    chromosomes = sp.unique(gf_chromosomes)
    chromosomes.sort()
    chr_dict = plinkfiles.get_chrom_dict(loci, chromosomes)

    tot_num_non_matching_nts = 0
    for chrom in chromosomes:
        chr_str = 'chrom_%d' % chrom
        print('Working on chromsome: %s' % chr_str)

        chrom_d = chr_dict[chr_str]
        try:
            ssg = ssf['chrom_%d' % chrom]
        except Exception as err_str:
            print(err_str)
            print('Did not find chromsome in SS dataset.')
            print('Continuing.')
            continue

        g_sids = chrom_d['sids']
        g_sid_set = set(g_sids)
        assert len(g_sid_set) == len(g_sids), 'Some SNPs appear to be duplicated?'
        ss_sids = ssg['sids'][...]
        ss_sid_set = set(ss_sids)
        assert len(ss_sid_set) == len(ss_sids), 'Some SNPs appear to be duplicated?'

        # Figure out filters:
        g_filter = sp.in1d(g_sids, ss_sids)
        ss_filter = sp.in1d(ss_sids, g_sids)

        # Order by SNP IDs
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
        assert not sp.any(sp.isnan(betas)), 'Some SNP effect estimates are NANs (not a number)'
        assert not sp.any(sp.isinf(betas)), 'Some SNP effect estimates are INFs (infinite numbers)'

        num_non_matching_nts = 0
        num_ambig_nts = 0
        ok_nts = []
        print('Found %d SNPs present in both datasets' % (len(g_indices)))

        if 'freqs' in list(ssg.keys()):
            ss_freqs = ssg['freqs'][...]

        ok_indices = {'g': [], 'ss': []}
        for g_i, ss_i in it.izip(g_indices, ss_indices):

            # Is the nucleotide ambiguous?
            g_nt = [g_nts[g_i][0], g_nts[g_i][1]]

            if not skip_coordination:
                if tuple(g_nt) in util.ambig_nts:
                    num_ambig_nts += 1
                    tot_num_non_matching_nts += 1
                    continue

                if (not g_nt[0] in util.valid_nts) or (not g_nt[1] in util.valid_nts):
                    num_non_matching_nts += 1
                    tot_num_non_matching_nts += 1
                    continue

                ss_nt = ss_nts[ss_i]

                # Are the nucleotides the same?
                flip_nts = False
                os_g_nt = sp.array(
                    [util.opp_strand_dict[g_nt[0]], util.opp_strand_dict[g_nt[1]]])
                if not (sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)):
                    # Opposite strand nucleotides
                    flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (
                        os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
                    if flip_nts:
                        betas[ss_i] = -betas[ss_i]
                        log_odds[ss_i] = -log_odds[ss_i]
                        if 'freqs' in list(ssg.keys()):
                            ss_freqs[ss_i] = 1 - ss_freqs[ss_i]
                    else:
                        num_non_matching_nts += 1
                        tot_num_non_matching_nts += 1

                        continue

            # everything seems ok.
            ok_indices['g'].append(g_i)
            ok_indices['ss'].append(ss_i)
            ok_nts.append(g_nt)

        print('%d SNPs were excluded due to ambiguous nucleotides.' % num_ambig_nts)
        print('%d SNPs were excluded due to non-matching nucleotides.' % num_non_matching_nts)

        # Resorting by position
        positions = sp.array(chrom_d['positions'])[ok_indices['g']]
        order = sp.argsort(positions)
        ok_indices['g'] = list(sp.array(ok_indices['g'])[order])
        ok_indices['ss'] = list(sp.array(ok_indices['ss'])[order])
        positions = positions[order]

        # Parse SNPs
        snp_indices = sp.array(chrom_d['snp_indices'])
        
        # Pinpoint where the SNPs are in the file.
        snp_indices = snp_indices[ok_indices['g']]
        raw_snps, freqs = plinkfiles.parse_plink_snps(
            genotype_file, snp_indices)
        print('raw_snps.shape=', raw_snps.shape)

        snp_stds = sp.sqrt(2 * freqs * (1 - freqs))  
        snp_means = freqs * 2  

        betas = betas[ok_indices['ss']]
        log_odds = log_odds[ok_indices['ss']]
        ps = ssg['ps'][...][ok_indices['ss']]
        nts = sp.array(ok_nts)[order]
        sids = ssg['sids'][...][ok_indices['ss']]

        # Check SNP frequencies..
        if check_mafs and 'freqs' in list(ssg.keys()):
            ss_freqs = ss_freqs[ok_indices['ss']]
            freq_discrepancy_snp = sp.absolute(ss_freqs - (1 - freqs)) > 0.15
            if sp.any(freq_discrepancy_snp):
                print('Warning: %d SNPs appear to have high frequency discrepancy between summary statistics and validation sample' % sp.sum(freq_discrepancy_snp))
                print(freqs[freq_discrepancy_snp])
                print(ss_freqs[freq_discrepancy_snp])

                # Filter freq_discrepancy_snps
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

        # Filter minor allele frequency SNPs.
        maf_filter = (freqs > min_maf) * (freqs < (1 - min_maf))
        maf_filter_sum = sp.sum(maf_filter)
        n_snps = len(maf_filter)
        assert maf_filter_sum <= n_snps, "Problems when filtering SNPs with low minor allele frequencies"
        if sp.sum(maf_filter) < n_snps:
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

            print('%d SNPs with MAF < %0.3f were filtered' % (n_snps - maf_filter_sum, min_maf))

        print('%d SNPs were retained on chromosome %d.' % (maf_filter_sum, chrom))

        rb_prs = sp.dot(sp.transpose(raw_snps), log_odds)
        if plinkf_dict['has_phenotype']:
            print('Normalizing SNPs')
            snp_means.shape = (len(raw_snps), 1)
            snp_stds.shape = (len(raw_snps), 1)
            snps = (raw_snps - snp_means) / snp_stds
            assert snps.shape == raw_snps.shape, 'Problems when normalizing SNPs (set to have variance 1 and 0 mean)'
            snp_stds = snp_stds.flatten()
            snp_means = snp_means.flatten()
            prs = sp.dot(sp.transpose(snps), betas)
            corr = sp.corrcoef(plinkf_dict['phenotypes'], prs)[0, 1]
            corr_list.append(corr)
            print('PRS correlation for chromosome %d was %0.4f' % (chrom, corr))
            rb_corr = sp.corrcoef(plinkf_dict['phenotypes'], rb_prs)[0, 1]
            rb_corr_list.append(rb_corr)
            print('Raw effect sizes PRS correlation for chromosome %d was %0.4f' % (chrom, rb_corr))

        sid_set = set(sids)
        if genetic_map_dir is not None:
            genetic_map = []
            with gzip.open(genetic_map_dir + 'chr%d.interpolated_genetic_map.gz' % chrom) as f:
                for line in f:
                    l = line.split()
                    if l[0] in sid_set:
                        genetic_map.append(l[0])

        print('Now storing coordinated data to HDF5 file.')
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
        ofg.create_dataset('betas', data=betas)
        ofg.create_dataset('log_odds', data=log_odds)
        ofg.create_dataset('log_odds_prs', data=rb_prs)
        if plinkf_dict['has_phenotype']:
            risk_scores += prs
        rb_risk_scores += rb_prs
        num_common_snps += len(betas)

    if plinkf_dict['has_phenotype']:
        
        # Now calculate the prediction R^2
        corr = sp.corrcoef(plinkf_dict['phenotypes'], risk_scores)[0, 1]
        rb_corr = sp.corrcoef(plinkf_dict['phenotypes'], rb_risk_scores)[0, 1]
        print('PRS R2 prediction accuracy for the whole genome was %0.4f (corr=%0.4f)' % (corr ** 2, corr))
        print('Log-odds (effects) PRS R2 prediction accuracy for the whole genome was %0.4f (corr=%0.4f)' % (rb_corr ** 2, rb_corr))
    print('There were %d SNPs in common' % num_common_snps)
    print('In all, %d SNPs were excluded due to nucleotide issues.' % tot_num_non_matching_nts)
    print('Done coordinating genotypes and summary statistics datasets.')


def coordinate_genotypes_ss_w_ld_ref(genotype_file=None,
                                     reference_genotype_file=None,
                                     hdf5_file=None,
                                     genetic_map_dir=None,
                                     check_mafs=False,
                                     min_maf=0.01,
                                     skip_coordination=False):
    print('Coordinating things w genotype file: %s \nref. genot. file: %s' % (genotype_file, reference_genotype_file))
    
    from plinkio import plinkfile
    plinkf = plinkfile.PlinkFile(genotype_file)

    # Loads only the individuals... (I think?)
    plinkf_dict = plinkfiles.get_phenotypes(plinkf)

    # Figure out chromosomes and positions.
    print('Parsing validation genotype bim file')
    loci = plinkf.get_loci()
    plinkf.close()
    gf_chromosomes = [l.chromosome for l in loci]

    chromosomes = sp.unique(gf_chromosomes)
    chromosomes.sort()

    chr_dict = plinkfiles.get_chrom_dict(loci, chromosomes)

    print('Parsing LD reference genotype bim file')
    plinkf_ref = plinkfile.PlinkFile(reference_genotype_file)
    loci_ref = plinkf_ref.get_loci()
    plinkf_ref.close()

    chr_dict_ref = plinkfiles.get_chrom_dict(loci_ref, chromosomes)

    # Open HDF5 file and prepare out data
    assert not 'iids' in list(hdf5_file.keys()), 'Something is wrong with the HDF5 file, no individuals IDs were found.'
    if plinkf_dict['has_phenotype']:
        hdf5_file.create_dataset('y', data=plinkf_dict['phenotypes'])

    hdf5_file.create_dataset('fids', data=plinkf_dict['fids'])
    hdf5_file.create_dataset('iids', data=plinkf_dict['iids'])
    ssf = hdf5_file['sum_stats']
    cord_data_g = hdf5_file.create_group('cord_data')

    maf_adj_risk_scores = sp.zeros(plinkf_dict['num_individs'])
    num_common_snps = 0
    # corr_list = []

    tot_g_ss_nt_concord_count = 0
    tot_rg_ss_nt_concord_count = 0
    tot_g_rg_nt_concord_count = 0
    tot_num_non_matching_nts = 0

    # Now iterate over chromosomes
    for chrom in chromosomes:
        ok_indices = {'g': [], 'rg': [], 'ss': []}

        chr_str = 'chrom_%d' % chrom
        print('Working on chromsome: %s' % chr_str)

        chrom_d = chr_dict[chr_str]
        chrom_d_ref = chr_dict_ref[chr_str]
        try:
            ssg = ssf['chrom_%d' % chrom]
        except Exception as err_str:
            print(err_str)
            print('Did not find chromsome in SS dataset.')
            print('Continuing.')
            continue

        ssg = ssf['chrom_%d' % chrom]
        g_sids = chrom_d['sids']
        rg_sids = chrom_d_ref['sids']
        ss_sids = ssg['sids'][...]
        print('Found %d SNPs in validation data, %d SNPs in LD reference data, and %d SNPs in summary statistics.' % (len(g_sids), len(rg_sids), len(ss_sids)))
        common_sids = sp.intersect1d(ss_sids, g_sids)
        common_sids = sp.intersect1d(common_sids, rg_sids)
        print('Found %d SNPs on chrom %d that were common across all datasets' % (len(common_sids), chrom))

        ss_snp_map = []
        g_snp_map = []
        rg_snp_map = []

        ss_sid_dict = {}
        for i, sid in enumerate(ss_sids):
            ss_sid_dict[sid] = i

        g_sid_dict = {}
        for i, sid in enumerate(g_sids):
            g_sid_dict[sid] = i

        rg_sid_dict = {}
        for i, sid in enumerate(rg_sids):
            rg_sid_dict[sid] = i

        for sid in common_sids:
            g_snp_map.append(g_sid_dict[sid])

        # order by positions
        g_positions = sp.array(chrom_d['positions'])[g_snp_map]
        order = sp.argsort(g_positions)
        # order = order.tolist()
        g_snp_map = sp.array(g_snp_map)[order]
        g_snp_map = g_snp_map.tolist()
        common_sids = sp.array(common_sids)[order]

        # Get the other two maps
        for sid in common_sids:
            rg_snp_map.append(rg_sid_dict[sid])

        for sid in common_sids:
            ss_snp_map.append(ss_sid_dict[sid])

        g_nts = sp.array(chrom_d['nts'])
        rg_nts = sp.array(chrom_d_ref['nts'])
        rg_nts_ok = sp.array(rg_nts)[rg_snp_map]
        ss_nts = ssg['nts'][...]
        betas = ssg['betas'][...]
        log_odds = ssg['log_odds'][...]

        if 'freqs' in list(ssg.keys()):
            ss_freqs = ssg['freqs'][...]

        g_ss_nt_concord_count = sp.sum(
            g_nts[g_snp_map] == ss_nts[ss_snp_map]) / 2.0
        rg_ss_nt_concord_count = sp.sum(rg_nts_ok == ss_nts[ss_snp_map]) / 2.0
        g_rg_nt_concord_count = sp.sum(g_nts[g_snp_map] == rg_nts_ok) / 2.0
        print('Nucleotide concordance counts out of %d genotypes: vg-g: %d, vg-ss: %d, g-ss: %d' % (len(g_snp_map), g_rg_nt_concord_count, g_ss_nt_concord_count, rg_ss_nt_concord_count))
        tot_g_ss_nt_concord_count += g_ss_nt_concord_count
        tot_rg_ss_nt_concord_count += rg_ss_nt_concord_count
        tot_g_rg_nt_concord_count += g_rg_nt_concord_count

        num_non_matching_nts = 0
        num_ambig_nts = 0

        # Identifying which SNPs have nucleotides that are ok..
        ok_nts = []
        for g_i, rg_i, ss_i in it.izip(g_snp_map, rg_snp_map, ss_snp_map):

            # To make sure, is the SNP id the same?
            assert g_sids[g_i] == rg_sids[rg_i] == ss_sids[ss_i], 'Some issues with coordinating the genotypes.'

            g_nt = g_nts[g_i]
            if not skip_coordination:

                rg_nt = rg_nts[rg_i]
                ss_nt = ss_nts[ss_i]

                # Is the nucleotide ambiguous.
                g_nt = [g_nts[g_i][0], g_nts[g_i][1]]
                if tuple(g_nt) in util.ambig_nts:
                    num_ambig_nts += 1
                    tot_num_non_matching_nts += 1
                    continue

                # First check if nucleotide is sane?
                if (not g_nt[0] in util.valid_nts) or (not g_nt[1] in util.valid_nts):
                    num_non_matching_nts += 1
                    tot_num_non_matching_nts += 1
                    continue

                os_g_nt = sp.array(
                    [util.opp_strand_dict[g_nt[0]], util.opp_strand_dict[g_nt[1]]])

                flip_nts = False
                if not ((sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)) and (sp.all(g_nt == rg_nt) or sp.all(os_g_nt == rg_nt))):
                    if sp.all(g_nt == rg_nt) or sp.all(os_g_nt == rg_nt):
                        flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (
                            os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
                        # Try flipping the SS nt
                        if flip_nts:
                            betas[ss_i] = -betas[ss_i]
                            log_odds[ss_i] = -log_odds[ss_i]
                            if 'freqs' in list(ssg.keys()):
                                ss_freqs[ss_i] = 1 - ss_freqs[ss_i]
                        else:
                            print("Nucleotides don't match after all?: g_sid=%s, ss_sid=%s, g_i=%d, ss_i=%d, g_nt=%s, ss_nt=%s" % \
                                (g_sids[g_i], ss_sids[ss_i], g_i,
                                 ss_i, str(g_nt), str(ss_nt)))
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

        print('%d SNPs had ambiguous nucleotides.' % num_ambig_nts)
        print('%d SNPs were excluded due to nucleotide issues.' % num_non_matching_nts)
        print('%d SNPs were retained on chromosome %d.' % (len(ok_indices['g']), chrom))

        # Resorting by position
        positions = sp.array(chrom_d['positions'])[ok_indices['g']]

        # Now parse SNPs ..
        snp_indices = sp.array(chrom_d['snp_indices'])
        # Pinpoint where the SNPs are in the file.
        snp_indices = snp_indices[ok_indices['g']]
        raw_snps, freqs = plinkfiles.parse_plink_snps(
            genotype_file, snp_indices)

        snp_indices_ref = sp.array(chrom_d_ref['snp_indices'])
        # Pinpoint where the SNPs are in the file.
        snp_indices_ref = snp_indices_ref[ok_indices['rg']]
        raw_ref_snps, freqs_ref = plinkfiles.parse_plink_snps(
            reference_genotype_file, snp_indices_ref)

        snp_stds_ref = sp.sqrt(2 * freqs_ref * (1 - freqs_ref))
        snp_means_ref = freqs_ref * 2

        snp_stds = sp.sqrt(2 * freqs * (1 - freqs))
        snp_means = freqs * 2

        betas = betas[ok_indices['ss']]  
        log_odds = log_odds[ok_indices['ss']]  

        ps = ssg['ps'][...][ok_indices['ss']]
        nts = sp.array(ok_nts)  # [order]
        sids = ssg['sids'][...][ok_indices['ss']]


        # Check SNP frequencies..
        if check_mafs and 'freqs' in list(ssg.keys()):
            ss_freqs = ss_freqs[ok_indices['ss']]
            freq_discrepancy_snp = sp.absolute(ss_freqs - (1 - freqs)) > 0.15
            if sp.any(freq_discrepancy_snp):
                print('Warning: %d SNPs were filtered due to high allele frequency discrepancy between summary statistics and validation sample' % sp.sum(freq_discrepancy_snp))

                # Filter freq_discrepancy_snps
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

        # Filter minor allele frequency SNPs.
        maf_filter = (freqs > min_maf) * (freqs < (1 - min_maf))
        maf_filter_sum = sp.sum(maf_filter)
        n_snps = len(maf_filter)
        assert maf_filter_sum <= n_snps, "Problems when filtering SNPs with low minor allele frequencies"
        if sp.sum(maf_filter) < n_snps:
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

        maf_adj_prs = sp.dot(log_odds, raw_snps)
        if plinkf_dict['has_phenotype']:
            maf_adj_corr = sp.corrcoef(
                plinkf_dict['phenotypes'], maf_adj_prs)[0, 1]
            print('Log odds, per genotype PRS correlation w phenotypes for chromosome %d was %0.4f' % (chrom, maf_adj_corr))

        genetic_map = []
        if genetic_map_dir is not None:
            with gzip.open(genetic_map_dir + 'chr%d.interpolated_genetic_map.gz' % chrom) as f:
                for line in f:
                    l = line.split()
#                     if l[0] in sid_set:
#                         genetic_map.append(l[0])

        print('Now storing coordinated data to HDF5 file.')
        ofg = cord_data_g.create_group('chrom_%d' % chrom)
        ofg.create_dataset('raw_snps_val', data=raw_snps, compression='lzf')
        ofg.create_dataset('snp_stds_val', data=snp_stds)
        ofg.create_dataset('snp_means_val', data=snp_means)
        ofg.create_dataset('freqs_val', data=freqs)
        ofg.create_dataset(
            'raw_snps_ref', data=raw_ref_snps, compression='lzf')
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

        # risk_scores += prs
        maf_adj_risk_scores += maf_adj_prs
        num_common_snps += len(betas)

    # Now calculate the prediction r^2
    if plinkf_dict['has_phenotype']:
        maf_adj_corr = sp.corrcoef(
            plinkf_dict['phenotypes'], maf_adj_risk_scores)[0, 1]
        print('Log odds, per PRS correlation for the whole genome was %0.4f (r^2=%0.4f)' % (maf_adj_corr, maf_adj_corr ** 2))
    print('Overall nucleotide concordance counts: g_rg: %d, g_ss: %d, rg_ss: %d' % (tot_g_rg_nt_concord_count, tot_g_ss_nt_concord_count, tot_rg_ss_nt_concord_count))
    print('There were %d SNPs in common' % num_common_snps)
    print('In all, %d SNPs were excluded due to nucleotide issues.' % tot_num_non_matching_nts)
    print('Done!')




def main():
    parameters = parser.parse_args()
    p_dict= vars(parameters)
 
    print """
    Note: For maximal accuracy all SNPs with LDpred weights should be included in the validation data set.
    If they are a subset of the validation data set, then we suggest recalculate LDpred for the overlapping SNPs.
    You can coordinate across the three data sets by either using the same LD reference and the validation data, or using
    the --vbim argument, and supply the validation data set PLINK formatted bim file.
    """
    
    bimfile = None
    if p_dict['N'] is None:
        print 'Please specify an integer value for the sample size used to calculate the GWAS summary statistics.'
    print 'Preparing to parse summary statistics'
    if p_dict['vbim'] is not None:
        bimfile = p_dict['vbim']
    elif p_dict['vgf'] is not None:
        bimfile = p_dict['vgf'] + '.bim'
    elif p_dict['gf'] is not None:
        bimfile = p_dict['gf'] + '.bim'
    else:
        print 'Set of validation SNPs is missing!  Please specify either a validation PLINK genotype file, ' \
              'or a PLINK BIM file with the SNPs of interest.'
    if os.path.isfile(p_dict['out']):
        print 'Output file (%s) already exists!  Delete, rename it, or use a different output file.'\
              % (p_dict['out'])
        raise Exception('Output file already exists!')

    h5f = h5py.File(p_dict['out'], 'w')
    
    ssp.parse_sum_stats(h5f, p_dict, bimfile)
    
    if not p_dict['vgf'] == None:
        assert p_dict['gf-format'] == 'PLINK', 'The validation genotype option currently only works with the PLINK format'
        coordinate_genotypes_ss_w_ld_ref(genotype_file=p_dict['vgf'], reference_genotype_file=p_dict['gf'],
                                         genetic_map_dir=p_dict['gmdir'], check_mafs=p_dict['check-maf'],
                                         hdf5_file=h5f, min_maf=p_dict['maf'], skip_coordination=p_dict['skip_coordination'])
    else:
        if p_dict['gf-format'] == 'PLINK':
            coordinate_genot_ss(genotype_file=p_dict['gf'], genetic_map_dir=p_dict['gmdir'], check_mafs=p_dict['check-maf'],
                                hdf5_file=h5f, min_maf=p_dict['maf'], skip_coordination=p_dict['skip-coordination'])
        else:
            raise Exception('Unknown genotype file format: %s' %
                            p_dict['gf-format'])

    h5f.close()


if __name__ == '__main__':
    main()
