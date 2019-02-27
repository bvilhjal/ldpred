"""
Code for handling plink files.

Uses plinkio.
"""
import scipy as sp
from plinkio import plinkfile

def get_chrom_dict(loci, chromosomes):
    chr_dict = {}
    for chrom in chromosomes:
        chr_str = 'chrom_%s' % chrom
        chr_dict[chr_str] = {'sids':[], 'snp_indices':[], 'positions':[], 'nts':[]}
      
    for i, l in enumerate(loci):
        chrom = l.chromosome
        pos = l.bp_position
        chr_str = 'chrom_%d' % chrom
        chr_dict[chr_str]['sids'].append(l.name)
        chr_dict[chr_str]['snp_indices'].append(i)
        chr_dict[chr_str]['positions'].append(pos)
        chr_dict[chr_str]['nts'].append([l.allele1, l.allele2])
      
    print('Genotype dictionary filled')
    return chr_dict


def parse_plink_snps(genotype_file, snp_indices):
    plinkf = plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()
    num_individs = len(samples)
    num_snps = len(snp_indices)
    raw_snps = sp.empty((num_snps, num_individs), dtype='int8')
    # If these indices are not in order then we place them in the right place while parsing SNPs.
    snp_order = sp.argsort(snp_indices)
    ordered_snp_indices = list(snp_indices[snp_order])
    ordered_snp_indices.reverse()
    # Iterating over file to load SNPs
    snp_i = 0
    next_i = ordered_snp_indices.pop()
    line_i = 0
    max_i = ordered_snp_indices[0]
    while line_i <= max_i: 
        if line_i < next_i:
            next(plinkf)
        elif line_i == next_i:
            line = next(plinkf)
            snp = sp.array(line, dtype='int8')
            bin_counts = line.allele_counts()
            if bin_counts[-1] > 0:
                mode_v = sp.argmax(bin_counts[:2])
                snp[snp == 3] = mode_v
            s_i = snp_order[snp_i]
            raw_snps[s_i] = snp
            if line_i < max_i:
                next_i = ordered_snp_indices.pop()
            snp_i += 1
        line_i += 1
    plinkf.close()
    assert snp_i == len(raw_snps), 'Parsing SNPs from plink file failed.'
    num_indivs = len(raw_snps[0])
    freqs = sp.sum(raw_snps, 1, dtype='float32') / (2 * float(num_indivs))
    return raw_snps, freqs


def get_phenotypes(plinkf):
    samples = plinkf.get_samples()
    num_individs = len(samples)
    Y = [s.phenotype for s in samples]
    fids = [s.fid for s in samples]
    iids = [s.iid for s in samples]
    unique_phens = sp.unique(Y)
    if len(unique_phens) == 1:
        print('Unable to find phenotype values.')
        has_phenotype = False
    elif len(unique_phens) == 2:
        cc_bins = sp.bincount(Y)
        assert len(cc_bins) == 2, 'Problems with loading phenotype'
        print('Loaded %d controls and %d cases' % (cc_bins[0], cc_bins[1]))
        has_phenotype = True
    else:
        print('Found quantitative phenotype values')
        has_phenotype = True
    return {'has_phenotype':has_phenotype, 'fids':fids, 'iids':iids, 'phenotypes':Y, 'num_individs':num_individs}