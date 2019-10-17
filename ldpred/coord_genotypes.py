#!/usr/bin/env python

import os
import scipy as sp
import gzip
import h5py
import sys
from ldpred import sum_stats_parsers
from ldpred import reporting
from ldpred import util
from ldpred import plinkfiles
from plinkio import plinkfile
import time


def _verify_coord_data_(data_dict):
    """
    Verify that merged data is ok
    """
    num_snps = len(data_dict['raw_snps_ref'])
    assert num_snps ==len(data_dict['snp_stds_ref']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['snp_means_ref']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['freqs_ref']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['ps']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['positions']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['nts']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['sids']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['betas']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['log_odds']), 'Inconsistencies in coordinated data sizes'
    assert num_snps ==len(data_dict['ns']), 'Inconsistencies in coordinated data sizes'
    if 'raw_snps_val' in data_dict:
        assert num_snps ==len(data_dict['raw_snps_val']), 'Inconsistencies in coordinated data sizes'
        assert num_snps ==len(data_dict['snp_stds_val']), 'Inconsistencies in coordinated data sizes'
        assert num_snps ==len(data_dict['snp_means_val']), 'Inconsistencies in coordinated data sizes'
        assert num_snps ==len(data_dict['freqs_val']), 'Inconsistencies in coordinated data sizes'
    


def write_coord_data(cord_data_g, coord_dict, debug=False):
    _verify_coord_data_(coord_dict)
    if debug:
        print('Storing coordinated data to HDF5 file.')
    ofg = cord_data_g.create_group(coord_dict['chrom'])
    ofg.create_dataset('raw_snps_ref', data=coord_dict['raw_snps_ref'], compression='lzf')
    ofg.create_dataset('snp_stds_ref', data=coord_dict['snp_stds_ref'])
    ofg.create_dataset('snp_means_ref', data=coord_dict['snp_means_ref'])
    ofg.create_dataset('freqs_ref', data=coord_dict['freqs_ref'])
    if 'raw_snps_val' in coord_dict:
        ofg.create_dataset('raw_snps_val', data=coord_dict['raw_snps_val'], compression='lzf')
        ofg.create_dataset('snp_stds_val', data=coord_dict['snp_stds_val'])
        ofg.create_dataset('snp_means_val', data=coord_dict['snp_means_val'])
        ofg.create_dataset('freqs_val', data=coord_dict['freqs_val'])
        ofg.create_dataset('log_odds_prs', data=coord_dict['log_odds_prs'])

    ofg.create_dataset('ps', data=coord_dict['ps'])
    ofg.create_dataset('positions', data=coord_dict['positions'])
    ofg.create_dataset('nts', data=sp.array(coord_dict['nts'],dtype=util.nts_dtype))
    ofg.create_dataset('sids', data=sp.array(coord_dict['sids'],dtype=util.sids_dtype))
    ofg.create_dataset('betas', data=coord_dict['betas'])
    ofg.create_dataset('log_odds', data=coord_dict['log_odds'])
    ofg.create_dataset('ns', data=coord_dict['ns'])

    
    if coord_dict['genetic_map'] is not None:
        ofg.create_dataset('genetic_map', data=coord_dict['genetic_map'])

                   
def get_snp_stds(raw_snps):
    return sp.std(raw_snps, axis=1, dtype='float32')


def get_mean_sample_size(n, cord_data_g):

    if n is None:
        all_ns = []
        for chrom_str in util.chromosomes_list:
            if chrom_str in cord_data_g:
                g = cord_data_g[chrom_str]
                all_ns.extend(g['ns'][...])
        assert all_ns is not None, 'Sample size missing. Please use --N flag, or ensure they are parsed as part of the summary statistics.'
        mean_n = sp.mean(all_ns)
    else:
        mean_n = n
    return mean_n



def coordinate_datasets(reference_genotype_file, hdf5_file, summary_dict,
                        validation_genotype_file=None,
                        genetic_map_dir=None,
                        min_maf=0.01,
                        skip_coordination=False, 
                        max_freq_discrep = 0.15,
                        debug=False):
    
    summary_dict[3.9]={'name':'dash', 'value':'Coordination'}
    t0 = time.time()
    if validation_genotype_file is not None:
        print('Coordinating datasets (Summary statistics, LD reference genotypes, and Validation genotypes).')
    else:
        print('Coordinating datasets (Summary statistics and LD reference genotypes).')
        
    plinkf = plinkfile.PlinkFile(reference_genotype_file)

    # Figure out chromosomes and positions.
    if debug:
        print('Parsing plinkf_dict_val reference genotypes')
    loci = plinkf.get_loci()
    plinkf.close()
    summary_dict[4]={'name':'Num individuals in LD Reference data:','value':plinkfiles.get_num_indivs(reference_genotype_file)}
    summary_dict[4.1]={'name':'SNPs in LD Reference data:','value':len(loci)}
    gf_chromosomes = [l.chromosome for l in loci]
    
    chromosomes = sp.unique(gf_chromosomes)
    chromosomes.sort()

    chr_dict = plinkfiles.get_chrom_dict(loci, chromosomes, debug)
    
    if validation_genotype_file is not None:
        if debug:
            print('Parsing LD validation bim file')
        plinkf_val = plinkfile.PlinkFile(validation_genotype_file)

        # Loads only the individuals... 
        plinkf_dict_val = plinkfiles.get_phenotypes(plinkf_val)
        
        loci_val = plinkf_val.get_loci()
        plinkf_val.close()
        summary_dict[5]={'name':'SNPs in Validation data:','value':len(loci_val)}

        chr_dict_val = plinkfiles.get_chrom_dict(loci_val, chromosomes, debug)

        # Open HDF5 file and prepare out data
        assert not 'iids' in hdf5_file, 'Something is wrong with the HDF5 file, no individuals IDs were found.'
        if plinkf_dict_val['has_phenotype']:
            hdf5_file.create_dataset('y', data=plinkf_dict_val['phenotypes'])
            summary_dict[6]={'name':'Num validation phenotypes:','value':plinkf_dict_val['num_individs']}
   
        hdf5_file.create_dataset('fids', data=sp.array(plinkf_dict_val['fids'], dtype=util.fids_dtype))
        hdf5_file.create_dataset('iids', data=sp.array(plinkf_dict_val['iids'], dtype=util.iids_dtype))

        maf_adj_risk_scores = sp.zeros(plinkf_dict_val['num_individs'])

    
    # Now summary statistics
    ssf = hdf5_file['sum_stats']
    cord_data_g = hdf5_file.create_group('cord_data')

    # corr_list = []



    chromosomes_found = set()
    num_snps_common_before_filtering =0
    num_snps_common_after_filtering =0
    tot_num_non_matching_nts = 0
    tot_num_non_supported_nts = 0
    tot_num_ambig_nts = 0
    tot_num_freq_discrep_filtered_snps = 0
    tot_num_maf_filtered_snps = 0
    tot_g_ss_nt_concord_count = 0
    if validation_genotype_file is not None:
        tot_g_vg_nt_concord_count = 0
        tot_vg_ss_nt_concord_count = 0
        
    # Now iterate over chromosomes
    chrom_i = 0
    for chrom in chromosomes:
        chrom_i +=1
        if not debug:
            sys.stdout.write('\r%0.2f%%' % (100.0 * (float(chrom_i) / (len(chromosomes)+1))))
            sys.stdout.flush()            
        try:
            chr_str = 'chrom_%d' % chrom
            ssg = ssf[chr_str]
                    
        except Exception as err_str:
                print(err_str)
                print('Did not find chromosome %d in SS dataset.'%chrom)
                print('Continuing.')
                continue
        
        if debug:
            print('Coordinating data for chromosome %s' % chr_str)

        chromosomes_found.add(chrom)
        
        #Get summary statistics chromosome group
        ssg = ssf['chrom_%d' % chrom]
        ss_sids = (ssg['sids'][...]).astype(util.sids_u_dtype)
        if validation_genotype_file is not None:
            chrom_d_val = chr_dict_val[chr_str]
            vg_sids = chrom_d_val['sids']
            common_sids = sp.intersect1d(ss_sids, vg_sids)
            
            # A map from sid to index for validation data        
            vg_sid_dict = {}
            for i, sid in enumerate(vg_sids):
                vg_sid_dict[sid] = i
        else:
            common_sids = ss_sids

        # A map from sid to index for summary stats        
        ss_sid_dict = {}
        for i, sid in enumerate(ss_sids):
            ss_sid_dict[sid] = i

        #The indices to retain for the LD reference genotypes
        chrom_d = chr_dict[chr_str]
        g_sids = chrom_d['sids']
        common_sids = sp.intersect1d(common_sids, g_sids)
        
        # A map from sid to index for LD reference data        
        g_sid_dict = {}
        for i, sid in enumerate(g_sids):
            g_sid_dict[sid] = i

        if debug:
            print('Found %d SNPs on chrom %d that were common across all datasets' % (len(common_sids), chrom))
            print('Ordering SNPs by genomic positions (based on LD reference genotypes).')
        
        g_snp_map = []
        for sid in common_sids:
            g_snp_map.append(g_sid_dict[sid])
        # order by positions (based on LD reference file)
        g_positions = sp.array(chrom_d['positions'])[g_snp_map]
        order = sp.argsort(g_positions)

        g_snp_map = sp.array(g_snp_map)[order]
        g_snp_map = g_snp_map.tolist()
        common_sids = sp.array(common_sids)[order]


        # Get the ordered sum stats SNPs indices.
        ss_snp_map = []
        for sid in common_sids:
            ss_snp_map.append(ss_sid_dict[sid])


        # Get the ordered validation SNPs indices
        if validation_genotype_file is not None:
            vg_snp_map = []
            for sid in common_sids:
                vg_snp_map.append(vg_sid_dict[sid])
            vg_nts = sp.array(chrom_d_val['nts'])
            vg_nts_ok = sp.array(vg_nts)[vg_snp_map]


        g_nts = sp.array(chrom_d['nts'])
        ss_nts = (ssg['nts'][...]).astype(util.nts_u_dtype)
        betas = ssg['betas'][...]
        log_odds = ssg['log_odds'][...]

        if 'freqs' in ssg:
            ss_freqs = ssg['freqs'][...]

        g_ss_nt_concord_count = sp.sum(
            g_nts[g_snp_map] == ss_nts[ss_snp_map]) / 2.0
        if validation_genotype_file is not None:
            vg_ss_nt_concord_count = sp.sum(vg_nts_ok == ss_nts[ss_snp_map]) / 2.0
            g_vg_nt_concord_count = sp.sum(g_nts[g_snp_map] == vg_nts_ok) / 2.0
            if debug:
                print('Nucleotide concordance counts out of %d genotypes, vg-rg: %d ; vg-ss: %d' % (len(g_snp_map), g_vg_nt_concord_count, vg_ss_nt_concord_count))
            tot_vg_ss_nt_concord_count += vg_ss_nt_concord_count
            tot_g_vg_nt_concord_count += g_vg_nt_concord_count
        tot_g_ss_nt_concord_count += g_ss_nt_concord_count
        if debug:
            print('Nucleotide concordance counts out of %d genotypes, rg-ss: %d' % (len(g_snp_map), g_ss_nt_concord_count))

        num_freq_discrep_filtered_snps = 0
        num_non_matching_nts = 0
        num_non_supported_nts = 0
        num_ambig_nts = 0

        # Identifying which SNPs have nucleotides that are ok..
        ok_nts = []
        ok_indices = {'g': [], 'ss': []}
        if validation_genotype_file is not None:
            ok_indices['vg']=[]

        #Now loop over SNPs to coordinate nucleotides.        
        if validation_genotype_file is not None:
            for g_i, vg_i, ss_i in zip(g_snp_map, vg_snp_map, ss_snp_map):
    
                # To make sure, is the SNP id the same?
                assert g_sids[g_i] == vg_sids[vg_i] == ss_sids[ss_i], 'Some issues with coordinating the genotypes.'
    
                g_nt = g_nts[g_i]
                if not skip_coordination:
    
                    vg_nt = vg_nts[vg_i]
                    ss_nt = ss_nts[ss_i]
    
                    # Is the nucleotide ambiguous.
                    g_nt = [g_nts[g_i][0], g_nts[g_i][1]]
                    if tuple(g_nt) in util.ambig_nts:
                        num_ambig_nts += 1
                        continue
    
                    # First check if nucleotide is sane?
                    if (not g_nt[0] in util.valid_nts) or (not g_nt[1] in util.valid_nts):
                        num_non_supported_nts += 1
                        continue
    
                    os_g_nt = sp.array(
                        [util.opp_strand_dict[g_nt[0]], util.opp_strand_dict[g_nt[1]]])
    
                    flip_nts = False
                    
                    #Coordination is a bit more complicate when validation genotypes are provided..
                    if not ((sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)) and (sp.all(g_nt == vg_nt) or sp.all(os_g_nt == vg_nt))):
                        if sp.all(g_nt == vg_nt) or sp.all(os_g_nt == vg_nt):
                            flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (
                                os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
                            # Try flipping the SS nt
                            if flip_nts:
                                betas[ss_i] = -betas[ss_i]
                                log_odds[ss_i] = -log_odds[ss_i]
                                if 'freqs' in ssg:
                                    ss_freqs[ss_i] = 1 - ss_freqs[ss_i]
                            else:
                                if debug:
                                    print("Nucleotides don't match after all?: g_sid=%s, ss_sid=%s, g_i=%d, ss_i=%d, g_nt=%s, ss_nt=%s" % \
                                          (g_sids[g_i], ss_sids[ss_i], g_i,
                                           ss_i, str(g_nt), str(ss_nt)))
                                num_non_matching_nts += 1
                                continue
    
                        else:
                            num_non_matching_nts += 1
                            continue
                            # Opposite strand nucleotides
    
                # everything seems ok.
                ok_indices['g'].append(g_i)
                ok_indices['vg'].append(vg_i)
                ok_indices['ss'].append(ss_i)
    
                ok_nts.append(g_nt)
        else:
            for g_i, ss_i in zip(g_snp_map, ss_snp_map):
    
                # To make sure, is the SNP id the same?
                assert g_sids[g_i] == ss_sids[ss_i], 'Some issues with coordinating the genotypes.'
    
                g_nt = g_nts[g_i]
                if not skip_coordination:
    
                    ss_nt = ss_nts[ss_i]
    
                    # Is the nucleotide ambiguous.
                    g_nt = [g_nts[g_i][0], g_nts[g_i][1]]
                    if tuple(g_nt) in util.ambig_nts:
                        num_ambig_nts += 1
                        continue
    
                    # First check if nucleotide is sane?
                    if (not g_nt[0] in util.valid_nts) or (not g_nt[1] in util.valid_nts):
                        num_non_matching_nts += 1
                        continue
    
                    os_g_nt = sp.array(
                        [util.opp_strand_dict[g_nt[0]], util.opp_strand_dict[g_nt[1]]])
    
                    flip_nts = False
                    
                    #Coordination is a bit more complicate when validation genotypes are provided..
                    if not (sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)):
                        flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (
                            os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
                        
                        # Try flipping the SS nt
                        if flip_nts:
                            betas[ss_i] = -betas[ss_i]
                            log_odds[ss_i] = -log_odds[ss_i]
                            if 'freqs' in ssg and ss_freqs[ss_i]>0:
                                ss_freqs[ss_i] = 1.0 - ss_freqs[ss_i]
                        else:
                            if debug:
                                print("Nucleotides don't match after all?: g_sid=%s, ss_sid=%s, g_i=%d, ss_i=%d, g_nt=%s, ss_nt=%s" % \
                                      (g_sids[g_i], ss_sids[ss_i], g_i,
                                       ss_i, str(g_nt), str(ss_nt)))
                            num_non_matching_nts += 1
                            continue
                   
                # everything seems ok.
                ok_indices['g'].append(g_i)
                ok_indices['ss'].append(ss_i)
                ok_nts.append(g_nt)
                
        if debug:
            print('%d SNPs had ambiguous nucleotides.' % num_ambig_nts)
            print('%d SNPs were excluded due to nucleotide issues.' % num_non_matching_nts)

        
        # Resorting by position
        positions = sp.array(chrom_d['positions'])[ok_indices['g']]

        # Now parse SNPs ..
        snp_indices = sp.array(chrom_d['snp_indices'])
        # Pinpoint where the SNPs are in the file.
        snp_indices = snp_indices[ok_indices['g']]
        raw_snps, freqs = plinkfiles.parse_plink_snps(
            reference_genotype_file, snp_indices)
        snp_stds = get_snp_stds(raw_snps)
        snp_means = sp.mean(raw_snps, axis=1, dtype='float32')

        betas = betas[ok_indices['ss']]  
        log_odds = log_odds[ok_indices['ss']]  

        ns = ssg['ns'][...][ok_indices['ss']]
        ps = ssg['ps'][...][ok_indices['ss']]
        nts = sp.array(ok_nts)  
        sids = (ssg['sids'][...]).astype(util.sids_u_dtype)
        sids = sids[ok_indices['ss']]

        #Parse validation genotypes, if available
        if validation_genotype_file is not None:
            snp_indices_val = sp.array(chrom_d_val['snp_indices'])
            # Pinpoint where the SNPs are in the file.
            snp_indices_val = snp_indices_val[ok_indices['vg']]
            raw_snps_val, freqs_val = plinkfiles.parse_plink_snps(
                validation_genotype_file, snp_indices_val)
    
            snp_stds_val = get_snp_stds(raw_snps_val)
            snp_means_val = freqs_val * 2

        # Check SNP frequencies, screen for possible problems..
        if max_freq_discrep<1 and 'freqs' in ssg:
            ss_freqs = ss_freqs[ok_indices['ss']]
            ok_freq_snps = sp.logical_or(sp.absolute(ss_freqs - freqs) < max_freq_discrep,sp.absolute(ss_freqs + freqs-1) < max_freq_discrep) #Array of np.bool values
            ok_freq_snps = sp.logical_or(ok_freq_snps,ss_freqs<=0) #Only consider SNPs that actually have frequencies
            num_freq_discrep_filtered_snps = len(ok_freq_snps)- sp.sum(ok_freq_snps)
            assert num_freq_discrep_filtered_snps>=0, "Problems when filtering SNPs with frequency discrepencies"
            if num_freq_discrep_filtered_snps>0:
                # Filter freq_discrepancy_snps
                raw_snps = raw_snps[ok_freq_snps]
                snp_stds = snp_stds[ok_freq_snps]
                snp_means = snp_means[ok_freq_snps]
                freqs = freqs[ok_freq_snps]
                ps = ps[ok_freq_snps]
                ns = ns[ok_freq_snps]
                positions = positions[ok_freq_snps]
                nts = nts[ok_freq_snps]
                sids = sids[ok_freq_snps]
                betas = betas[ok_freq_snps]
                log_odds = log_odds[ok_freq_snps]
                if validation_genotype_file is not None:
                    raw_snps_val = raw_snps_val[ok_freq_snps]
                    snp_stds_val = snp_stds_val[ok_freq_snps]
                    snp_means_val = snp_means_val[ok_freq_snps]
                    freqs_val = freqs_val[ok_freq_snps]
            if debug:
                print('Filtered %d SNPs due to frequency discrepancies'%num_freq_discrep_filtered_snps)

        # Filter minor allele frequency SNPs.
        maf_filter = (freqs > min_maf) * (freqs < (1 - min_maf))
        num_maf_filtered_snps = len(maf_filter)-sp.sum(maf_filter)
        assert num_maf_filtered_snps>=0, "Problems when filtering SNPs with low minor allele frequencies"
        if num_maf_filtered_snps>0:
            raw_snps = raw_snps[maf_filter]
            snp_stds = snp_stds[maf_filter]
            snp_means = snp_means[maf_filter]
            freqs = freqs[maf_filter]
            ps = ps[maf_filter]
            ns = ns[maf_filter]
            positions = positions[maf_filter]
            nts = nts[maf_filter]
            sids = sids[maf_filter]
            betas = betas[maf_filter]
            log_odds = log_odds[maf_filter]
            if validation_genotype_file is not None:
                raw_snps_val = raw_snps_val[maf_filter]
                snp_stds_val = snp_stds_val[maf_filter]
                snp_means_val = snp_means_val[maf_filter]
                freqs_val = freqs_val[maf_filter]
            if debug:
                print('Filtered %d SNPs due to low MAF'%num_maf_filtered_snps)

        genetic_map = []
        if genetic_map_dir is not None:
            with gzip.open(genetic_map_dir + 'chr%d.interpolated_genetic_map.gz' % chrom) as f:
                for line in f:
                    l = line.split()
#                     if l[0] in sid_set:
#                         genetic_map.append(l[0])
        else:
            genetic_map = None

        coord_data_dict = {'chrom': 'chrom_%d' % chrom, 
                           'raw_snps_ref': raw_snps, 
                           'snp_stds_ref': snp_stds, 
                           'snp_means_ref': snp_means, 
                           'freqs_ref': freqs,
                           'ps': ps,
                           'ns': ns,
                           'positions': positions,
                           'nts': nts,
                           'sids': sids,
                           'genetic_map': genetic_map,
                           'betas': betas,
                           'log_odds': log_odds}
        if validation_genotype_file is not None:
            maf_adj_prs = sp.dot(log_odds, raw_snps_val)
            if debug and plinkf_dict_val['has_phenotype']:
                maf_adj_corr = sp.corrcoef(plinkf_dict_val['phenotypes'], maf_adj_prs)[0, 1]
                print('Log odds, per genotype PRS correlation w phenotypes for chromosome %d was %0.4f' % (chrom, maf_adj_corr))
            coord_data_dict['raw_snps_val']=raw_snps_val
            coord_data_dict['snp_stds_val']=snp_stds_val
            coord_data_dict['snp_means_val']=snp_means_val
            coord_data_dict['freqs_val']=freqs_val
            coord_data_dict['log_odds_prs']=maf_adj_prs
            maf_adj_risk_scores += maf_adj_prs
         
         
        write_coord_data(cord_data_g, coord_data_dict, debug=debug)
        if debug:
            print('%d SNPs were retained on chromosome %d.' % (len(sids), chrom))
        
        
        num_snps_common_before_filtering += len(common_sids)
        num_snps_common_after_filtering += len(sids)
        tot_num_ambig_nts += num_ambig_nts
        tot_num_non_supported_nts += num_non_supported_nts
        tot_num_non_matching_nts += num_non_matching_nts
        tot_num_freq_discrep_filtered_snps += num_freq_discrep_filtered_snps
        tot_num_maf_filtered_snps += num_maf_filtered_snps

    if not debug:
        sys.stdout.write('\r%0.2f%%\n' % (100.0))
        sys.stdout.flush()                        


    # Now calculate the prediction r^2
    if validation_genotype_file:
        if debug and plinkf_dict_val['has_phenotype']:
            maf_adj_corr = sp.corrcoef(
                plinkf_dict_val['phenotypes'], maf_adj_risk_scores)[0, 1]
            print('Log odds, per PRS correlation for the whole genome was %0.4f (r^2=%0.4f)' % (maf_adj_corr, maf_adj_corr ** 2))
            print('Overall nucleotide concordance counts: rg_vg: %d, rg_ss: %d, vg_ss: %d' % (tot_g_vg_nt_concord_count, tot_g_ss_nt_concord_count, tot_vg_ss_nt_concord_count))
    else:
        if debug:
            print('Overall nucleotide concordance counts, rg_ss: %d' % (tot_g_ss_nt_concord_count))        
    
    summary_dict[7]={'name':'Num chromosomes used:','value':len(chromosomes_found)}
    summary_dict[8]={'name':'SNPs common across datasets:','value':num_snps_common_before_filtering}
    summary_dict[9]={'name':'SNPs retained after filtering:','value':num_snps_common_after_filtering}
    if tot_num_ambig_nts>0:
        summary_dict[10]={'name':'SNPs w ambiguous nucleotides filtered:','value':tot_num_ambig_nts}
    if tot_num_non_supported_nts>0:
        summary_dict[10.1]={'name':'SNPs w unknown/unsupported nucleotides filtered:','value':tot_num_non_supported_nts}
    if tot_num_non_matching_nts>0:
        summary_dict[11]={'name':'SNPs w other nucleotide discrepancies filtered:','value':tot_num_non_matching_nts}
    if min_maf>0:
        summary_dict[12]={'name':'SNPs w MAF<%0.3f filtered:'%min_maf,'value':tot_num_maf_filtered_snps}
    if max_freq_discrep<0.5:
        summary_dict[13]={'name':'SNPs w allele freq discrepancy > %0.3f filtered:'%max_freq_discrep,'value':tot_num_freq_discrep_filtered_snps}

    t1 = time.time()
    t = (t1 - t0)
    summary_dict[13.9]={'name':'dash', 'value':'Running times'}
    summary_dict[15]={'name':'Run time for coordinating datasets:','value': '%d min and %0.2f sec'%(t / 60, t % 60)}



def main(p_dict):

    bimfile = None
    if p_dict['vbim'] is not None:
        bimfile = p_dict['vbim']
    elif p_dict['vgf'] is not None:
        bimfile = p_dict['vgf'] + '.bim'
    elif p_dict['gf'] is not None:
        bimfile = p_dict['gf'] + '.bim'
    else:
        print('Set of validation SNPs is missing!  Please specify either a validation PLINK genotype file, ' \
              'or a PLINK BIM file with the SNPs of interest.')
    if os.path.isfile(p_dict['out']):
        print('Output file (%s) already exists!  Delete, rename it, or use a different output file.'\
              % (p_dict['out']))
        raise Exception('Output file already exists!')

    h5f = h5py.File(p_dict['out'], 'w')
    
    summary_dict = {}
    summary_dict[0]={'name':'Summary statistics filename:','value':p_dict['ssf']}
    summary_dict[1]={'name':'LD reference genotypes filename:','value':p_dict['gf']}
    summary_dict[3]={'name':'Coordinated data output filename:','value':p_dict['out']}
    if p_dict['vgf'] is not None:
        summary_dict[2]={'name':'Validation genotypes filename:','value':p_dict['vgf']}

    sum_stats_parsers.parse_sum_stats(h5f, p_dict, bimfile, summary_dict)
    coordinate_datasets(p_dict['gf'], h5f,summary_dict,
                        validation_genotype_file=p_dict['vgf'], 
                        max_freq_discrep=p_dict['max_freq_discrep'],
                        min_maf=p_dict['maf'], 
                        skip_coordination=p_dict['skip_coordination'], 
                        debug=p_dict['debug'])
    h5f.close()
    reporting.print_summary(summary_dict, 'Summary of coordination step')
    return summary_dict
