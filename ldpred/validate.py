import os
import scipy as sp
from scipy import linalg
import h5py
from ldpred import plinkfiles
from ldpred import util
import time
from ldpred import reporting

def get_prs(genotype_file, rs_id_map, phen_map=None, only_score = False, verbose=False):
    plinkf = plinkfiles.plinkfile.PlinkFile(genotype_file)
    samples = plinkf.get_samples()

    if not only_score:
        # 1. Figure out indiv filter and get true phenotypes
        indiv_filter = sp.zeros(len(samples), dtype='bool8')
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
                    if 'pcs' in list(phen_map[sample.iid].keys()):
                        pcs.append(phen_map[sample.iid]['pcs'])
                    if 'sex' in list(phen_map[sample.iid].keys()):
                        sex.append(phen_map[sample.iid]['sex'])
                    if 'covariates' in list(phen_map[sample.iid].keys()):
                        covariates.append(phen_map[sample.iid]['covariates'])
            if len(pcs) > 0:
                assert len(pcs) == len(
                    true_phens), 'PC information missing for some individuals with phenotypes'
            if len(sex) > 0:
                assert len(sex) == len(
                    true_phens), 'Sex information missing for some individuals with phenotypes'
            if len(covariates) > 0:
                assert len(covariates) == len(
                    true_phens), 'Covariates missing for some individuals with phenotypes'
        else:
            for samp_i, sample in enumerate(samples):
                if sample.affection != 2:
                    indiv_filter[samp_i] = True
                    true_phens.append(sample.affection)
                    iids.append(sample.iid)
    
        num_individs = sp.sum(indiv_filter)
        assert num_individs > 0, 'Issues in parsing the phenotypes and/or PCs?'
    
        assert not sp.any(sp.isnan(
            true_phens)), 'Phenotypes appear to have some NaNs, or parsing failed.'
    
        if verbose:
            print('%d individuals have phenotype and genotype information.' % num_individs)
    else:
        iids = [sample.iid for sample in samples]
        num_individs = len(samples)
    num_non_matching_nts = 0
    num_flipped_nts = 0

    pval_derived_effects_prs = sp.zeros(num_individs)
    # If these indices are not in order then we place them in the right place
    # while parsing SNPs.
    if verbose:
        print('Iterating over BED file to calculate risk scores.')
    locus_list = plinkf.get_loci()
    snp_i = 0

    for locus, row in zip(locus_list, plinkf):
        upd_pval_beta = 0
        try:
            # Check rs-ID
            sid = locus.name
            rs_info = rs_id_map[sid]
        except Exception:  # Move on if rsID not found.
            continue

        if rs_info['upd_pval_beta'] == 0:
            continue

        # Check whether the nucleotides are OK, and potentially flip it.
        ss_nt = rs_info['nts']
        g_nt = [locus.allele1, locus.allele2]
        flip_nts = False
        os_g_nt = sp.array(
            [util.opp_strand_dict[g_nt[0]], util.opp_strand_dict[g_nt[1]]])
        if not (sp.all(g_nt == ss_nt) or sp.all(os_g_nt == ss_nt)):
            # Opposite strand nucleotides
            flip_nts = (g_nt[1] == ss_nt[0] and g_nt[0] == ss_nt[1]) or (
                os_g_nt[1] == ss_nt[0] and os_g_nt[0] == ss_nt[1])
            if flip_nts:
                upd_pval_beta = -rs_info['upd_pval_beta']
                num_flipped_nts += 1
            else:
                num_non_matching_nts += 1
                continue
        else:
            upd_pval_beta = rs_info['upd_pval_beta']

        # Parse SNP, and fill in the blanks if necessary.
        if only_score:
            snp = sp.array(row, dtype='int8')
        else:
            snp = sp.array(row, dtype='int8')[indiv_filter]
        bin_counts = row.allele_counts()
        if bin_counts[-1] > 0:
            mode_v = sp.argmax(bin_counts[:2])
            snp[snp == 3] = mode_v


        # Update scores and move on.
        pval_derived_effects_prs += snp * upd_pval_beta
        assert not sp.any(sp.isnan(pval_derived_effects_prs)
                          ), 'Some individual weighted effects risk scores are NANs (not a number).  They are corrupted.'

        if snp_i > 0 and snp_i % 100000 == 0:
            if verbose:
                print('%d SNPs parsed'%snp_i)
                print('Number of non-matching NTs: %d' % num_non_matching_nts)
            if not only_score:
                pval_eff_r2 = (sp.corrcoef(
                    pval_derived_effects_prs, true_phens)[0, 1]) ** 2
                if verbose:
                    print('Current PRS r2: %0.4f' % pval_eff_r2)

        snp_i += 1

    plinkf.close()

    if verbose:
        print("DONE!")
        print('Number of non-matching NTs: %d' % num_non_matching_nts)
        print('Number of flipped NTs: %d' % num_flipped_nts)
    if not only_score:
        pval_eff_corr = sp.corrcoef(pval_derived_effects_prs, true_phens)[0, 1]
        pval_eff_r2 = pval_eff_corr ** 2

        if verbose:
            print('Current PRS correlation: %0.4f' % pval_eff_corr)
            print('Current PRS r2: %0.4f' % pval_eff_r2)
    
        ret_dict = {'pval_derived_effects_prs': pval_derived_effects_prs,
                    'true_phens': true_phens[:], 'iids': iids, 'prs_r2':pval_eff_r2, 'prs_corr':pval_eff_corr}
    
        if len(pcs) > 0:
            ret_dict['pcs'] = pcs
        if len(sex) > 0:
            ret_dict['sex'] = sex
        if len(covariates) > 0:
            ret_dict['covariates'] = covariates

    else:    
        ret_dict = {'pval_derived_effects_prs': pval_derived_effects_prs,'iids': iids}

    
    return ret_dict


def parse_phen_file(pf, pf_format, verbose=False, summary_dict=None):
    print(pf)
    phen_map = {}
    num_phens_found = 0
    if pf != None:
        if verbose:
            print('Parsing Phenotypes')
        if pf_format == 'FAM':
            """
            Individual's family ID ('FID')
            Individual's within-family ID ('IID'; cannot be '0')
            Within-family ID of father ('0' if father isn't in dataset)
            Within-family ID of mother ('0' if mother isn't in dataset)
            Sex code ('1' = male, '2' = female, '0' = unknown)
            Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)
            """
            with open(pf, 'r') as f:
                for line in f:
                    l = line.split()
                    iid = l[1]
                    sex = int(l[4])
                    phen = float(l[5])
                    if sex != 0 and phen != -9:
                        phen_map[iid] = {'phen': phen, 'sex': sex}
                        num_phens_found +=1

            summary_dict[1]={'name':'Phenotype file (plink format):','value':pf}
            
        if pf_format == 'STANDARD':
            """
            IID   PHE
            """
            with open(pf, 'r') as f:
                for line in f:
                    l = line.split()
                    iid = l[0]
                    phen = float(l[1])
                    phen_map[iid] = {'phen': phen}
                    num_phens_found +=1

            summary_dict[1]={'name':'Phenotype file (STANDARD format):','value':pf}

        if pf_format == 'LSTANDARD':
            """
            FID IID PHE
            """
            with open(pf, 'r') as f:
                # Read first line and see if it is header
                line = f.readline()
                l = line.split()
                try:
                    iid=l[1]
                    phen=float(l[2])
                    phen_map[iid] = {'phen':phen}
                    num_phens_found+=1
                except:
                    # Do nothing
                    pass
                for line in f:
                    l = line.split()
                    iid=l[1]
                    phen=float(l[2])
                    phen_map[iid] = {'phen':phen}
                    num_phens_found+=1
            summary_dict[1]={'name':'Phenotype file (LSTANDARD format):','value':pf}
 
        elif pf_format == 'S2':
            """
            IID Age Sex Height_Inches
            """
            with open(pf, 'r') as f:
                for line in f:
                    l = line.split()
                    iid = l[0]
                    age = float(l[1])
                    if l[2] == 'Male':
                        sex = 1
                    elif l[2] == 'Female':
                        sex = 2
                    else:
                        raise Exception('Sex missing')
                    phen = float(l[3])
                    phen_map[iid] = {'phen': phen, 'age': age, 'sex': sex}
                    num_phens_found +=1
            summary_dict[1]={'name':'Phenotype file (S2 format):','value':pf}

    print("Parsed %d phenotypes successfully"%num_phens_found)
    return phen_map


def parse_ldpred_res(file_name):
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    ldpred_inf_beta    ldpred_beta
    1, 798959, rs11240777, C, T, -1.1901e-02, 3.2443e-03, 2.6821e-04
    """
    with open(file_name, 'r') as f:
        next(f)
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1, nt2]
            raw_beta = float(l[5])
            upd_pval_beta = float(l[6])
            rs_id_map[rs_id] = {'chrom': chrom, 'pos': pos, 'nts': nts, 'raw_beta': raw_beta,
                                'upd_pval_beta': upd_pval_beta}
    return rs_id_map


def parse_pt_res(file_name):
    non_zero_chromosomes = set()
    rs_id_map = {}
    """
    chrom    pos    sid    nt1    nt2    raw_beta    raw_pval_beta    upd_beta    upd_pval_beta
    1    798959    rs11240777    C    T    -1.1901e-02    -1.1901e-02    2.6821e-04    2.6821e-04
    """
    with open(file_name, 'r') as f:
        next(f)
        for line in f:
            l = line.split()
            chrom_str = l[0]
            chrom = int(chrom_str[6:])
            pos = int(l[1])
            rs_id = l[2].strip()
            nt1 = l[3].strip()
            nt2 = l[4].strip()
            nts = [nt1, nt2]
            raw_beta = float(l[5])
            upd_pval_beta = float(l[8])
            if raw_beta != 0:
                rs_id_map[rs_id] = {'chrom': chrom, 'pos': pos, 'nts': nts, 'raw_beta': raw_beta,
                                    'upd_pval_beta': upd_pval_beta}
                non_zero_chromosomes.add(chrom)

    return rs_id_map, non_zero_chromosomes

def write_scores_file(out_file, prs_dict, pval_derived_effects_prs, adj_pred_dict, 
                      output_regression_betas=False, weights_dict=None):
    num_individs = len(prs_dict['iids'])
    with open(out_file, 'w') as f:
        print ('Writing polygenic scores to file %s'%out_file)
        out_str = 'IID, true_phens, PRS'
        if 'sex' in prs_dict:
            out_str = out_str + ', sex'
        if 'pcs' in prs_dict:
            pcs_str = ', '.join(['PC%d' % (1 + pc_i)
                                 for pc_i in range(len(prs_dict['pcs'][0]))])
            out_str = out_str + ', ' + pcs_str
        out_str += '\n'
        f.write(out_str)
        for i in range(num_individs):
            out_str = '%s, %0.6e, %0.6e' % (prs_dict['iids'][i], prs_dict['true_phens'][i],
                                                     pval_derived_effects_prs[i])
            if 'sex' in prs_dict:
                out_str = out_str + ', %d' % prs_dict['sex'][i]
            if 'pcs' in prs_dict:
                pcs_str = ', '.join(map(str, prs_dict['pcs'][i]))
                out_str = out_str +', '+ pcs_str
            out_str += '\n'
            f.write(out_str)

    if len(list(adj_pred_dict.keys())) > 0:
        with open(out_file + '.adj', 'w') as f:
            adj_prs_labels = list(adj_pred_dict.keys())
            out_str = 'IID, true_phens, PRS, ' + \
                ', '.join(adj_prs_labels)
            out_str += '\n'
            f.write(out_str)
            for i in range(num_individs):
                out_str = '%s, %0.6e, %0.6e' % (prs_dict['iids'][i], prs_dict['true_phens'][i],
                                                       pval_derived_effects_prs[i])
                for adj_prs in adj_prs_labels:
                    out_str += ', %0.4f' % adj_pred_dict[adj_prs][i]
                out_str += '\n'
                f.write(out_str)
    if output_regression_betas and weights_dict != None:
        hdf5file = out_file + '.weights.hdf5'
        print ('Writing PRS regression weights to file %s'%hdf5file)
        oh5f = h5py.File(hdf5file, 'w')
        for k1 in list(weights_dict.keys()):
            kg = oh5f.create_group(k1)
            for k2 in weights_dict[k1]:
                kg.create_dataset(k2, data=sp.array(weights_dict[k1][k2]))
        oh5f.close()


def write_only_scores_file(out_file, prs_dict, pval_derived_effects_prs):
    num_individs = len(prs_dict['iids'])
    with open(out_file, 'w') as f:
        print ('Writing polygenic scores to file %s'%out_file)
        out_str = 'IID, PRS \n'
        f.write(out_str)
        for i in range(num_individs):
            out_str = '%s, %0.6e\n' % (prs_dict['iids'][i],
                                                     pval_derived_effects_prs[i])
            f.write(out_str)


def calc_risk_scores(bed_file, rs_id_map, phen_map, out_file=None, 
                     split_by_chrom=False, adjust_for_sex=False,
                     adjust_for_covariates=False, adjust_for_pcs=False, 
                     non_zero_chromosomes=None, only_score = False, 
                     verbose=False, summary_dict=None):
    print('Parsing PLINK bed file: %s' % bed_file)

    if split_by_chrom:
        num_individs = len(phen_map)
        assert num_individs > 0, 'No individuals found.  Problems parsing the phenotype file?'
        pval_derived_effects_prs = sp.zeros(num_individs)

        for i in range(1, 23):
            if non_zero_chromosomes is None or i in non_zero_chromosomes:
                genotype_file = bed_file + '_%i_keep' % i
                if os.path.isfile(genotype_file + '.bed'):
                    if verbose:
                        print('Working on chromosome %d' % i)
                    prs_dict = get_prs(genotype_file, rs_id_map, phen_map, only_score=only_score, verbose=verbose)

                    pval_derived_effects_prs += prs_dict['pval_derived_effects_prs']
            elif verbose:
                    print('Skipping chromosome')

    else:
        prs_dict = get_prs(bed_file, rs_id_map, phen_map, only_score=only_score, verbose=verbose)
        num_individs = len(prs_dict['iids'])
        pval_derived_effects_prs = prs_dict['pval_derived_effects_prs']

    if only_score:
        write_only_scores_file(out_file, prs_dict, pval_derived_effects_prs)
        res_dict = {}
    else:
        # Report prediction accuracy
        assert len(phen_map) > 0, 'No individuals found.  Problems parsing the phenotype file?'

        pval_eff_corr = sp.corrcoef(pval_derived_effects_prs, prs_dict['true_phens'])[0, 1]
        pval_eff_r2 = pval_eff_corr ** 2
    
        res_dict = {'pred_r2': pval_eff_r2}
    
        pval_derived_effects_prs.shape = (len(pval_derived_effects_prs), 1)
        true_phens = sp.array(prs_dict['true_phens'])
        true_phens.shape = (len(true_phens), 1)
    
        # Store covariate weights, slope, etc.
        weights_dict = {}
    
        # Store Adjusted predictions
        adj_pred_dict = {}
    
        # Direct effect
        Xs = sp.hstack([pval_derived_effects_prs, sp.ones((len(true_phens), 1))])
        (betas, rss00, r, s) = linalg.lstsq(
            sp.ones((len(true_phens), 1)), true_phens)
        (betas, rss, r, s) = linalg.lstsq(Xs, true_phens)
        pred_r2 = 1 - rss / rss00
        weights_dict['unadjusted'] = {
            'Intercept': betas[1][0], 'ldpred_prs_effect': betas[0][0]}
    
        if verbose:
            print('PRS correlation: %0.4f' % pval_eff_corr)
        print('Variance explained (Pearson R2) by PRS: %0.4f' % pred_r2)
    
        # Adjust for sex
        if adjust_for_sex and 'sex' in prs_dict and len(prs_dict['sex']) > 0:
            sex = sp.array(prs_dict['sex'])
            sex.shape = (len(sex), 1)
            (betas, rss0, r, s) = linalg.lstsq(
                sp.hstack([sex, sp.ones((len(true_phens), 1))]), true_phens)
            Xs = sp.hstack([pval_derived_effects_prs, sex,
                            sp.ones((len(true_phens), 1))])
            (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
            weights_dict['sex_adj'] = {
                'Intercept': betas[2][0], 'ldpred_prs_effect': betas[0][0], 'sex': betas[1][0]}
            if verbose:
                print('Fitted effects (betas) for PRS, sex, and intercept on true phenotype:', betas)
            adj_pred_dict['sex_adj'] = sp.dot(Xs, betas)
            pred_r2 = 1 - rss_pd / rss0
            print('Variance explained (Pearson R2) by PRS adjusted for Sex: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
            res_dict['Sex_adj_pred_r2'] = pred_r2
            pred_r2 = 1 - rss_pd / rss00
            print('Variance explained (Pearson R2) by PRS + Sex : %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
            res_dict['Sex_adj_pred_r2+Sex'] = pred_r2
    
        # Adjust for PCs
        if adjust_for_pcs and 'pcs' in prs_dict and len(prs_dict['pcs']) > 0:
            pcs = prs_dict['pcs']
            (betas, rss0, r, s) = linalg.lstsq(
                sp.hstack([pcs, sp.ones((len(true_phens), 1))]), true_phens)
            Xs = sp.hstack([pval_derived_effects_prs,
                            sp.ones((len(true_phens), 1)), pcs])
            (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
            weights_dict['pc_adj'] = {
                'Intercept': betas[1][0], 'ldpred_prs_effect': betas[0][0], 'pcs': betas[2][0]}
            adj_pred_dict['pc_adj'] = sp.dot(Xs, betas)
            pred_r2 = 1 - rss_pd / rss0
            print('Variance explained (Pearson R2) by PRS adjusted for PCs: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
            res_dict['PC_adj_pred_r2'] = pred_r2
            pred_r2 = 1 - rss_pd / rss00
            print('Variance explained (Pearson R2) by PRS + PCs: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
            res_dict['PC_adj_pred_r2+PC'] = pred_r2
    
            # Adjust for both PCs and Sex
            if adjust_for_sex and 'sex' in prs_dict and len(prs_dict['sex']) > 0:
                sex = sp.array(prs_dict['sex'])
                sex.shape = (len(sex), 1)
                (betas, rss0, r, s) = linalg.lstsq(
                    sp.hstack([sex, pcs, sp.ones((len(true_phens), 1))]), true_phens)
                Xs = sp.hstack([pval_derived_effects_prs, sex,
                                sp.ones((len(true_phens), 1)), pcs])
                (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
                weights_dict['sex_pc_adj'] = {
                    'Intercept': betas[2][0], 'ldpred_prs_effect': betas[0][0], 'sex': betas[1][0], 'pcs': betas[3][0]}
                adj_pred_dict['sex_pc_adj'] = sp.dot(Xs, betas)
                pred_r2 = 1 - rss_pd / rss0
                print('Variance explained (Pearson R2) by PRS adjusted for PCs and Sex: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
                res_dict['PC_Sex_adj_pred_r2'] = pred_r2
                pred_r2 = 1 - rss_pd / rss00
                print('Variance explained (Pearson R2) by PRS+PCs+Sex: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
                res_dict['PC_Sex_adj_pred_r2+PC_Sex'] = pred_r2
    
        # Adjust for covariates
        if adjust_for_covariates and 'covariates' in prs_dict and len(prs_dict['covariates']) > 0:
            covariates = prs_dict['covariates']
            (betas, rss0, r, s) = linalg.lstsq(
                sp.hstack([covariates, sp.ones((len(true_phens), 1))]), true_phens)
            Xs = sp.hstack([pval_derived_effects_prs, covariates,
                            sp.ones((len(true_phens), 1))])
            (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
            adj_pred_dict['cov_adj'] = sp.dot(Xs, betas)
            pred_r2 = 1 - rss_pd / rss0
            print('Variance explained (Pearson R2) by PRS adjusted for Covariates: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
            res_dict['Cov_adj_pred_r2'] = pred_r2
            pred_r2 = 1 - rss_pd / rss00
            print('Variance explained (Pearson R2) by PRS + Cov: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
            res_dict['Cov_adj_pred_r2+Cov'] = pred_r2
    
            if adjust_for_pcs and 'pcs' in prs_dict and len(prs_dict['pcs']) and 'sex' in prs_dict and len(prs_dict['sex']) > 0:
                pcs = prs_dict['pcs']
                sex = sp.array(prs_dict['sex'])
                sex.shape = (len(sex), 1)
                (betas, rss0, r, s) = linalg.lstsq(
                    sp.hstack([covariates, sex, pcs, sp.ones((len(true_phens), 1))]), true_phens)
                Xs = sp.hstack([pval_derived_effects_prs, covariates,
                                sex, pcs, sp.ones((len(true_phens), 1))])
                (betas, rss_pd, r, s) = linalg.lstsq(Xs, true_phens)
                adj_pred_dict['cov_sex_pc_adj'] = sp.dot(Xs, betas)
                pred_r2 = 1 - rss_pd / rss0
                print('Variance explained (Pearson R2) by PRS adjusted for Cov+PCs+Sex: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
                res_dict['Cov_PC_Sex_adj_pred_r2'] = pred_r2
                pred_r2 = 1 - rss_pd / rss00
                print('Variance explained (Pearson R2) by PRS+Cov+PCs+Sex: %0.4f (%0.6f)' % (pred_r2, (1 - pred_r2) / sp.sqrt(num_individs)))
                res_dict['Cov_PC_Sex_adj_pred_r2+Cov_PC_Sex'] = pred_r2
    
    
        # Now calibration
        y_norm = (true_phens - sp.mean(true_phens)) / sp.std(true_phens)
        denominator = sp.dot(pval_derived_effects_prs.T, pval_derived_effects_prs)
        numerator = sp.dot(pval_derived_effects_prs.T, y_norm)
        regression_slope = (numerator / denominator)[0][0]
        if verbose:
            print('The slope for predictions with weighted effects is: %0.4f'% regression_slope)
    
    
        num_individs = len(prs_dict['pval_derived_effects_prs'])

        # Write PRS out to file.
        if out_file != None:
            write_scores_file(out_file, prs_dict, pval_derived_effects_prs, adj_pred_dict, weights_dict=weights_dict)

    return res_dict



def main(p_dict):
    summary_dict = {}
    non_zero_chromosomes = set()
    verbose = p_dict['debug']

    t0 = time.time()

    summary_dict[0]={'name':'Validation genotype file (prefix):','value':p_dict['gf']}
    summary_dict[0.1]={'name':'Output scores file(s) (prefix):','value':p_dict['out']}

    adjust_for_pcs=False
    adjust_for_covs=False

    if not p_dict['only_score']:
        print('Parsing phenotypes')
        if p_dict['pf'] is None:
            if p_dict['gf'] is not None:
                phen_map = parse_phen_file(p_dict['gf'] + '.fam', 'FAM', verbose=verbose, summary_dict=summary_dict)
            else:
                raise Exception('Validation phenotypes were not found.')
        else:
            phen_map = parse_phen_file(p_dict['pf'], p_dict['pf_format'], verbose=verbose, summary_dict=summary_dict)
        t1 = time.time()
        t = (t1 - t0)
        summary_dict[1.1]={'name':'Individuals with phenotype information:','value':len(phen_map)}
        summary_dict[1.2]={'name':'Running time for parsing phenotypes:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}
    
        if p_dict['cov_file'] != None:
            adjust_for_covs=True
            if verbose:
                print('Parsing additional covariates')
    
            with open(p_dict['cov_file'], 'r') as f:
                num_missing = 0
                for line in f:
                    l = line.split()
                    iid = l[0]
                    if iid in phen_map:
                        covariates = list(map(float, l[1:]))
                        phen_map[iid]['covariates'] = covariates
                    else:
                        num_missing += 1
                if num_missing > 0:
                    summary_dict[2.1]={'name':'Individuals w missing covariate information:','value':num_missing}
                    if verbose:
                        print('Unable to find %d iids in phen file!' % num_missing)
            summary_dict[2]={'name':'Parsed covariates file:','value':p_dict['cov_file']}
    
        if p_dict['pcs_file']:
            adjust_for_pcs=True
            if verbose:
                print('Parsing PCs')
    
            with open(p_dict['pcs_file'], 'r') as f:
                num_missing = 0
                for line in f:
                    l = line.split()
                    iid = l[1]
                    if iid in phen_map:
                        pcs = list(map(float, l[2:]))
                        phen_map[iid]['pcs'] = pcs
                    else:
                        num_missing += 1
                if num_missing > 0:
                    summary_dict[3.1]={'name':'Individuals w missing PCs:','value':num_missing}
                    if verbose:
                        print('Unable to find %d iids in phen file!' % num_missing)
            summary_dict[3]={'name':'Parsed PCs file:','value':p_dict['pcs_file']}
    
        num_individs = len(phen_map)
        assert num_individs > 0, 'No phenotypes were found!'
    else:
        phen_map = None

    t0 = time.time()
    prs_file_is_missing = True
    res_dict = {}
    if p_dict['rf_format'] == 'LDPRED' or p_dict['rf_format']=='ANY':
        weights_file = '%s_LDpred-inf.txt' % (p_dict['rf'])
        
        if os.path.isfile(weights_file):
            print('')
            print('Calculating LDpred-inf risk scores')
            rs_id_map = parse_ldpred_res(weights_file)
            out_file = '%s_LDpred-inf.txt' % (p_dict['out'])
            res_dict['LDpred_inf'] = calc_risk_scores(p_dict['gf'], rs_id_map, phen_map, out_file=out_file, 
                             split_by_chrom=p_dict['split_by_chrom'],
                             adjust_for_pcs=adjust_for_pcs,
                             adjust_for_covariates=adjust_for_covs,
                             only_score=p_dict['only_score'],
                             verbose=verbose, summary_dict=summary_dict)
            prs_file_is_missing = False
        


        for p in p_dict['f']:
            weights_file = '%s_LDpred_p%0.4e.txt' % (p_dict['rf'], p)
            if os.path.isfile(weights_file):
                print('')
                print('Calculating LDpred risk scores using f=%0.3e' % p)
                rs_id_map = parse_ldpred_res(weights_file)
                out_file = '%s_LDpred_p%0.4e.txt' % (p_dict['out'], p)
                method_str = 'LDpred_p%0.4e' % (p)
                res_dict[method_str] = calc_risk_scores(p_dict['gf'], rs_id_map, phen_map, out_file=out_file,
                                                        split_by_chrom=p_dict['split_by_chrom'],
                                                        adjust_for_pcs=adjust_for_pcs,
                                                        adjust_for_covariates=adjust_for_covs,
                                                        only_score=p_dict['only_score'],
                                                        verbose=verbose, summary_dict=summary_dict)
                prs_file_is_missing=False

        # Plot results?

    if p_dict['rf_format'] == 'P+T' or p_dict['rf_format']=='ANY':

        for max_r2 in p_dict['r2']:
            for p_thres in p_dict['p']:
                weights_file = '%s_P+T_r%0.2f_p%0.4e.txt' % (p_dict['rf'], max_r2, p_thres)
                if os.path.isfile(weights_file):
                    print('')
                    print('Calculating P+T risk scores using p-value threshold of %0.3e, and r2 threshold of %0.2f' % (p_thres, max_r2))
                    rs_id_map, non_zero_chromosomes = parse_pt_res(weights_file)
                    out_file = '%s_P+T_p%0.4e.txt' % (p_dict['out'], p_thres)
                    method_str = 'P+T_p%0.4e' % (p_thres)
                    res_dict[method_str] = calc_risk_scores(p_dict['gf'], rs_id_map, phen_map, out_file=out_file,
                                                            split_by_chrom=p_dict['split_by_chrom'],
                                                            non_zero_chromosomes=non_zero_chromosomes, 
                                                            adjust_for_pcs=adjust_for_pcs,
                                                            adjust_for_covariates=adjust_for_covs,
                                                            only_score=p_dict['only_score'],
                                                            verbose=verbose, summary_dict=summary_dict)
                    prs_file_is_missing=False

        # Plot results?
    else:
        raise NotImplementedError(
            'Results file format missing or unknown: %s' % p_dict['rf_format'])
    
    #Identifying the best prediction
    best_pred_r2 = 0
    best_method_str = None
    for method_str in res_dict:
        if len(res_dict[method_str]) and (res_dict[method_str]['pred_r2']) >best_pred_r2:
            best_pred_r2 = res_dict[method_str]['pred_r2']
            best_method_str = method_str
    if best_method_str is not None:
        print('The highest (unadjusted) Pearson R2 was %0.4f, and provided by %s'%(best_pred_r2,best_method_str))
        summary_dict[5]={'name':'Method with highest (unadjusted) Pearson R2:','value':best_method_str}
        summary_dict[5.1]={'name':'Best (unadjusted) Pearson R2:','value':'%0.4f'%best_pred_r2}
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[4]={'name':'Running time for calculating scores:','value':'%d min and %0.2f secs'% (t / 60, t % 60)}

    if prs_file_is_missing:
        print('SNP weights files were not found.  This could be due to a mis-specified --rf flag, or other issues.')
    
    reporting.print_summary(summary_dict,'Scoring Summary')
