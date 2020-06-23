import scipy as sp
from scipy import stats
from scipy import isfinite
import gzip
from ldpred import util
import time
import sys

def parse_sum_stats(h5f, p_dict, bimfile, summary_dict):
    t0 = time.time()
   
    ss_args = {'filename':p_dict['ssf'], 'bimfile':bimfile, 'hdf5_file':h5f, 
               'only_hm3':p_dict['only_hm3'], 'n':p_dict['N'], 
               'debug':p_dict['debug'], 'summary_dict':summary_dict, 
               'match_genomic_pos':p_dict['match_genomic_pos'], 
               'z_from_se':p_dict['z_from_se']}
    if p_dict['ssf_format'] == 'LDPRED':
        """
        CHR     POS     SNP_ID    REF     ALT     REF_FRQ    PVAL    BETA    SE    N
        chr1    1020428    rs6687776    C       T       0.85083    0.0587  -0.0100048507289348    0.0100    8000
        """
        parse_sum_stats_custom(ch='CHR', pos='POS', A1='REF', A2='ALT', reffreq='REF_FRQ', info='info',
                    rs='SNP_ID', pval='PVAL', eff='BETA', se='SE', ncol='N', eff_type='LINREG', **ss_args)   
    elif p_dict['ssf_format'] == 'STANDARD':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        """
        chr     pos     ref     alt     reffrq  info    rs       pval    effalt
        chr1    1020428 C       T       0.85083 0.98732 rs6687776    0.0587  -0.0100048507289348
        """        
        parse_sum_stats_custom(ch='chr', pos='pos', A1='ref', A2='alt', reffreq='reffrq', info='info',
                    rs='rs', pval='pval', eff='effalt', eff_type='LINREG', **ss_args)      
    elif p_dict['ssf_format'] == 'PGC_OLD':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        """
        hg19chrc        snpid   a1      a2      bp      info    or      se      p       ngt
        chr1    rs4951859       C       G       729679  0.631   0.97853 0.0173  0.2083  0    
        """
        parse_sum_stats_custom(ch='hg19chrc', pos='bp', A1='a1', A2='a2', info='info', rs='snpid', pval='p', eff='or', se='se',
                    eff_type='OR', **ss_args)
    elif p_dict['ssf_format'] == 'BASIC':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        """
        hg19chrc    snpid    a1    a2    bp    or    p       
        chr1    rs4951859    C    G    729679    0.97853    0.2083      
        """
        parse_sum_stats_custom(ch='hg19chrc', pos='bp', A1='a1', A2='a2', info='info', rs='snpid', pval='p',
                    eff='or', eff_type='OR', **ss_args)
    elif p_dict['ssf_format'] == 'PGC':
        """    
        CHR    SNP    BP    A1    A2    FRQ_A_30232    FRQ_U_40578    INFO    OR    SE    P    ngt    Direction    HetISqt    HetChiSq    HetDf    HetPVa
        ...    
        """
        parse_sum_stats_custom(ch='CHR', pos='BP', A1='A1', A2='A2', reffreq='Freq.Hapmap.Ceu',
                    case_freq='FRQ_A_34129', control_freq='FRQ_U_45512', case_n='Nca',
                    control_n='Nco', info='INFO', rs='SNP', pval='P', eff='OR', se='SE',
                    eff_type='OR', **ss_args)
    elif p_dict['ssf_format'] == 'GIANT':
        """
        MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
        rs10 a c 0.0333 0.8826 78380    
        """
        parse_sum_stats_custom(A1='Allele1', A2='Allele2', reffreq='Freq.Allele1.HapMapCEU', rs='MarkerName',
                    pval='p', eff='b', ncol='N', eff_type='LINREG', **ss_args)
    elif p_dict['ssf_format'] == 'GIANT2':
        """
        MarkerName A1 A2 Freq.Hapmap.Ceu BETA SE.2gc P.2gc N
        rs4747841 a g 0.55 0.0025 0.0061 0.68 60558.2
        """
        parse_sum_stats_custom(ncol='N', A1='A1', A2='A2', reffreq='Freq.Hapmap.Ceu', rs='MarkerName',
                    pval='P.2gc', eff='BETA', se='SE.2gc', eff_type='LINREG', **ss_args)
    elif p_dict['ssf_format'] == 'CUSTOM':
        parse_sum_stats_custom(ch=p_dict['chr'],
                    A1=p_dict['A1'], A2=p_dict['A2'], reffreq=p_dict['reffreq'], info=p_dict['info'],
                    rs=p_dict['rs'], pval=p_dict['pval'], eff=p_dict['eff'], ncol=p_dict['ncol'],
                    pos=p_dict['pos'], eff_type=p_dict['eff_type'], case_freq=p_dict['case_freq'], 
                    control_freq=p_dict['control_freq'], case_n=p_dict['case_n'], se=p_dict['se'],
                    control_n=p_dict['control_n'], **ss_args)
    else:
        raise Exception('Unknown Summary Statistics Format.')
    t1 = time.time()
    t = (t1 - t0)
    summary_dict[14]={'name':'Run time for parsing summary stats:','value': '%d min and %0.2f sec'%(t / 60, t % 60)}


def parse_sample_size(line_dict,n,ncol,case_n,control_n,header_dict,):
    if n is None:
        if ncol not in header_dict:
            case_N = float(line_dict[header_dict[case_n]])
            control_N = float(line_dict[header_dict[control_n]])
            N = case_N + control_N
        else:
            N = float(line_dict[header_dict[ncol]])
    else:
        N = n
    return N


def parse_nucleotides(line_dict, header_dict, A1, A2):
    return [line_dict[header_dict[A1]].upper(), line_dict[header_dict[A2]].upper()]
                
def parse_freq(line_dict, header_dict, reffreq, control_freq, case_freq, control_n, case_n):
    
    if reffreq is not None and reffreq in header_dict:
        if line_dict[header_dict[reffreq]] == '.' or line_dict[header_dict[reffreq]] == 'NA':
            return -1
        else:
            return float(line_dict[header_dict[reffreq]])
    elif (case_freq is not None and control_freq is not None 
          and case_freq in header_dict and control_freq in header_dict):
        if (case_n is not None and control_n is not None 
              and case_n in header_dict and control_n in header_dict) :
            if (line_dict[header_dict[control_n]] == '.' or line_dict[header_dict[control_n]] == 'NA' 
                or line_dict[header_dict[case_n]] == '.' or line_dict[header_dict[case_n]] == 'NA' 
                or line_dict[header_dict[control_freq]] == '.' or line_dict[header_dict[control_freq]] == 'NA' 
                or line_dict[header_dict[case_freq]] == '.' or line_dict[header_dict[case_freq]] == 'NA'):
                return -1
            else:
                case_N = float(line_dict[header_dict[case_n]])
                control_N = float(line_dict[header_dict[control_n]])
                tot_N = case_N + control_N
                a_scalar = case_N / float(tot_N)
                u_scalar = control_N / float(tot_N)
                freq = float(line_dict[header_dict[case_freq]]) * a_scalar + float(line_dict[header_dict[control_freq]]) * u_scalar
                return freq
        else:
            if (line_dict[header_dict[case_freq]] == '.' or line_dict[header_dict[case_freq]] == 'NA' 
                or line_dict[header_dict[control_freq]] == '.' or line_dict[header_dict[control_freq]] == 'NA'):
                return -1
            else:
                return (float(line_dict[header_dict[case_freq]]) + float(line_dict[header_dict[control_freq]]))/2.0
    else:  
        return -1

def get_raw_beta(beta_read, eff_type):
    if eff_type=='OR':
        raw_beta = sp.log(beta_read)
    elif eff_type=='LINREG' or eff_type=='LOGOR' or eff_type=='BLUP':
        raw_beta = beta_read
    else: 
        raise Exception('Unknown effect type')
    return raw_beta

def get_beta_from_se(beta_read, se_read, eff_type, raw_beta, N):
    if se_read==0:
        return 
    if eff_type=='LINREG' or eff_type=='LOGOR':
        abs_beta = sp.absolute(beta_read)/se_read
    elif eff_type=='OR':
        abs_beta = sp.absolute(1-beta_read)/se_read
    else: 
        raise Exception('Unknown effect type')      
    return sp.sign(raw_beta) * abs_beta/ sp.sqrt(N)


def get_beta(pval_read, raw_beta, beta_read, line_dict, header_dict, se, 
             z_from_se, N, eff_type, se_inferred_zscores):
    if eff_type=='BLUP':
        return raw_beta
    if z_from_se:
        assert se in header_dict, "SE was not specified in summary statistics provided, which is necessary when the 'z_from_se' flag is used."
        se_read = float(line_dict[header_dict[se]])
        if se_read==0 or not isfinite(se_read):
            return None
        se_inferred_zscores += 1
        return get_beta_from_se(beta_read, se_read, eff_type, raw_beta, N)
    else:
        if pval_read==0 or not isfinite(stats.norm.ppf(pval_read)):
            #Attempt to Parse SEs to infer Z-score
            if not se in header_dict:
                return None
            else:
                se_read = float(line_dict[header_dict[se]])
                if se_read==0 or not isfinite(se_read):
                    return None

                se_inferred_zscores += 1
                return get_beta_from_se(beta_read, se_read, eff_type, raw_beta, N)
        else:               
            return -1 * sp.sign(raw_beta) * stats.norm.ppf(pval_read / 2.0)/ sp.sqrt(N) 



def parse_sum_stats_custom(filename=None, bimfile=None, only_hm3=False, hdf5_file=None, n=None, ch=None, pos=None,
                    A1=None, A2=None, reffreq=None, case_freq=None, control_freq=None, case_n=None,
                    control_n=None, info=None, rs=None, pval=None, eff=None, ncol=None, se=None,
                    eff_type='OR', match_genomic_pos=False, debug=False, 
                    z_from_se=False, summary_dict = None):
    # Check required fields are here
    assert not A2 is None, 'Require header for non-effective allele'
    assert not A1 is None, 'Require header for effective allele'
    assert not rs is None, 'Require header for RS ID'
    assert not eff is None, 'Require header for Statistics'
    assert not pval is None, 'Require header for pval'
    assert not ncol is None or not n is None or (control_n is not None and case_n is not None), 'Require either N or NCOL information'

    if ch is None:
        assert not bimfile is None, 'Require bimfile when chromosome header not provided'
        print("Chromosome Header not provided, will use info from bim file")
    if pos is None:
        assert not bimfile is None, 'Require bimfile when position header not provided'
        print("Position Header not provided, will use info from bim file")

    num_lines = util.count_lines(filename)
    snps_pos_map = {}
    if only_hm3:
        if debug:
            print('Loading HapMap3 SNPs')
        hm3_sids = util.load_hapmap_SNPs()
  
    if bimfile is not None:
        valid_sids = set()
        if debug:
            print('Parsing bim file: %s' % bimfile)
        
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                chrom = util.get_chrom_num(l[0])
                if chrom not in util.ok_chromosomes:
                    continue
                sid = l[1]
                if only_hm3:
                    if sid in hm3_sids:
                        valid_sids.add(sid)
                        snps_pos_map[sid] = {'pos':int(l[3]), 'chrom':chrom}
                else:
                    valid_sids.add(sid)
                    snps_pos_map[sid] = {'pos':int(l[3]), 'chrom':chrom}

        if len(valid_sids)==0:
            raise Exception('Unable to parse BIM file')
    else:
        raise Exception('BIM file missing. Please check genotype paths provided.')
        
    invalid_chr = 0
    invalid_pos = 0
    invalid_p = 0
    invalid_beta = 0
    se_inferred_zscores = 0
    chrom_dict = {}
    opener = open
    if util.is_gz(filename):
        opener = gzip.open
    print('Parsing summary statistics file: %s' % filename)
    with opener(filename) as f:
        header = f.readline()
        if util.is_gz(filename):
            header = header.decode('utf-8')
        if debug:
            print('File header:')
            print(header)
        header_dict = {}
        columns = (header.strip()).split()
        index = 0
        for col in columns:
            header_dict[col] = index
            index += 1
                
        assert ch is None or ch in header_dict, 'Chromosome header cannot be found in summary statistic file'
        assert A2 in header_dict, 'Non-effective allele column cannot be found in summary statistic file'
        assert A1 in header_dict, 'Effective allele column cannot be found in summary statistic file'
        assert eff in header_dict, 'Effect size column not found in summary statistic file'
        assert rs in header_dict, 'SNP ID column not found in summary statistic file'
        assert pos is None or pos in header_dict, 'Position column not found in summary statistic file'
        assert pval in header_dict, 'P Value column not found in summary statistic file'
        assert not n is None or ncol in header_dict or (control_n in header_dict and case_n in header_dict), 'Sample size column not found in summary statistic ' \
                                                     'file and N not provided'
                                                     
        if z_from_se:
            assert se is not None, 'SE column must be specified to infer z-scores from SEs'
            assert se in header_dict, 'SE column not found in summary stats file, this is required to infer z-scores from SEs'

        # header_dict now contains the header column name for each corresponding input
        bad_chromosomes = set()
        line_i = 1
        for line in f:
            line_i +=1
            if line_i%1000==0 and num_lines>0:
                sys.stdout.write('\r%0.2f%%' % (100.0 * (float(line_i) / (num_lines))))
                sys.stdout.flush()            
            if util.is_gz(filename):
                line = line.decode('utf-8')
            l = (line.strip()).split()
            # get the SNP ID first
            sid = l[header_dict[rs]]
            # check the SNP ID
            if sid in valid_sids:
                # Get the chromosome information
                chrom = 0
                if not ch is None and ch in header_dict:
                    chrom = util.get_chrom_num(l[header_dict[ch]])
                    # Check if the chromosome of the SNP is correct
                    if not chrom == snps_pos_map[sid]['chrom']:
                        invalid_chr += 1
                        continue
                else:
                    chrom = snps_pos_map[sid]['chrom']
                
    
                #Parse position
                pos_read = 0
                if not pos is None and pos in header_dict:
                    pos_read = int(l[header_dict[pos]])
                    if not pos_read == snps_pos_map[sid]['pos']:
                        invalid_pos += 1
                        if match_genomic_pos:
                            continue
                else:
                    pos_read = snps_pos_map[sid]['pos']

                #Get the sample size
                N = parse_sample_size(l,n,ncol,case_n,control_n,header_dict)

                #Parse raw beta      
                beta_read = float(l[header_dict[eff]])
                if not isfinite(beta_read):
                    invalid_beta += 1
                    continue

                raw_beta = get_raw_beta(beta_read, eff_type)

                #Parse p-value and effect size
                pval_read = float(l[header_dict[pval]])
                if pval_read==0 or not isfinite(stats.norm.ppf(pval_read)):
                    invalid_p += 1


                beta = get_beta(pval_read, raw_beta, beta_read, l, header_dict, se, 
                                z_from_se, N, eff_type, se_inferred_zscores)
                    
                if beta==None:
                    continue
                
                #All necessary information was found, so we should store things
                
                if not chrom in chrom_dict:
                    chrom_dict[chrom] = {'ps':[], 'infos':[], 'freqs':[],
                                         'nts': [], 'sids': [], 'positions': [], 'ns':[],
                                          'log_odds':[], 'betas':[]}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos_read)
                chrom_dict[chrom]['ps'].append(pval_read)
                chrom_dict[chrom]['ns'].append(N)
                chrom_dict[chrom]['log_odds'].append(raw_beta)
                chrom_dict[chrom]['betas'].append(beta)
                
                # Get the INFO score
                info_sc = -1
                if info is not None and info in header_dict:
                    info_sc = float(l[header_dict[info]])
                chrom_dict[chrom]['infos'].append(info_sc)
                
                #Parse nucleotides
                nt = parse_nucleotides(l, header_dict, A1, A2)
                chrom_dict[chrom]['nts'].append(nt)
                
                # Parse the frequency
                freq = parse_freq(l, header_dict, reffreq, control_freq, case_freq, control_n, case_n)  
                chrom_dict[chrom]['freqs'].append(freq)

        if len(bad_chromosomes) > 0:
            if debug:
                print('Ignored chromosomes: %s' % (','.join(list(bad_chromosomes))))
                print('Please note that only data on chromosomes 1-23, and X are parsed.')

    if num_lines>0:
        sys.stdout.write('\r%0.2f%%\n' % (100.0))
        sys.stdout.flush()            
    print('SS file loaded, now sorting and storing in HDF5 file.')
    assert not 'sum_stats' in hdf5_file, 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    num_non_finite = 0
    for chrom in chrom_dict:
        if debug:
            print ('%d SNPs on chromosome %s' % (len(chrom_dict[chrom]['positions']), chrom))
        assert len(chrom_dict[chrom]['positions'])==len(chrom_dict[chrom]['betas'])==len(chrom_dict[chrom]['ps'])==len(chrom_dict[chrom]['nts']), 'Problems with parsing summary stats'
        sl = list(zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['infos'],
                 chrom_dict[chrom]['freqs'], chrom_dict[chrom]['ps'], chrom_dict[chrom]['ns']))
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        infos = []
        freqs = []
        ns = []
        prev_pos = -1
        for pos, sid, nt, beta, lo, info, frq, p, num_ind in sl:
            if pos == prev_pos:
                if debug:
                    print('duplicated position %d' % pos)
                continue
            else:
                prev_pos = pos
            if not sp.isfinite(beta):
                num_non_finite += 1
                continue
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            infos.append(info)
            freqs.append(frq)
            ns.append(num_ind)
        nts = sp.array(nts, dtype=util.nts_dtype)
        sids = sp.array(sids, dtype=util.sids_dtype)
        if debug:
            if not num_non_finite == 0:
                print('%d SNPs have non-finite statistics on chromosome %s' % (num_non_finite, chrom))
            print ('Still %d SNPs on chromosome %s' % (len(ps), chrom))
        g = ssg.create_group('chrom_%s' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('freqs', data=freqs)
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps += len(log_odds)
        g.create_dataset('infos', data=infos)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        g.create_dataset('ns', data=ns)
        hdf5_file.flush()
    
    summary_dict[3.09]={'name':'dash', 'value':'Summary statistics'}
    summary_dict[3.1]={'name':'Num SNPs parsed from sum stats file','value':num_snps}
    if invalid_p>0 or debug:
        summary_dict[3.2]={'name':'Num invalid P-values in sum stats','value':invalid_p}
    if invalid_beta>0 or debug:
        summary_dict[3.21]={'name':'Num invalid betas in sum stats','value':invalid_beta}
    if se_inferred_zscores>0 or debug:
        summary_dict[3.22]={'name':'Num z-scores inferred from SEs and effects','value':se_inferred_zscores}
    if invalid_chr>0 or debug:
            summary_dict[3.4]={'name':'SNPs w non-matching chromosomes excluded','value':invalid_chr}
    if invalid_pos>0 or debug:
        if match_genomic_pos:
            summary_dict[3.3]={'name':'SNPs w non-matching positions excluded','value':invalid_pos}
        else:
            summary_dict[3.3]={'name':'SNPs w non-matching positions (not excluded)','value':invalid_pos}

