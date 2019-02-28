import scipy as sp
from scipy import stats
from scipy import isinf
import re
import gzip
from ldpred import util


def parse_sum_stats(h5f, p_dict, bimfile):
    ss_args = {'filename':p_dict['ssf'], 'bimfile':bimfile, 'hdf5_file':h5f, 'only_hm3':p_dict['only_hm3'], 'n':p_dict['N'], 'debug':p_dict['debug']}
    if p_dict['ssf_format'] == 'STANDARD':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        """
        chr     pos     ref     alt     reffrq  info    rs       pval    effalt
        chr1    1020428 C       T       0.85083 0.98732 rs6687776    0.0587  -0.0100048507289348
        """        
        parse_sum_stats_custom(ch='chr', pos='pos', A1='ref', A2='alt', reffreq='reffrq', info='info',
                    rs='rs', pval='pval', eff='effalt', input_is_beta=True, **ss_args)      
    elif p_dict['ssf_format'] == 'PGC_OLD':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        """
        hg19chrc        snpid   a1      a2      bp      info    or      se      p       ngt
        chr1    rs4951859       C       G       729679  0.631   0.97853 0.0173  0.2083  0    
        """
        parse_sum_stats_custom(ch='hg19chrc', pos='bp', A1='a1', A2='a2', info='info', rs='snpid', pval='p', eff='or',
                    input_is_beta=False, **ss_args)
    elif p_dict['ssf_format'] == 'BASIC':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        """
        hg19chrc    snpid    a1    a2    bp    or    p       
        chr1    rs4951859    C    G    729679    0.97853    0.2083      
        """
        parse_sum_stats_custom(ch='hg19chrc', pos='bp', A1='a1', A2='a2', info='info', rs='snpid', pval='p',
                    eff='or', input_is_beta=False, **ss_args)
    elif p_dict['ssf_format'] == 'PGC':
        """    
        CHR    SNP    BP    A1    A2    FRQ_A_30232    FRQ_U_40578    INFO    OR    SE    P    ngt    Direction    HetISqt    HetChiSq    HetDf    HetPVa
        ...    
        """
        parse_sum_stats_custom(ch='CHR', pos='BP', A1='A1', A2='A2', reffreq='Freq.Hapmap.Ceu',
                    case_freq='FRQ_A_30232', control_freq='FRQ_U_40578', case_n='FRQ_A_30232',
                    control_n='FRQ_U_40578', info='INFO', rs='SNP', pval='P', eff='OR',
                    ncol='ngt', input_is_beta=False, **ss_args)
    elif p_dict['ssf_format'] == 'GIANT':
        """
        MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
        rs10 a c 0.0333 0.8826 78380    
        """
        parse_sum_stats_custom(A1='Allele1', A2='Allele2', reffreq='Freq.Allele1.HapMapCEU', rs='MarkerName',
                    pval='p', eff='b', ncol='N', input_is_beta=True, **ss_args)
    elif p_dict['ssf_format'] == 'GIANT2':
        """
        MarkerName A1 A2 Freq.Hapmap.Ceu BETA SE.2gc P.2gc N
        rs4747841 a g 0.55 0.0025 0.0061 0.68 60558.2
        """
        parse_sum_stats_custom(ncol='N', A1='A1', A2='A2', reffreq='Freq.Hapmap.Ceu', rs='MarkerName',
                    pval='P.2gc', eff='BETA', input_is_beta=True, **ss_args)
    elif p_dict['ssf_format'] == 'CUSTOM':
        parse_sum_stats_custom(ch=p_dict['chr'],
                    A1=p_dict['A1'], A2=p_dict['A2'], reffreq=p_dict['reffreq'], info=p_dict['info'],
                    rs=p_dict['rs'], pval=p_dict['pval'], eff=p_dict['eff'], ncol=p_dict['ncol'],
                    pos=p_dict['pos'], input_is_beta=p_dict['beta'], **ss_args)
    else:
        raise Exception('Unknown Summary Statistics Format.')


def is_gz(name):
    return name.lower().endswith(('.gz', '.gzip'))


def parse_sum_stats_custom(filename=None, bimfile=None, only_hm3=False, hdf5_file=None, n=None, ch=None, pos=None,
                    A1=None, A2=None, reffreq=None, case_freq=None, control_freq=None, case_n=None,
                    control_n=None, info=None, rs=None, pval=None, eff=None, ncol=None,
                    input_is_beta=False, debug=False):
    # Check required fields are here
    assert not A2 is None, 'Require header for non-effective allele'
    assert not A1 is None, 'Require header for effective allele'
    assert not rs is None, 'Require header for RS ID'
    assert not eff is None, 'Require header for Statistics'
    assert not pval is None, 'Require header for pval'
    assert not ncol is None or not n is None, 'Require either N or NCOL information'

    if ch is None:
        assert not bimfile is None, 'Require bimfile when chromosome header not provided'
        print("Chromosome Header not provided, will use info from bim file")
    if pos is None:
        assert not bimfile is None, 'Require bimfile when position header not provided'
        print("Position Header not provided, will use info from bim file")

    snps_pos_map = {}
    if only_hm3:
        hm3_sids = util.load_hapmap_SNPs()
    
    if bimfile is not None:
        valid_sids = set()
        if debug:
            print('Parsing bim file: %s' % bimfile)
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                # Bim file format is CHR SNP BP
                sid = l[1]
                if only_hm3:
                    if sid in hm3_sids:
                        valid_sids.add()
                        snps_pos_map[sid] = {'pos':int(l[3]), 'chrom':l[0]}
                else:
                    valid_sids.add(sid)
                    snps_pos_map[sid] = {'pos':int(l[3]), 'chrom':l[0]}

        if len(valid_sids)==0:
            raise Exception('Unable to parse BIM file')
    else:
        raise Exception('Unable to load LD reference or validation genotypes bim file')
        
    chr_filter = 0
    pos_filter = 0
    invalid_p = 0
    chrom_dict = {}
    opener = open
    if is_gz(filename):
        opener = gzip.open
    print('Parsing summary statistics file: %s' % filename)
    with opener(filename) as f:
        header = f.readline()
        if is_gz(filename):
            header = header.decode('utf-8')
        if debug:
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
        assert not n is None or ncol in header_dict, 'Sample size column not found in summary statistic ' \
                                                     'file and N not provided'
        # header_dict now contains the header column name for each corresponding input
        bad_chromosomes = set()
        for line in f:
            if is_gz(filename):
                line = line.decode('utf-8')
            l = (line.strip()).split()
            # get the SNP ID first
            sid = l[header_dict[rs]]
            # check the SNP ID
            if sid in valid_sids:
                # Get the chromosome information
                chrom = 0
                if not chr is None and chr in header_dict:
                    chrom = l[header_dict[ch]]
                    chrom = re.sub("chr", "", chrom)
                    if not chrom == snps_pos_map[sid]['chrom']:
                        chr_filter += 1
                else:
                    chrom = snps_pos_map[sid]['chrom']
                if not chrom in util.ok_chromosomes:
                    bad_chromosomes.add(chrom)
                    continue
                # Check if the chromosome of the SNP is correct
    
                pos_read = 0
                if not pos is None and pos in header_dict:
                    pos_read = int(l[header_dict[pos]])
                    if not pos_read == snps_pos_map[sid]['pos']:
                        pos_filter += 1
                        continue
                else:
                    pos_read = snps_pos_map[sid]['pos']

                if not chrom in chrom_dict:
                    chrom_dict[chrom] = {'ps':[], 'log_odds':[], 'infos':[], 'freqs':[],
                             'betas':[], 'nts': [], 'sids': [], 'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos_read)
                # Check the frequency
                if reffreq is not None and reffreq in header_dict:
                    if l[header_dict[reffreq]] == '.' or l[header_dict[reffreq]] == 'NA':
                        chrom_dict[chrom]['freqs'].append(-1)
                    else:
                        chrom_dict[chrom]['freqs'].append(float(l[header_dict[reffreq]]))
                elif (case_n is not None and control_n is not None 
                      and case_n in header_dict and control_n in header_dict 
                      and case_freq is not None and control_freq is not None 
                      and case_freq in header_dict and control_freq in header_dict) :
                    if (l[header_dict[control_n]] == '.' or l[header_dict[control_n]] == 'NA' 
                        or l[header_dict[case_n]] == '.' or l[header_dict[case_n]] == 'NA' 
                        or l[header_dict[control_freq]] == '.' or l[header_dict[control_freq]] == 'NA' 
                        or l[header_dict[case_freq]] == '.' or l[header_dict[case_freq]] == 'NA'):
                        chrom_dict[chrom]['freqs'].append(-1)
                    else:
                        case_N = float(l[header_dict[case_n]])
                        control_N = float(l[header_dict[control_n]])
                        N = case_N + control_N
                        a_scalar = case_N / N
                        u_scalar = control_N / N
                        freq = float(l[header_dict[case_freq]]) * a_scalar + float(l[header_dict[control_freq]]) * u_scalar
                        chrom_dict[chrom]['freqs'].append(freq)
                else:  
                    chrom_dict[chrom]['freqs'].append(-1)
                # Get the INFO score
                info_sc = -1
                if info is not None and info in header_dict:
                    info_sc = float(l[header_dict[info]])
                chrom_dict[chrom]['infos'].append(info_sc)
                pval_read = float(l[header_dict[pval]])
                chrom_dict[chrom]['ps'].append(pval_read)
                if isinf(stats.norm.ppf(pval_read)):
                    invalid_p += 1
                    continue

                nt = [l[header_dict[A1]].upper(), l[header_dict[A2]].upper()]
                chrom_dict[chrom]['nts'].append(nt)
                raw_beta = float(l[header_dict[eff]])
                if not input_is_beta:
                    raw_beta = sp.log(raw_beta)
                    chrom_dict[chrom]['log_odds'].append(raw_beta)
                    beta = sp.sign(raw_beta) * stats.norm.ppf(pval_read / 2.0)
                    if n is None:
                        chrom_dict[chrom]['betas'].append(beta / sp.sqrt(int(header_dict[ncol])))
                    else:
                        chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
                else:
                    beta = sp.sign(raw_beta) * stats.norm.ppf(pval_read / 2.0)
                    if n is None:
                        chrom_dict[chrom]['log_odds'].append(beta / sp.sqrt(int(header_dict[ncol])))
                        chrom_dict[chrom]['betas'].append(beta / sp.sqrt(int(header_dict[ncol])))
                    else:
                        chrom_dict[chrom]['log_odds'].append(beta / sp.sqrt(n))
                        chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))

        if len(bad_chromosomes) > 0:
            print('Ignored chromosomes: %s' % (','.join(list(bad_chromosomes))))
            print('Please note that only data on chromosomes 1-23, and X are parsed.')

    print('SS file loaded, now sorting and storing in HDF5 file.')
    assert not 'sum_stats' in hdf5_file, 'Something is wrong with HDF5 file?'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    num_non_finite = 0
    for chrom in chrom_dict:
        if debug:
            print ('%d SNPs on chromosome %s' % (len(chrom_dict[chrom]['positions']), chrom))
        sl = list(zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['infos'],
                 chrom_dict[chrom]['freqs'], chrom_dict[chrom]['ps']))
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        infos = []
        freqs = []
        prev_pos = -1
        for pos, sid, nt, beta, lo, info, frq, p in sl:
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
        hdf5_file.flush()
    print('%d SNPs excluded due to invalid chromosome ID.' % chr_filter)
    print('%d SNPs excluded due to invalid chromosome position' % pos_filter)
    print('%d SNPs excluded due to invalid P value' % invalid_p)
    print('%d SNPs parsed from summary statistics file.' % num_snps)

