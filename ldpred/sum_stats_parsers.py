import scipy as sp
from scipy import stats
from scipy import isinf
import random
import re
from ldpred import util


sids_dtype = "S30" #30 character limit
nts_dtype = "S1"


def parse_sum_stats(h5f, p_dict, bimfile):
    if p_dict['ssf_format'] == 'STANDARD':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        parse_sum_stats_standard(
            filename=p_dict['ssf'], bimfile=bimfile, hdf5_file=h5f, n=p_dict['N'], debug=p_dict['debug'])
    elif p_dict['ssf_format'] == 'PGC_OLD':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        parse_sum_stats_pgc_small(
            filename=p_dict['ssf'], bimfile=bimfile, hdf5_file=h5f, n=p_dict['N'], debug=p_dict['debug'])
    elif p_dict['ssf_format'] == 'BASIC':
        if p_dict['N'] is None: 
            raise Exception('This summary statistics format requires summary statistics sample size to be given, i.e. set --N flag.')        
        parse_sum_stats_basic(
            filename=p_dict['ssf'], bimfile=bimfile, hdf5_file=h5f, n=p_dict['N'], debug=p_dict['debug'])
    elif p_dict['ssf_format'] == 'PGC':
        parse_sum_stats_pgc(
            filename=p_dict['ssf'], bimfile=bimfile, hdf5_file=h5f, debug=p_dict['debug'])
    elif p_dict['ssf_format'] == 'GIANT':
        parse_sum_stats_giant(
            filename=p_dict['ssf'], bimfile=bimfile, hdf5_file=h5f, debug=p_dict['debug'])
    elif p_dict['ssf_format'] == 'GIANT2':
        parse_sum_stats_giant2(
            filename=p_dict['ssf'], bimfile=bimfile, hdf5_file=h5f, debug=p_dict['debug'])
    elif p_dict['ssf_format'] == 'CUSTOM':
        parse_sum_stats_custom(filename=p_dict['ssf'], bimfile=bimfile, hdf5_file=h5f, n=p_dict['N'], chr=p_dict['chr'],
                    A1=p_dict['A1'], A2=p_dict['A2'], reffreq=p_dict['reffreq'], info=p_dict['info'],
                    rs=p_dict['rs'], pval=p_dict['pval'], eff=p_dict['eff'], ncol=p_dict['ncol'],
                    pos=p_dict['pos'], input_is_beta=p_dict['beta'], debug=p_dict['debug'])


def parse_sum_stats_custom(filename=None,
                    bimfile=None,
                    hdf5_file=None,
                    n=None,
                    chr=None,
                    pos=None,
                    A1=None,
                    A2=None,
                    reffreq=None,
                    info=None,
                    rs=None,
                    pval=None, 
                    eff=None,
                    ncol=None,
                    input_is_beta=False, 
                    debug=False):
    # Check required fields are here
    assert not A2 is None, 'Require header for non-effective allele'
    assert not A1 is None, 'Require header for effective allele'
    assert not rs is None, 'Require header for RS ID'
    assert not eff is None, 'Require header for Statistics'
    assert not pval is None, 'Require header for pval'

    assert not ncol is None or not n is None, 'Require either N or NCOL information'


    if chr is None:
        assert not bimfile is None, 'Require bimfile when chromosome header not provided'
        print("Chromosome Header not provided, will use info from bim file")
    if pos is None:
        assert not bimfile is None, 'Require bimfile when position header not provided'
        print("Position Header not provided, will use info from bim file")

    snps_pos_map = {}
    if bimfile is not None:
        valid_sids = set()
        if debug:
            print('Parsing bim file: %s' % bimfile)
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                # Bim file format is CHR SNP BP
                valid_sids.add(l[1])
                snps_pos_map[l[1]] = {'pos':int(l[3]), 'chrom':l[0]}
    chr_filter = 0
    pos_filter = 0
    invalid_p = 0
    chrom_dict = {}
    print('Parsing summary statistics file: %s' % filename)
    with open(filename) as f:

        header = next(f)
        if debug:
            print(header)
        header_dict={}
        columns = (header.strip()).split()
        index = 0
        for col in columns:
            header_dict[col] = index
            index+=1
        assert chr is None or chr in header_dict, 'Chromosome header cannot be found in summary statistic file'
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
            l = (line.strip()).split()
            # get the SNP ID first
            sid = l[header_dict[rs]]
            # Get the chromosome information
            chrom = 0
            if not chr is None and chr in header_dict:
                chrom = l[header_dict[chr]]
                chrom = re.sub("chr", "", chrom)
                if not chrom == snps_pos_map[sid]['chrom']:
                    chr_filter += 1
                    continue
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
            # Now check the SNP ID
            if sid in valid_sids:
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
                        chrom_dict[chrom]['betas'].append(beta/ sp.sqrt(int(header_dict[ncol])))
                    else:
                        chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
                else:
                    beta = sp.sign(raw_beta) * stats.norm.ppf(pval_read / 2.0)
                    if n is None:
                        chrom_dict[chrom]['log_odds'].append(beta/ sp.sqrt(int(header_dict[ncol])))
                        chrom_dict[chrom]['betas'].append(beta/ sp.sqrt(int(header_dict[ncol])))
                    else:
                        chrom_dict[chrom]['log_odds'].append(beta / sp.sqrt(n))
                        chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))


        if len(bad_chromosomes) > 0:
            print('Ignored chromosomes: %s'%(','.join(list(bad_chromosomes))))
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
                num_non_finite+=1
                continue
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            infos.append(info)
            freqs.append(frq)
        nts = sp.array(nts,dtype=nts_dtype)
        sids = sp.array(sids,dtype=sids_dtype)
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


def parse_sum_stats_standard(filename=None,
                             bimfile=None,
                             hdf5_file=None,
                             n=None,
                             debug=False):
    """
    Input format:

     chr     pos     ref     alt     reffrq  info    rs       pval    effalt
    chr1    1020428 C       T       0.85083 0.98732 rs6687776    0.0587  -0.0100048507289348
    chr1    1020496 G       A       0.85073 0.98751 rs6678318    0.1287  -0.00826075392985992

    """

    if bimfile is not None:
        valid_sids = set()
        if debug:
            print('Parsing bim file: %s' % bimfile)
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])

    chrom_dict = {}

    print('Parsing summary statistics file: %s' % filename)
    with open(filename) as f:
        if debug:
            print(next(f))
        missing_chromosomes = set()
        for line in f:
            l = (line.strip()).split()
            chrom = l[0][3:]
            if not chrom in util.valid_chromosomes:
                missing_chromosomes.add(chrom)
                continue
            pos = int(l[1])
            sid = l[6]
            if sid in valid_sids:
                if not chrom in chrom_dict:
                    chrom_dict[chrom] = {'ps': [], 'log_odds': [], 'infos': [], 'freqs': [],
                                         'betas': [], 'nts': [], 'sids': [], 'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                chrom_dict[chrom]['freqs'].append(float(l[4]))
                info = float(l[5])
                chrom_dict[chrom]['infos'].append(info)
                pval = float(l[7])
                chrom_dict[chrom]['ps'].append(pval)
                nt = [l[2], l[3]]
                chrom_dict[chrom]['nts'].append(nt)
                raw_beta = float(l[8])
                chrom_dict[chrom]['log_odds'].append(raw_beta)
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval / 2.0)

                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))

        util.check_chromosomes(missing_chromosomes)

    print('SS file loaded, now sorting and storing in HDF5 file.')
    assert not 'sum_stats' in hdf5_file, 'Something is wrong with HDF5 file.  Summary stats already found.'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict:
        if debug:
            print('%d SNPs on chromosome %s' % (len(chrom_dict[chrom]['positions']), chrom))
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
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
            infos.append(info)
            freqs.append(frq)
        nts = sp.array(nts,dtype=nts_dtype)
        sids = sp.array(sids,dtype=sids_dtype)
        if debug:
            print('Still %d SNPs on chromosome %s' % (len(ps), chrom))
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
    print('%d SNPs parsed from summary statistics file.' % num_snps)


def parse_sum_stats_giant(filename=None,
                          bimfile=None,
                          hdf5_file=None,
                          debug=False):
    """
    Input format:

    MarkerName      Allele1 Allele2 Freq.Allele1.HapMapCEU  b       SE      p       N
    MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
    rs10 a c 0.0333 0.8826 78380
    rs1000000 a g 0.3667 0.1858 133822

    """

    snps_pos_map = {}
    assert bimfile is not None, 'BIM file is required'
    valid_sids = set()
    if debug:
        print('Parsing bim file: %s' % bimfile)
    with open(bimfile) as f:
        for line in f:
            l = line.split()
            valid_sids.add(l[1])
            snps_pos_map[l[1]] = {'pos': int(l[3]), 'chrom': l[0]}

    chrom_dict = {}

    print('Parsing the file: %s' % filename)
    with open(filename) as f:
        if debug:
            print(next(f))
        missing_chromosomes = set()
        for line in f:
            l = (line.strip()).split()
            sid = l[0]
            if debug and random.random() > 0.1:
                continue
            if sid in valid_sids:
                pos = snps_pos_map[sid]['pos']
                chrom = snps_pos_map[sid]['chrom']
                if not chrom in util.valid_chromosomes:
                    missing_chromosomes.add(chrom)
                    continue
                if not chrom in chrom_dict:
                    chrom_dict[chrom] = {'ps': [], 'log_odds': [], 'freqs': [],
                                         'betas': [], 'nts': [], 'sids': [], 'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                if l[3] == '.' or l[3] == 'NA':
                    freq = -1
                else:
                    freq = float(l[3])
                chrom_dict[chrom]['freqs'].append(freq)
                raw_beta = float(l[4])
                pval = float(l[6])
                n = float(l[7]) #Sample size
                chrom_dict[chrom]['ps'].append(pval)
                nt = [l[1], l[2]]
                chrom_dict[chrom]['nts'].append(nt)
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval / 2.0)
                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
                chrom_dict[chrom]['log_odds'].append(beta / sp.sqrt(n))

        util.check_chromosomes(missing_chromosomes)

    print('SS file loaded, now sorting and storing in HDF5 file.')
    assert not 'sum_stats' in hdf5_file, 'Something is wrong with HDF5 file?  Summary stats already found.'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict:
        print('%d SNPs on chromosome %d' % (len(chrom_dict[chrom]['positions']), chrom))
        sl = list(zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'],
                 chrom_dict[chrom]['freqs'], chrom_dict[chrom]['ps']))
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        freqs = []
        prev_pos = -1
        for pos, sid, nt, beta, lo, frq, p in sl:
            if pos == prev_pos:
                print('duplicated position %d' % pos)
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
        print('Still %d SNPs on chromosome %s' % (len(ps), chrom))
        g = ssg.create_group('chrom_%s' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('freqs', data=freqs)
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps += len(log_odds)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print('%d SNPs parsed from summary statistics file.' % num_snps)


def parse_sum_stats_giant2(filename=None,
                           bimfile=None,
                           hdf5_file=None,
                           debug=False):
    """
    Input format:

    MarkerName A1 A2 Freq.Hapmap.Ceu BETA SE.2gc P.2gc N
    rs4747841 a g 0.55 0.0025 0.0061 0.68 60558.2
    rs4749917 t c 0.45 -0.0025 0.0061 0.68 60558.1
    rs737656 a g 0.3667 -0.0073 0.0064 0.25 60529.2

    """
    lc_CAPs_dict = {'a': 'A', 'c': 'C', 'g': 'G', 't': 'T'}

    snps_pos_map = {}
    assert bimfile is not None, 'BIM file is required'
    valid_sids = set()
    print('Parsing bim file: %s' % bimfile)
    with open(bimfile) as f:
        for line in f:
            l = line.split()
            valid_sids.add(l[1])
            snps_pos_map[l[1]] = {'pos': int(l[3]), 'chrom': l[0]}

    chrom_dict = {}

    print('Parsing the file: %s' % filename)
    with open(filename) as f:
        if debug:
            print(next(f))
        missing_chromosomes = set()
        for line in f:
            l = (line.strip()).split()
            sid = l[0]
            if debug and random.random() > 0.1:
                continue
            if sid in valid_sids:
                pos = snps_pos_map[sid]['pos']
                chrom = snps_pos_map[sid]['chrom']
                if not chrom in util.valid_chromosomes:
                    missing_chromosomes.add(chrom)
                    continue

                if not chrom in chrom_dict:
                    chrom_dict[chrom] = {'ps': [], 'log_odds': [], 'freqs': [],
                                         'betas': [], 'nts': [], 'sids': [], 'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                if l[3] == '.' or l[3] == 'NA':
                    freq = -1
                else:
                    freq = float(l[3])
                chrom_dict[chrom]['freqs'].append(freq)
                raw_beta = float(l[4])
                pval = float(l[6])
                n = float(l[7])
                chrom_dict[chrom]['ps'].append(pval)
                nt = [lc_CAPs_dict[l[1]], lc_CAPs_dict[l[2]]]
                chrom_dict[chrom]['nts'].append(nt)
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval / 2.0)
                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
                chrom_dict[chrom]['log_odds'].append(beta / sp.sqrt(n))

        
        util.check_chromosomes(missing_chromosomes)

    print('SS file loaded, now sorting and storing in HDF5 file.')
    assert not 'sum_stats' in hdf5_file, 'Something is wrong with HDF5 file?  Summary stats already found.'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict:
        print('%d SNPs on chromosome %d' % (len(chrom_dict[chrom]['positions']), chrom))
        sl = list(zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'],
                 chrom_dict[chrom]['freqs'], chrom_dict[chrom]['ps']))
        sl.sort()
        ps = []
        betas = []
        nts = []
        sids = []
        positions = []
        log_odds = []
        freqs = []
        prev_pos = -1
        for pos, sid, nt, beta, lo, frq, p in sl:
            if pos == prev_pos:
                print('duplicated position %d' % pos)
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
        print('Still %d SNPs on chromosome %s' % (len(ps), chrom))
        g = ssg.create_group('chrom_%s' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('freqs', data=freqs)
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps += len(log_odds)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print('%d SNPs parsed from summary statistics file.' % num_snps)



def parse_sum_stats_pgc(filename=None,
                        bimfile=None,
                        hdf5_file=None):
    """
    Input format:

    CHR    SNP    BP    A1    A2    FRQ_A_30232    FRQ_U_40578    INFO    OR    SE    P    ngt    Direction    HetISqt    HetChiSq    HetDf    HetPVa
    ...

    """

    if bimfile is not None:
        print('Parsing SNP list')
        valid_sids = set()
        print('Parsing bim file: %s' % bimfile)
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])
        print('Parsed %d SNP IDs from bim file.'%len(valid_sids))

    chrom_dict = {}

    print('Parsing the file: %s' % filename)
    with open(filename) as f:
        print(next(f))
        missing_chromosomes = set()
        line_i = 0
        for line in f:
            line_i +=1
            l = (line.strip()).split()
            chrom = l[0]
            if not chrom in util.valid_chromosomes:
                missing_chromosomes.add(chrom)
                continue
            pos = int(l[2])
            sid = l[1]
            if sid in valid_sids:
                if not chrom in chrom_dict:
                    chrom_dict[chrom] = {'ps': [], 'log_odds': [], 'infos': [], 'freqs': [],
                                         'betas': [], 'nts': [], 'sids': [], 'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                freq_a = float(l[5])
                freq_u = float(l[6])
                n_cases = float(l[15])
                n_controls = float(l[16])
                n = n_cases+n_controls
                a_scalar =  n_cases/ n
                u_scalar = n_controls / n
                freq = freq_a * a_scalar + freq_u * u_scalar
                chrom_dict[chrom]['freqs'].append(freq)
                info = float(l[7])
                chrom_dict[chrom]['infos'].append(info)
                pval = float(l[10])
                chrom_dict[chrom]['ps'].append(pval)
                nt = [l[3], l[4]]
                chrom_dict[chrom]['nts'].append(nt)
                raw_beta = sp.log(float(l[8]))
                chrom_dict[chrom]['log_odds'].append(raw_beta)
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval / 2.0)

                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
        util.check_chromosomes(missing_chromosomes)

    print('SS file loaded, now sorting and storing in HDF5 file.')
    assert not 'sum_stats' in hdf5_file, 'Something is wrong with HDF5 file? Summary stats already found.'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict:
        print('Parsed summary stats for %d SNPs on chromosome %s' % (len(chrom_dict[chrom]['positions']), chrom))
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
                print('duplicated position %d' % pos)
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
        nts = sp.array(nts,dtype=nts_dtype)
        sids = sp.array(sids,dtype=sids_dtype)
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
    print('In all, %d SNPs parsed from summary statistics file.' % num_snps)


def parse_sum_stats_pgc_small(filename=None,
                              bimfile=None,
                              hdf5_file=None,
                              n=None, 
                              debug=False):
    """
    Input format:

    hg19chrc        snpid   a1      a2      bp      info    or      se      p       ngt
    chr1    rs4951859       C       G       729679  0.631   0.97853 0.0173  0.2083  0
    chr1    rs142557973     T       C       731718  0.665   1.01949 0.0198  0.3298  0
    ...

    """

    if bimfile is not None:
        print('Parsing SNP list')
        valid_sids = set()
        print('Parsing bim file: %s' % bimfile)
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])
        print(len(valid_sids))

    chrom_dict = {}

    print('Parsing the file: %s' % filename)
    with open(filename) as f:
        if debug:
            print (next(f))
        missing_chromosomes = set()
        for line in f:
            l = (line.strip()).split()
            chrom_str = l[0]
            chrom = chrom_str[3:]
            if not chrom in util.valid_chromosomes:
                missing_chromosomes.add(chrom)
                continue
            pos = int(l[4])
            sid = l[1]
            if sid in valid_sids:
                if not chrom in chrom_dict:
                    chrom_dict[chrom] = {'ps': [], 'log_odds': [], 'infos': [],
                                         'betas': [], 'nts': [], 'sids': [],
                                         'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                info = float(l[5])
                chrom_dict[chrom]['infos'].append(info)
                pval = float(l[8])
                chrom_dict[chrom]['ps'].append(pval)
                nt = [l[2], l[3]]
                chrom_dict[chrom]['nts'].append(nt)
                raw_beta = sp.log(float(l[6]))
                chrom_dict[chrom]['log_odds'].append(raw_beta)
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval / 2.0)

                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))
        util.check_chromosomes(missing_chromosomes)

    print('SS file loaded, now sorting and storing in HDF5 file.')
    assert not 'sum_stats' in hdf5_file, 'Something is wrong with HDF5 file? Summary stats already found.'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict:
        print('Parsed summary stats for %d SNPs on chromosome %s' % (len(chrom_dict[chrom]['positions']), chrom))
        sl = list(zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['infos'],
                 chrom_dict[chrom]['ps']))
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
                print('duplicated position %d' % pos)
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
        nts = sp.array(nts,dtype=nts_dtype)
        sids = sp.array(sids,dtype=sids_dtype)
        g = ssg.create_group('chrom_%s' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps += len(log_odds)
        g.create_dataset('infos', data=infos)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print('In all, %d SNPs parsed from summary statistics file.' % num_snps)


def parse_sum_stats_basic(filename=None,
                          bimfile=None,
                          hdf5_file=None,
                          n=None,
                          debug=False):
    """
    Input format:

    hg19chrc    snpid    a1    a2    bp    or    p       
    chr1    rs4951859    C    G    729679    0.97853    0.2083  
    chr1    rs142557973    T    C    731718    1.01949    0.3298  
    ...

    """

    if bimfile is not None:
        valid_sids = set()
        print('Parsing bim file: %s' % bimfile)
        with open(bimfile) as f:
            for line in f:
                l = line.split()
                valid_sids.add(l[1])
    chrom_dict = {}

    print('Parsing the file: %s' % filename)
    with open(filename) as f:
        missing_chromosomes = set()
        if debug:
            print(next(f))
        for line in f:
            l = (line.strip()).split()
            chrom_str = l[0]
            chrom = chrom_str[3:]
            if not chrom in util.valid_chromosomes:
                missing_chromosomes.add(chrom)
                continue
            pos = int(l[4])
            sid = l[1]
            if sid in valid_sids:
                if not chrom in chrom_dict:
                    chrom_dict[chrom] = {'ps': [], 'log_odds': [], 'infos': [],
                                         'betas': [], 'nts': [], 'sids': [],
                                         'positions': []}
                chrom_dict[chrom]['sids'].append(sid)
                chrom_dict[chrom]['positions'].append(pos)
                pval = float(l[6])
                chrom_dict[chrom]['ps'].append(pval)
                nt = [l[2], l[3]]
                chrom_dict[chrom]['nts'].append(nt)
                odds_ratio = float(l[5])
                assert odds_ratio > 0, 'The odds ratios appear to have negative values.  Please use an appropriate data format.'
                raw_beta = sp.log(odds_ratio)
                chrom_dict[chrom]['log_odds'].append(raw_beta)
                beta = sp.sign(raw_beta) * stats.norm.ppf(pval / 2.0)
                chrom_dict[chrom]['betas'].append(beta / sp.sqrt(n))

        util.check_chromosomes(missing_chromosomes)

    print('SS file loaded, now sorting and storing in HDF5 file.')
    assert not 'sum_stats' in hdf5_file, 'Something is wrong with HDF5 file? Summary stats already found.'
    ssg = hdf5_file.create_group('sum_stats')
    num_snps = 0
    for chrom in chrom_dict:
        print('Parsed summary stats for %d SNPs on chromosome %s' % (len(chrom_dict[chrom]['positions']), chrom))
        sl = list(zip(chrom_dict[chrom]['positions'], chrom_dict[chrom]['sids'], chrom_dict[chrom]['nts'],
                 chrom_dict[chrom]['betas'], chrom_dict[chrom]['log_odds'], chrom_dict[chrom]['ps']))
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
                print('duplicated position %d' % pos)
                continue
            else:
                prev_pos = pos
            ps.append(p)
            betas.append(beta)
            nts.append(nt)
            sids.append(sid)
            positions.append(pos)
            log_odds.append(lo)
        nts = sp.array(nts,dtype=nts_dtype)
        sids = sp.array(sids,dtype=sids_dtype)
        g = ssg.create_group('chrom_%s' % chrom)
        g.create_dataset('ps', data=sp.array(ps))
        g.create_dataset('betas', data=betas)
        g.create_dataset('log_odds', data=log_odds)
        num_snps += len(log_odds)
        g.create_dataset('nts', data=nts)
        g.create_dataset('sids', data=sids)
        g.create_dataset('positions', data=positions)
        hdf5_file.flush()
    print('In all, %d SNPs parsed from summary statistics file.' % num_snps)


