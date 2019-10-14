#!/usr/bin/env python

import argparse
from ldpred import coord_genotypes
from ldpred import LDpred_gibbs
from ldpred import LDpred_inf
from ldpred import LD_pruning_thres
from ldpred import validate
from ldpred import reporting

import sys
import textwrap
import random
import scipy

__version__ = '1.0.7'
__date__ = '11 October 2019'
 
#Set a random seed to avoid varying results.
random_seed=42
random.seed(random_seed)
scipy.random.seed(random_seed)

title = '\033[96mLDpred v. %s\033[0m'%__version__
title_len= len(title)-9
num_dashes = reporting.window_length-title_len-2
left_dashes = '='*int(num_dashes/2)
right_dashes = '='*int(num_dashes/2+num_dashes%2)
title_string = '\n'+left_dashes  + ' '+title+' '+right_dashes +'\n'

description_string = """
Typicial workflow:
    1. Use\033[1m LDpred coord\033[0m to parse summary statistics and genotypes 
       and coordinate data. See LDpred coord --help for further usage 
       description and options.
       
    2. Use\033[1m LDpred gibbs|inf|p+t\033[0m to obtain SNP weights for polygenic
       scoring using one or more methods. See LDpred gibbs|inf|p+t --help for 
       further usage description and options.
    
    3. Use\033[1m LDpred score\033[0m to calculate polygenic scores using SNP
       weights from previous step. See LDpred score --help for further usage 
       description and options.
         
Thank you for using\033[1m LDpred\033[0m!

%s 
(c) Bjarni J Vilhjalmsson: bjarni.vilhjalmsson@gmail.com

"""%(__date__)
description_string += '='*reporting.window_length

parser = argparse.ArgumentParser(prog='LDpred',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent(description_string))
                                                             

subparsers = parser.add_subparsers(help='Select ldpred action among the following options: ', dest='ldpred_action')
parser_coord = subparsers.add_parser('coord', help='Parse and coordinate summary statistics and genotypes.')
parser_gibbs = subparsers.add_parser('gibbs', help='Estimate LDpred (Gibbs sampler) SNP weights. (Requires a coordinated dataset.)')
parser_inf = subparsers.add_parser('inf', help='Estimate LDpred-inf SNP weights. (Requires a coordinated dataset.)')
parser_pt = subparsers.add_parser('p+t', help='Obtain pruning+thresholding SNP weights. (Requires a coordinated dataset.)')
parser_score = subparsers.add_parser('score', help='Calculate polygenic scores using given SNP weights.')

#General arguments
parser.add_argument('--debug', default=False, action='store_true',
                    help="Activate debugging mode (more verbose)")

#coord arguments 
parser_coord.add_argument('--gf', type=str, required=True,
                    help='LD Reference Genotype File. '
                         'Should be a (full path) filename prefix to a standard PLINK bed file (without .bed). '
                         'Make sure that the fam and bim files with same names are in the same directory. ')
parser_coord.add_argument('--ssf', type=str, required=True,
                    help='Summary Statistic File. '
                         'Filename for a text file with the GWAS summary statistics')
parser_coord.add_argument('--N', type=int, default=None,
                    help='Number of Individuals in Summary Statistic File.  Required for most summary '
                         'statistics formats.')
parser_coord.add_argument('--out', type=str, required=True,
                    help='Output Prefix')
parser_coord.add_argument('--vbim', type=str, default=None,
                    help='Validation SNP file. '
                         'This is a PLINK BIM file (with .bim) which can be used to filter the set of SNPs down '
                         'to the set of validation SNPs. To maximize accuracy, we recommend calculating LDpred '
                         'weights for the subset of SNPs that are used to calculate the risk scores in the '
                         'target (validation) sample.')
parser_coord.add_argument('--vgf', type=str, default=None,
                    help='Validation genotype file. '
                         'This is a PLINK BIM file which can be used to filter the set of SNPs down to the '
                         'set of validation SNPs. This is a filename prefix (without .bed). To maximize '
                         'accuracy, we recommend calculating LDpred weights for the subset of SNPs that are '
                         'used to calculate the risk scores in the '
                         'target (validation) sample.')
parser_coord.add_argument('--only-hm3', default=False, action='store_true',
                    help='Restrict analysis to 1.2M HapMap 3 SNPs.')
parser_coord.add_argument('--ilist', type=str,
                    help='List of individuals to include in the analysis. ', default=None)
parser_coord.add_argument('--skip-coordination', default=False, action='store_true',
                    help="Assumes that the alleles have already been coordinated between LD reference, "
                         "validation samples, and the summary statistics files")
parser_coord.add_argument('--beta', default=False, action='store_true',
                    help='Assumes the summary statistics are BETA (linear regression) instead of OR (logistic '
                         'regression)')
parser_coord.add_argument('--match-genomic-pos', default=False, action='store_true',
                    help='Excludes SNPs from sum stats file if their genomics positions differ from validation '
                    'data.  These will also differ if reference genome differ, e.g. hg19 vs hg38.')
parser_coord.add_argument('--maf', type=float, default=0.01,
                    help='MAF filtering threshold.  Set to 0 to disable MAF filtering.')
parser_coord.add_argument('--max-freq-discrep', type=float, default=0.1,
                    help='Max frequency discrepancy allowed between reported sum stats frequency and frequency '
                    'in the LD reference data.  To disable check, set to 1.')
parser_coord.add_argument('--ssf-format', type=str, default="CUSTOM", choices={'LDPRED','CUSTOM','STANDARD','GIANT', 'BASIC', 'PGC'},
                    help='This is the format type of the summary statistics file. '
                    'By default the CUSTOM format requires the user to specify the file format using additional '
                    'arguments.')
parser_coord.add_argument('--rs', type=str, default="SNP",
                    help="Column header of SNP ID")
parser_coord.add_argument('--A1', type=str, default="A1",
                    help="Column header containing the effective allele. "
                         "There isn't any standardized label for the effective allele, "
                         "therefore extra care must be taken to ensure the correct label is provided, "
                         "otherwise, the effect will be flipped.")
parser_coord.add_argument('--A2', type=str, default="A2",
                    help="Column header containing non-effective allele.")
parser_coord.add_argument('--pos', type=str, default="BP",
                    help="Column header containing the coordinate of SNPs.")
parser_coord.add_argument('--info', type=str, default="INFO",
                    help="Column header containing the INFO score.")
parser_coord.add_argument('--chr', type=str, default="CHR",
                    help="Column header containing the chromosome information.")
parser_coord.add_argument('--reffreq', type=str, default="MAF",
                    help="Column header containing the reference MAF")
parser_coord.add_argument('--pval', type=str, default="P",
                    help="Column header containing the P-value information")
parser_coord.add_argument('--eff', type=str, default="OR",
                    help="Column header containing effect size information")
parser_coord.add_argument('--se', type=str, default="SE",
                    help='Column header containing standard error. Note: these are only used if corresponding p-values'
                    ' are smaller than approximately 1E-325, or if --z-from-se is set.')
parser_coord.add_argument('--ncol', type=str, default="N",
                    help="Column header containing sample size information")
parser_coord.add_argument('--case-freq', type=str, default=None,
                    help="Column header containing case frequency information")
parser_coord.add_argument('--control-freq', type=str, default=None,
                    help="Column header containing control frequency information")
parser_coord.add_argument('--case-n', type=str, default=None,
                    help="Column header containing case sample size information")
parser_coord.add_argument('--control-n', type=str, default=None,
                    help="Column header containing control sample size information")
parser_coord.add_argument('--z-from-se', default=False, action='store_true',
                    help='If this option is set, then the effects use by LDpred are derived using the effect estimates '
                    'and their standard errors. (Only possible if standard errors are available in the summary stats)')
# parser_coord.add_argument('--gmdir', type=str,
#                     help='The directory of genetic map.', default=None)



#gibbs arguments 
parser_gibbs.add_argument('--cf', type=str, required=True,
                    help='Coordinated file (generated using ldpred coord). '
                         'Should be a (full path) filename. ')
parser_gibbs.add_argument('--ldr', type=int, required=True,
                    help='LD radius.  An integer number which denotes the number of SNPs on each side '  
                    'of the focal SNP for which LD should be adjusted. A value corresponding M/3000, '
                    'where M is the number of SNPs in the genome is recommended')
parser_gibbs.add_argument('--ldf', type=str, required=True,
                    help='LD file (prefix). A path and filename prefix for the LD file. If it does not '
                    'exist, it will be generated.  This can take up to several hours, depending on '
                    'LD radius used.')
parser_gibbs.add_argument('--out', type=str, required=True,
                    help='Output Prefix for SNP weights')
parser_gibbs.add_argument('--f', default=[1,0.3,0.1,0.03,0.01,0.003,0.001], nargs='+', type=float,
                    help="Fraction of causal variants used in the Gibbs sampler")

parser_gibbs.add_argument('--N', type=int, default=None,
                    help='Number of individuals on which the summary statistics are assumed to be based on.')
parser_gibbs.add_argument('--n-iter', type=int, default=60,
                    help='The number of iterations used by the Gibbs sampler. The default is 60.')
parser_gibbs.add_argument('--n-burn-in', type=int, default=5,
                    help='The number of burn-in iterations used by the Gibbs sampler. The default is 5.')
parser_gibbs.add_argument('--h2', type=float, default=None, 
                    help='The genome-wide heritability assumed by LDpred, which is then partitioned '
                    'proportional to the number of SNPs on each chromosome.  By default it estimates the '
                    'heritability for each chromosome from the GWAS summary statistics using LD score '
                    'regression (Bulik-Sullivan et al., Nat Genet 2015).')
parser_gibbs.add_argument('--use-gw-h2', default=False, action='store_true',
                    help='Estimate heritability genome-wide and partition it proportional to the number of SNPs on each chromsome '
                    'instead of estimating it for each chromosome separately.  This is like a likely good choice if '
                    'the summary statistics used are based on small sample sizes (approx <50K), or if the trait is '
                    'not very heritable.')
parser_gibbs.add_argument('--no-ld-compression', default=False, action='store_true',
                    help='Do not compress LD information.  Saves storing and loading time of LD information, '
                    'but takes more space on disk.')
parser_gibbs.add_argument('--hickle-ld', default=False, action='store_true',
                    help='Use hickle instead of pickle for storing LD files.  This saves memory, but generally'
                    'takes more time to write and load. Requires hickle to be installed on your system, '
                    'see http://telegraphic.github.io/hickle/ for how to install.')

# parser_gibbs.add_argument('--gm-ldr', type=float, default=None, 
#                     help='If this option is set, then a genetic map will be used to calculate LD-radius. '
#                     'A value around 1 is arguably reasonable.')


#inf arguments 
parser_inf.add_argument('--cf', type=str, required=True,
                    help='Coordinated file (generated using ldpred coord). '
                         'Should be a (full path) filename. ')
parser_inf.add_argument('--ldr', type=int, required=True,
                    help='LD radius.  An integer number which denotes the number of SNPs on each side '  
                    'of the focal SNP for which LD should be adjusted. A value corresponding M/3000, '
                    'where M is the number of SNPs in the genome is recommended')
parser_inf.add_argument('--ldf', type=str, required=True,
                    help='LD file (prefix). A path and filename prefix for the LD file. If it does not '
                    'exist, it will be generated.  This can take up to several hours, depending on '
                    'LD radius used.')
parser_inf.add_argument('--out', type=str, required=True,
                    help='Output Prefix for SNP weights')
parser_inf.add_argument('--N', type=int, default=100000, required=True,
                    help='Number of individuals on which the summary statistics are assumed to be based on.')
parser_inf.add_argument('--h2', type=float, default=None, 
                    help='The genome-wide heritability assumed by LDpred, which is then partitioned '
                    'proportional to the number of SNPs on each chromosome.  By default it estimates the '
                    'heritability for each chromosome from the GWAS summary statistics using LD score '
                    'regression (Bulik-Sullivan et al., Nat Genet 2015).')
parser_inf.add_argument('--use-gw-h2', default=False, action='store_true', 
                    help='Estimate heritability genome-wide and partition it proportional to the number of SNPs on each chromsome '
                    'instead of estimating it for each chromosome separately.  This is like a likely good choice if '
                    'the summary statistics used are based on small sample sizes (approx <50K), or if the trait is '
                    'not very heritable.')
parser_inf.add_argument('--no-ld-compression', default=False, action='store_true',
                    help='Do not compress LD information.  Saves storing and loading time of LD information, '
                    'but takes more space on disk.')
parser_inf.add_argument('--hickle-ld', default=False, action='store_true',
                    help='Use hickle instead of pickle for storing LD files.  This saves memory, but generally'
                    'takes more time to write and load. Requires hickle to be installed on your system, '
                    'see http://telegraphic.github.io/hickle/ for how to install.')
# parser_inf.add_argument('--gm-ldr', type=float, default=None, 
#                     help='If this option is set, then a genetic map will be used to calculate LD-radius. '
#                     'A value around 1 is arguably reasonable.')


parser_pt.add_argument('--cf', type=str, required=True,
                    help='Coordinated file (generated using ldpred coord). '
                         'Should be a (full path) filename. ')
parser_pt.add_argument('--ldr', type=int, required=True,
                    help='LD radius.  An integer number which denotes the number of SNPs on each side '  
                    'of the focal SNP for which LD should be adjusted. A value corresponding M/3000, '
                    'where M is the number of SNPs in the genome is recommended')
parser_pt.add_argument('--out', type=str, required=True,
                    help='Output Prefix for SNP weights')
parser_pt.add_argument('--p', default=[1,0.3,0.1,0.03,0.01,0.003,0.001,3*1E-4,
                   1E-4,3*1E-5,1E-5,1E-6,1E-7,1E-8], nargs='+', type=float,
                    help="P value thresholds")
parser_pt.add_argument('--r2', default=[0.2], nargs='+', type=float,
                    help='LD R2 thresholds. The maximum LD squared correlation allowed between any SNPs '
                    'within the given LD-radius.  Default max R2 value is set to 0.2')
# parser_pt.add_argument('--gm-ldr', type=float, default=None, 
#                     help='If this option is set, then a genetic map will be used to calculate LD-radius. '
#                     'A value around 1 is arguably reasonable.')


parser_score.add_argument('--gf', type=str, default=None,
                    help='Validation genotype file. '
                         'PLINK formatted genotypes for which we want to calculate risk scores.')
parser_score.add_argument('--rf', type=str, required=True,
                    help='SNP weights file (prefix used for output in preceding step), e.g. LDpred SNP weights. '
                    'Note that the prefix is the same as given as output in the previous ste, and should exclude '
                    'the _LDpred_.. extension that are added to the file.')
parser_score.add_argument('--out', type=str, required=True,
                    help='The prefix of risk score output file.')
parser_score.add_argument('--pf', type=str, default=None,
                    help='A file with individual IDs and phenotypes.')
parser_score.add_argument('--pf-format', type=str, default='STANDARD', choices = {'STANDARD','FAM','LSTANDARD'},
                    help='Phenotype file format. Two formats are supported a (PLINK) FAM format, and '
                    'STANDARD format (default), which is a whitespace/tab delimited file with two '
                    'columns IID and PHEN.')
parser_score.add_argument('--rf-format', type=str, default='ANY', choices = {'ANY','LDPRED', 'P+T'},
                    help='The format to expect the results to be in.')
parser_score.add_argument('--cov-file', type=str, default=None,
                    help='Covariate file, format is whitespace-delimted file with columns IID, COV1, COV2, etc.')
parser_score.add_argument('--pcs-file', type=str, default=None,
                    help='PCs file, format is whitespace-delimted file with columns FID, IID, PC1, PC2, etc.')
parser_score.add_argument('--split-by-chrom', default=False, action='store_true',
                    help='Split analysis by chromosome')
parser_score.add_argument('--only-score', default=False, action='store_true',
                    help='Only output scores, and do not validate predictions. No phenotypes or covariates '
                    'will be parsed.')
parser_score.add_argument('--f', default=[1,0.3,0.1,0.03,0.01,0.003,0.001], nargs='+', type=float,
                    help="Fraction of causal variants used in the Gibbs sampler")
parser_score.add_argument('--p', default=[1,0.3,0.1,0.03,0.01,0.003,0.001,3*1E-4,
                                          1E-4,3*1E-5,1E-5,1E-6,1E-7,1E-8], nargs='+', type=float,
                    help="P value thresholds (P+T)")
parser_score.add_argument('--r2', default=[1,0.2], nargs='+', type=float,
                    help="LD R2 thresholds (P+T)")



def main_with_args(args):
    print(title_string)
    if len(args)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parameters = parser.parse_args(args)
    p_dict= vars(parameters)
    if p_dict['debug']:
        print ('Parsed parameters:')
        print(p_dict)
    
    action = p_dict['ldpred_action']
    if action=='coord':
        coord_genotypes.main(p_dict)
    elif action=='gibbs':
        LDpred_gibbs.main(p_dict)
    elif action=='inf':
        LDpred_inf.main(p_dict)
    elif action=='p+t':
        LD_pruning_thres.main(p_dict)
    elif action=='score':
        validate.main(p_dict)
    elif action=='all':
        pass
    
    
    
def main():
  main_with_args(sys.argv[1:])


if __name__ == '__main__':
    main()
