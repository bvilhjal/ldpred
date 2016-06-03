Copyright (C) 2015, Bjarni J. Vilhjalmsson (bjarni.vilhjalmsson@gmail.com)

# LDpred #


LDpred is a Python based software package that adjusts GWAS summary statistics 
for the effects of linkage disequilibrium (LD).  The details of the method is 
described in Vilhjalmsson et al. (AJHG 2015) [http://www.cell.com/ajhg/abstract/S0002-9297(15)00365-1]

* The current version is 0.6.


## Getting Started ##

### Requirements ###
LDpred currently requires three Python packages to be installed and in path.  These 
are **h5py** [http://www.h5py.org/](http://www.h5py.org/), **scipy** [http://www.scipy.org/](http://www.scipy.org/) 
and **libplinkio** [https://github.com/mfranberg/libplinkio](https://github.com/mfranberg/libplinkio).  Lastly, LDpred
has currently only been tested with **Python 2.x**.  

The first two packages **h5py** and **scipy** are commonly used Python packages, and pre-installed on many computer systems. The last **libplinkio** package can be installed using **pip** (https://pip.pypa.io/en/latest/quickstart.html), which is also pre-installed on many systems.

With **pip**, one can install **libplinkio** using the following command:

`pip install plinkio`

or if you need to install it locally you can try

`pip install --user plinkio`

With these three packages in place, you should be all set to install and use LDpred.

### Installing LDpred ###

As with most Python packages, configurating LDpred is simple.  You can either use **git** (which is installed on most systems) and clone this repository using the following git command:

`git clone git@bitbucket.org:bjarni_vilhjalmsson/ldpred.git`

Alternatively, you can simply download the source files and place them somewhere.

With the Python source code in place and the three packages **h5py**, **scipy**, and **libplinkio** installed, then you should be ready to use LDpred.


### How to run tests ###
A couple of simulated data examples can be found in the **test_data** directory.  These datasets were simulated using two different values of *p* (fraction of causal markers) and with heritability set to 0.1. The sample size used when simulating the summary statistics is 10,000.


### Code Contributions ###
We encourage users to extend the code, and adapt it too their needs.  Currently there are no formal guidelines set for 
contributions, and pull requests will be reviewed on a case by case basis.  

### Who do I talk to? ###
If you have any questions or trouble getting the method to work, please contact Bjarni Vilhjalmsson (bjarni.vilhjalmsson@gmail.com).

## Using LDpred ##
To run LDpred, at least two steps are required:

1. The first step is a data synchronization step, where two or three data sets, genotypes and summary statistics are synchronized.  This generates a HDF5 file which contains the synchronized genotypes.  This step is implemented in the file ** coord_genotypes.py **.  This step requires at least one genotype file (the LD reference genotypes), where we recommend at least 1000 unrelated individuals with the same ancestry make-up as the individuals for which summary statistics datasets are obtained from.  Another genotype file can also be given if the user intends to validate the predictions using a separate set of genotypes.

2. After generating the coordinated data file then the one can apply LDpred and run it on the synchronized dataset.  This step is implemented in ** LDpred.py **.  This step generates two files, a LD file with LD information for the given LD radius, and the re-weighted effect estimates.  The LD file enables the user to not have to generate the LD file again when trying, e.g., different values of *p* (the fraction of causal variants). However, it is re-generated if a different LD radius is given.  The other file that LDpred generates contains the LDpred-adjusted effect estimates. 

## Generating individual risk scores ##
Individual risk scores can be generated using the **validate.py** script.  It calculates polygenic risk scores for the individuals in the validation data if given, otherwise it treats the LD reference genotypes as validation genotypes.  A phenotype file can be provided, covariate file, as well as plink-formatted principal components file.  



## LD-pruning + Thresholding ##
In addition to the LDpred.py script, which implements the LDpred algorithm, the package also implements LD-pruning + Thresholding, as an alternative method.  This method often yields better predictions than LDpred when the LD reference panel is small.  This method is implemented in the ** LD_pruning_thres.py ** script.


## Individual scripts and their usage information ##


### coord_genotypes.py ###
Coordinate genotypes and summary statistics datasets for calculating polygenic risk scores.  Only SNPs that overlap between the two (or three) different datasets are retained.  

Usage: 

`python coord_genotypes.py --gf=PLINK_LD_REF_GENOTYPE_FILE --ssf=SUM_STATS_FILE --N=SS_SAMPLE_SIZE  --out=OUT_COORD_FILE [--vgf=PLINK_VAL_GENOTYPE_FILE --vbim=VAL_PLINK_BIM_FILE  --ssf_format=SSF_FORMAT --gmdir=GENETIC_MAP_DIR --gf_format=GENOTYPE_FILE_FORMAT --indiv_list=INDIV_LIST_FILE  --skip_coordination]`

 * PLINK_LD_REF_GENOTYPE_FILE (and PLINK_VAL_GENOTYPE_FILE) should be a (full path) filename prefix to a standard PLINK bed file (without .bed) Make sure that the fam and bim files with same names are in the same directory. PLINK_LD_REF_GENOTYPE_FILE refers LD reference genotypes, and PLINK_VAL_GENOTYPE_FILE refers to validation genotypes. It is not necessary to have LD validation genotypes at this stage.

 * SUM_STATS_FILE should be a (full path) filename prefix for a text file with the GWAS summary statistics.  Several formats are supported, the STANDARD format is as follows:
```
    chr     pos     ref     alt     reffrq  info    rs       pval    effalt
    chr1    1020428 C       T       0.85083 0.98732 rs6687776    0.0587  -0.0100048507289348
    chr1    1020496 G       A       0.85073 0.98751 rs6678318    0.1287  -0.00826075392985992
    ..
    ..
```
    
 * SS_SAMPLE_SIZE should be the approximate number of individuals used for calculating the GWAS summary statistics.

 * OUT_COORD_FILE is the output file.  This file will follow a HDF5 format and contain both LD-reference genotypes and summary statistics.
 
 * VAL_PLINK_BIM_FILE (optional): This is a PLINK BIM file which can be used to filter the set of SNPs down to the set of validation SNPs.  To maximize accuracy, it is best to calculate the LDpred weights for the SNPs that are used to calculate the risk scores.
 
 * SSF_FORMAT (optional): This is the format type of the summary statistics file.  Currently there are two implementations, "STANDARD", "BASIC", "PGC", and "PGC_large".  The standard format is described above.   
```
BASIC" format, which contains of the basic required information, is as follows:
    hg19chrc    snpid    a1    a2    bp    or    p       
    chr1    rs4951859    C    G    729679    0.97853    0.2083  
    chr1    rs142557973    T    C    731718    1.01949    0.3298  
    ..
    ..
```

 * GENOTYPE_FILE_FORMAT (optional): The expected genotype format.  The standard format is PLINK.  Other formats implemented is DECODE format.  If the DECODE format is used, then the program assumes that the data directory is supplied in the --gf flag.
   
 * INDIV_LIST_FILE (optional): List of individuals to include in the analysis.  Currently required for the DECODE format.
 
 * The *--skip_cordination* flag assumes that the alleles have already been coordinated between LD reference, validation samples,
   and the summary statistics files.

 * GENETIC_MAP_DIR (optional): The directory of genetic a genetic map. 


### LDpred.py ###
Implements LDpred, an approximate Gibbs sampler that calculate posterior means of effects, conditional on LD information.  The method requires the user to have generated a coordinated dataset using coord_genotypes.py

Usage: 

`python LDpred.py --coord=COORD_DATA_FILE  --ld_radius=LD_RADIUS   --local_ld_file_prefix=LD_FILE_NAME  --PS=FRACTIONS_CAUSAL  --N=SAMPLE_SIZE  --out=OUTPUT_FILE_PREFIX  [ --num_iter=NUM_ITER  --H2=HERTIABILITY  --gm_ld_radius=GEN_MAP_RADIUS]`
    
 * COORD_DATA_FILE: The HDF5 file obtained by running the coord_genotypes.py
 
 * LD_RADIUS: An integer number which denotes the number of SNPs on each side of the focal SNP for which LD should be adjusted. A value corresponding M/3000, where M is the number of SNPs used for the analysis is reasonable for genome length of 3000Mb.  This should result in a LD-radius of about 1Mb on average.
 
 * LD_FILE_NAME: A path and filename prefix for the LD file.  If it doesn't exist, it will be generated.  This can take up to several hours, depending on LD radius number of SNPs, etc.  If it does exits, that file will be used.
    
 * FRACTION_CAUSAL: A list of comma separated (without space) values between 1 and 0, excluding 0.  1 corresponds to the infinitesimal model and will yield results similar to LDpred-inf.  Default is --PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001
 
 * N: This is the sample size which LDpred assumes was used to calculate the GWAS summary statistics.

 * OUTPUT_FILE_PREFIX:  The prefix of output file.  

 * NUM_ITER (optional): The number of iterations used by the Gibbs sampler. The default is 60, and burn-in is fixed to 5.
 
 * HERTIABILITY (optional): The heritability assumed by LDpred.  By default it estimates the heritability from the GWAS summary statistics.

* GEN_MAP_RADIUS (optional):  If this option is set, then a genetic map will be used to calculate LD-radius.  A value around 1 is arguably reasonable.


### validate.py ###
Takes **LDpred.py** (or **LD_pruning_thres.py**) effect estimates , and (validation) genotypes in PLINK bed format as input.  The script then works out overlap and outputs predictions or risk scores as well as some prediction accuracy statistics.
    
Note that for maximal accuracy all SNPs with LDpred weights should be included in the validation dataset. If they are a subset of the validation dataset, then we suggest recalculate LDpred for the overlapping SNPs.

Usage: 

`python validate.py --vgf=PLINK_VAL_GENOTYPE_FILE  --rf=RESULT_FILE_PREFIX  --out=OUTPUT_FILE_PREFIX  [--res_format=LDPRED --split_by_chrom --pf=PHEN_FILE --pf_format=STANDARD --cov_file=COVARIATE_FILE --pcs_file=PCS_FILE --PS=FRACTIONS_CAUSAL  --TS=PVAL_THRESHOLDS]`
    
 * PLINK_VAL_GENOTYPE_FILE: PLINK formatted genotypes for which we want to calculate risk scores.
 
 * RESULT_FILE_PREFIX: SNP weights file, e.g. LDpred SNP weights.
 
 * OUTPUT_FILE_PREFIX:  The prefix of output file.  
    
 * RESULT_FILE_FORMAT: The format to expect the results to be in.  The default format is LDPRED, which refers to the format which running LDpred output. LDPRED-INF and P+T (LD-pruning + p-value thresholding) are also implemented.
 
 * PHEN_FILE: Is a file with individual IDs and phenotypes

 * PVAL_THRESHOLDS: This option is only valid if a P+T result file prefix is supplied.  It's a list of p-value thresholds, separated by a comma (without space), to be used for LDpred. Default values are --TS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001,3E-5,1E-5,1E-6,1E-7,1E-8

 * FRACTIONS_CAUSAL: This option is only valid if a LDPRED result file prefix is supplied.  A list of comma separated (without space) values between 1 and 0, excluding 0.  1 corresponds to the infinitesimal model and will yield results similar to LDpred-inf.  Default values are --PS=1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001