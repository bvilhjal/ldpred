
# LDpred #


LDpred is a Python based software package that adjusts GWAS summary statistics
for the effects of linkage disequilibrium (LD).  The details of the method is
described in Vilhjalmsson et al. (AJHG 2015) [http://www.cell.com/ajhg/abstract/S0002-9297(15)00365-1]

* The current version is 1.0.9

## Getting Started ##
### News ###

Oct 17th, 2019: Version 1.0.9 released and is available on pip using 

`pip install ldpred`

And it should actually work now :)

Also, tests have been markedly improved.


### Requirements ###
LDpred currently requires three Python packages to be installed and in path.  These
are **h5py** [http://www.h5py.org/](http://www.h5py.org/), **scipy** [http://www.scipy.org/](http://www.scipy.org/)
and **libplinkio** [https://github.com/mfranberg/libplinkio](https://github.com/mfranberg/libplinkio).  Lastly, LDpred
has currently only been tested with **Python 3.6+**.

The first two packages **h5py** and **scipy** are commonly used Python packages, and pre-installed on many computer systems. The last **libplinkio** package can be installed using **pip** (https://pip.pypa.io/en/latest/quickstart.html), which is also pre-installed on many systems.

With **pip**, one can install **libplinkio** using the following command:

`pip install plinkio`

or if you need to install it locally you can try

`pip install --user plinkio`

With these three packages in place, you should be all set to install and use LDpred.

### Installing LDpred ###

As with most Python packages, configurating LDpred is simple.  You can use **pip** to install it by typing

`pip install ldpred`

This should automatically take care of dependencies.  The examples below assume ldpred has been installed using pip.

Alternatively you can use **git** (which is installed on most systems) and clone this repository using the following git command:

`git clone https://github.com/bvilhjal/ldpred.git`

Finally, you can also download the source files and place them somewhere.

With the Python source code in place and the three packages **h5py**, **scipy**, and **libplinkio** installed, then you should be ready to use LDpred.


### How to run tests ###
A couple of simulated data examples can be found in the **test_data** directory.  These datasets were simulated using two different values of *p* (fraction of causal markers) and with heritability set to 0.1. The sample size used when simulating the summary statistics is 10,000.


### Code Contributions ###
I encourage users to extend the code, and adapt it too their needs.  Currently there are no formal guidelines set for
contributions, and pull requests will be reviewed on a case by case basis.

### Who do I talk to? ###
If you have any questions or trouble getting the method to work, try first to look at issues, to see if it is reported there.  Also, you can check if some of the cloned LDpred repos have addressed your issue.

In emergencies, please contact Bjarni Vilhjalmsson (bjarni.vilhjalmsson@gmail.com), but expect slow replies.  

## Using LDpred ##
A typical LDpred workflow consists of 3 steps:

### Step 1: Coordinate data ###
The first step is a data synchronization step, where two or three data sets, genotypes and summary statistics are synchronized.  This generates a HDF5 file which contains the synchronized genotypes.  This step can be done by running 

`ldpred coord`

use --help for detailed options.  This step requires at least one genotype file (the LD reference genotypes), where we recommend at least 1000 unrelated individuals with the same ancestry make-up as the individuals for which summary statistics datasets are obtained from.  Another genotype file can also be given if the user intends to validate the predictions using a separate set of genotypes.

### Step 2: Generate LDpred SNP weights ###
After generating the coordinated data file then the one can apply LDpred and run it on the synchronized dataset.  This step can be done by running 

`ldpred gibbs`

use --help for detailed options.  This step generates two files, a LD file with LD information for the given LD radius, and the re-weighted effect estimates.  The LD file enables the user to not have to generate the LD file again when trying, e.g., different values of **p** (the fraction of causal variants). However, it is re-generated if a different LD radius is given.  The other file that LDpred generates contains the LDpred-adjusted effect estimates.

### Step 3: Generating individual risk scores ###
Individual risk scores can be generated using the following command

`ldpred score`

use --help for detailed options.  It calculates polygenic risk scores for the individuals in the validation data if given, otherwise it treats the LD reference genotypes as validation genotypes.  A phenotype file can be provided, covariate file, as well as plink-formatted principal components file.


### Additional methods: LD-pruning + Thresholding ###
In addition to the LDpred gibbs sampler and infinitesimal model methods, the package also implements LD-pruning + Thresholding as an alternative method. You can run this using the following command

`ldpred p+t`

This method often yields better predictions than LDpred when the LD reference panel is small, or when the training data is very large (due to problems with gibbs sampler convergence).

### Tests###
You can run two tests to see if LDpred work on your system by running the following tests

`ldpred-unittest`

Note that passing these tests may not guarantee that LDpred work in all situations.

### Citation ###
Please cite [this paper](https://doi.org/10.1016/j.ajhg.2015.09.001)

### Acknowledges ###
Thanks to all who provided bug reports and contributed code.
