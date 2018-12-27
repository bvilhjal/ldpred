"""
Copyright (C) 2016, Bjarni J. Vilhjalmsson (bjarni.vilhjalmsson@gmail.com)

LDpred is a Python based software package that adjusts GWAS summary statistics 
for the effects of linkage disequilibrium (LD).  The details of the method is 
described in Vilhjalmsson et al. (AJHG 2015) [http://www.cell.com/ajhg/abstract/S0002-9297(15)00365-1]

"""

#Load ldpred packages..
from . import LDpred
from . import validate
from . import LD_pruning_thres
from . import LDpred_inf
from . import coord_genotypes