"""
Copyright (C) 2016, Bjarni J. Vilhjalmsson (bjarni.vilhjalmsson@gmail.com)

LDpred is a Python based software package that adjusts GWAS summary statistics 
for the effects of linkage disequilibrium (LD).  The details of the method is 
described in Vilhjalmsson et al. (AJHG 2015) [http://www.cell.com/ajhg/abstract/S0002-9297(15)00365-1]

"""

#Load ldpred packages..
import LDpred
import validate
import LD_pruning_thres
import LDpred_inf
import coord_genotypes