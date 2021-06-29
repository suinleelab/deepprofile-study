###############################
#Author: Ayse Dincer
#Date: May 7 2019
#Script for training DAE models
###############################
import sys
cancer_type = sys.argv[1]

for run in range(10):
    get_ipython().magic(u"run -i 'DAE_2Layers_Model.py' " + cancer_type + " " + str(run))
