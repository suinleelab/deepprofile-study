###############################
#Author: Ayse Dincer
#Script for running multiple FETs
###############################

import sys

cancer_type = sys.argv[1]
pathway = sys.argv[2]

get_ipython().magic(u"run -i 'Fishers_Test.py' '" +  cancer_type + "' " + pathway + " " + "DeepProfile" + " " + str(0) + " " + str(1))
get_ipython().magic(u"run -i 'Fishers_Test.py' '" +  cancer_type + "' " + pathway + " " + "PCA" + " " + str(0) + " " + str(1))
get_ipython().magic(u"run -i 'Fishers_Test.py' '" +  cancer_type + "' " + pathway + " " + "ICA" + " " + str(0) + " " + str(10))
get_ipython().magic(u"run -i 'Fishers_Test.py' '" +  cancer_type + "' " + pathway + " " + "RP" + " " + str(0) + " " + str(10))
get_ipython().magic(u"run -i 'Fishers_Test.py' '" +  cancer_type + "' " + pathway + " " + "AE" + " " + str(-1) + " " + str(9))
get_ipython().magic(u"run -i 'Fishers_Test.py' '" +  cancer_type + "' " + pathway + " " + "DAE" + " " + str(-1) + " " + str(9))
