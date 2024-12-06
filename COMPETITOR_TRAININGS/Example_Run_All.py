###############################
#Example for training competitor models for a cancer type
###############################

get_ipython().magic(u"run -i Create_PCA_Data.py BRCA")

get_ipython().magic(u"run -i Create_ICA_Data.py BRCA")

get_ipython().magic(u"run -i Create_RP_Data.py BRCA")

get_ipython().magic(u"run -i Train_AE_Models.py BRCA")
get_ipython().magic(u"run -i Get_AE_IG_Attributions.py BRCA 0")

get_ipython().magic(u"run -i Train_DAE_Models.py BRCA")
get_ipython().magic(u"run -i Get_DAE_IG_Attributions.py BRCA 0")
