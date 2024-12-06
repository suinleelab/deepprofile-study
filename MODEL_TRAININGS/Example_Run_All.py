###############################
#Example for training model for a cancer type
###############################
import sys

##STEP 1: Creating PCs
get_ipython().magic(u"run -i Create_PCs_for_DeepLearning_Models.py BRCA 1000")

##STEP 2: Training VAE models
get_ipython().magic(u"run -i Run_VAE_Models.py BRCA 5 0 100")
get_ipython().magic(u"run -i Run_VAE_Models.py BRCA 10 0 100")
get_ipython().magic(u"run -i Run_VAE_Models.py BRCA 25 0 100")
get_ipython().magic(u"run -i Run_VAE_Models.py BRCA 50 0 100")
get_ipython().magic(u"run -i Run_VAE_Models.py BRCA 75 0 100")
get_ipython().magic(u"run -i Run_VAE_Models.py BRCA 100 0 100")

##STEP 3: Running IG for VAE models
get_ipython().magic(u"run -i Get_VAE_IG_Attributions.py BRCA 5 0 100")
get_ipython().magic(u"run -i Get_VAE_IG_Attributions.py BRCA 10 0 100")
get_ipython().magic(u"run -i Get_VAE_IG_Attributions.py BRCA 25 0 100")
get_ipython().magic(u"run -i Get_VAE_IG_Attributions.py BRCA 50 0 100")
get_ipython().magic(u"run -i Get_VAE_IG_Attributions.py BRCA 75 0 100")
get_ipython().magic(u"run -i Get_VAE_IG_Attributions.py BRCA 100 0 100")

##STEP 4: Learning ensemble labels 
get_ipython().magic(u"run -i Create_Ensemble_Labels.py BRCA 150")

##STEP 5: Creating DeepProfile ensemble training embedding
get_ipython().magic(u"run -i Create_DeepProfile_Training_Embeddings.py BRCA")

##STEP 6: Creating DeepProfile ensemble gene attribution matrices
get_ipython().magic(u"run -i Create_DeepProfile_Ensemble_Weights.py BRCA")

