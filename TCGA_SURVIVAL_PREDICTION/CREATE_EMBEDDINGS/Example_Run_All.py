###############################
#Example for training TCGA models for a cancer type
###############################
import sys

##STEP 1: Preprocessing Data
get_ipython().magic(u"run -i Preprocess_TCGA_Rnaseq_Expression.py BRCA BRCA")
get_ipython().magic(u"run -i Create_TCGA_Rnaseq_PCs.py BRCA BRCA")

##STEP 2: Encoding Expression with DeepProfile
get_ipython().magic(u"run -i Create_All_VAE_Embeddings.py BRCA BRCA")
get_ipython().magic(u"run -i Create_DeepProfile_TCGA_Embeddings.py BRCA BRCA")

##STEP 3: Encoding Expression with Competitor Models
get_ipython().magic(u"run -i Encode_TCGA_Data_with_PCA.py BRCA BRCA")
get_ipython().magic(u"run -i Encode_TCGA_Data_with_ICA.py BRCA BRCA")
get_ipython().magic(u"run -i Encode_TCGA_Data_with_RP.py BRCA BRCA")
get_ipython().magic(u"run -i Encode_TCGA_Data_with_AE.py BRCA BRCA")
get_ipython().magic(u"run -i Encode_TCGA_Data_with_DAE.py BRCA BRCA")