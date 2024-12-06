###############################
#Example for generating healthy tissue embeddings for a cancer type
###############################

#Preprocess data
get_ipython().magic(u"run -i Preprocess_Gtex_Rnaseq_Expressions.py BRCA")
get_ipython().magic(u"run -i Create_Gtex_Rnaseq_PCs.py BRCA")

#Create DeepProfile embeddings
get_ipython().magic(u"run -i Encode_GTEX_Data_with_VAE.py BRCA 5 0 100")
get_ipython().magic(u"run -i Encode_GTEX_Data_with_VAE.py BRCA 10 0 100")
get_ipython().magic(u"run -i Encode_GTEX_Data_with_VAE.py BRCA 25 0 100")
get_ipython().magic(u"run -i Encode_GTEX_Data_with_VAE.py BRCA 50 0 100")
get_ipython().magic(u"run -i Encode_GTEX_Data_with_VAE.py BRCA 75 0 100")
get_ipython().magic(u"run -i Encode_GTEX_Data_with_VAE.py BRCA 100 0 100")

get_ipython().magic(u"run -i Create_DeepProfile_GTEX_Embeddings.py BRCA")

#Train healthy tissue classifiers
get_ipython().magic(u"run -i Normal_Tissue_Classifier.py BRCA")