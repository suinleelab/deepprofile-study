###############################
#Author: Ayse Dincer
#Script for creating TCGA RNA-Seq PCA embeddings
###############################

import numpy as np
import pandas as pd
import csv
from sklearn.decomposition import PCA
import sklearn.preprocessing
from scipy.stats.mstats import winsorize
import sys

#Read cancer type and TCGA type
cancer_type = sys.argv[1]
tcga_type = sys.argv[2]
print("CANCER NAME: " + cancer_type)
print("TEST NAME: " + tcga_type)

input_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/TCGA_FILES/' 

#Read training data
data_df = pd.read_table( input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
print("Training data ", data_df.shape)
training_data = data_df.values
training_data = np.nan_to_num(training_data)

#Train PCA model
pca = PCA(n_components = 150)
pca.fit(training_data)
components = pca.components_
print("PCA components ", components.shape)

#Read TCGA RNA-Seq expression data
tcga_df = pd.read_table(output_folder+ '/TCGA_' + tcga_type + '_PREPROCESSED_RNASEQ_EXPRESSION.tsv', index_col= 0)
print("TCGA data ", tcga_df.shape)

#Encode TCGA data with PCA model
test_data = tcga_df.values
encoded_data = pca.transform(test_data)
print("Encoded TCGA data ", encoded_data.shape)

encoded_df = pd.DataFrame(encoded_data, index = tcga_df.index)
encoded_df.to_csv(output_folder + '/TCGA_RNASEQ_' + tcga_type + '_PCA_150L.tsv', sep = '\t')
