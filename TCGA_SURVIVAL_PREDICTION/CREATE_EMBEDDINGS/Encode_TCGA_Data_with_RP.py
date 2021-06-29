###############################
#Author: Ayse Dincer
#Script for creating TCGA RNA-Seq RP embeddings
###############################

import numpy as np
import pandas as pd
import csv
from sklearn.random_projection import GaussianRandomProjection
import sklearn.preprocessing
from scipy.stats.mstats import winsorize

#Read cancer type and TCGA type
import sys
cancer_type = sys.argv[1]
tcga_type = sys.argv[2]
print("CANCER NAME: " + cancer_type)
print("TEST NAME: " + tcga_type)

input_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/TCGA_FILES/' 

#Read training data
data_df = pd.read_table(input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
print("Training data ", data_df.shape)
training_data = data_df.values
training_data = np.nan_to_num(training_data)

#Read TCGA RNA-Seq expression data
tcga_df = pd.read_table(input_folder + 'TCGA_FILES/TCGA_' + tcga_type + '_PREPROCESSED_RNASEQ_EXPRESSION.tsv', index_col= 0)
print("TCGA data ", tcga_df.shape)
test_data = tcga_df.values

#Train all models
for run in range(10):
    #Train model
    transformer = GaussianRandomProjection(n_components = 150, random_state = run * 12345)
    transformer.fit(training_data)
    components = transformer.components_
    print("RP components ", components.shape)

    #Save the encoded data
    encoded_data = transformer.transform(test_data)
    print("Encoded TCGA data ", encoded_data.shape)
    encoded_df = pd.DataFrame(encoded_data, index = tcga_df.index)
    encoded_df.to_csv(output_folder + 'TCGA_RNASEQ_' + tcga_type + '_RP_150L_' + str(run + 1) + '.tsv', sep = '\t')
    
