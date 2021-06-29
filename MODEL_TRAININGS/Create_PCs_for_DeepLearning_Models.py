###############################
#Author: Ayse Dincer
#This script is for PCA transforming the input data to pass to deep learning models
###############################

import numpy as np
import pandas as pd
import csv
from sklearn.decomposition import PCA
import sys

#Read cancer type input
cancer_type = sys.argv[1]
#Read number of components
component_count = int(sys.argv[2])

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'

#Method for creating PCs
def createPCs(cancer_type):
    
    print("************************* " + cancer_type)
    
    #Read training data
    data_df = pd.read_table(input_folder  + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
    print("Training data ", data_df.shape)
    training_data = data_df.values
    training_data = np.nan_to_num(training_data)

    #Transform training data to top principal components
    pca = PCA(n_components = component_count)
    pca.fit(training_data)
    components = pca.components_
    print("PCA components ", components.shape)

    #Save the learned components
    component_df = pd.DataFrame(components.T, index = data_df.columns)
    component_df.to_csv(output_folder + cancer_type + '_DATA_TOP2_JOINED_PCA_' + str(component_count) + 'L_COMPONENTS.tsv', sep = '\t')

    #Save the encoded data
    encoded_data = pca.transform(training_data)
    print("PCA encoded data ", encoded_data.shape)
    encoded_df = pd.DataFrame(encoded_data, index = data_df.index)
    encoded_df.to_csv(output_folder + cancer_type + '_DATA_TOP2_JOINED_PCA_' + str(component_count) + 'L.tsv', sep = '\t')

createPCs(cancer_type)