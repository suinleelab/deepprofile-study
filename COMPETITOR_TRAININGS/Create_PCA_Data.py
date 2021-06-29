###############################
#Author: Ayse Dincer
#Script for training PCA models
###############################

import numpy as np
import pandas as pd
import csv
from sklearn.decomposition import PCA
import sys

#Read cancer type
cancer_type = sys.argv[1]

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/PCA_FILES/'

#Method for defining PCs for training data
def createData(cancer_type):
    
    L = 150
    print("Number of latent nodes " + str(L))
    
    data_df = pd.read_table(input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
    print("Training data ", data_df.shape)

    training_data = data_df.values
    training_data = np.nan_to_num(training_data)

    #Fit PCA model
    pca = PCA(n_components = L)
    pca.fit(training_data)
    components = pca.components_
    print("PCA components ", components.shape)

    #Save the learned components
    component_df = pd.DataFrame(components.T, index = data_df.columns)
    component_df.to_csv(output_folder + cancer_type + '_DATA_TOP2_JOINED_PCA_COMPONENTS_' + str(L) + 'L.tsv', sep = '\t')


createData(cancer_type)