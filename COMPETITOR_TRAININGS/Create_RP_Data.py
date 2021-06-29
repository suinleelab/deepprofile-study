###############################
#Author: Ayse Dincer
#Script for training Random Projection models
###############################

import numpy as np
import pandas as pd
import csv
from sklearn.random_projection import GaussianRandomProjection
import sys

#Read cancer type
cancer_type = sys.argv[1]

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/RP_FILES/'

#Method for defining ICA for training data
def createData(cancer_type):
    
    L = 150
    print("Number of latent nodes " + str(L))
    
    data_df = pd.read_table(input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
    print("Training data ", data_df.shape)

    training_data = data_df.values
    training_data = np.nan_to_num(training_data)

    #Fit 10 different RP models with different random seeds
    for run in range(10):
        transformer = GaussianRandomProjection(n_components = L, random_state = run * 12345)
        transformer.fit(training_data)
        
        components = transformer.components_
        print(components.shape)

        #Save the learned components
        component_df = pd.DataFrame(components.T, index = data_df.columns)
        component_df.to_csv(output_folder + cancer_type + '_DATA_TOP2_JOINED_RP_COMPONENTS_fold' + str(run + 1) + '.tsv', sep = '\t')
        

createData(cancer_type)