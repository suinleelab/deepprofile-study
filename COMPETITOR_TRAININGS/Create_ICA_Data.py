###############################
#Author: Ayse Dincer
#Script for training ICA models
###############################

import numpy as np
import pandas as pd
import math
import csv
from sklearn.decomposition import FastICA
import sys

#Read cancer type
cancer_type = sys.argv[1]

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/ICA_FILES/'

L = 150
print("Number of latent nodes " + str(L))
    
data_df = pd.read_table(input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
print("Training data ", data_df.shape)

training_data = data_df.values
training_data = np.nan_to_num(training_data)

#Fit 10 different ICA models with different random seeds
for run in range(10):
    ica = FastICA(n_components = L, random_state = 12345 * run, max_iter = 100000)
    print(ica)

    ica.fit(training_data) 

    components = ica.components_
    print(components.shape)

    #Save the learned components
    component_df = pd.DataFrame(components.T, index = data_df.columns)
    component_df.to_csv(output_folder + cancer_type + '_DATA_TOP2_JOINED_ICA_COMPONENTS_' + str(L) + 'L_fold' + str(run + 1) + '.tsv', sep = '\t')
