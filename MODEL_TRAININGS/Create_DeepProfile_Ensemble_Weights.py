###############################
#This script is for creating gene attribution matrices for DeepProfile
###############################

import numpy as np
import pandas as pd
import csv
import sys

#Read user input
cancer_type = sys.argv[1]

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/' 

#Read all VAE model gene attributions
L = 150
data_df = pd.read_table(input_folder + 'VAE_WEIGHTS/' + cancer_type + '_DATA_VAE_Cluster_Weights_TRAINING_' + str(100) + 'L_fold' + str(1) + '.tsv', index_col = 0)
print(data_df.shape)
basic_length = data_df.shape[0]

weight_list = []
dims = [5, 10, 25, 50, 75, 100]
run_count = 100
for dim in dims:
    VAE_weights = np.zeros((run_count * dim, basic_length))
    for i in range(run_count):
        data_df = pd.read_table(input_folder + 'VAE_WEIGHTS/' + cancer_type + '_DATA_VAE_Cluster_Weights_TRAINING_' + str(dim) + 'L_fold' + str(i) + '.tsv', index_col = 0)
        data_df = data_df.T
        #print(data_df.shape)
        start = dim * i
        end = dim * (i + 1)
        VAE_weights[start:end, :] = data_df.values
    weight_list.append(VAE_weights)

#Read the ensemble labels
labels_df = pd.read_table(input_folder + cancer_type + '_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_' + str(L) + 'L.txt', header= None)
labels = labels_df.values
print("Ensemble labels ", len(labels))

#Concatenate all the gene attributions
joined_weights = np.concatenate(weight_list)
print("Joined weights ", joined_weights.shape)

#Create ensemble weights
ensemble_weights = np.zeros((L, joined_weights.shape[1]))
for label in range(L):
    indices = np.where(labels == label)[0]
    average_weights = np.mean(joined_weights[indices, :], axis = 0)
    ensemble_weights[label, :] = average_weights

print("Ensemble weights ", ensemble_weights.shape)

#Record ensemble weights
ensemble_weights_df = pd.DataFrame(ensemble_weights, index = np.arange(L), columns = data_df.columns)
ensemble_weights_df.to_csv(output_folder + cancer_type + '_DeepProfile_Ensemble_Gene_Importance_Weights_' + str(L) + 'L.tsv', sep = '\t')