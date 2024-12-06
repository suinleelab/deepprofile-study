###############################
#This script is for creating training embeddings
###############################

import pandas as pd
import numpy as np
import csv
import sys

#Read user input
cancer_type = sys.argv[1]

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/' 

#Read all training embeddings
dims  = [5, 10, 25, 50, 75, 100]
data_list = []
for dim in dims:
    run = 100
    for i in range(run):
        data_df = pd.read_table(input_folder + 'VAE_FILES/' + cancer_type + '_DATA_TOP2_JOINED_encoded_' + str(dim) + 'L_TRAINING_fold' + str(i) + '.tsv', index_col = 0)      
        print(data_df.shape)
        data_list.append(data_df.values)

joined_data = np.concatenate(data_list, axis=1)
print("Joined training embeddings" , joined_data.shape)

#Read the ensemble labels
L = 150
labels_df = pd.read_table(input_folder + cancer_type + '_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_' + str(L) + 'L.txt', header= None)
labels = labels_df.values
print("Ensemble labels ", len(labels))

#Create ensemble training embeddings
ensemble_embeddings = np.zeros((joined_data.shape[0], L))
for label in range(L):
    indices = np.where(labels == label)[0]
    average_values = np.mean(joined_data[:, indices], axis = 1)
    ensemble_embeddings[:, label] = average_values

print("Training ensemble embedding ", ensemble_embeddings.shape)

#Save the training embedding
ensemble_embeddings_df = pd.DataFrame(ensemble_embeddings, index = data_df.index, columns = np.arange(L))
ensemble_embeddings_df.to_csv(output_folder + cancer_type + '_DeepProfile_Training_Embedding_' + str(L) + 'L.tsv', sep = '\t')