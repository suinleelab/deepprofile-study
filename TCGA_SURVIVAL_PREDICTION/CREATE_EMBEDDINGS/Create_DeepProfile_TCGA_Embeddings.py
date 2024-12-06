###############################
#Script for creating TCGA RNA-Seq DeepProfile embeddings
###############################

import pandas as pd
import numpy as np
import csv
import sys

#Read cancer type from user
cancer_type = sys.argv[1]
tcga_type = sys.argv[2]

#Read all VAE embeddings
dims  = [5, 10, 25, 50, 75, 100]
run = 100

input_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/TCGA_FILES/' 

data_list = []
for dim in dims:
    for i in range(run):
        data_df = pd.read_table(input_folder + 'TCGA_FILES/TCGA_' + tcga_type + '_RNASeq_Expression_VAE_encoded_' + str(dim) + 'L_' + str(i) + '.tsv', index_col = 0)   
        print("TCGA VAE embedding ", data_df.shape)
        data_list.append(data_df.values)

#Concatenate all embeddings     
joined_data = np.concatenate(data_list, axis=1)
print("Joined VAE embedding ",joined_data.shape)

#Read DeepProfile ensemble labels
L = 150
labels_df = pd.read_table( input_folder + cancer_type + '_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_' + str(L) + 'L.txt', header= None)
labels = labels_df.values
print("DeepProfile ensemble labels ", len(labels))

#Create ensemble embedding
ensemble_embeddings = np.zeros((joined_data.shape[0], L))
for label in range(L):
    indices = np.where(labels == label)[0]
    average_values = np.mean(joined_data[:, indices], axis = 1)
    ensemble_embeddings[:, label] = average_values

#Record the ensemble embeddings
print(ensemble_embeddings.shape)
ensemble_embeddings_df = pd.DataFrame(ensemble_embeddings, index = data_df.index, columns = np.arange(L))
ensemble_embeddings_df.to_csv(output_folder + tcga_type + '_DeepProfile_TCGA_RNASeq_Embedding_' + str(L) + 'L.tsv', sep = '\t')