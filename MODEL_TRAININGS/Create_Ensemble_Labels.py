###############################
#This script is for learning ensemble labels for VAE models
###############################

import pandas as pd
import numpy as np
import csv
import sys
from sklearn.cluster import KMeans

#Read user inputs
cancer_type = sys.argv[1]
final_dim = int(sys.argv[2])
print("FINAL DIM " + str(final_dim))

#Read all training embeddings
dims  = [5, 10, 25, 50, 75, 100]
data_list = []

for dim in dims:
    run = 100
    for i in range(run):
        print(i)
        data_df = pd.read_table('../ALL_CANCER_FILES/' + cancer_type + '/VAE_FILES/' + cancer_type + '_DATA_TOP2_JOINED_encoded_' + str(dim) + 'L_TRAINING_fold' + str(i) + '.tsv', index_col = 0)      
        print(data_df.shape)
        data_list.append(data_df.values)

joined_df = np.concatenate(data_list, axis=1)
print("Joined training embeddings" , joined_df.shape)

#Apply kmeans clustering to this data
X = joined_df

kmeans = KMeans(n_clusters= final_dim, random_state=123).fit(X.transpose())
print("K-means labels ", kmeans.labels_)

#Save labels
np.savetxt('../ALL_CANCER_FILES/' + cancer_type + '/' + cancer_type + '_TRAINING_DATA_kmeans_ENSEMBLE_LABELS_' + str(final_dim) + 'L.txt' , kmeans.labels_, delimiter=',')
