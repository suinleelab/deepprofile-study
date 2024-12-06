###############################
#Script for creating PCs for expression matrices
###############################

import numpy as np
import pandas as pd
import csv
from sklearn.decomposition import PCA
import sklearn.preprocessing
import statsmodels.api as sm
from sklearn.preprocessing import scale

def createData(cancer_type):
    
    input_folder ='../ALL_CANCER_FILES/' + cancer_type + '/'
    
    #Read training data
    data_df = pd.read_table(input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
    print("Training expression dataframe ", data_df.shape)

    training_data = data_df.values
    training_data = np.nan_to_num(training_data)

    #Train PCA models
    pca = PCA(n_components = 1000)
    pca.fit(training_data)
    components = pca.components_
    print("PCA Components ", components.shape)

    #Read GTEX expression dataframe
    test_df = pd.read_table(input_folder + 'HEALTHY_TISSUE_FILES/' + 'GTEX_' + cancer_type + '_PREPROCESSED_RNASEQ_EXPRESSION.tsv', sep = '\t', index_col=0)
    print("Gtex expression dataframe ", test_df.shape)
    
    #Get genes available in training dataset
    joined_df = pd.concat([data_df, test_df], sort=False, join = 'outer')
    joined_df = joined_df[data_df.columns]
    joined_df = joined_df.iloc[-1 * test_df.shape[0]:, :]
    test_df = joined_df
    
    print("Gtex expression dataframe ", test_df.shape)
    
    #Encode test data using trained PCA model
    test_df = test_df.fillna(test_df.mean().fillna(0))
    test_data = test_df.values

    #Save the encoded data
    encoded_data = pca.transform(test_data)
    encoded_df = pd.DataFrame(encoded_data, index = test_df.index)
    print("GTEX PCA data ", encoded_df.shape)
    print("GTEX PCA data ", encoded_df.head)
    encoded_df.to_csv(input_folder + '/HEALTHY_TISSUE_FILES/GTEX_' + cancer_type + '_DATA_1K_PCs.tsv', sep = '\t')

import sys

cancer_type = sys.argv[1]
createData(cancer_type)
