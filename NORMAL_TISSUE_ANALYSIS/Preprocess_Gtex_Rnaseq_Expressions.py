###############################
#Script for creating expression matrices for GTEX healthy samples
###############################

import numpy as np
import pandas as pd
import csv
from sklearn.decomposition import PCA
import sklearn.preprocessing
import statsmodels.api as sm
from sklearn.preprocessing import scale

#Read all GTEX expression file
MAIN_df = pd.read_table('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct', sep = '\t', index_col=0)
print("Gtex expression dataframe ", MAIN_df.shape)
all_samples = np.asarray(MAIN_df.columns)

#Method for creating tissue-specific GTEX 
def save_tissue_expression(cancer):
    
    data_folder = '../ALL_CANCER_FILES/' + cancer + '/HEALTHY_TISSUE_FILES/'
    
    #Read sample names of tissue-specific samples
    index_df = pd.read_table(data_folder + 'GTEX_' + cancer + '_SAMPLES.txt', sep = '\n', index_col=0)
    cancer_specific_samples = np.asarray(index_df.index)
    print("Samples ", cancer_specific_samples)

    #Find list of matching samples
    matching_samples = np.intersect1d(cancer_specific_samples, all_samples)
    print("MATCHING SAMPLES COUNT ", len(matching_samples))

    #Get the expression for these patients
    cancer_df = MAIN_df[matching_samples]
    gene_names = MAIN_df['Description'].values
    cancer_df = pd.DataFrame(cancer_df.values.T, index = cancer_df.columns, columns = gene_names)
    print("Samples ", cancer_df.shape)
    print('Range ', (np.max(cancer_df.values) - np.min(cancer_df.values) ))
    
    #Mean impute the missing values
    cancer_df = cancer_df.fillna(cancer_df.mean().fillna(0))
    
    #Log scale the data and make 0-mean univariate
    scaled_expression_values = np.log(cancer_df.values)
    scaled_expression_values[scaled_expression_values == np.NINF] = 0
    normalized_data = sklearn.preprocessing.scale(scaled_expression_values)
    print("Mean values ", np.mean(normalized_data, axis = 0))
    print("Mean values ", len(np.mean(normalized_data, axis = 0)))
    print("Std values ", np.std(normalized_data, axis = 0))
    print("Std values ", len(np.std(normalized_data, axis = 0)))
    
    #Save the final expressiom matrix
    cancer_df = pd.DataFrame(normalized_data, index = cancer_df.index, columns = cancer_df.columns)
    print("Final dataframe ", cancer_df.shape)
    print("Final dataframe ", cancer_df.head())
    print('Final dataframe range: ', (np.max(cancer_df.values) - np.min(cancer_df.values) ))
    
    cancer_df.to_csv(data_folder + 'GTEX_' + cancer + '_PREPROCESSED_RNASEQ_EXPRESSION.tsv', sep = '\t')

import sys

cancer_type = sys.argv[1]
save_tissue_expression(cancer_type)