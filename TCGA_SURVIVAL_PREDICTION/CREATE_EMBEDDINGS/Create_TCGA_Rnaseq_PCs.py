###############################
#Author: Ayse Dincer
#Script for recording top PCs or TCGA RNA-Seq data
###############################

import numpy as np
import pandas as pd
import csv
from sklearn.decomposition import PCA
import sklearn.preprocessing

#Define method for preprocessing data
def create_data(cancer_type, tcga_type):

    input_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/'
    output_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/' + 'TCGA_FILES/'
    
    #Read training data
    data_df = pd.read_table(input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
    print("Training data ", data_df.shape)

    #Apply PCA to training data
    training_data = data_df.values
    training_data = np.nan_to_num(training_data)
    
    pca = PCA(n_components = 1000)
    pca.fit(training_data)
    components = pca.components_
    print("PCA components ", components.shape)
           
    #Read TCGA RNA-Seq expression
    tcga_df = pd.read_table(output_folder + 'TCGA_' + tcga_type + '_PREPROCESSED_RNASEQ_EXPRESSION.tsv', index_col= 0)
    print("TCGA expression ", tcga_df.shape)
    print('RANGE: ', (np.max(tcga_df.values) - np.min(tcga_df.values) ))
    
    #Encode test data using trained PCA model
    test_data = tcga_df.values
    encoded_data = pca.transform(test_data)
    print("Encoded TCGA data ", encoded_data.shape)
        
    #Record expression data
    encoded_df = pd.DataFrame(encoded_data, index = tcga_df.index)
    encoded_df.to_csv(output_folder + 'TCGA_RNASEQ_' + tcga_type + '_PCA_1000L.tsv', sep = '\t')

    
import sys
cancer_type = sys.argv[1]
tcga_type = sys.argv[2]
create_data(cancer_type, tcga_type)
