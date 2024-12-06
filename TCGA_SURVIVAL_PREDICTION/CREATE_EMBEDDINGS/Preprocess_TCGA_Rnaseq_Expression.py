###############################
#Script for preprocessing TCGA RNA-Seq expression
###############################

import numpy as np
import pandas as pd
import sklearn.preprocessing

#Define method for preprocessing data
def create_data(cancer_type, tcga_type):

    input_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/'
    output_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/' + 'TCGA_FILES/'
    
    #Read TCGA RNA-Seq expression
    print("Reading expression...")
    tcga_df = pd.read_table('../TCGA_DATA/TCGA_RNASEQ/' + tcga_type + '.uncv2.mRNAseq_RSEM_normalized_log2.txt', sep = '\t', index_col= 0)
    tcga_df = tcga_df.transpose()
    print("TCGA expression ", tcga_df.shape)
    print("TCGA expression ", tcga_df.head())
    
    #Map to gene names and eliminate unknown genes
    print("Correcting gene names...")
    gene_names = tcga_df.columns
    gene_names = [n[:n.index('|')] for n in gene_names]
    tcga_df = pd.DataFrame(tcga_df.values, index = tcga_df.index, columns = gene_names)
 
    #Eliminate unknown genes
    tcga_df = tcga_df.iloc[:, tcga_df.columns.values != '?']
    tcga_df = tcga_df.loc[:,~tcga_df.columns.duplicated()]
    print("TCGA expression ", tcga_df.shape)
    print("TCGA expression ", tcga_df)
   
    #Eliminate normal samples
    print("Eliminating normal samples..")
    sample_codes = [s[-2:] for s in  tcga_df.index]
    print("Sample codes ", np.unique(sample_codes))
    normal_codes = [s[-2] for s in  tcga_df.index]
    cancer_samples = np.where(np.asarray(normal_codes) == '0')[0]
    print("Total number of samples ", len(tcga_df.index))
    print("Total number of cancer samples ", len(cancer_samples))
    tcga_df = tcga_df.iloc[cancer_samples, :]
    print("TCGA expression ", tcga_df.shape)
    print("TCGA expression cancer samples ", tcga_df.index)
    
    #Read training data
    print("Combining with training data...")
    data_df = pd.read_table(input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', sep = '\t', index_col=0)
    print("Training data ", data_df.shape)

    #Get only training genes from the expression data
    joined_df = pd.concat([data_df, tcga_df], join = 'outer')
    joined_df = joined_df[data_df.columns]
    joined_df = joined_df.iloc[-1 * tcga_df.shape[0]:, :]
    joined_df = joined_df.fillna(joined_df.mean().fillna(0))
    print("TCGA expression ", joined_df.shape)
    print("TCGA expression ", joined_df.head())
    
    #Standardize data to make 0 mean univariate
    print("Standardizing the data...")
    scaled_expression_values = joined_df.values
    normalized_data = sklearn.preprocessing.scale(scaled_expression_values)
    print("TCGA expression mean ", np.mean(normalized_data, axis = 0))
    print("TCGA expression mean ", len(np.mean(normalized_data, axis = 0)))
    print("TCGA expression std ", np.std(normalized_data, axis = 0))
    print("TCGA expression std ", len(np.std(normalized_data, axis = 0)))
    
    #Record joined dataframe
    joined_df = pd.DataFrame(normalized_data, index = joined_df.index, columns = joined_df.columns)
    print("Final dataframe ", joined_df.shape)
    print('RANGE: ', (np.max(joined_df.values) - np.min(joined_df.values) ))
    
    #Record expression data
    joined_df.to_csv(output_folder + 'TCGA_' + tcga_type + '_PREPROCESSED_RNASEQ_EXPRESSION.tsv', sep = '\t')

    
import sys
cancer_type = sys.argv[1]
tcga_type = sys.argv[2]
create_data(cancer_type, tcga_type)
