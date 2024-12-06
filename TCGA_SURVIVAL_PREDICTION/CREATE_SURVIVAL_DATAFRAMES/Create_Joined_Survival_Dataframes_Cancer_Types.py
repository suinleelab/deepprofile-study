###############################
#Script for creating joined cancer types TCGA survival dataframes 
###############################

import pandas as pd
import numpy as np

#Method for combining datasets
def create_Data(cancer):
    
    input_folder = '../../ALL_CANCER_FILES/' + cancer + '/' + 'TCGA_FILES/'
    output_folder = '../../ALL_CANCER_FILES/' + cancer + '/' + 'TCGA_FILES/'
    
    c = np.where(np.asarray(cancer_types) == cancer)[0][0]
    df_list = []
    for test in test_cases[c]:
        print("TCGA type ", test)
        surv_df = pd.read_table(input_folder + '/DeepProfile_Embedding_and_' + test + '_Survival_df.tsv', sep = '\t', index_col = 0)
        print("Survival dataframe ", surv_df.shape)
        df_list.append(surv_df)

    #Combine dataframes
    joined_df = pd.concat(df_list)
    print("Joined survival dataframe ", joined_df.shape)
    joined_df.to_csv(output_folder + '/DeepProfile_Embedding_and_' + cancer + '_Survival_df.tsv', sep = '\t')

cancer_types = ['LUNG']
test_cases = [ ['LUAD', 'LUSC']]
  
for i in range(len(cancer_types)):
    print("Cancer type ", cancer_types[i])
    create_Data(cancer_types[i])

#Method for combining survival dataframes
def create_Data(cancer):
    
    input_folder = '../../ALL_CANCER_FILES/' + cancer + '/' + 'TCGA_FILES/'
    output_folder = '../../ALL_CANCER_FILES/' + cancer + '/' + 'TCGA_FILES/'
    
    c = np.where(np.asarray(cancer_types) == cancer)[0][0]
    df_list = []
    for test in test_cases[c]:
        print("TCGA type ", test)
        surv_df = pd.read_table(input_folder + 'TCGA_' + test + '_Survival_df.tsv', sep = ',', index_col = 0)
        print("Survival dataframe ", surv_df.shape)
        df_list.append(surv_df)

    #Combine dataframes
    joined_df = pd.concat(df_list)
    print("Joined survival dataframe ", joined_df.shape)
    joined_df.to_csv(output_folder + 'TCGA_' + cancer + '_Survival_df.tsv', sep = '\t')

cancer_types = ['LUNG']
test_cases = [ ['LUAD', 'LUSC']]

for i in range(len(cancer_types)):
    print("Cancer type ", cancer_types[i])
    create_Data(cancer_types[i])

   