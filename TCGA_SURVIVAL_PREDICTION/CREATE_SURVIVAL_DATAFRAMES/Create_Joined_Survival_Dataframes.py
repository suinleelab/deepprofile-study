###############################
#Script for creating joined TCGA survival dataframes and DeepProfile embeddings
###############################

import pandas as pd
import numpy as np
import sys

def createJoinedDf(tcga_type, cancer_type):
    
    input_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/' + 'TCGA_FILES/'
    output_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/' + 'TCGA_FILES/'
    
    #Read survival dataframe
    surv_df = pd.read_table(input_folder + 'TCGA_' + tcga_type + '_Survival_df.tsv', index_col = 0, sep = ',')
    surv_df = surv_df.astype(float)
    print("Survival dataframe ", surv_df.shape)
    
    #Drop nan samples
    indices_to_drop1 = np.where(np.isnan(surv_df.values))[0]
    indices_to_drop2 = np.where(surv_df['Survival_in_days'].values <= 0)[0]
    indices_to_drop = np.unique(np.concatenate((indices_to_drop1, indices_to_drop2)))
    surv_df = surv_df.drop(surv_df.index[indices_to_drop])
    surv_df = pd.DataFrame(surv_df.values, index = surv_df.index, columns = ['fustat', 'futime'])
    print("Survival dataframe ", surv_df.shape)
    
    #Read DeepProfile embedding
    data_df = pd.read_table(input_folder  + tcga_type + '_DeepProfile_TCGA_RNASeq_Embedding_150L.tsv', index_col = 0)
    print("DeepProfile embedding ", data_df.shape)
    
    #Match sample indices
    surv_df_sample_names = surv_df.index
    data_df_sample_names = data_df.index
    print("Surv samples ", surv_df_sample_names)
    print("Data samples ", data_df_sample_names)
    
    new_indices = [s.upper() for s in surv_df.index]
    surv_df = pd.DataFrame(surv_df.values, index = new_indices, columns = surv_df.columns)
    
    new_columns = ['Node ' + str(i) for i in range(1, 151)]
    new_indices = [s[:12] for s in data_df.index]
    data_df = pd.DataFrame(data_df.values, index = new_indices, columns = new_columns)
    
    surv_df_sample_names = surv_df.index
    data_df_sample_names = data_df.index
    print("Surv samples ", surv_df_sample_names)
    print("Data samples ", data_df_sample_names)
    
    #Take the samples available in both datasets
#     intersect_indices = np.intersect1d(data_df.index, surv_df.index)
#     print("Common indices ", intersect_indices)
    
    #Create joined dataframe
    joined_df = data_df.merge(surv_df, left_index=True, right_index=True)
    joined_df = joined_df.sort_index()
    joined_df = joined_df.loc[~joined_df.index.duplicated(keep='first')]
    print("Joined dataframe ", joined_df.shape)
    print(joined_df)
    joined_df.to_csv(output_folder + '/DeepProfile_Embedding_and_' + tcga_type + '_Survival_df.tsv', sep = '\t')
    
#Read cancer types
cancer_type = sys.argv[1]
tcga_type = sys.argv[2]

createJoinedDf(tcga_type, cancer_type)
