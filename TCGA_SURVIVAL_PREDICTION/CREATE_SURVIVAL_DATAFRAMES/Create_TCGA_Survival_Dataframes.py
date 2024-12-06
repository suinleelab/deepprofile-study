###############################
#Script for creating TCGA survival dataframes
###############################

import numpy as np
import pandas as pd
import math

#Method for defining the survival dataframe
def createSurvivalDF(cancer_type, tcga_type):
    
    input_folder = '../TCGA_DATA/TCGA_CLINICAL_DATA/'
    output_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/' + 'TCGA_FILES/'
    
    #Read clinical data
    survival_df = pd.read_table( input_folder + tcga_type + '.clin.merged.picked.txt', index_col = 0)
    survival_df = survival_df.transpose()
    print("TCGA clinical dataframe ", survival_df.shape)
    print("TCGA clinical dataframe ", survival_df.columns)
    
    #Extract vital status, days to death, and days to follow up
    vital_df = survival_df['vital_status']
    dead_df = survival_df['days_to_death']
    alive_df = survival_df['days_to_last_followup']
    
    #Create joined arrays
    vital_status_array = []
    days_status_array = []
    for i in range(vital_df.shape[0]):
        if int(vital_df.values[i])== 0:
            vital_status_array.append(False)
            days_status_array.append(alive_df.values[i])
        else:
            vital_status_array.append(True)
            days_status_array.append(dead_df.values[i])


    #Create joined dataframe
    vital_status_df = pd.DataFrame(vital_status_array, index = survival_df.index, columns = ['Status'])
    days_status_df = pd.DataFrame(days_status_array, index = survival_df.index, columns = ['Survival_in_days'])
    joined_df = pd.concat([vital_status_df, days_status_df], axis = 1)
    print("TCGA survival dataframe ", joined_df)
    joined_df.to_csv(output_folder + '/TCGA_' + tcga_type + '_Survival_df.tsv')
    
    
import sys
cancer_type = sys.argv[1]
tcga_type = sys.argv[2]

createSurvivalDF(cancer_type, tcga_type)
