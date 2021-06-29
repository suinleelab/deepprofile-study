###############################
#Author: Ayse Dincer
#Script for predicting survival status of patients
###############################

import numpy as np
import pandas as pd
from scipy import stats
from sklearn import metrics
import random

from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import KFold
from sklearn import linear_model
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.metrics import accuracy_score
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc

#Define method for training models
def trait_classification_accuracy(X, Y):
    
    #Do cross validation
    loo = KFold(20, shuffle = True, random_state = 123456)
    
    predictions = np.zeros(X.shape[0])
    probabilities = np.zeros(X.shape[0])
    
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        Y_train, Y_test = Y[train_index], Y[test_index]

        #Normalize training data
        scaler = StandardScaler()
        scaler.fit(X_train)

        X_std = scaler.transform(X_train)
        X_test_std = scaler.transform(X_test)
     
        #Tune parameters
        tuned_parameters = [{'C': [0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 100]}]
        
        model = LogisticRegression(random_state=12345, penalty = 'l1', max_iter=1000, 
                                   solver = 'liblinear')
        clf = GridSearchCV(model, tuned_parameters, cv = 5, scoring='roc_auc', n_jobs = -1)
        clf.fit(X_std, Y_train)
        
        #Record predictions and probabilities
        predicted_Y = clf.predict(X_test_std)
        predictions[test_index] = predicted_Y
        
        probs = clf.predict_proba(X_test_std)

        probabilities[test_index] = probs[:, 1]
        

    #Calculate accuracy and ROC-AUC
    accuracy = accuracy_score(Y, predictions)
    score = roc_auc_score(Y, probabilities)
    
    return [accuracy, score]


#Define method for predicting survival
def predict_survival(cancer_type, tcga_types, method, run_index, seed):
    print(tcga_types)
    accuracies = []
    aucs = []
    
    df_list = []
    for tcga_type in tcga_types:
        if method == 'PCA':
            X_df  = pd.read_table(input_folder + 'TCGA_RNASEQ_' + tcga_type + '_PCA_150L.tsv', index_col = 0)

        if method == 'ICA' or method == 'RP':
            X_df  = pd.read_table(input_folder + 'TCGA_RNASEQ_' + tcga_type + '_' + method + '_150L_' + str(run_index) + '.tsv', index_col = 0)

        if method == 'AE' or method == 'DAE':
            X_df  = pd.read_table(input_folder + 'TCGA_' + tcga_type + '_RNASeq_Expression_' + method + '_encoded_150L_' + str(run_index) + '.tsv', index_col = 0)
         
        if method == 'VAE':
            X_df  = pd.read_table(input_folder + 'TCGA_' + tcga_type + '_RNASeq_Expression_' + method + '_encoded_' + VAE_dim + 'L_' + str(run_index) + '.tsv', index_col = 0)
            
        if method == 'DeepProfile':
            X_df  = pd.read_table(input_folder + tcga_type + '_DeepProfile_TCGA_RNASeq_Embedding_150L.tsv', index_col = 0, sep = '\t')


        df_list.append(X_df)

    X_df = pd.concat(df_list, axis = 0)
    print("Expression data joined ", X_df.shape)

    #Now, replace X_df index to match with Y_df index
    mapper = lambda t: t[:12]
    vfunc = np.vectorize(mapper)
    newX_index = vfunc( X_df.index)

    X_df = pd.DataFrame(X_df.values, index = newX_index, columns = X_df.columns)

    #Take intersecting samples in datasets
    X_samples = X_df.index
    Y_samples = Y_df.index
    intersecting_samples = np.intersect1d(X_samples, Y_samples)

    subX_df = X_df.T[intersecting_samples].T
    subY_df = Y_df.T[intersecting_samples].T

    print("X dataframe ", subX_df.shape)
    print("Y dataframe ", subY_df.shape)

    print("X dataframe ", subX_df.index)
    print("Y dataframe ", subY_df.index)

    sample_indices = [np.where(subY_df.values == False)[0], np.where(subY_df.values == True)[0]]
    sample_counts = [len(sample_indices[0]), len(sample_indices[1])]
    print("SAMPLES LABEL 0: ", sample_counts[0], " SAMPLES LABEL 1: ", sample_counts[1])

    #Now select the class with highest number of samples and subsample
    low_class = np.argmin(sample_counts)
    high_class = np.argmax(sample_counts)
    print("Lower class size ", sample_counts[low_class], "samples subsampled from class ", high_class)
    random.seed(12345 * seed)
    random_indices = random.sample(list(np.arange(0, sample_counts[high_class])), sample_counts[low_class])
    selected_indices = np.sort(sample_indices[high_class][random_indices])

    subX_df = pd.concat([subX_df.iloc[sample_indices[low_class]], subX_df.iloc[selected_indices]])
    subY_df = pd.concat([subY_df.iloc[sample_indices[low_class]], subY_df.iloc[selected_indices]])                           
    subX_df = subX_df.sort_index()
    subY_df = subY_df.sort_index()

    results = trait_classification_accuracy(subX_df.values, np.ravel(subY_df.values))
    return results

#Read user inputs
import sys
cancer_type = sys.argv[1] #main cancer type
method = sys.argv[2] #name of the method
run_index = int(sys.argv[3]) #run index
if len(sys.argv) > 4:
    VAE_dim = sys.argv[4]

if cancer_type == 'LUNG':
    tcga_types = ['LUAD', 'LUSC']

input_folder = '../../ALL_CANCER_FILES/' + cancer_type + '/TCGA_FILES/'
output_folder = 'Prediction_Results/'

#Join data for cancer subtypes
df_list = []
for tcga_type in tcga_types:
    Y_df  = pd.read_table(input_folder + 'DeepProfile_Embedding_and_' + tcga_type + '_Survival_df.tsv', index_col = 0, sep = '\t')
    print("Survival dataframe ", Y_df.shape)
    df_list.append(Y_df)
    
Y_df = pd.concat(df_list, axis = 0)
print("JOINED survival dataframe ", Y_df.shape)

print("ALIVE..")
print( Y_df[Y_df['fustat']  == 0]['futime'])
print( np.mean(Y_df[Y_df['fustat']  == 0]['futime'].values))

#Select all dead patients, only if they died within a year
Y_df_dead = Y_df[Y_df['fustat']  == 1]
indices_dead  = np.where(Y_df_dead['futime'] < 5 * 365)[0]
print("Dead within 5 year ", Y_df_dead.iloc[indices_dead][['fustat', 'futime']])
print("Dead within 5 year  ", np.max(Y_df_dead.iloc[indices_dead]['futime']))

#Select all alive patients, only if they lived more than a year
indices_alive  = np.where(Y_df['futime'] > 5 * 365)[0]
print("Alive after 5 year ", Y_df.iloc[indices_alive][['fustat', 'futime']])
print("Alive after 5 year ", np.min(Y_df.iloc[indices_alive]['futime']))

indices = list(indices_dead) + list(indices_alive)
indices = np.unique(indices)
Y_df = Y_df['fustat']
Y_df = Y_df.iloc[indices]
Y_df = Y_df.dropna()
print("Survival dataframe \n ", Y_df)

class0_count = len(np.where(Y_df.values == 0)[0])
class1_count = len(np.where(Y_df.values == 1)[0])

all_accuracies = [] 
all_aucs = [] 

for sampling_index in range(50):
    result = predict_survival(cancer_type, tcga_types, method, run_index, sampling_index)
    print("Accuracy: ", result) 
    print("ROC-AUC: ", result[1]) 
    all_accuracies.append(result[0])
    all_aucs.append(result[1])

print("FINAL RESULTS...")
print("MEAN ACCURACY ", np.mean(np.asarray(all_accuracies)))
print("MEAN ROC-AUC ", np.mean(np.asarray(all_aucs)))

#Save results to a file
if method == 'VAE':
    np.savetxt(output_folder + cancer_type + '/TCGA_Survival_5year_LR_Balanced_Subsample_20FOLD_50Runs_' + cancer_type + '_' + method + '_' + VAE_dim + 'L_' + str(run_index) + '_ACCs.txt', np.asarray(all_accuracies), delimiter='\n')
    np.savetxt(output_folder + cancer_type + '/TCGA_Survival_5year_LR_Balanced_Subsample_20FOLD_50Runs_' + cancer_type + '_' + method + '_' + VAE_dim + 'L_' + str(run_index) + '_AUCs.txt', np.asarray(all_aucs), delimiter='\n')

else:
    np.savetxt(output_folder + cancer_type + '/TCGA_Survival_5year_LR_Balanced_Subsample_20FOLD_50Runs_' + cancer_type + '_' + method + '_' + str(run_index) + '_ACCs.txt', np.asarray(all_accuracies), delimiter='\n')
    np.savetxt(output_folder + cancer_type + '/TCGA_Survival_5year_LR_Balanced_Subsample_20FOLD_50Runs_' + cancer_type + '_' + method + '_' + str(run_index) + '_AUCs.txt', np.asarray(all_aucs), delimiter='\n')
    

    
    
