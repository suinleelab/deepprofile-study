###############################
#Author: Ayse Dincer
#Script for training classifiers for separating healthy and cancer tissue embeddings
###############################

import pandas as pd
import seaborn as sb
import numpy as np
import pickle
import random
from tqdm import *
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.utils import resample
from sklearn.preprocessing import StandardScaler


def trainClassifier(cancer_type):

    input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
    output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/HEALTHY_TISSUE_FILES/'

    #Read cancer embedding
    cancer_data = pd.read_csv(input_folder + cancer_type + '_DeepProfile_Training_Embedding_150L.tsv',sep='\t',index_col=0)
    print("Cancer embedding ", cancer_data.shape)
    
    #Read GTEX embedding
    healthy_data = pd.read_csv(input_folder + 'HEALTHY_TISSUE_FILES/' + cancer_type + '_DeepProfile_GTEX_Healthy_Tissue_Embedding_150L.tsv',sep='\t',index_col=0)
    print("GTEX embedding ", healthy_data.shape)
    
    #Combine datasets
    FULL_FRAME = pd.concat([cancer_data, healthy_data],axis=0)

    #Define healthy tissue labels
    healthy_label = [x < cancer_data.shape[0] for x in range(FULL_FRAME.shape[0])]

    #Train 100 L2 models with bootstrapping
    bootstrap_weights = []
    for i in tqdm(range(500)):
        X_re,y_re = resample(FULL_FRAME,healthy_label, random_state = 1234 * i)
        clf = LogisticRegression(penalty = 'l2', solver = 'liblinear')
        clf.fit(X_re,y_re)

        bootstrap_weights.append(clf.coef_)
    
    #Save the results
    pickle.dump(bootstrap_weights,open(output_folder + 'bootstrap_' + cancer_type + '_weights.p','wb'))

import sys

cancer_type = sys.argv[1]
trainClassifier(cancer_type)
