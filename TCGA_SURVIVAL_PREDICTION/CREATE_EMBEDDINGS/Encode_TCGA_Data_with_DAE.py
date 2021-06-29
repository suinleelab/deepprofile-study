###############################
#Author: Ayse Dincer
#Script for creating TCGA RNA-Seq AE embeddings
###############################

import os
import numpy as np
import pandas as pd

import math 
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
from keras.models import model_from_json
from sklearn import preprocessing

import tensorflow as tf
from keras.layers import Input, Dense, Lambda, Layer, Activation
from keras.layers.normalization import BatchNormalization
from keras.models import Model
from keras import backend as K
from keras import metrics, optimizers
from keras.callbacks import Callback
import keras

import csv
import sys

#Prevent tensorflow from using all the memory
config = tf.ConfigProto()
config.gpu_options.allow_growth=True
sess = tf.Session(config=config)

#Define reconstruction loss
def reconstruction_loss(x_input, x_decoded):
    return metrics.mse(x_input, x_decoded)

#Read user inputs
import sys
cancer = sys.argv[1]
tcga_name = sys.argv[2]
print("CANCER NAME: " + cancer)
print("TEST NAME: " + tcga_name)

input_folder = '../../ALL_CANCER_FILES/' + cancer + '/'
output_folder = '../../ALL_CANCER_FILES/' + cancer + '/TCGA_FILES/' 

start = 0
end = 10
dimension = 150

#Read TCGA RNA-Seq input data
input_df_test = pd.read_table(input_folder + 'TCGA_FILES/TCGA_RNASEQ_' + tcga_name + '_PCA_1000L.tsv', index_col = 0)      
print("RNA-Seq expression dataframe ", input_df_test.shape)

#Encode test data with all 10 DAE models
for fold in range(start, end):
    print("DAE model with " + str(dimension) + " nodes and fold " + str(fold))
    
    #Load DAE models
    json_file = open(input_folder + 'DAE_FILES/DAE_' + cancer + '_encoder_' + str(dimension) + 'L_' + str(fold) + '.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    encoder = model_from_json(loaded_model_json)
    
    encoder.load_weights(input_folder + 'DAE_FILES/DAE_' + cancer + '_encoder_' + str(dimension) + 'L_' + str(fold) + '.h5')
    print("Loaded model from disk")

    # Encode test data using the DAE model
    test_encoded = encoder.predict(input_df_test)
    test_encoded_df = pd.DataFrame(test_encoded, index = input_df_test.index)
    print("Encoded data ", test_encoded_df.shape)
    test_encoded_df.to_csv(output_folder + 'TCGA_' + tcga_name + '_RNASeq_Expression_DAE_encoded_' + str(dimension) + 'L_' + str(fold) + '.tsv', sep = '\t')



