###############################
#Script for encoding TCGA microarray expression using VAE models
###############################

import os
import numpy as np
import pandas as pd

import math 
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

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
from keras.models import model_from_json
from sklearn import preprocessing

#Prevent tensorflow from using all the memory
config = tf.ConfigProto()
config.gpu_options.allow_growth=True
sess = tf.Session(config=config)

#Method for defining the VAE loss
def vae_loss(x_input, x_decoded):
    reconstruction_loss = original_dim * metrics.mse(x_input, x_decoded)
    kl_loss = - 0.5 * K.sum(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1)
    return K.mean(reconstruction_loss + (K.get_value(beta) * kl_loss))


#Read user inputs
import sys
cancer = sys.argv[1]
tcga_name = sys.argv[2]
dimension = int(sys.argv[3])
start = int(sys.argv[4])
end = int(sys.argv[5])

print("CANCER NAME: " + cancer)
print("TEST NAME: " + tcga_name)

input_folder = '../../ALL_CANCER_FILES/' + cancer + '/'
output_folder = '../../ALL_CANCER_FILES/' + cancer + '/TCGA_FILES/' 

#Read input data
input_df_test = pd.read_table(input_folder + 'TCGA_FILES/TCGA_MICROARRAY_' + tcga_name + '_PCA_1000L.tsv', index_col = 0)    
print("TCGA expression dataframe ", input_df_test.shape)
    
#Read GTEX expression
for fold in range(start, end):
    print("VAE model with " + str(dimension) + " nodes and fold " + str(fold))

    #Load VAE models
    json_file = open( input_folder + 'VAE_FILES/VAE_' + cancer + '_encoder_' + str(dimension) + 'L_' + str(fold) + '.json', 'r')   
    loaded_model_json = json_file.read()
    json_file.close()
    encoder = model_from_json(loaded_model_json)

    encoder.load_weights(input_folder + 'VAE_FILES/VAE_' + cancer + '_encoder_' + str(dimension) + 'L_' + str(fold) + '.h5')
    print("Loaded model from disk")

    #Define placeholder VAE model
    original_dim = input_df_test.shape[1]
    intermediate1_dim = 100
    intermediate2_dim = 25
    latent_dim = dimension

    batch_size = 50
    epochs = 50
    learning_rate = 0.0005
    beta = K.variable(1)
    kappa = 0
    init_mode = 'glorot_uniform'

    x = Input(shape=(original_dim, ))

    net = Dense(intermediate1_dim, kernel_initializer=init_mode)(x)
    net2 = BatchNormalization()(net)
    net3 = Activation('relu')(net2)

    net4 = Dense(intermediate2_dim, kernel_initializer=init_mode)(net3)
    net5 = BatchNormalization()(net4)
    net6 = Activation('relu')(net5)

    z_mean = Dense(latent_dim, kernel_initializer=init_mode)(net6)
    z_log_var = Dense(latent_dim, kernel_initializer=init_mode)(net6)

    adam = optimizers.Adam(lr=learning_rate)

    #Encode test data using the VAE model
    test_encoded = encoder.predict(input_df_test, batch_size = batch_size)
    test_encoded_df = pd.DataFrame(test_encoded, index = input_df_test.index)
    test_encoded_df.to_csv(output_folder + 'TCGA_' + tcga_name + '_MICROARRAY_Expression_VAE_encoded_' + str(dimension) + 'L_' + str(fold) + '.tsv', sep = '\t')

