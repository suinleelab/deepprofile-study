###############################
#Author: Ayse Dincer
#Script for running integrated gradients to get gene-level attributions of each node
###############################

import os
import numpy as np
import pandas as pd
import math 
from sklearn.metrics import mean_squared_error
import tensorflow as tf
from keras.layers import Input, Dense, Lambda, Layer, Activation
from keras.layers.normalization import BatchNormalization
from keras.models import Model
from keras import backend as K
from keras import metrics, optimizers
from keras.callbacks import Callback
import keras
import csv
from keras.models import model_from_json
import sys

#Prevent tensorflow from using all the memory
config = tf.ConfigProto()
config.gpu_options.allow_growth=True
sess = tf.Session(config=config)

#Read all user inputs
cancer = sys.argv[1]
dimension = int(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])

print("CANCER " + str(cancer))
print("DIM " + str(dimension))
print("START " + str(start)) 
print("END " + str(end)) 

input_folder = '../ALL_CANCER_FILES/' + cancer + '/' 
output_folder = '../ALL_CANCER_FILES/' + cancer + '/VAE_WEIGHTS/' 

#Load PCA weights
pca_df = pd.read_table(input_folder + cancer + '_DATA_TOP2_JOINED_PCA_1000L_COMPONENTS.tsv', index_col = 0)
print("PCA COMPONENTS ",  pca_df.shape)
pca_components = pca_df.values

 #Read input data
input_df = pd.read_table(input_folder + cancer + '_DATA_TOP2_JOINED_PCA_1000L.tsv', index_col=0)
print("INPUT FILE ", input_df.shape)

#VAE loss definition
def vae_loss(x_input, x_decoded):
    reconstruction_loss = original_dim * metrics.mse(x_input, x_decoded)
    kl_loss = - 0.5 * K.sum(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1)
    return K.mean(reconstruction_loss + (K.get_value(beta) * kl_loss))

#Save the weight for each run
for vae_run in range(start, end):
    
    print("MODEL " + str(vae_run))
    
    #Load model
    json_file = open(input_folder + 'VAE_FILES/VAE_' + cancer + '_encoder_' + str(dimension) + 'L_' + str(vae_run) + '.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    encoder = model_from_json(loaded_model_json)
    
    #Load weights
    encoder.load_weights(input_folder + 'VAE_FILES/VAE_' + cancer + '_encoder_' + str(dimension) + 'L_' + str(vae_run) + '.h5')
    print("Loaded model from disk")

    #Define hyperparameters
    input_df_training = input_df
    original_dim = input_df_training.shape[1]
    intermediate1_dim = 100
    intermediate2_dim = 25
    latent_dim = dimension

    batch_size = 50
    epochs = 50
    learning_rate = 0.0005
    beta = K.variable(1)
    kappa = 0

    #Encoder network
    x = Input(shape=(original_dim, ))

    net = Dense(intermediate1_dim)(x)
    net2 = BatchNormalization()(net)
    net3 = Activation('relu')(net2)

    net4 = Dense(intermediate2_dim)(net3)
    net5 = BatchNormalization()(net4)
    net6 = Activation('relu')(net5)

    z_mean = Dense(latent_dim)(net6)
    z_log_var = Dense(latent_dim)(net6)

    adam = optimizers.Adam(lr=learning_rate)
    encoder.compile(optimizer=adam, loss = vae_loss)
    encoder.summary()

    #Encode training data using the model
    training_encoded = encoder.predict(input_df_training, batch_size = batch_size)
    print("ENCODED TRAINING DATA ", training_encoded.shape)


    #Measure weights and save absolute value of importance, averaged over samples
    from IntegratedGradients import *

    ig = integrated_gradients(encoder)

    overall_weights = np.zeros((pca_components.shape[0], dimension))

    #Go over each node
    for latent in range(dimension):
        print("Node " + str(latent + 1))
        weights = np.zeros((pca_components.shape[0], input_df_training.shape[0]))
        
        #Go over each sample
        for i in range(input_df_training.shape[0]):
            #print("Sample " + str(i + 1))
            vals = ig.explain(input_df_training.values[i, :], latent)    
            new_vals = np.matmul(vals, pca_components.T)
            weights[:, i] = new_vals
            
        #Take absolute values avg over all samples 
        overall_weights[:, latent] = np.mean(np.abs(weights), axis = 1)

    ig_df = pd.DataFrame(overall_weights, index = pca_df.index)
    print("EXPLANATIONS DF ", ig_df.shape)
    
    ig_df.to_csv(output_folder + cancer + '_DATA_VAE_Cluster_Weights_TRAINING_' + str(dimension) + 'L_fold' + str(vae_run) + '.tsv', sep='\t', quoting = csv.QUOTE_NONE)
    
