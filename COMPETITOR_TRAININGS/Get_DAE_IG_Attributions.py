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
vae_run = int(sys.argv[2])
dimension = 150

input_folder = '../ALL_CANCER_FILES/' + cancer + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer + '/DAE_FILES/'

#Load PCA weights
pca_df = pd.read_table(input_folder + cancer + '_DATA_TOP2_JOINED_PCA_1000L_COMPONENTS.tsv', index_col = 0)
print("PCA COMPONENTS ",  pca_df.shape)
pca_components = pca_df.values

#Define reconstruction loss
def reconstruction_loss(x_input, x_decoded):
    return metrics.mse(x_input, x_decoded)

#Save the weight for each 100 runs
print("MODEL " + str(vae_run))

#Load model
json_file = open(input_folder + 'DAE_FILES/DAE_' + cancer + '_encoder_' + str(dimension) + 'L_' + str(vae_run) + '.json', 'r')
loaded_model_json = json_file.read()
json_file.close()

encoder = model_from_json(loaded_model_json)
encoder.load_weights(input_folder + 'DAE_FILES/DAE_' + cancer  + '_encoder_' + str(dimension) + 'L_' + str(vae_run) + '.h5')
print("Loaded model from disk")

#Read input data
input_df = pd.read_table(input_folder + cancer + '_DATA_TOP2_JOINED_PCA_1000L.tsv', index_col=0)
print("INPUT FILE ", input_df.shape)

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

#Define encoder
x = Input(shape=(original_dim, ))

net = Dense(intermediate1_dim)(x)
net2 = BatchNormalization()(net)
net3 = Activation('relu')(net2)

net4 = Dense(intermediate2_dim)(net3)
net5 = BatchNormalization()(net4)
net6 = Activation('relu')(net5)

adam = optimizers.Adam(lr=learning_rate)
encoder.compile(optimizer=adam, loss = reconstruction_loss)
encoder.summary()

#Encode training data using the model
training_encoded = encoder.predict(input_df_training, batch_size = batch_size)
print("ENCODED TRAINING DATA ", training_encoded.shape)

#Measure weights and save absolute value of importance, averaged over samples
from IntegratedGradients import *

ig = integrated_gradients(encoder)

overall_weights = np.zeros((pca_components.shape[0], dimension))

for latent in range(dimension):
    print("Node " + str(latent + 1))
    weights = np.zeros((pca_components.shape[0], input_df_training.shape[0]))

    for i in range(input_df_training.shape[0]):
        vals = ig.explain(input_df_training.values[i, :], latent)    
        new_vals = np.matmul(vals, pca_components.T)
        weights[:, i] = new_vals
        
    #Take absolute values avg over all samples 
    overall_weights[:, latent] = np.mean(np.abs(weights), axis = 1)

ig_df = pd.DataFrame(overall_weights, index = pca_df.index)
print("EXPLANATIONS DF ", ig_df.shape)

ig_df.to_csv(output_folder + cancer + '_DATA_DAE_Weights_TRAINING_' + str(dimension) + 'L_fold' + str(vae_run) + '.tsv', sep='\t', quoting = csv.QUOTE_NONE)
print(ig_df.shape)
