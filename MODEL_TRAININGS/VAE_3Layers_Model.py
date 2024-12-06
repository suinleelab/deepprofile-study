###############################
#VAE model

#Code is modified from https://github.com/keras-team/keras/blob/master/examples/variational_autoencoder.py
###############################

import os
import numpy as np
import pandas as pd
import math 
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import tensorflow as tf
from keras.layers import Input, Dense, Lambda, Layer, Activation, Dropout
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


# Method for reparameterization trick to make model differentiable
def sampling(args):
    
    # Function with args required for Keras Lambda function
    z_mean, z_log_var = args

    # Draw epsilon of the same shape from a standard normal distribution
    epsilon = K.random_normal(shape=K.shape(z_mean), mean=0., stddev=1.0)
    
    # The latent vector is non-deterministic and differentiable
    # in respect to z_mean and z_log_var
    z = z_mean + K.exp(z_log_var / 2) * epsilon
    return z

#Method for defining the VAE loss
def vae_loss(x_input, x_decoded):
    
    reconstruction_loss = original_dim * metrics.mse(x_input, x_decoded)
    kl_loss = - 0.5 * K.sum(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1)
    
    return K.mean(reconstruction_loss + (K.get_value(beta) * kl_loss))

#Method for calculating the reconstruction loss
def reconstruction_loss(x_input, x_decoded):
    
    return metrics.mse(x_input, x_decoded)

#Method for calculating the KL-divergence loss
def kl_loss(x_input, x_decoded):
    return - 0.5 * K.sum(1 + z_log_var - K.square(z_mean) - K.exp(z_log_var), axis=-1)

class WarmUpCallback(Callback):
    def __init__(self, beta, kappa):
        self.beta = beta
        self.kappa = kappa
    
    # Behavior on each epoch
    def on_epoch_end(self, epoch, logs={}):
        if K.get_value(self.beta) <= 1:
            K.set_value(self.beta, K.get_value(self.beta) + self.kappa)

#Read input file
cancer_type = sys.argv[1]

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/VAE_FILES/'

input_filename = input_folder + cancer_type + '_DATA_TOP2_JOINED_PCA_1000L.tsv'
output_filename = cancer_type + '_DATA_TOP2_JOINED_encoded_'

input_df = pd.read_table(input_filename, index_col=0)
print("INPUT FILE", input_df.shape)
print(input_df.head(5))

# Set hyperparameters
original_dim = input_df.shape[1]
intermediate1_dim = int(sys.argv[2])
intermediate2_dim = int(sys.argv[3])
latent_dim = int(sys.argv[4])
fold = int(sys.argv[5])

#SET RANDOM SEEDS
from numpy.random import seed
seed(123456 * fold)
from tensorflow import set_random_seed
set_random_seed(123456 * fold)


init_mode = 'glorot_uniform'
batch_size = 50
epochs = 50
learning_rate = 0.0005
beta = K.variable(1)
kappa = 0

input_df_training = input_df

#Define encoder
x = Input(shape=(original_dim, ))

net = Dense(intermediate1_dim, kernel_initializer=init_mode)(x)
net2 = BatchNormalization()(net)
net3 = Activation('relu')(net2)

net4 = Dense(intermediate2_dim, kernel_initializer=init_mode)(net3)
net5 = BatchNormalization()(net4)
net6 = Activation('relu')(net5)

z_mean = Dense(latent_dim, kernel_initializer=init_mode)(net6)
z_log_var = Dense(latent_dim, kernel_initializer=init_mode)(net6)

# Sample from mean and var
z = Lambda(sampling, output_shape=(latent_dim,))([z_mean, z_log_var])

#Define decoder
decoder_h = Dense(intermediate2_dim, activation='relu', kernel_initializer=init_mode)
decoder_h2 = Dense(intermediate1_dim, activation='relu', kernel_initializer=init_mode)
decoder_mean = Dense(original_dim, kernel_initializer=init_mode)

h_decoded = decoder_h(z)
h_decoded2 = decoder_h2(h_decoded)
x_decoded_mean = decoder_mean(h_decoded2)

#VAE model
vae = Model(x, x_decoded_mean)

adam = optimizers.Adam(lr=learning_rate)
vae.compile(optimizer=adam, loss = vae_loss, metrics = [reconstruction_loss, kl_loss])
vae.summary()

#Train model
history  = vae.fit(np.array(input_df_training), np.array(input_df_training),
               shuffle=True,
               epochs=epochs,
               batch_size=batch_size,
               verbose = 2,
               callbacks=[WarmUpCallback(beta, kappa)])

# DEFINE ENCODER
encoder = Model(x, z_mean)

#DEFINE DECODER
decoder_input = Input(shape=(latent_dim, )) 
_h_decoded = decoder_h(decoder_input)
_h_decoded2 = decoder_h2(_h_decoded)
_x_decoded_mean = decoder_mean(_h_decoded2)
decoder = Model(decoder_input, _x_decoded_mean)


training_encoded = encoder.predict(input_df_training, batch_size = batch_size)
training_encoded_df = pd.DataFrame(training_encoded, index = input_df_training.index)

# How well does the model reconstruct the input data
training_reconstructed = decoder.predict(np.array(training_encoded_df))
training_reconstructed_df = pd.DataFrame(training_reconstructed, index = input_df_training.index, columns = input_df_training.columns)

recons_error = mean_squared_error(np.array(input_df_training), np.array(training_reconstructed_df))

print("TRAINING RECONSTRUCTION ERROR: " + str(recons_error))

#Save encoded test data
training_encoded_df.to_csv(output_folder + output_filename + str(latent_dim) + "L_TRAINING_fold" + str(fold) + ".tsv", sep = '\t')


#SAVE ENCODER MODEL
from keras.models import model_from_json

model_json = encoder.to_json()
with open(output_folder + "VAE_" + cancer_type + "_encoder_" + str(latent_dim) + "L_"+ str(fold) + ".json", "w") as json_file:
    json_file.write(model_json)

encoder.save_weights(output_folder + "VAE_" + cancer_type + "_encoder_" + str(latent_dim) + "L_"+ str(fold) + ".h5")
print("Saved model to disk")


model_json = decoder.to_json()
with open(output_folder + "VAE_" + cancer_type + "_decoder_" + str(latent_dim) + "L_"+ str(fold) + ".json", "w") as json_file:
    json_file.write(model_json)

decoder.save_weights(output_folder + "VAE_" + cancer_type + "_decoder_" + str(latent_dim) + "L_"+ str(fold) + ".h5")
print("Saved model to disk")


#Calculate training r squared
from sklearn.metrics import r2_score

training_r2_vals = np.zeros(input_df_training.shape[0])
for i in range(input_df_training.shape[0]):
    training_r2 = r2_score(input_df_training.values[i, :], training_reconstructed_df.values[i, :])
    training_r2_vals[i] = training_r2

print("TRAINING R2 " + str(np.mean(training_r2_vals)))
