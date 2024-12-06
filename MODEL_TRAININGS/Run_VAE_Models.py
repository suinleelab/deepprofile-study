###############################
#Script for training VAE models
###############################
import sys

cancer_type = sys.argv[1]
latent = int(sys.argv[2])
start = int(sys.argv[3])
end = int(sys.argv[4])

if latent == 5:  
    dim1 = 100
    dim2 = 25
if latent == 10:   
    dim1 = 250
    dim2 = 50
if latent == 25:  
    dim1 = 250
    dim2 = 100
if latent == 50:  
    dim1 = 250
    dim2 = 100
if latent == 75:  
    dim1 = 250
    dim2 = 100
if latent == 100:   
    dim1 = 250
    dim2 = 100

for run in range(start, end):
    get_ipython().magic(u"run -i 'VAE_3Layers_Model.py' '" +  cancer_type + "' " + str(dim1) + " " + str(dim2) + " " + str(latent) + " " + str(run))
