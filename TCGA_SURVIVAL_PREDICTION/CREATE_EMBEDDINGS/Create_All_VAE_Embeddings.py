import sys

cancer_type = sys.argv[1]
tcga_type = sys.argv[2]

dims = [5, 10, 25, 50, 75, 100]
for dim in dims:
    get_ipython().magic(u"run -i 'Encode_TCGA_Data_with_VAE.py' '" +  cancer_type + "' " + tcga_type + " " + str(dim) + " " + str(0) + " " + str(100))