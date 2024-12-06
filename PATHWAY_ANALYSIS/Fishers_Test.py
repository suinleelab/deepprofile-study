###############################
#Script for running fisher's test for pathway enrichment analysis
###############################

import numpy as np
import pandas as pd
import csv
import statsmodels.api as sm
import sys

#Read user inputs
cancer_type = sys.argv[1]
pathway_type = sys.argv[2]
method = sys.argv[3]
start = int(sys.argv[4])
end = int(sys.argv[5])
if len(sys.argv) > 6:
    dimension = int(sys.argv[6])
    L = dimension
else:
    L = 150

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/PATHWAY_FILES/'

pathway_matrix = pd.read_table(input_folder + 'PATHWAY_FILES/PATHWAY_' + pathway_type + '_MATRIX_INTERSECTION_GENES.tsv', index_col = 0)
print(pathway_matrix.shape)
pathway_df = pathway_matrix
pathway_matrix = pathway_matrix.values

if pathway_type == 'C2':
    N = 51 #average number of pathways
if pathway_type == 'C4_CM':
    N = 113 #average number of pathways
if pathway_type == 'C4_CGN':
    N = 99 #average number of pathways
if pathway_type == 'C6':
    N = 166 #average number of pathways
if pathway_type == 'C5_BP':
    N = 114 #average number of pathways
if pathway_type == 'C5_CC':
    N = 151 #average number of pathways
if pathway_type == 'C5_MF':
    N = 106 #average number of pathways
if pathway_type == 'H':
    N = 146 #average number of pathways

#Run test for each run
for run in range(start, end):
    if method == 'PCA':
        data_df = pd.read_table(input_folder + 'PCA_FILES/' + cancer_type + '_DATA_TOP2_JOINED_PCA_COMPONENTS_150L.tsv', index_col = 0)
        print(data_df.shape)

        ensemble_weights = np.abs(data_df.values.T)
        print(ensemble_weights.shape)
        
    if method == 'ICA':
        data_df = pd.read_table(input_folder + 'ICA_FILES/' + cancer_type + '_DATA_TOP2_JOINED_ICA_COMPONENTS_150L_fold' + str(run + 1) + '.tsv', index_col = 0)
        print(data_df.shape)

        ensemble_weights = np.abs(data_df.values.T)
        print(ensemble_weights.shape)

    if method == 'RP':
        data_df = pd.read_table(input_folder + 'RP_FILES/' + cancer_type + '_DATA_TOP2_JOINED_RP_COMPONENTS_fold' + str(run + 1) + '.tsv', index_col = 0)
        print(data_df.shape)

        ensemble_weights = np.abs(data_df.values.T)
        print(ensemble_weights.shape)
        
    if method == 'AE':
        data_df = pd.read_table(input_folder + 'AE_FILES/' + cancer_type + '_DATA_AE_Weights_TRAINING_150L_fold' + str(run + 1) + '.tsv', index_col = 0)
        print(data_df.shape)

        ensemble_weights = data_df.values.T
        print(ensemble_weights.shape)
        
    if method == 'DAE':
        data_df = pd.read_table(input_folder + 'DAE_FILES/' + cancer_type + '_DATA_DAE_Weights_TRAINING_150L_fold' + str(run + 1) + '.tsv', index_col = 0)
        print(data_df.shape)

        ensemble_weights = data_df.values.T
        print(ensemble_weights.shape)
    
    if method == 'DeepProfile':
        data_df = pd.read_table(input_folder + cancer_type + '_DeepProfile_Ensemble_Gene_Importance_Weights_150L.tsv', index_col = 0)
        print(data_df.shape)

        ensemble_weights = data_df.values
        print(ensemble_weights.shape)
        
    if method == 'VAE':
        data_df = pd.read_table(input_folder + 'VAE_WEIGHTS/' + cancer_type + '_DATA_VAE_Cluster_Weights_TRAINING_' + str(dimension) + 'L_fold' + str(run) + '.tsv', index_col = 0)
        print(data_df.shape)

        ensemble_weights = data_df.values.T
        print(ensemble_weights.shape)
        

    #Apply fisher test
    p_vals = np.zeros((ensemble_weights.shape[0], pathway_matrix.shape[1]))

    print("Running for top ", N, " genes")
    
    for i in range(p_vals.shape[0]):
        print(i)
        for j in range(p_vals.shape[1]):

            #Create contingency matrix
            matrix = np.zeros((2, 2))

            pathway_indices = np.where(pathway_matrix[:, j] == 1)[0]
            #print(pathway_df.index[pathway_indices])

            gene_indices = ensemble_weights[i, :].argsort()[-N:][::-1]
            #print(len(gene_indices))
            #print(pathway_df.index[gene_indices])

            in_pathway_firstN = len(np.intersect1d(pathway_indices ,gene_indices))
            #print(pathway_df.index[np.intersect1d(pathway_indices ,gene_indices)])

            out_pathway_firstN = N - in_pathway_firstN
            #print(out_pathway_firstN)

            in_pathway_other = len(pathway_indices) - in_pathway_firstN
            #print(in_pathway_other)

            out_pathway_other = pathway_matrix.shape[0] - in_pathway_other
            #print(out_pathway_other)

            matrix[0, 0] = in_pathway_firstN
            matrix[0, 1] = in_pathway_other
            matrix[1, 0] = out_pathway_firstN
            matrix[1, 1] = out_pathway_other

            import scipy.stats as stats
            oddsratio, pvalue = stats.fisher_exact(matrix)

            p_vals[i, j] = pvalue


    #Record uncorrected p-values
    if method == 'VAE':
        p_vals_df = pd.DataFrame(p_vals, index = np.arange(L) + 1, columns = pathway_df.columns)
        p_vals_df.to_csv(output_folder + cancer_type + '_FISHER_UNCORRECTED_PVALS_' + pathway_type + '_' + method + '_' + str(dimension) + 'L_' + str(run + 1) + '.tsv', sep = '\t')
    else:
        p_vals_df = pd.DataFrame(p_vals, index = np.arange(L) + 1, columns = pathway_df.columns)
        p_vals_df.to_csv(output_folder + cancer_type + '_FISHER_UNCORRECTED_PVALS_' + pathway_type + '_' + method + '_' + str(run + 1) + '.tsv', sep = '\t')

    new_p_values = np.zeros(((p_vals.shape[0], p_vals.shape[1])))

    #Record corrected p-values
    for i in range(pathway_matrix.shape[1]):
        corrected_pval = sm.stats.multipletests( p_vals[:, i], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
        new_p_values[:, i]  = corrected_pval                

    x = np.where([new_p_values < 0.05])[2]
    unique_count = len(np.unique(x))
    print("UNIQUE PATHWAY COUNT: " + str(unique_count))

    p_vals_df = pd.DataFrame(new_p_values, index = np.arange(L) + 1, columns = pathway_df.columns)
    #print(p_vals_df)
    if method == 'VAE':
        p_vals_df.to_csv(output_folder + cancer_type + '_FISHER_FDR_CORRECTED_PVALS_' + pathway_type + '_' + method + '_' + str(dimension) + 'L_' + str(run + 1) + '.tsv', sep = '\t')
    else:
        p_vals_df.to_csv(output_folder + cancer_type + '_FISHER_FDR_CORRECTED_PVALS_' + pathway_type + '_' + method + '_' + str(run + 1) + '.tsv', sep = '\t')

    x = np.where([p_vals_df.values < 0.05])[2]
    unique_count = len((x))
    print("AVERAGE PATHWAY COUNT: ", unique_count / 150)
    
