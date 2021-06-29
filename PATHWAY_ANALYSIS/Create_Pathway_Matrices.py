###############################
#Author: Ayse Dincer
#Script for creating pathway matrices for cancer type
###############################

import numpy as np
import pandas as pd
import csv
import sys

#Read cancer name and pathway file
cancer_type = sys.argv[1]
pathway_name = sys.argv[2]

input_folder = '../ALL_CANCER_FILES/' + cancer_type + '/'
output_folder = '../ALL_CANCER_FILES/' + cancer_type + '/PATHWAY_FILES/'

def create_pathway_matrix(cancer_type, pathway_file):
   
    #1) Read input data
    data_df = pd.read_table(input_folder + cancer_type + '_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv', index_col=0)
    print("Input data ", data_df.shape)
    gene_names = data_df.columns
    
    #2) Read pathway data
    if pathway_file == 'C2':
        filename = 'MSIGDB_PATHWAYS/c2.v6.2.symbols.gmt'
    if pathway_file == 'H':
        filename = 'MSIGDB_PATHWAYS/h.all.v6.2.symbols.gmt'
    if pathway_file == 'C4_CGN':
        filename = 'MSIGDB_PATHWAYS/c4.cgn.v6.2.symbols.gmt'
    if pathway_file == 'C4_CM':
        filename = 'MSIGDB_PATHWAYS/c4.cm.v6.2.symbols.gmt'
    if pathway_file == 'C5_BP':
        filename = 'MSIGDB_PATHWAYS/c5.bp.v6.2.symbols.gmt'
    if pathway_file == 'C5_CC':
        filename = 'MSIGDB_PATHWAYS/c5.cc.v6.2.symbols.gmt'
    if pathway_file == 'C5_MF':
        filename = 'MSIGDB_PATHWAYS/c5.mf.v6.2.symbols.gmt'
    if pathway_file == 'C6':
        filename = 'MSIGDB_PATHWAYS/c6.all.v6.2.symbols.gmt'
    if pathway_file == 'C7':
        filename = 'MSIGDB_PATHWAYS/c7.all.v6.2.symbols.gmt'
    
           
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content] 

    pathway_count = len(content)
    print("Pathway count " + str(pathway_count))
    
    pathway = np.zeros((len(gene_names), pathway_count), dtype = np.int)
    pathway_names = []
    pathway_lens = []
    
    for i in range(pathway_count):
        data = content[i].split("\t")
        genes = data[2:]
        pathway_name = data[0]
        pathway_names.append(pathway_name)

        pathway_lens.append(len(genes))
        
        #Loop through all genes
        for j in range(len(genes)):

            index = np.where(gene_names == genes[j])[0]
            if len(index) != 0:
                pathway[index[0], i] = 1

    #3) Save matrix
    new_df = pd.DataFrame(pathway, index = gene_names, columns = pathway_names)
    print("Pathway matrix ", new_df.shape)
    print("Average pathway length ", np.mean(pathway_lens))
    print("Average pathway length ", pathway_lens)
    new_df.to_csv(output_folder + 'PATHWAY_' + pathway_file + '_MATRIX_INTERSECTION_GENES.tsv', sep='\t', quoting = csv.QUOTE_NONE)

    
    #Also record gene symbols
    with open(output_folder + 'Gene_Symbols.txt', 'w') as f:
        for item in gene_names:
            f.write("%s\n" % item)
    

create_pathway_matrix(cancer_type, pathway_name)