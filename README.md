# deepprofile-study
Repository with scripts for all model training and analysis for paper "A deep profile of gene expression across 18 human cancers"

The folder **MODEL_TRAININGS** includes all scripts and notebooks for training VAE models and obtaining attributions.

The script **Example_Run_All.py** includes all commands for training DeepProfile model for one cancer type. 

<font color=blue>**STEP 1: Creating PCs for each data**</font>

**Create_PCs_for_DeepLearning_Models.py** takes a cancer type and component_count as input and applies PCA on the training data to train deep learning models.

<font color=blue>**STEP 2: Training VAE models**</font>

**VAE_3Layers_Model.py** is the Keras implementation of VAE model.
**Run_VAE_Models.py** takes the cancer type, number of hidden nodes, and start-end folds to train VAE models for the given cancer type.

<font color=blue>**STEP 3: Running IG for VAE models**</font>

**IntegratedGradients.py** is the Keras implementation for Integrated Gradients feature attribution method.
**Get_VAE_IG_Attributions.py** is the script for running IG and get gene-level explanations for each of the nodes. It takes the cancer type, number of hidden nodes, and start-end folds to get explanations for the VAE models for the given cancer type.

<font color=blue>**STEP 4: Learning ensemble labels**</font>

**Create_Ensemble_Labels.py** is the script for running k-means clustering to learn ensemble weights. It takes the cancer type and number of final latent nodes as the input and saves the ensemble labels. 
**Select_Latent_Dimension_with_Gmeans** is the notebook for running g-means clustering to select  the ensemble latent dimension size.

<font color=blue>**STEP 5: Creating DeepProfile ensemble training embedding**</font>

**Create_DeepProfile_Training_Embeddings.py** is the script for joining all the training data VAE embeddings and ensembling them using the learned ensemble labels. It takes the cancer type as the input and creates training DeepProfile ensemble embedding. 

<font color=blue>**STEP 6: Creating DeepProfile ensemble gene attribution matrices**</font>

**Create_DeepProfile_Ensemble_Weights.py** is the script for joining all the VAE gene attributions and ensembling them using the learned ensemble labels. It takes the cancer type as the input and creates DeepProfile gene attribution matrix. 

### <font color=red>  PART 2: TRAINING COMPETITOR MODELS </font>

The script **Example_Run_All.py** includes all commands for training competitor models for one cancer type. 

In **COMPETITOR_TRAININGS**, all the scripts for comparing DeepProfile to other methods is included

<font color=blue>**STEP 1: Training PCA Models**</font>

**Create_PCA_Data.py** takes a cancer type and creates PCA components for the training data.

<font color=blue>**STEP 2: Training ICA Models**</font>

**Create_ICA_Data.py** takes a cancer type and creates ICA components for the training data, repeating 10 times.

<font color=blue>**STEP 3: Training RP Models**</font>

**Create_RP_Data.py** takes a cancer type and creates RP components for the training data, repeating 10 times.

<font color=blue>**STEP 4: Training AE Models**</font>

**AE_2Layers_Model.py** is the Keras implementation of AE model.
**Train_AE_Models.py** takes a cancer type as input and trains 10 AE models with different random seeds.
**Get_AE_IG_Attributions.py** is the script for running IG and get gene-level explanations for each of the nodes. It takes the cancer type and fold to get explanations for the AE models for the given cancer type.

<font color=blue>**STEP 5: Training DAE Models**</font>

**DAE_2Layers_Model.py** is the Keras implementation of DAE model.
**Train_DAE_Models.py** takes a cancer type as input and trains 10 DAE models with different random seeds.
**Get_DAE_IG_Attributions.py** is the script for running IG and get gene-level explanations for each of the nodes. It takes the cancer type and fold to get explanations for the DAE models for the given cancer type.



### <font color=red>  PART 3: TCGA SURVIVAL PREDICTIONS </font>

In **TCGA_SURVIVAL_PREDICTION** folder, all files and scripts are included for predicting TCGA expression survival.

In folder **TCGA_DATA**, **TCGA_CLINICAL_DATA** folder includes clinical data for TCGA samples. **TCGA_MICROARRAY** folder includes microarray expression data and **TCGA_RNASEQ** folder includes RNA-Seq expression data.

<font color=blue>**STEP 1: Preprocessing Data**</font>

**CREATE_EMBEDDINGS** folder includes all scripts to generate TCGA RNA-Seq embeddings.

The script **Example_Run_All.py** includes all commands for generating TCGA expression embeddings for one cancer type. 

**Preprocess_TCGA_Rnaseq_Expression.py** script takes the cancer type and TCGA cancer type as input and preprocesses the expression data to train models.

**Create_TCGA_Rnaseq_PCs.py** script again takes the cancer type and TCGA cancer type as input and applies PCA to preprocessed expression to record top PCs to train deep learning models.

<font color=blue>**STEP 2: Encoding Expression with DeepProfile**</font>

**Encode_TCGA_Data_with_VAE.py** takes the preprocessed PCAed expression and encodes it using the already trained VAE models. The script takes cancer type, TCGA type, VAE dimension, start and end runs to encode the expression.

**Create_All_VAE_Embeddings.py** takes cancer type and TCGA type as input and encoder TCGA expression with all trained VAE models.

**Create_DeepProfile_TCGA_Embeddings.py** takes the cancer type and TCGA type as input and generates the DeepProfile embedding. The script loads in all the VAE embeddings and ensemble labels to generate an ensemble DeepProfile embedding.

<font color=blue>**STEP 3: Encoding Expression with Competitor Models**</font>

**Encode_TCGA_Data_with_PCA.py** takes the cancer type and TCGA type as input and generated PCA embedding for TCGA RNA-Seq expressions.

**Encode_TCGA_Data_with_ICA.py** takes the cancer type and TCGA type as input and generated ICA embedding for TCGA RNA-Seq expressions.

**Encode_TCGA_Data_with_RP.py** takes the cancer type and TCGA type as input and generated RP embedding for TCGA RNA-Seq expressions.

**Encode_TCGA_Data_with_AE.py** takes the cancer type and TCGA type as input and generated AE embedding for TCGA RNA-Seq expressions.

**Encode_TCGA_Data_with_DAE.py** takes the cancer type and TCGA type as input and generated DAE embedding for TCGA RNA-Seq expressions.

<font color=blue>**STEP 4: Generating Survival DataFrames**</font>

Folder **CREATE_SURVIVAL_DATAFRAMES** includes all scripts for generating survival data frames. 

**Create_TCGA_Survival_Dataframes.py** takes the cancer type and TCGA type as input and extract the necssary fields from clinical data to define the survival dataframe.

**Create_Joined_Survival_Dataframes.py** takes the cancer type and TCGA type as input and comnbines the DeepProfile RNA-Seq embeddings with survival data frames.

**Create_Joined_Survival_Dataframes_Cancer_Types.py** combines data frames for cancer subtypes under the main cancer type. 

<font color=blue>**STEP 5: Predicting Survival**</font>

Folder **PREDICT_SURVIVAL** contains scripts for predicting survival.

**Predict_Survival.py** trains lasso regression models with subsampling taking the TCGA RNA-Seq embeddings as the input. 

**Predict_Survival_Subtypes_Joined.py** trains lasso regression models with subsampling taking the TCGA RNA-Seq embeddings as the input while joining multiple TCGA cancer types when there are multiple TCGA cancer subtypes corresponding to one major cancer type we have.

**Run_Models.py** trains all prediction models for all models and cancer types.

**Plots_of_Survival_Prediction.ipynb** and **Plots_of_Survival_Prediction_VAEs.ipynb** are notebooks for generating plots of comparing survival predictions of models. 

<font color=blue>**STEP 6: Comparing RNA-Seq and microarray DeepProfile embeddings**</font>

**COMPARING_RNASEQ_and_MICROARRAY** folder includes all scripts to generate TCGA microarray embeddings and to compare the embeddings with RNA-Seq embeddings.

**Preprocess_TCGA_Rnaseq_Expression.py** script takes the cancer type and TCGA cancer type as input and preprocesses the expression data to train models.

**Create_TCGA_Microarray_PCs.py** script again takes the cancer type and TCGA cancer type as input and applies PCA to preprocessed expression to record top PCs to train deep learning models.

**Encode_TCGA_Microarray_Data_with_VAE.py** takes the preprocessed PCAed expression and encodes it using the already trained VAE models. The script takes cancer type, TCGA type, VAE dimension, start and end runs to encode the expression.

**Create_DeepProfile_TCGA_Microarray_Embeddings.py** takes the cancer type and TCGA type as input and generates the DeepProfile embedding. The script loads in all the VAE embeddings and ensemble labels to generate an ensemble DeepProfile embedding.

**Rnaseq_and_Microarray_Embedding_Correlation_Plots** notebook calculates correlation between RNA-Seq and microarray embeddings and generates plots.



### <font color=red>  PART 4: PATHWAY ENRICHMENT TESTS </font>

In **PATHWAY_ANALYSIS** folder, the scripts and files for pathway analysis are included.

**MSIGDB_PATHWAYS** folder, the files for Molecular Signature Database pathways are included.

<font color=blue>**STEP 1: Running enrichment tests**</font>

**Create_Pathway_Matrices.py** is the script for creating binary pathway matrices for the genes that are present in the training datasets. It takes a cancer type and pathway type as input and creates an binary matrix of pathway overlaps.

**Fishers_Test.py** is the script for running Fisher's test. It takes the cancer type, pathway type, the method name, and the range of runs and records uncorrected and FDR-corrected p-values. 

**Run_Multiple_Fishers_Test.py** is the script for running multiple tests consecutively. It takes the cancer type and pathway name as input and carries enrichment tests for all methods. 

<font color=blue>**STEP 2: Comparing pathway coverages**</font>

**PATHWAY COVERAGE ANALYSIS** folder includes all scripts for comparing pathway coverage of models.

**Plot_of_Average_Pathway_Coverages** genereates plots of average pathway coverage to compare DeepProfile and other dimensionality reduction methods.

**Plot_of_Pathway_Coverage_Distributions** generates plots of distribution of pathway count of each node of DeepProfile and other dimensionality reduction methods.

**Plot_of_Node_Level_Pathway_Annotations** generates plots of percent of nodes annotated by at least one pathway across multiple thresholds.

**Plot_of_Pathway_Detection_Comparison_VAEs_vs_DeepProfile** creates plots for comparing pathways captured by DeepProfile vs individual VAE models.

**Plot_of_Pathway_Percent_Comparison_VAEs_vs_DeepProfile**  creates plots for comparing pathways captured by DeepProfile vs individual VAE models based on percentages.



### <font color=red>  PART 5: NORMAL TISSUE ANALYSIS </font>

In **NORMAL_TISSUE_ANALYSIS** folder, the scripts for normal tissue analysis are included.

The script **Example_Run_All.py** includes all commands for generating normal tissue expression embeddings for one cancer type. 

**Gtex_Tissue_Name_Mappings** is the notebook for mapping GTEX tissue names to cancer types we have. The GTEX expression data includes samples from many different tissues. We extract the GTEX sample names for each cancer type we have.

**Preprocess_Gtex_Rnaseq_Expressions.py** is the script for creating preprocessed GTEX gene expression. It takes the cancer type as input and preprocesses the GTEX RNA-Seq expression using the same preprocessing steps applied to our training data.

**Create_Gtex_Rnaseq_PCs.py** is the script for taking top PCs of the GTEX expression profiles to train DeepProfile model. It takes the cancer type as input and records the top PCs of GTEX expression.

**Encode_GTEX_Data_with_VAE.py** is the script for encoding GTEX expression using trained VAE models. The inputs to the model are the cancer type, the number of latent nodes, and start and end runs. 

**Create_DeepProfile_GTEX_Embeddings.py** is the script for creating DeepProfile embedding using generated VAE embeddings. It takes the cancer type as input and records the final DeepProfile embedding for GTEX normal tissue samples.

**Normal_Tissue_Classifier.py** is the script for training the classifier to separate normal vs cancer tissue embeddings. It takes the cancer type as input and records the bootstrapped classifier weights.



### <font color=red>  PART 6: PANCANCER ANALYSIS </font>

In **PANCANCER_ANALYSIS** folder, the scripts for various pancancer analysis are included.

<font color=blue>**1. COMMON GENES ANALYSIS**</font>

**Create_Joined_Gene_Ranking_Matrix** is the script for combining gene importance matrices for all cancer types and recording the joined matrix.

**Create_Plot_of_Common_Genes** script creates a plot of top cancer common genes, showing their average percentile scores.

**Common_Top_Genes_Pathway_Enrichment_Plots** script runs Fisher's Exact Test for top cancer common genes. It then generates plots for showing the top enriched pathways. 

**PCA_Common_Top_Genes_and_Pathway_Enrichments** script records the top universally important genes for PCA and carries enrichment tests for the top common PCA genes. 

**DeepProfile_Immune_Cell_Type_Enrichment_Tests** script runs enrichment tests for the top universally important DeepProfile genes using immune cell type gene lists.

**Receptor_Files** folder includes all receptor gene lists.

**DeepProfile_Receptor_Enrichment_Tests** script runs enrichment tests for the top universally important DeepProfile genes using all the receptor gene lists.

**PCA_Receptor_Enrichment_Tests** script runs enrichment tests for the top universally important PCA genes using all the receptor gene lists.

**Create_Plots_of_Receptor Scores** script generates plots of receptor enrichment scores.


<font color=blue>**2. COMMON PATHWAYS ANALYSIS**</font>

**Create_Joined_PValue_Matrix** script reads p-value scores for each cancers and joins all matrices to define the p-value matrix.

**Record_Pancancer_Cancer_Character_Scores** script calculates the average cancer character scores for all cancers.

**Plot_of_Cancer_Common_Pathways** script generates the plots for cancer common pathways.

**Record_Cancer_Common_Pathways_Table** script records the cancer common pathway scores.


<font color=blue>**3.CANCER SPECIFIC GENES ANALYSIS**</font>

**Plot_of_Cancer_Specific_Genes** script creates plot of top cancer specific genes for each cancer.

**Plot_of_Cancer_Top_Genes** script creates plot of top genes for each cancer.

**Signature_Gene_Lists** folder includes all cancer type specific gene lists.

**Cancer_Specific_Genes_Enrichment_Tests.ipynb** script runs enrichment tests for cancer-specific genes using cancer subtype signatures.


<font color=blue>**4.CANCER SPECIFIC PATHWAYS ANALYSIS**</font>

**Plot_of_Cancer_Specific_Pathways** script creates the plots for top cancer specific pathways.

**Plot_of_Cancer_Specific_Pathways_with_Cancer_and_Character** script creates plots of cancer specific pathways along with cancer character vectors.



### <font color=red>  PART 7: SURVIVAL AND MUTATION ANALYSIS </font>

SURVIVAL_AND_MUTATION_ANALYSIS folder includes all scripts for carrying survival and mutation analysis


<font color=blue>**STEP 1:  Pan-cancer Survival Analysis**</font>

The folder **SURVIVAL_ANALYSIS** contains scripts for pancancer survival analysis. 

**Record_DeepProfile_Node_Survival_Scores.R** fits Cox regression model to each DeepProfile node and records results. It takes the main cancer type and the TCGA cancer type as inputs and records the survival coefficient p-values and z-scores. 

**Calculate_Universal_Pathway_Enrichment_and_Survival_Scores** records the top cancer common survival related pathways by calculating universal enrichment and survival scores.

**Create_Cancer_Specific_Survival_and_Mutation_Plots** creates plots of cancer specific survival related pathways and mutation related pathways.


<font color=blue>**STEP 2:  Pan-cancer Mutation Analysis**</font>

The folder **MUTATION_ANALYSIS** contains scripts for pancancer mutation analysis and **TCGA_MUTATION_DATA** contains the downloaded TCGA mutation profile datasets.

**Generate_Mutation_Counts_Data** script converts the sample-level mutation profiles to total number of mutations for each sample.

**Combine_DeepProfile_Embeddings_and_Mutation_Data** script maps the DeepProfile embeddings to mutation counts and calculates the correlation between DeepProfile nodes and mutation.

**Calculate_Universal_Pathway_Enrichment_and_Mutation_Scores** records the top cancer common mutation related pathways by calculating universal enrichment and mutation scores.


<font color=blue> **STEP 3:  Downstream Analysis** </font>

The folder **DOWNSTREAM_ANALYSIS** contains scripts for all the survival and mutation downstream analysis.

The folder **SURVIVAL** contains survival related downstream analysis scripts.

**Record_MMR_Genes_Survival_Scores.R** fits Cox regression model to each gene in mismatch repair pathway and records results. It takes the main cancer type and the TCGA cancer type as inputs and records the survival coefficient p-values and z-scores. 

**Create_MMR_Genes_Clustermaps** creates clustermap of gene-level MMR survival scores.

**Record_MHC_II_Genes_Survival_Scores.R** fits Cox regression model to each gene in mhc II pathway and records results. It takes the main cancer type and the TCGA cancer type as inputs and records the survival coefficient p-values and z-scores. 

**Create_MHC_II_Genes_Clustermap** creates clustermap of gene-level mhc II survival scores.

**Calculate_Average_Pathway_Expression** script records average expression of each pathway and saves the survival dataframes.

**Create_KMPlots** creates Kaplan Meier plots for average expression of selected pathways.

The folder **IMMUNE_CELL_SIGNATURES** contains immune cell signature analysis scripts.

**Immune_Cell_Signatures** folder contains all immune cell signature files. 

**Immune_Cell_Types_Gene_Ranking_Analysis** script calculates the average gene ranking of the gene signatures for different immune cell types.

**Immune_Cell_Types_Correlation_Analysis** script calculates the correlation between HLA-D expression and cell type signatures for different immune cell types.

**Macrophages_M1_and_M2_Analysis** script carries gene ranking analysis to compare M1 and M2 macrophage signatures.

**Create_Immune_Cell_Type_Plots** generates all plots for the immune cell signatures. 






