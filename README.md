# DeepProfile

Repository with scripts for all model training and analysis for paper "A deep profile of gene expression across 18 human cancers"

All fully pre-processed input data for training the models can be found on our Mendeley Data repository. For each cancer, the basic data we used is **'CANCER_DATA_TOP2_JOINED_BATCH_CORRECTED_CLEANED.tsv'** where CANCER is the name of the cancer type. This data is GEO datasets collected from top 2 platforms, intersecting genes taken, and batch correction applied.

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




