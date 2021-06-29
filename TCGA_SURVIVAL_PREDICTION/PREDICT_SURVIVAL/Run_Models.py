import sys
run = int(sys.argv[1])

cancer_types = ['BRCA', 'AML', 
                'COLON', 
                'BRAIN', 'OV', 
                'SARCOMA', 'KIDNEY', 
                'LIVER', 'STOMACH', 
                'SKIN', 'UCEC', 
                'HEAD_NECK', 'PANCREAS',
                'CERVICAL', 'BLADDER', 'LUNG']

tcga_types = ['BRCA', 'LAML', 
                'COADREAD', 
                'GBMLGG', 'OV', 
                'SARC', 'KIPAN', 
                'LIHC', 'STAD', 
                'SKCM', 'UTERINE',
                'HNSC', 'PAAD',
                'CESC', 'BLCA', 'LUNG']

for c in range(len(cancer_types)):
    cancer_type = cancer_types[c]
    tcga_type = tcga_types[c]
    print("------------")
    print(cancer_type)
    print(tcga_type)
    
    if cancer_type == 'LUNG':
     
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "PCA " + str(run))
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "DeepProfile " + str(run))
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "ICA " + str(run + 1))
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "RP " + str(run + 1))
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "AE " + str(run))
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "DAE " + str(run))
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "VAE " + str(run) + " 5")
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "VAE " + str(run) + " 10")
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "VAE " + str(run) + " 25")
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "VAE " + str(run) + " 50")
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "VAE " + str(run) + " 75")
        get_ipython().magic(u"run -i     Predict_Survival_Subtypes_Joined.py " +  cancer_type + " " + "VAE " + str(run) + " 100")

    else:
        
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "PCA " + str(run))
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "DeepProfile " + str(run))
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "ICA " + str(run + 1))
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "RP " + str(run + 1))
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "AE " + str(run))
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "DAE " + str(run))
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "VAE " + str(run) + " 5")
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "VAE " + str(run) + " 10")
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "VAE " + str(run) + " 25")
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "VAE " + str(run) + " 50")
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "VAE " + str(run) + " 75")
        get_ipython().magic(u"run -i     Predict_Survival.py " +  cancer_type + " "  +  tcga_type + " " + "VAE " + str(run) + " 100")
