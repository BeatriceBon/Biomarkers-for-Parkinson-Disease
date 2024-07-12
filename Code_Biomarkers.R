### HOMEWORK: BIOMARKERS

#setwd('C:/Users/matti/OneDrive/Desktop/Bioingegneria/PrimoAnno/SecondoSemestre/Biomarkers/Homework/Lavoro')
setwd('C:/Profiles/Desktop/MAGISTRALE/A_Biomarkers Precision Medicine and Drug Development/Project_Biomarkers')

#Adding the libraries further used in the code
library(readxl)     # reading excel files
library(caret)        # classification tool
library(glmnet)       # Logistic regression with lasso
library(psych)        # Correlation visualization
library(dplyr)        # Principal Component Analysis
library(plyr)         # Necessary for ordering
library(gbm)          # Gradient Boosting
library(randomForest) # Random Forest
library(e1071)        # SVM
library(MASS)         # Linear Discriminant Analysis
library(lattice)      # Matrices visualization
library(pROC)         # ROC 
library(corrplot)     # Correlation Plot
library(ggplot2)      # Linear regression visualization

######################### DATA PREPROCESSING ####################################
### Loading data relative to the radiomic features and to demographics
# Before the loading, the data sheets present in the data given by the professor 
# were separated into different files. For what concerns the classification part, 
# the radiomics data sheet was saved as 'Homework_radiomics.xlsx', while the 
# sheet concerning the demographics data as 'Homework_Demographics_clinical.xlsx'.
# Here after they are loaded into the R environment.
rad_data <- read_excel("Homework_radiomics.xlsx")
demo_data <- read_excel("Homework_Demographics_clinical.xlsx")
# Changing into factors the "Group" and "Gender" column of the demographics data
demo_data$`Group (1=PD, 0 = Controls)`<- as.factor(demo_data$`Group (1=PD, 0 = Controls)`)
demo_data$`Gender (M=0, F=1)`<- as.factor(demo_data$`Gender (M=0, F=1)`)
summary(demo_data)

### CHECKING OF THE PRESENCE OF MISSING DATA
nrow(rad_data) #53 rows
rad_data <- na.omit(rad_data)
nrow(rad_data) #53 rows
# There are no missing values since the number of rows of the dataset remains unchanged. 
# If there were any, the function na.omit() would have removed the entire line.


### CREATING A NEW DATASET CALLED cl_data
# It is created by adding to the radiomics data the column relative to the class 
# of belonging of every subject in the dataset
# !!! NOTE: for the classification using ALL the data in the demographic data sheet, 
#           see the dedicated section. 
method='radiomics'
cl_data <- rad_data
cl_data$group <- demo_data$`Group (1=PD, 0 = Controls)`

### Z-SCORING
# Creating a new dataset with the z-scored data (important for some classifiers)
cl_data_zsc <- cl_data
cl_data_zsc[, c((2:173),(175:178))] <- scale(cl_data[,c((2:173),(175:178))], center=TRUE, scale=TRUE)
# NOTE: for the z-scoring the column 174 was omitted. This was done because the 
# column is a list of 1s for every subject in the dataset, so the variance was 
# NULL and it originated a list of NAs. It's important to note that, since the 
# values are the same for every subject, this column would be absolutely useless
# for classification.

### DIVIDING INTO TRAINING AND TEST SET
# paying attention to the train-test-all_dataset balance of data partition
set.seed(2)
idx <- createDataPartition(cl_data$group, p = .7, list = FALSE, times = 1) # creates balanced trn/tst
trn=cl_data[idx,c(2:179)]
tst=cl_data[-idx,c(2:179)]
trn_zsc=cl_data_zsc[idx,c(2:179)]
tst_zsc=cl_data_zsc[-idx,c(2:179)]
# Checking that the dataset is homogeneous
table(cl_data$group)/nrow(cl_data) 
table(trn$group)/nrow(trn) 
table(tst$group)/nrow(tst) # --> it is homogeneous! 
# Percentages of HC/Parkinson groups remain homogeneous across training and test set


### ANALYSIS TAKING INTO CONSIDERATION ALSO THE INFORMATION RELATIVE TO THE DEMOGRAPHICS DATA
# In this section we consider for the analysis also the information relative to the demographics data, 
# to see if the addition of such information will change significantly the features selected thereafter. 
# First of all, to make the data loading easier, we created a excel file called 'Homework_dati_uniti.xlsx'
# that contains all the information related to the radiomics alongside the infos on the demographic data 
# (we just added the corresponding columns at the end of the radiomics ones)
# method='all'
# rad_data <- read_excel('Homework_dati_uniti.xlsx')
# ### CHECKING OF THE PRESENCE OF MISSING DATA
# nrow(rad_data) #392 rows
# rad_data <- na.omit(rad_data)
# nrow(rad_data) #392 <-- no missing data!
# # Changing the variables group and gender into factors!
# rad_data$`Group (1=PD, 0 = Controls)`<- as.factor(rad_data$`Group (1=PD, 0 = Controls)`)
# rad_data$`Gender (M=0, F=1)`<- as.factor(rad_data$`Gender (M=0, F=1)`)
# # Creating the cl_data with all the infos
# cl_data <- rad_data
# colnames(cl_data)[c(181)] = c('gender')
# cl_data$group = cl_data$`Group (1=PD, 0 = Controls)`
# cl_data[179] = NULL
# # Creating balanced test and training set, both taking into consideration the group and the gender factors
# set.seed(2)
# idx <- createDataPartition(c(cl_data$group, cl_data$gender), p = .7, list = FALSE, times = 1)
# trn=cl_data[idx,c(2:185)]
# tst=cl_data[-idx,c(2:185)]
# trn = na.omit(trn)
# # Scaling the features
# pre_proc_val <- preProcess(cl_data[,-174], smethod = c("center", "scale"))
# trn_zsc = predict(pre_proc_val, trn)
# tst_zsc = predict(pre_proc_val, tst)
# Runnig these datasets on the following lines of code, it has been possible to notice
# that there's no significant difference in the radiomics features selected or in the 
# performances of classifiers. Therefore, we will from now on use only the smaller 
# features dataset (the one containing only the radiomics features) to proceed with the analysis.





########################## FEATURE SELECTION ###################################
# We will try now different methods to see if the selected features are the same

##### FIRST METHOD: LASSO LOGISTIC REGRESSION

y <- trn_zsc$group # defining response variable
x <- data.matrix(trn_zsc[,c(1:177)]) # defining matrix of predictor variables
# Performing k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(x, y, alpha = 1, family="binomial")
# Finding optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
best_lambda #0.006700195
# Producing plot of test MSE by lambda value
par(mfrow=c(1,1))
plot(cv_model) 
#finding coefficients of best model
best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda, family="binomial")
coef_lasso <- coef(best_model)
Lasso_idx <- which(best_model$beta != 0) # indexes of the selected features by Lasso
# Extracting the data 
if (method=='all') {
  Lasso_training <- trn_zsc[,c(Lasso_idx, 184)]} else {
    Lasso_training <- trn_zsc[,c(Lasso_idx, 178)]}
Lasso_vars <-names(Lasso_training)
Lasso_vars


### SECOND METHOD: BOTTOMUP - ADD FEATURES 1 by 1

# The idea here is to create a logistic regression model adding 1 by 1 every time 
# a new feature. Then perform an analysis with Anova to see if the new feature 
# added is informative for the model or not. At every iteration just the features 
# that give an Anova p-value smaller than alpha=0.05 are kept, the other ones 
# removed from the model.

nfeat_tot=length(trn)-1
x = colnames(trn)
counts=vector(, length=nfeat_tot) 
feat_sel_groups=list()
idx_features=list()
c_idx_feat=vector()
AIC_feat <- vector(, length=nfeat_tot)
alpha=0.05 #alpha threshold
for (fId in 1:nfeat_tot) {
  if (method=='all') {
    c_idx_feat<- c(c_idx_feat, fId, 184)} else {
      c_idx_feat<- c(c_idx_feat, fId, 178)}
  c_trn <- trn[,c_idx_feat] # current training set
  c_lm <- glm(group ~., family=binomial, data=c_trn) #current logit model
  c_statistics <- anova(c_lm, test='Chisq') # anova analysis
  c_indici<- unlist(c_statistics$`Pr(>Chi)`)<alpha #indexes of TRUE-FALSE with respect to the condition
  c_righe <- rownames(c_statistics) #extracting the names of the overall rows considered
  c_idx_true <- which(c_indici == TRUE) # extracting indexes for which p-value less than alpha
  # if there is no feature which satisfies the condition (it happens with the 
  # initial iterations, then the feature with lower p-value is taken)
  if (length( c_idx_true)==0) {
    c_cond <- unlist(c_statistics$`Pr(>Chi)`)==min(c_statistics$`Pr(>Chi)`, na.rm=TRUE)
    c_idx_true <- which(c_cond==TRUE)}
  idx_features[fId] <- list(c_idx_true-1) #selected features at the current iteration
  c_idx_feat <- which(x %in% c_righe[c_idx_true]) # extracting the overall right indexes
  feat_sel_groups[fId]=list(c_righe[c_idx_true]) # extracting the right name of the features
  AIC_feat[fId] = AIC(c_lm) # performing AIC
  counts[fId] <- length(unlist(feat_sel_groups[fId])) # just to control
}
# It is possible to see, showing the feat_sel_groups voice, that the features 
# that are selected at every iteration are stat_mean and stat_skew, so what it 
# happens is that the prediction doesn't improve after 5th iteration.
AIC_feat
# The AIC_feat vector confirms what has been previously found. The lower value 
# in fact is present at the 5th iteration. 
Bottomup_th0.05_vars <- unlist(feat_sel_groups[nfeat_tot])
Bottomup_th0.05_vars
Bottomup_th0.05_idx <- c_idx_feat #since they are the ones surviving at the last iteration





### THIRD METHOD: SELECTING THE MOST IMPORTANT FEATURES, GROUP BY GROUP

# The idea here is to group all the features into the groups described in the 
# homework, and to use the most important features of each group to try 
# predicting the presence of Parkinson's Disease.
var_names <- names(cl_data)
# Dividing features in groups by extracting indeces
# All the indeces have -1 at the end since in training and test set we don't have 
# the initial column of ids.
idx_groups <- structure(list("idx_local_intensity_feat"=(2:3)-1, "idx_stat_feat"=(4:21)-1
                             , "idx_ivh_feat"=(22:35)-1, "idx_morph_feat"=(36:60)-1, "idx_ih_feat"=(61:83)-1
                             , "idx_cm_feat"=(84:108)-1 , "idx_rlm_feat"=(109:124)-1, "idx_szm_feat"=(125:140)-1 
                             , "idx_dzm_feat"=(141:156)-1, "idx_ngt_feat"=(157:161)-1, "idx_ngl_feat"=(162:178)-1))
#Checking the correlation among features of the same group
colori=colorRampPalette(c("yellow", "red"))
for (gId in 1:length(idx_groups)) {
  current_features = unlist(idx_groups[gId])
  if (gId==length(idx_groups)) {
    current_features = current_features[-13]
  }
  cfeat_corr <- corr.test(trn[,current_features], method='kendall', adjust='bonferroni')
  # The 2 following level plots have been explicitly run for every iteration
  levelplot(cfeat_corr$r, pretty=TRUE, col.regions=colori,  scales=list(x=list(rot=90)))
  levelplot(cfeat_corr$p, pretty=TRUE, col.regions=colori,  scales=list(x=list(rot=90)))
}
par(mfrow=c(1,1))

ngruppi=length(idx_groups) # number of groups
counts=vector(, length=ngruppi) # counts will contain, at the end, the number 
# of important features selected for each group
group_length=vector(, length=ngruppi) # here I will have the lengths of each group
ind_fin=0;
feat_sel_groups=list() # list of selected features at each group iteration
idx_features=list()    # list of indeces of the selected features at each iteration
for (gId in 1:ngruppi){
  if (method=='all') {
    c_trn <- trn[c(unlist(idx_groups[gId]), 184)]} else {
      c_trn <- trn[c(unlist(idx_groups[gId]), 178)]}
  c_lm <- glm(group ~., family=binomial, data=c_trn)  # creating a logistic regression model
  c_statistics <- anova(c_lm, test='Chisq')   # using Anova to evaluate feature importance
  c_indici<-unlist(c_statistics$`Pr(>Chi)`)<0.05 # vector of TRUE/FALSE
  c_righe=rownames(c_statistics) # names of the rows (aka names of the variables)
  group_length[gId]=length(unlist(idx_groups[gId])) # number of elements in the group
  c_idx_true <- which(c_indici==TRUE) #extracting indices with p-value < 0.05
  idx_features[gId] <- list((c_idx_true-1) + ind_fin) # original index in the dataset
  print(idx_features[gId])
  ind_fin=ind_fin+group_length[gId]; 
  feat_sel_groups[gId] = list(c_righe[c_idx_true]) # selected features
  counts[gId] <- length(unlist(feat_sel_groups[gId])) # just to control
}

Group_vars <- c(unlist(feat_sel_groups))
Group_vars
Group_idx <- c(unlist(idx_features))
formula_group = group ~ names(trn[,Group_idx])



############################# CLASSIFICATION ##################################
# We tried different classification methods for the data with the features already selected


train.control <- trainControl(method='LOOCV')

### CLASSIFICATION METHOD 1: LOGISTIC REGRESSION

# LASSO VARS
glm_LASSO.train <- train(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75 + ivh_diff_i25_i75 +
                           morph_comp_1 + morph_sphericity + morph_area_dens_aee + ih_skew_fbs_w0.0125 + 
                           ih_iqr_fbs_w0.0125 + ih_max_grad_fbs_w0.0125 + ngl_dcnu_norm_d1_a0.0_3d_fbs_w0.0125, 
                         data=trn, method = 'glm', trControl = train.control)
glm_LASSO.train #Accuracy with LOOCV : 97%
pred_glm_LASSO<-predict(glm_LASSO.train, tst)
confusionMatrix(pred_glm_LASSO, tst$group, positive = "1")
# Accuracy on test set: 100%

# BOTTOM UP VARS
glm_bu.train <- train(group~stat_mean+stat_skew, data=trn, method = 'glm', 
                      trControl = train.control)
glm_bu.train #Accuracy on LOOCV : 92%
pred_glm_bu <- predict(glm_bu.train, tst)
confusionMatrix(pred_glm_bu, as.factor(tst$group), positive = "1") 
# Accuracy = 93.33%

# GROUP VARS
glm_group.train <- train(group~., data=trn[,c(Group_idx,178)], method = 'glm', trControl = train.control)
glm_group.train #Accuracy with LOOCV: 68.42%
pred_glm_group<-predict(glm_group.train, tst)
confusionMatrix(pred_glm_group, tst$group, positive = "1")
# Accuracy : 53.33% 



### CLASSIFICATION METHOD 2: K NEAREST NEIGHBORS - KNN

# LASSO VARS
knn_LASSO.train <- train(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75 + ivh_diff_i25_i75 +
                           morph_comp_1 + morph_sphericity + morph_area_dens_aee + ih_skew_fbs_w0.0125 + 
                           ih_iqr_fbs_w0.0125 + ih_max_grad_fbs_w0.0125 + ngl_dcnu_norm_d1_a0.0_3d_fbs_w0.0125,
                         data=trn, method = 'knn', trControl = train.control)
knn_LASSO.train
pred_knn_LASSO<-predict(knn_LASSO.train, tst)
confusionMatrix(pred_knn_LASSO, tst$group, positive = "1")
#Accuracy=100%


# BOTTOM UP VARS
knn_bu.train <- train(group~stat_mean+stat_skew, data=trn, method = 'knn', 
                      trControl = train.control)
knn_bu.train
pred_knn_bu <- predict(knn_bu.train, tst)
confusionMatrix(pred_knn_bu, as.factor(tst$group), positive = "1")
# Accuracy = 100%


# GROUP VARS
knn_group.train <- train(group~., data=trn[,c(Group_idx,178)], method = 'knn', trControl = train.control)
knn_group.train #k=9
pred_knn_group<-predict(knn_group.train, tst)
confusionMatrix(pred_knn_group, tst$group, positive = "1")
# Accuracy: 93%



# ### CLASSIFICATION METHOD 3: GRADIENT BOOSTING
#
# grid <- expand.grid(list(n.trees=seq(5, 105, by=10),
#                          interaction.depth=1:10,
#                          shrinkage=seq(0.01, 0.1, by=0.01),
#                          n.minobsinnode=1:5))
# 
# # LASSO VARS
# gbm_LASSO.train <- train(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75 + ivh_diff_i25_i75 +
#                            morph_comp_1 + morph_sphericity + morph_area_dens_aee + ih_skew_fbs_w0.0125 +
#                            ih_iqr_fbs_w0.0125 + ih_max_grad_fbs_w0.0125 + ngl_dcnu_norm_d1_a0.0_3d_fbs_w0.0125,
#                           data=trn, method = 'gbm', trControl = train.control, tuneGrid=grid, verbose=F)
# gbm_LASSO.train
# # The final values used for the model were n.trees = 5, interaction.depth = 8,
# # shrinkage = 0.09 and n.minobsinnode = 2.
# summary(gbm_LASSO.train) # plot of the feature importance
# # only stat_skew, stat_qcod and ivh_diff_v25_v75 have been used in the prediction
# pred_gbm_LASSO<-predict(gbm_LASSO.train, tst)
# confusionMatrix(pred_gbm_LASSO, tst$group, positive = "1")
# # Accuracy=93%
# # Let's see if the accuracy improves selecting only the 3 most important features
# gbm_LASSO2.train <- train(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75,
#                          data=trn, method = 'gbm', trControl = train.control, tuneGrid=grid, verbose=F)
# gbm_LASSO2.train
# # The final values used for the model were n.trees = 5, interaction.depth = 1,
# # shrinkage = 0.07 and n.minobsinnode = 5.
# summary(gbm_LASSO2.train) # plot of the feature importance
# # only stat_skew and ivh_diff_v25_v75 have been used in the prediction
# pred_gbm_LASSO2<-predict(gbm_LASSO2.train, tst)
# confusionMatrix(pred_gbm_LASSO2, tst$group, positive = "1")
# # Accuracy: 86.67% DECREASED!!!
# 
# 
# # BOTTOM UP VARS
# gbm_bu.train <- train(group~stat_mean+stat_skew, data=trn, method = 'gbm',
#                        trControl = train.control, tuneGrid=grid, verbose=F)
# gbm_bu.train
# # The final values used for the model were n.trees = 5, interaction.depth = 1,
# # shrinkage = 0.1 and n.minobsinnode = 5.
# summary(gbm_bu.train) # plot of the features importance
# pred_gbm_bu<-predict(gbm_bu.train, tst)
# confusionMatrix(pred_gbm_bu, tst$group, positive = "1")
# # Accurcacy= 93%
# 
# # GROUP VARS
# gbm_group.train <- train(group~., data=trn[,c(Group_idx,178)], method = 'gbm',
#                          trControl = train.control, tuneGrid=grid, verbose=F)
# gbm_group.train
# # Accuracy was used to select the optimal model using the largest value.
# # The final values used for the model were n.trees = 5, interaction.depth = 2, shrinkage = 0.1 and n.minobsinnode = 2.
# summary(gbm_group.train)
# # Only rlm_sre_3d_avg_fbs_w0.0125 was used
# pred_gbm_group<-predict(gbm_group.train, tst)
# confusionMatrix(pred_gbm_group, tst$group, positive = "1")
# # Acccuracy=93.33%


### CLASSIFICATION METHOD 4: RANDOM FOREST

# # LASSO VARS
rf_LASSO.train <- train(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75 + ivh_diff_i25_i75 +
                          morph_comp_1 + morph_sphericity + morph_area_dens_aee + ih_skew_fbs_w0.0125 +
                          ih_iqr_fbs_w0.0125 + ih_max_grad_fbs_w0.0125 + ngl_dcnu_norm_d1_a0.0_3d_fbs_w0.0125,
                        data=trn, method = 'rf', trControl = train.control, tuneLength=7, verbose=F)
rf_LASSO.train #mtry=2
# Implementing the training with the function randomForest cause it allow us
# to see the importance of the features
rf_LASSO.train <- randomForest(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75 + ivh_diff_i25_i75 +
                                 morph_comp_1 + morph_sphericity + morph_area_dens_aee + ih_skew_fbs_w0.0125 +
                                 ih_iqr_fbs_w0.0125 + ih_max_grad_fbs_w0.0125 + ngl_dcnu_norm_d1_a0.0_3d_fbs_w0.0125,
                               data=trn, importance=T, mtry=2)
varImpPlot(rf_LASSO.train)
# stat_skew, ih_skew_fbs_w0.0125, ivh_diff_v25_v75 and stat_qcod are the most important features
pred_rf_LASSO<-predict(rf_LASSO.train, tst)
confusionMatrix(pred_rf_LASSO, tst$group, positive = "1")
# Accuracy=93.33%
#Let's see if the accuracy improves using only the 4 most important features
rf_LASSO2.train <- train(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75 + ih_skew_fbs_w0.0125,
                         data=trn, method = 'rf', trControl = train.control, tuneLength=7, verbose=F)
rf_LASSO2.train # mtry=2
rf_LASSO2.train <- randomForest(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75 + ih_skew_fbs_w0.0125,
                                data=trn, importance=T, mtry=2)
varImpPlot(rf_LASSO2.train)
# stat_skew, ih_skew_fbs_w0.0125, ivh_diff_v25_v75 are the most important features
pred_rf_LASSO2<-predict(rf_LASSO2.train, tst)
confusionMatrix(pred_rf_LASSO2, tst$group, positive = "1")
# Accuracy= 93.33% non è migliorata


# BOTTOM UP VARS
rf_bu.train <- train(group~stat_mean+stat_skew, data=trn, method = 'rf',
                     trControl = train.control, tuneLength=7, verbose=F)
rf_bu.train
#output: mtry=2 (also because we have 2 features)
# Implementing the training with the function randomForest cause it allow us
# to see the importance of the features
rf_bu.train <- randomForest(group~stat_mean+stat_skew, data=trn, importance=T, mtry=2)
pred_rf_bu<-predict(rf_bu.train, tst)
varImpPlot(rf_bu.train)
rf_bu.train$importance #stat_skew more important than stat_mean
confusionMatrix(pred_rf_bu, tst$group, positive = "1")
# Accuracy=93.33%

#
# # GROUP VARS
rf_group.train <- train(group~., data=trn[,c(Group_idx,178)], method = 'rf', trControl = train.control, tuneLength=7, verbose=F)
rf_group.train #mtry=2
rf_group.train <- randomForest(group~., data=trn[,c(Group_idx,178)], importance=T, mtry=2)
# varImpPlot(rf_group.train)
pred_rf_group<-predict(rf_group.train, tst)
confusionMatrix(pred_rf_group, tst$group, positive = "1")
# Accuracy=100%


### CLASSIFICATION METHOD 5: SUPPORT VECTOR MACHINE - SVM

param_grid <- expand.grid(cost = c(0.1, 1, 10, 20),  kernel = c("linear", "polynomial", "radial", "sigmoid"))

# LASSO VARS
svm_lasso <- tune(svm, group ~  stat_skew + stat_qcod + ivh_diff_v25_v75 + morph_comp_2 + morph_com + 
                    morph_area_dens_aee + ih_skew_fbs_w0.0125 + cm_inv_diff_d1_3d_avg_fbs_w0.0125 + 
                    ngl_dcnu_d1_a0.0_3d_fbs_w0.0125, data = trn_zsc, ranges = param_grid)
prediction_svm_class_lasso <- predict(svm_lasso$best.model, tst_zsc)
confusionMatrix(prediction_svm_class_lasso, as.factor(tst_zsc$group), positive = "1")  
#Accuracy: 100%

# BOTTOM UP VARS
svm_bottomup_th0.05 <- tune(svm, group ~ stat_mean + stat_skew, data = trn_zsc, ranges = param_grid)
prediction_class_svm_bottomup_th0.05 <- predict(svm_bottomup_th0.05$best.model, tst_zsc)
confusionMatrix(prediction_class_svm_bottomup_th0.05, as.factor(tst_zsc$group), positive = "1")  
#Accuracy: 93.33%

# GROUP VARS
svm_group <- tune(svm,group~., data=trn_zsc[,c(Group_idx,178)], ranges = param_grid)
prediction_svm_class_group <- predict(svm_group$best.model, tst_zsc)
confusionMatrix(prediction_svm_class_group, as.factor(tst_zsc$group), positive = "1") 
# Accuracy=100%


### CLASSIFICATION METHOD 6: LINEAR DISCRIMINANT ANALYSIS

# LASSO VARS
lda_lasso <- lda(group ~ stat_skew + stat_qcod + ivh_diff_v25_v75 + ivh_diff_i25_i75 +
                   morph_comp_1 + morph_sphericity + morph_area_dens_aee + ih_skew_fbs_w0.0125 + 
                   ih_iqr_fbs_w0.0125 + ih_max_grad_fbs_w0.0125 + ngl_dcnu_norm_d1_a0.0_3d_fbs_w0.0125, data = trn_zsc)
prediction_lda_lasso <- predict(lda_lasso, tst_zsc)
prediction_lda_class_lasso <- as.factor(ifelse(prediction_lda_lasso$x<0.5, 0, 1))
confusionMatrix(prediction_lda_class_lasso, as.factor(tst_zsc$group), positive = "1")  
# Accuracy: 100%

# BOTTOM UP VARS
lda_bottomup_th0.05 <- lda(group ~ stat_mean + stat_skew, data = trn_zsc)
prediction_lda_bottomup_th0.05 <- predict(lda_bottomup_th0.05, tst_zsc)
prediction_class_lda_bottomup_th0.05 <- as.factor(ifelse(prediction_lda_bottomup_th0.05$x<0.5, 0, 1))
confusionMatrix(prediction_class_lda_bottomup_th0.05, as.factor(tst_zsc$group), positive = "1")  
#Accuracy: 93.33%

# GROUP VARS
lda_group <- lda(group~., data=trn_zsc[,c(Group_idx,178)])
#Warning message:
#In lda.default(x, grouping, ...) : presence of collinearity among variables
prediction_lda_group <- predict(lda_group, tst_zsc)
prediction_lda_class_group <- as.factor(ifelse(prediction_lda_group$x<0.5, 0, 1))
confusionMatrix(prediction_lda_class_group, as.factor(tst_zsc$group), positive = "1")  
# Accuracy: 86.67%






############################ STATISTICAL ANALYSIS #################################

##################### ANALYSIS OF FEATURES' NORMALITY #############################
# The Shapiro-Wilk test, often referred to as the Shapiro test, is a statistical test
# used to assess the normality of a dataset. 
# The null hypothesis for the Shapiro-Wilk test is that the data follow a normal 
# distribution. The alternative hypothesis is that the data do not follow a normal
# distribution. 
selected_b <- c()  # Empty vector to store the values of b
for (i in c(2:173, 175:178)) {
  selected_cols <- as.numeric(unlist(cl_data[, i]))
  result <- shapiro.test(selected_cols)
  b <- names(cl_data[, i])
  par(mfrow = c(1,2))
  if (result$p.value < 0.05) {
    selected_b <- c(selected_b, b)  # Adding b to selected_b vector
    # Setting the position and text content for the statistical test result
    print(c('Null hypothesis rejected with a pvalue of ', result$p.value))
  }
  qqnorm(selected_cols, main = names(cl_data[, i]))
  hist(selected_cols, breaks = "FD", probability = TRUE, main = b)
  subtitle <- paste('H-rejected, p-value',format(result$p.value, digits = 5))
  mtext(subtitle, side = 3, line=0, cex = 1, col = 'black')
}
#print(selected_b)
num_NONgaussian <- length(selected_b)+1 #+1 because counting also for the column of ones
print(paste('There are ', num_NONgaussian, ' features among 177 that are not gaussian.'))



######################## WILCOXON RANK SUM TEST ################################
# Wilcoxon test
group <- unlist(cl_data[,179])
radiomic_data <- scale(cl_data[,-c(1,174,179)], scale=T, center = T) # removing outcome from data
features_data <- scale(cl_data[,-c(1,174,179)], scale=T, center = T) # Needing this version to boxplot

# Separating data by outcome
parkinson_data <- radiomic_data[group == 1,]
control_data <- radiomic_data[group == 0,]

i = 0
acc = list()
for (i in c(1:dim(parkinson_data)[2])){
  par(mfrow = c(1,1))
  wilcoxon_test = wilcox.test(parkinson_data[,i], control_data[,i])
  #print(ifelse(wilcoxon_test$p.value<0.05,paste('Null Hypothesys is rejected with a p-value of ',wilcoxon_test$p.value, ' for the feature called ',names(cl_data[,i])), paste('Null Hypothesys is accepted for the feature called ',names(cl_data[,i]))))
  
  # Add significance level to plot
  if (wilcoxon_test$p.value < 0.05) {
    boxplot(features_data[,i] ~ group, ylab = 'Features', main = c('Distribution of radiomic feature',names(cl_data[,i+1])),names=c("Parkinson's", "Controls"))
    abline(h=median(features_data[,i]), col="red", lty=2)
    text(x=1.5, y=median(features_data[,i]), labels="*", cex=2)
    subtitle <- paste('H-rejected, p-value',format(wilcoxon_test$p.value, digits = 5))
    mtext(subtitle, side = 3, line=0, cex = 1, col = 'black')
  }
  else{
    acc = c(acc, names(cl_data[,i+1]))
    boxplot(features_data[,i] ~ group, ylab = 'Features', main = c('Distribution of radiomic feature',names(cl_data[,i+1])),names=c("Parkinson's", "Controls"))
    abline(h=median(features_data[,i]), col="red", lty=2)
    text(x=1.5, y=median(features_data[,i]), labels="*", cex=2)
    subtitle <- paste('H-accepted, p-value',format(wilcoxon_test$p.value, digits = 5))
    mtext(subtitle, side = 3, line=0, cex = 1, col = 'black')
  }
}


################## ROC ANALYSIS (AUC - SENSITIVITY - SPECIFICITY) ################
# Defining a function to compute all the Confusion Matrix Information 
confmatinfo <- function(confusion, pred) {
  # Extract true positive, true negative, false positive, and false negative values
  tp <- confusion[2, 2]  # True Positives
  tn <- confusion[1, 1]  # True Negatives
  fp <- confusion[1, 2]  # False Positives
  fn <- confusion[2, 1]  # False Negatives
  # Calculating sensitivity (true positive rate)
  sensitivity <- tp / (tp + fn)
  # Calculating specificity (true negative rate)
  specificity <- tn / (tn + fp)
  #AUC on test set: 
  roc_curve <- roc(response = tst$group, predictor = as.numeric(pred))
  auc_r = auc(roc_curve)
  plot(roc_curve, main = 'ROC curve',print.auc=TRUE)
  confmatinfo <- data.frame(tp, tn, fp, fn, sensitivity, specificity, auc_r)
}


# LASSO VARS 
# Logistic Regression 
confusion <- table(tst$group,pred_glm_LASSO)
LassoGlm <- confmatinfo(confusion, pred_glm_LASSO)
LassoGlm
Lasso_roc_curve_glm <- roc(response = tst$group, predictor = as.numeric(pred_glm_LASSO))
plot(Lasso_roc_curve_glm , main = 'ROC curve',print.auc=TRUE)
# K-nearest-neighbours
confusion <- table(tst$group,pred_knn_LASSO)
LassoKnn <- confmatinfo(confusion, pred_knn_LASSO)
LassoKnn
Lasso_roc_curve_knn <- roc(response = tst$group, predictor = as.numeric(pred_knn_LASSO))
plot(Lasso_roc_curve_knn, main = 'ROC curve',print.auc=TRUE)
#Gradient Boosting
# confusion <- table(tst$group,pred_gbm_LASSO)
# LassoGbm <- confmatinfo(confusion, pred_gbm_LASSO)
# LassoGbm
# Lasso_roc_curve_gbm <- roc(response = tst$group, predictor = as.numeric(pred_gbm_LASSO))
# plot(Lasso_roc_curve_gbm, main = 'ROC curve',print.auc=TRUE)
# Random Forest
confusion <- table(tst$group,pred_rf_LASSO)
LassoRf <- confmatinfo(confusion, pred_rf_LASSO)
LassoRf
Lasso_roc_curve_rf <- roc(response = tst$group, predictor = as.numeric(pred_rf_LASSO))
plot(Lasso_roc_curve_rf, main = 'ROC curve',print.auc=TRUE)
#Support Vector Machine
confusion <- table(tst$group,prediction_svm_class_lasso)
LassoSvm <- confmatinfo(confusion, prediction_svm_class_lasso)
LassoSvm
Lasso_roc_curve_svm <- roc(response = tst$group, predictor = as.numeric(prediction_svm_class_lasso))
plot(Lasso_roc_curve_svm, main = 'ROC curve',print.auc=TRUE)
# Linear Discriminant Analysis
confusion <- table(tst$group,prediction_lda_class_lasso)
LassoLda <- confmatinfo(confusion, prediction_lda_class_lasso)
LassoLda
Lasso_roc_curve_lda <- roc(response = tst$group, predictor = as.numeric(prediction_lda_class_lasso))
plot(Lasso_roc_curve_lda, main = 'ROC curve',print.auc=TRUE)


# BOTTOM UP VARS
# Logistic Regression
confusion <- table(tst$group,pred_glm_bu)
BottomUpGlm <- confmatinfo(confusion, pred_glm_bu)
BottomUpGlm
Bu_roc_curve_glm <- roc(response = tst$group, predictor = as.numeric(pred_glm_bu))
plot(Bu_roc_curve_glm, main = 'ROC curve',print.auc=TRUE)
# K-nearest-neighbours
confusion <- table(tst$group,pred_knn_bu)
BottomUpKnn <- confmatinfo(confusion, pred_knn_bu)
BottomUpKnn
Bu_roc_curve_knn <- roc(response = tst$group, predictor = as.numeric(pred_knn_bu))
plot(Bu_roc_curve_knn, main = 'ROC curve',print.auc=TRUE)
# Gradient Boosting
# confusion <- table(tst$group,pred_gbm_bu)
# BottomUpGbm <- confmatinfo(confusion, pred_gbm_bu)
# BottomUpGbm
# Bu_roc_curve_gbm <- roc(response = tst$group, predictor = as.numeric(pred_gbm_bu))
# plot(Bu_roc_curve_gbm, main = 'ROC curve',print.auc=TRUE)
# Random Forest
confusion <- table(tst$group,pred_rf_bu)
BottomUpRf <- confmatinfo(confusion, pred_rf_bu)
BottomUpRf
Bu_roc_curve_rf <- roc(response = tst$group, predictor = as.numeric(pred_rf_bu))
plot(Bu_roc_curve_rf, main = 'ROC curve',print.auc=TRUE)
# Support Vector Machine
confusion <- table(tst$group,prediction_class_svm_bottomup_th0.05)
BottomUpSvm <- confmatinfo(confusion, prediction_class_svm_bottomup_th0.05)
BottomUpSvm
Bu_roc_curve_svm <- roc(response = tst$group, predictor = as.numeric(prediction_class_svm_bottomup_th0.05))
plot(Bu_roc_curve_svm, main = 'ROC curve',print.auc=TRUE)
# Linear Discriminant Analysis
confusion <- table(tst$group,prediction_class_lda_bottomup_th0.05)
BottomUpLda <- confmatinfo(confusion, prediction_class_lda_bottomup_th0.05)
BottomUpLda
Bu_roc_curve_lda <- roc(response = tst$group, predictor = as.numeric(prediction_class_lda_bottomup_th0.05))
plot(Bu_roc_curve_lda, main = 'ROC curve',print.auc=TRUE)


#GROUP VARS
confusion <- table(tst$group,pred_glm_group)
GroupGlm <- confmatinfo(confusion, pred_glm_group)
GroupGlm
Group_roc_curve_glm <- roc(response = tst$group, predictor = as.numeric(pred_glm_group))
plot(Group_roc_curve_glm, main = 'ROC curve',print.auc=TRUE)
# K-nearest-neighbours
confusion <- table(tst$group,pred_knn_group)
GroupKnn <- confmatinfo(confusion, pred_knn_group)
GroupKnn
Group_roc_curve_knn <- roc(response = tst$group, predictor = as.numeric(pred_knn_group))
plot(Group_roc_curve_knn, main = 'ROC curve',print.auc=TRUE)
# # Gradient Boosting
# confusion <- table(tst$group,pred_gbm_group)
# GroupGbm <- confmatinfo(confusion, pred_gbm_group)
# GroupGbm
# Group_roc_curve_gbm <- roc(response = tst$group, predictor = as.numeric(pred_gbm_group))
# plot(Group_roc_curve_gbm, main = 'ROC curve',print.auc=TRUE)
# Random Forest
confusion <- table(tst$group,pred_rf_group)
GroupRf <- confmatinfo(confusion, pred_rf_group)
GroupRf
Group_roc_curve_rf <- roc(response = tst$group, predictor = as.numeric(pred_rf_group))
plot(Group_roc_curve_rf, main = 'ROC curve',print.auc=TRUE)
# Support Vector Machine
confusion <- table(tst$group,prediction_svm_class_group)
GroupSvm <- confmatinfo(confusion, prediction_svm_class_group)
GroupSvm
Group_roc_curve_svm <- roc(response = tst$group, predictor = as.numeric(prediction_svm_class_group))
plot(Group_roc_curve_svm, main = 'ROC curve',print.auc=TRUE)
# Linear Discriminant Analysis
confusion <- table(tst$group,prediction_lda_class_group)
GroupLda <- confmatinfo(confusion, prediction_lda_class_group)
GroupLda
Group_roc_curve_lda <- roc(response = tst$group, predictor = as.numeric(prediction_lda_class_group))
plot(Group_roc_curve_lda, main = 'ROC curve',print.auc=TRUE)




############################## DE LONG TEST ####################################
## Calculate the ROC curves for each model and compare them in order to say if there are significant differences in AUC and so if there are some discriminatory abilities.
full_variables <- glm(group ~., data=trn_zsc, family=binomial)
prediction_full_variables <- predict(full_variables, tst_zsc, type='response')
prediction_class_full_variables <- as.factor(ifelse(prediction_full_variables<0.5, 0, 1))

model1_roc <- Bu_roc_curve_glm
auc1 = BottomUpGlm$auc_r
model2_roc <- Lasso_roc_curve_glm
auc2= LassoGlm$auc_r
model3_roc <- Group_roc_curve_glm
auc3= GroupGlm$auc_r
model4_roc <- Bu_roc_curve_svm
auc4= BottomUpSvm$auc_r
model5_roc <- Lasso_roc_curve_svm
auc5= LassoSvm$auc_r
model6_roc <- Group_roc_curve_svm
auc6= GroupSvm$auc_r
model7_roc <- Bu_roc_curve_lda
auc7= BottomUpLda$auc_r
model8_roc <- Lasso_roc_curve_lda
auc8= LassoLda$auc_r
model9_roc <- Group_roc_curve_lda
auc9= GroupLda$auc_r
model10_roc <- Bu_roc_curve_knn
auc10= BottomUpKnn$auc_r
model11_roc <- Lasso_roc_curve_knn
auc11= LassoKnn$auc_r
model12_roc <- Group_roc_curve_knn
auc12= GroupKnn$auc_r
model13_roc <- Bu_roc_curve_rf
auc13 <- BottomUpRf$auc_r
model14_roc <- Lasso_roc_curve_rf
auc14 <- LassoRf$auc_r
model15_roc <- Group_roc_curve_rf
auc_15 <- GroupRf$auc_r
# Prediction with all the variables
model16_roc<- roc(response = tst_zsc$group, predictor = as.numeric(prediction_class_full_variables))
auc16= auc(model16_roc)


# Creating a list for all the models
model_list = list(model1_roc,model2_roc,model3_roc,model4_roc,model5_roc,model6_roc,model7_roc,model8_roc,model9_roc,model10_roc,model11_roc,model12_roc,model13_roc, model14_roc, model15_roc, model16_roc)
delong = list()

for (i in 1:length(model_list)){
  a = i*length(model_list)-length(model_list)
  for(j in 1:length(model_list)){
    b = a + j
    delong[[b]] = roc.test(model_list[[i]], model_list[[j]], method = 'delong', p.adjust.methods= 'bonferroni')
    p_value <- delong[[b]]$p.value
    if (p_value < 0.05) {
      # reject null hypothesis
      cat("Model", i, "is significantly different from model", j, "\n")
    } else {
      # do not reject null hypothesis
      cat("Model", i, "is not significantly different from model", j, "\n")
    }
  }
}

par(mfrow = c(3,3))

plot(model1_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 1st model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model1_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model2_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 2nd model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model2_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model3_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 3rd model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model3_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model4_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 4th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model4_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model5_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 5th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model5_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model6_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 6th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model6_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model7_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 7th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model7_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model8_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 8th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model8_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model9_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 9th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model9_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model10_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 10th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model10_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model11_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 11th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model11_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model12_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 12th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model12_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model13_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 13th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model13_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model14_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 14th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model14_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))

plot(model15_roc, col = 'red',legacy.axes = TRUE, main = 'ROC curve of the 15th model')
lines(model16_roc, col = 'black', lwd = 1, lty = 2)
text(0.3, 0.45, paste("AUC =", round(auc(model15_roc), 3)), adj = c(0, 1), col = 'red')
text(0.3, 0.3, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1))


par(mfrow = c(1,1))

plot(model2_roc, col = '#006400', main = 'ROC curves', legacy.axes = T)
lines(model7_roc, col = '#00B900')
lines(model1_roc, col = '#00FF00')
lines(model9_roc, col = '#FFCC00')
lines(model16_roc, col = 'orange')
lines(model3_roc, col = 'red')


text(0.55, 0.32, paste("AUC =", round(auc(model2_roc), 3)), adj = c(0, 1), col = '#006400')
text(0.55, 0.24, paste('Lasso-LR / Lasso-KNN / Lasso-SVM'),cex = 0.6, adj = c(0, 1), col = '#006400')
text(0.55, 0.17, paste('Lasso-LDA / Bottomup-KNN / Group-RF / Group-SVM'), cex = 0.6, adj = c(0, 1), col = '#006400')
text(0.55, 0.32, paste("AUC =", round(auc(model7_roc), 3)), adj = c(0, 1), col = '#00B900')
text(0.55, 0.24, paste('Bottomup_LDA / Group-GB'), cex = 0.6, adj = c(0, 1), col = '#00B900')
text(0.65, 0.32, paste("AUC =", round(auc(model1_roc), 3)), adj = c(0, 1), col = '#00FF00')
text(0.65, 0.24, paste('Lasso-GB / Lasso-RF / Bottomup_LR'), cex = 0.6, adj = c(0, 1), col = '#00FF46')
text(0.65, 0.17, paste('Bottomup-GB / Bottomup-RF / Bottomup-SVM / Group-KNN'), cex = 0.6, adj = c(0, 1), col = '#00FF00')
text(0.35, 0.32, paste("AUC =", round(auc(model9_roc), 3)), adj = c(0, 1), col = '#FFCC00')
text(0.35, 0.24, paste('Group-LDA'), cex = 0.6, adj = c(0, 1), col = '#FFCC00')
text(0.35, 0.32, paste("AUC =", round(auc(model16_roc), 3)), adj = c(0, 1), col = 'orange')
text(0.35, 0.24, paste('All features - LR'), cex = 0.6, adj = c(0, 1), col = 'orange')
text(0.35, 0.32, paste("AUC =", round(auc(model3_roc), 3)), adj = c(0, 1), col = 'red')
text(0.35, 0.24, paste('Group-LR'), cex = 0.6, adj = c(0, 1),col='red')



############### SECOND PART: ANALYSIS OF THE DISEASE'S SEVERITY #################

########################### CORRELATION ANALYSIS ###############################
# Loading the data which refer to the gravity of Parkison's Disease for subjects with Parkinson
clin_data <- read_excel("Homework_clinical.xlsx")
summary(clin_data) # there are no values that should be changed into factors

### CHECKING OF THE PRESENCE OF MISSING DATA
nrow(clin_data) #33 rows
clin_data <- na.omit(clin_data)
nrow(clin_data) #33 rows
# There are no missing values since the number of rows of the dataset remains unchanged. 
# If there were any, the function na.omit() would have removed the entire line.

set.seed(2)
idx<-sample(nrow(rad_data[1:33,]), floor(0.7*nrow(rad_data[1:33,])))

### Creating some dataplots to visualize the data 
# Checking the correlation among all the metrics in the data
eval_data <- clin_data[idx,-1]
clin_corr <- corr.test(eval_data, method='kendall', adjust='bonferroni')
levelplot(as.matrix(clin_corr$r), col.regions=colori, scales=list(x=list(rot=90)))
levelplot(as.matrix(abs(clin_corr$r)), col.regions=colori, scales=list(x=list(rot=90)))
levelplot(as.matrix(clin_corr$p<0.05/(length(eval_data)^2)), col.regions=colori, scales=list(x=list(rot=90)))


# From literature we know that for the MMSE scale, values >26 represent a 
# person who is normal from the cognitive point of view. From visual evaluation, 
# it doesn't seem that there are any subjects with a score lower than 26. 
# It seems different the situation for the MoCA score.
# I check this hypothesis:
count_MMSE=0
id_cog_imp_MMSE=vector()
count_MoCA=0
id_cog_imp_MoCA=vector()
nsubs=nrow(eval_data)
for (i in seq(1, nsubs)) {
  if (eval_data$MMSE[i] < 26){
    count_MMSE=count_MMSE+1
    id_cog_imp_MMSE= c(id_cog_imp_MMSE, i)
  } else if (eval_data$MoCA[i] < 26) {
    count_MoCA=count_MoCA+1
    id_cog_imp_MoCA= c(id_cog_imp_MoCA, i)}
}
count_MMSE
id_cog_imp_MMSE
count_MoCA
id_cog_imp_MoCA
# So every subject, if with Parkinson's Disease, is cognitively normal for MMSE. 
# This does not happen for the MoCA scale, that recognizes 10 subjects which are not normal.

# # A first division can be done, for subject cognitevely normal and not normal.
# id_cogn_norm <- id_cog_imp_MoCA
# id_cogn_imp <- assex[-id_cogn_norm]
# # We will use this division later to see if there is a way to distinguish subjects that 
# # are cognitevely normal for the ones that are not


# The idea now is: we have got many methods that assess the severity of Parkinson's Disease.
# To see if we can assess something on the Parkinson's gravity thanks to some of 
# our radiomics features, we are checking the correlation among the covariates and
# the given severity measures.

ind1 = c(2:178)
ind2 = c(177:185)
ind3 = c(1:176)
ind4 = 177
ind5 = c((2:173),(175:178))
ind6 = 1
ind7 = c(178:186)
ind8 = ind6
ind9 = 178
ind10 = 179
n_feat <- 177

#method = 'all' #set the method to all to see also the results with the demographic data

if (method == 'all') {
  rad_data <- read_excel('Homework_dati_uniti.xlsx')
  ### CHECKING OF THE PRESENCE OF MISSING DATA
  nrow(rad_data) #392 rows
  rad_data <- na.omit(rad_data)
  nrow(rad_data) #392 <-- no missing data!
  # Changing the variables group and gender into factors!
  rad_data$`Group (1=PD, 0 = Controls)`<- as.factor(rad_data$`Group (1=PD, 0 = Controls)`)
  rad_data$`Gender (M=0, F=1)`<- as.factor(rad_data$`Gender (M=0, F=1)`)
  # Creating the cl_data with all the infos
  cl_data <- rad_data
  colnames(cl_data)[c(181)] = c('gender')
  cl_data$group = cl_data$`Group (1=PD, 0 = Controls)`
  cl_data[179] = NULL
  
  ind1 = c(2:179, 181:184)
  ind2 = c(182:190)
  ind3 = c(1:181)
  ind4 = 182
  ind5 = c((2:173),(175:179),(181:184))
  ind6 = 1
  ind7 = 183:191
  ind8 = c(1,185)
  ind9 = 183
  ind10 = c(180,185)
}



All_data=cbind(cl_data[idx,ind1], eval_data)
trn_all = All_data

# Computing overall matrix of correlation (with all features)
cor_matrix <- corr.test(trn_all[, -173], method='kendall', adjust='bonferroni')
image(as.matrix(cor_matrix$r))
image(as.matrix(cor_matrix$p<0.05/(length(cor_matrix$r))))
levelplot(as.matrix(cor_matrix$r), col.regions = colori, pretty=TRUE,  scales=list(x=list(rot=90)))
# Adding +1 to all the indices above 173 
# Visualizing in 2 colors only the correlations among the 0.35 threshold. 
# This value was chosen because among the features and the severity measures there
# are no correlations above the 0.5 threshold, and only one above the 0.4 threshold
cor_matrix_mod <- ifelse(abs(cor_matrix$r)<0.4, 0, 1)
cor_matrix_mod <- as.matrix(cor_matrix_mod)
image(cor_matrix_mod, axes=FALSE)


# Selecting the rows from the correlation matrix we are interested in. 
# Since I am only interested to the correlation with the severity measures, 
# which are the columns from 177 to 185 (since we have removed the first column)
# we are cutting the correlation matrix for those values.
cor_matrix <- corr.test(trn_all[, -173], method='kendall', adjust='bonferroni')
cor_matrix <- as.matrix(cor_matrix$r)
cor_intr <- cor_matrix[ind2,ind3]
cor_intr_mod <- cor_matrix_mod[ind2,ind3]
levelplot(as.matrix(cor_intr_mod), pretty=TRUE, col.regions = colori, scales=list(x=list(rot=90)))
# extracting the names of the most correlated features

n_metval <- nrow(cor_intr_mod)
idx_corr <- list()
features_corr <- list()
nomi_features <- names(All_data[,-173])
corr_values <- list()
if (method=='all') {n_feat=length(All_data)-n_metval-1} else {n_feat=n_feat-1}#-1 since the column 173 was not considered
for (m in 1:n_metval) {
  cidx=vector()
  cvals=vector()
  for (f in 1:n_feat) {
    if (cor_intr_mod[m,f]==1) {
      cvals = c(cvals, cor_intr[m,f])
      if (f<173) {cidx=c(cidx,f)
      } else {cidx=c(cidx,(f+1))}
    }
  }
  features_corr[m] <- list(nomi_features[cidx])
  idx_corr[m]=list(which(colnames(trn_zsc) %in% unlist(features_corr[m])))
  corr_values[m]=list(cvals)
}
corr_values
idx_corr
features_corr

# Now we are plotting our findings
for (fId in seq(1:length(features_corr))) {
  current_metric=trn_all[, ind4+fId]
  nome_metrica <- names(eval_data)[fId]
  current_features=as.matrix(trn_all[, unlist(idx_corr[fId])])
  curr_corr_values=unlist(corr_values[fId])
  if (length(unlist(idx_corr[fId]))!=0) {
    par(mfrow=(c(1, length(features_corr[fId]))))
    for (pId in 1:length(unlist(idx_corr[fId]))){
      # We don't need to correct by Bonferroni cause we are doing a pairwise correlation
      test_cor <- cor.test(current_features[, pId], current_metric, method='kendall')
      # WARNING MESSAGE:
      # In cor.test.default(current_metric, current_features[, pId], method = "kendall") :
      # Impossibile calcolare p-value esatti in presenza di ties
      # Il p-value è un po' "incerto" poiché ti porta ad accettare/rifutare
      # l'ipotesi cambiando di poco il livello di significatività scelto.
      if (length(unlist(idx_corr[fId])) == 1) {
        nome_var <- unlist(features_corr[fId])
      } else {nome_var <- names(All_data[, unlist(idx_corr[fId])])[pId]}
      plot(current_metric~unlist(current_features[, pId]), main=c('Correlation_value: ', as.character(round(test_cor$estimate, digits=4)), 
                                                                  'p-value: ', as.character(round(test_cor$p.value, digits=4))),
           xlab=nome_var, ylab= nome_metrica)
      abline(lm(unlist(current_metric~current_features[, pId]) ), col='red')
    }
  }
}
# Luckily all the p-values are significant





##################### FEATURE SELECTION TO ANALYSE DISEASE'S SEVERITY  ###########################

# In this section we are selecting the features with the aim of predicting the disease's severity
PD_data<-cl_data[1:33,]
#Z-scoring the data
PD_data_zsc<-PD_data
PD_data_zsc[,ind5] <- scale(PD_data[, ind5], center=TRUE, scale=TRUE)
PD_data<-PD_data_zsc
# NOTE: for the z-scoring the column 174 was omitted. This was done because the 
# column is a list of 1s for every subject in the dataset, so the variance was 
# NULL and it originated a list of NAs. It's important to note that, since the 
# values are the same for every subject, this column would be absolutely useless
# for the analysis.
merged_data <- cbind(PD_data[,-ind10],clin_data[, -1])
summary(merged_data)
merged_data<-merged_data[,-ind6] # removing the Id group at the beginning
colnames(merged_data)[ind7]<-c("LEDD_TOTAL","UPDRS_I","UPDRS_II","UPDRS_III","UPDRS_IV","UPDRS_TOTAL","NMSQ","MMSE","MoCA")
formula_llr = list(LEDD_TOTAL ~ ., UPDRS_I~ .,UPDRS_II~ .,UPDRS_III~ .,UPDRS_IV~ .,UPDRS_TOTAL~ .,NMSQ~ .,MMSE~ .,MoCA~ .)
set.seed(2)
idx<-sample(nrow(merged_data), floor(0.7*nrow(merged_data)))
trn_PD<-merged_data[idx,]
tst_PD<-merged_data[-idx,]


### LASSO REGRESSION FOR DISEASE'S SEVERITY
FEATURES_DS_LASSO = list()
IDX_FEATURES_DS_LASSO = list()
for (evId in 1:(length(clin_data)-1)) {
  set.seed(2)
  y <- trn_PD[, ind4+evId] # defining response variable
  x <- data.matrix(trn_PD[,c(1:ind4)]) # defining matrix of predictor variables (already z-scored)
  x_test <- data.matrix(tst_PD[,c(1:ind4)])
  y_test <- tst_PD[, ind4+evId]
  # Setting alpha = 1 implements lasso regression
  lambdas <- 10^seq(2, -3, by = -.001)
  lasso_reg <- cv.glmnet(x, y, alpha = 1, lambda = lambdas, nfolds=7)
  #lasso_reg <- cv.glmnet(x, y, alpha = 1, nfolds=7)
  # Best 
  lambda_best <- lasso_reg$lambda.min 
  lambda_best
  # Compute R^2 from true and predicted values
  lasso_model <- glmnet(x, y, alpha = 1, lambda = lambda_best)
  coef_lasso <- coef(lasso_model)
  Lasso_idx_UPDRSIII <- which(lasso_model$beta != 0) # indexes of the selected features by Lasso
  if (length(Lasso_idx_UPDRSIII) == 0){
    print(c(names(merged_data)[ind4 + evId], ' does not work'))
    next
  }
  # Extracting the data 
  Lasso_vars_UPDRSIII <-names(trn_PD[,c(Lasso_idx_UPDRSIII)])
  Lasso_vars_UPDRSIII
  # Searching for multicollinearity among the Lasso selected features
  predictions <- predict(lasso_model, newx = x_test)
  test_error <- sqrt(mean(predictions - data.matrix(tst_PD[, ind4+evId]))^2)
  cat("RMSE on test set:", test_error, "\n")  
  Lasso_training<- trn_PD[,c(Lasso_idx_UPDRSIII, ind4+evId)]
  Lg_lasso <- train(formula_llr[[evId]], data=Lasso_training, method='lm', trControl=train.control)
  Lg_lasso <- Lg_lasso$finalModel
  print(Lg_lasso)
  # Computing R^2 adjust
  r_squared <- summary(Lg_lasso)$r.squared
  print(paste('R squared:', summary(Lg_lasso)$r.squared))
  # Get the number of predictors
  num_predictors <- length(Lg_lasso$coefficients) - 1
  # Get the number of observations
  num_observations <- length(Lg_lasso$residuals)
  # Compute the adjusted R-squared value
  adjusted_r_squared <- 1 - (1 - r_squared) * (num_observations - 1) / (num_observations - num_predictors - 1)
  print(paste('R squared adjusted:' , adjusted_r_squared))
  # Computing AIC
  print(paste('AIC: ', AIC(Lg_lasso)))
  print(c(names(merged_data)[ind4 + evId],' works'))
  FEATURES_DS_LASSO[evId] = list(Lasso_vars_UPDRSIII)
  IDX_FEATURES_DS_LASSO[evId] = list(Lasso_idx_UPDRSIII)
}
FEATURES_DS_LASSO
IDX_FEATURES_DS_LASSO


### BOTTOM UP FOR DISEASE'S SEVERITY
nfeat_tot = ind4
nomi = colnames(trn_PD)
alpha=0.05 #alpha threshold
FEATURES_DS_BU = list()
IDX_FEATURES_DS_BU = list()
for (evId in 1:(length(clin_data)-1)) {
  set.seed(2)
  counts_DS=vector(, length=nfeat_tot) 
  feat_sel_groups_DS=list()
  idx_features_DS=list()
  c_idx_feat_DS=vector()
  AIC_feat_DS <- vector(, length=nfeat_tot)
  for (fId in 1:nfeat_tot) {
    c_idx_feat_DS<- c(c_idx_feat_DS, fId)
    c_lm <- glm(formula_llr[[evId]], data=trn_PD[, c(c_idx_feat_DS, (ind4+evId))]) #current logit model
    c_statistics <- anova(c_lm, test='Chisq') # anova analysis
    c_indici<- unlist(c_statistics$`Pr(>Chi)`)<alpha #indexes of TRUE-FALSE with respect to the condition
    c_righe <- rownames(c_statistics) #extracting the names of the overall rows considered
    c_idx_true <- which(c_indici == TRUE) # extracting indexes for which p-value less than alpha
    # if there is no feature which satisfies the condition (it happens with the
    # initial iterations, then the feature with lower p-value is taken)
    if (length( c_idx_true)==0) {
      c_cond <- unlist(c_statistics$`Pr(>Chi)`)==min(c_statistics$`Pr(>Chi)`, na.rm=TRUE)
      c_idx_true <- which(c_cond==TRUE)}
    idx_features_DS[fId] <- list(c_idx_true-1) #selected features at the current iteration
    c_idx_feat_DS <- which(nomi %in% c_righe[c_idx_true]) # extracting the overall right indexes
    feat_sel_groups_DS[fId]=list(c_righe[c_idx_true]) # extracting the right name of the features
    AIC_feat_DS[fId] = AIC(c_lm) # performing AIC
    counts_DS[fId] <- length(unlist(feat_sel_groups_DS[fId])) # just to control
  }
  print(paste('Predicting: ',names(clin_data)[evId+1], ' with BottomUp'))
  # The AIC_feat vector confirms what has been previously found. The lower value 
  # in fact is present at the 5th iteration. 
  Bottomup_vars_DS <- unlist(feat_sel_groups_DS[nfeat_tot])
  print(Bottomup_vars_DS)
  Bottomup_idx_DS <- c_idx_feat_DS #since they are the ones surviving at the last iteration}
  print(Bottomup_idx_DS)
  # It is possible to see, showing the feat_sel_groups voice, that the features 
  # that are selected at every iteration are stat_mean and stat_skew, so what it 
  # happens is that the prediction doesn't improve after 5th iteration.
  best_model <- train(formula_llr[[evId]], data=trn_PD[, c(Bottomup_idx_DS, ind4+evId)], method='lm')
  predictions <- predict(best_model, newdata = tst_PD[,1:ind4])
  test_error <- sqrt(mean(predictions - tst_PD[, ind4+evId])^2)
  cat("RMSE on test set:", test_error, "\n")  #RMSE on test set: 31.39389
  BottomUp_training<- trn_PD[,c(Bottomup_idx_DS, ind4+evId)]
  Bu_model <- train(formula_llr[[evId]], data=BottomUp_training, method='lm', trControl=train.control)
  Bu_model <- Bu_model$finalModel
  print(Bu_model)
  # Computing R^2 adjust
  r_squared <- summary(Bu_model)$r.squared
  print(paste('R squared:', summary(Bu_model)$r.squared))
  # Get the number of predictors
  num_predictors <- length(Bu_model$coefficients) - 1
  # Get the number of observations
  num_observations <- length(Bu_model$residuals)
  # Compute the adjusted R-squared value
  adjusted_r_squared <- 1 - (1 - r_squared) * (num_observations - 1) / (num_observations - num_predictors - 1)
  print(paste('R squared adjusted:' , adjusted_r_squared))
  # Computing AIC
  print(paste('AIC: ', AIC(Bu_model)))
  FEATURES_DS_BU[evId] = list(Bottomup_vars_DS)
  IDX_FEATURES_DS_BU[evId] = list(Bottomup_idx_DS)
}
FEATURES_DS_BU
IDX_FEATURES_DS_BU



###################### LINEAR REGRESSION FROM CORRELATION ANALYSIS##############################
# For this analysis we are selecting the following metrics: 
# UPDRS III (column 5), MMSE (column 8)
CLINICAL_data<-clin_data[,c(1,3,9)]
head(CLINICAL_data)
merged_data2 <- cbind(PD_data[,-ind10],CLINICAL_data[,-c(1,3)])
summary(merged_data2)
merged_data2<-merged_data2[,-ind8] # removing the Id group at the beginning
colnames(merged_data2)[ind9]<-c("UPDRS_I")
print(colnames(merged_data2)[ind9])
# 4: Regression for UPDRS III with features:morph_pca_flatness, morph_vol_dens_aabb, morph_area_dens_conv_hull
idx_corr<-25
trn_PD<-merged_data2[idx,]
tst_PD<-merged_data2[-idx,]
#head(trn_UPDRS_III)
#head(tst_UPDRS_III)
trn_UPDRS_I_idx_corr<-trn_PD[,c(idx_corr,ind9)]
tst_UPDRS_I_idx_corr<-tst_PD[,c(idx_corr,ind9)]
dependent_variable <- "UPDRS_I"
independent_variables <- setdiff(colnames(trn_UPDRS_I_idx_corr), dependent_variable)
lm_model_UPDRS_I_4 <- train(UPDRS_I ~ ., data = trn_UPDRS_I_idx_corr,
                              method = "lm", trControl = train.control)
lm_model_UPDRS_I_4 #RMSE:4.263847   Rsquared: 0.1885762  
lm_model_UPDRS_I_4 <- lm_model_UPDRS_I_4$finalModel 
AIC(lm_model_UPDRS_I_4) #134.5254
# anova(lm_model_UPDRS_I_4)
# drop1(lm_model_UPDRS_I_4,test="F")
predictions <- predict(lm_model_UPDRS_I_4, newdata = tst_UPDRS_I_idx_corr)
test_error <- sqrt(mean((predictions - tst_UPDRS_I_idx_corr[, dependent_variable])^2))
print(lm_model_UPDRS_I_4)
cat("RMSE on test set:", test_error, "\n")  #RMSE on test set:  6.10712


# Analysis for MMSE
merged_data3 <- cbind(PD_data[-ind10],CLINICAL_data[,-c(1,2)])
summary(merged_data3)
merged_data3<-merged_data3[,-ind8] # removing the Id group at the beginning
trn_PD_ivh <- merged_data3[idx,c(15,18,25,73,ind9)]
tst_PD_ivh <- merged_data3[-idx,c(15,18,25,73,ind9)]
lm_ivh_MMSE <- train(MMSE ~ ., data = trn_PD_ivh, method = "lm", trControl = train.control)
lm_ivh_MMSE<- lm_ivh_MMSE$finalModel
AIC(lm_ivh_MMSE) #61.21795
pred_ivh_MMSE <- predict(lm_ivh_MMSE, newdata = tst_PD_ivh)
RMSE_ivh_MMSE <- RMSE(pred_ivh_MMSE, tst_PD_ivh[,2])
RMSE_ivh_MMSE #29.65321

