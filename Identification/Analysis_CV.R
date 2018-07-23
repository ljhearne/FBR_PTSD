# Analysis script for FBR identification

# this script creates a training and testing set of data then does a series 
# of logistics regressions and a LASSO regression and tests each model on the test data set.
# I'm not very well versed in R, so the code is very procedural and less efficient than it could be.

library(caret)
library(e1071)
library(glmnet)
set.seed(17072018)
rm(list=ls())

# INPUTS, PARAMETERS ----------------------------------------------------------------------

Path = "C:\\Users\\HearneL\\AnacondaProjects\\FBR_Screen\\" #path to data and t-test function
Data = read.csv(paste(Path,"Raw_data.csv",sep=""), stringsAsFactors = FALSE)
threshold = 0.5                                #threshold for reg cutoff (I wouldn't change)
perms = 1:500                                #number of permutations
N = nrow(Data)                                 #Sample size
S = round(104*.80,0)                           #train/test split - current set to 70%
#LambdaType = "lambda.1se" #lamba that gives most regularized model such that error is within one SE of minimum
LambdaType = "lambda.min" #lamba that gives minimum mean CV error

Sens_perm  = matrix(data=NA,nrow=max(perms),ncol=4)
Spec_perm  = matrix(data=NA,nrow=max(perms),ncol=4)
Accu_perm  = matrix(data=NA,nrow=max(perms),ncol=4)
LASSO_coef = matrix(data=NA,nrow=83,ncol=max(perms)) #

#organize data
K6 <- Data[,grepl("K6", names(Data))]
GEM <- Data[,grepl("GEM_", names(Data))]
TSCL <- Data[,grepl("TSCL", names(Data))]
PTCI <- Data[,grepl("PTCI", names(Data))]
BSSS <- Data[,grepl("ptsd_screen", names(Data))]    #Individual BSSS regressors
BSSStot <- Data[,grepl("PTSD_SCREEN", names(Data))] #BSSS clinical cutoff
PTSD <- Data[,grepl("ICD_PTSD12m", names(Data))]

# PERMUTATION ANALYSIS --------------------------------------------------------------------

df = data.frame(Data$CODE)
df$Y <- PTSD #dependent variable is PTSD diagnosis

for (i in perms){
  print(i)
  trainIDX = sample(1:N,S) #the same training IDX is used across all tests.
  
  t <-sum(df$Y[-trainIDX])  #ensure a certain level of PTSD balance
  while((t< 5 || t > 15)){
    trainIDX = sample(1:N,S)
    t <-sum(df$Y[-trainIDX])
  }
  
  #BSSS clincal
  df$X <- as.matrix(BSSStot)    #the predictor(s)
  testData    = df[-trainIDX, ] #test predictors(x) and outcomes (y)
  trainData   = df[trainIDX, ]  #training "" ""
  model       = glm(Y~X,family="binomial",data=trainData)           #do logistic regression on trainData
  result = (predict(model, testData,type = "response")>threshold)*1 #apply to new data
  cm = confusionMatrix(as.factor(result), as.factor(testData$Y), positive="1") #create Confusion matrix
  Sens_perm[i,1] = cm$byClass['Sensitivity']   #save results
  Spec_perm[i,1] = cm$byClass['Specificity']
  Accu_perm[i,1] = cm$byClass['Balanced Accuracy']
  
  #BSSS item regression
  df$X = as.matrix(BSSS)
  testData    = df[-trainIDX, ]
  trainData   = df[trainIDX, ]
  model       = glm(Y~X,family="binomial",data=trainData)
  result = (predict(model, testData,type = "response")>threshold)*1
  cm = confusionMatrix(as.factor(result), as.factor(testData$Y), positive="1") 
  Sens_perm[i,2] = cm$byClass['Sensitivity']
  Spec_perm[i,2] = cm$byClass['Specificity']
  Accu_perm[i,2] = cm$byClass['Balanced Accuracy']
  
  #K6/GEM
  df$X = as.matrix(cbind(K6,GEM))
  testData    = df[-trainIDX, ]
  trainData   = df[trainIDX, ]
  model       = glm(Y~X,family="binomial",data=trainData)
  result = (predict(model, testData,type = "response")>threshold)*1
  cm = confusionMatrix(as.factor(result), as.factor(testData$Y), positive="1") 
  Sens_perm[i,3] = cm$byClass['Sensitivity']
  Spec_perm[i,3] = cm$byClass['Specificity']
  Accu_perm[i,3] = cm$byClass['Balanced Accuracy']
  
  #LASSO regression with all items
  df$X <- as.matrix(cbind(BSSS,TSCL,PTCI,K6,GEM))
  testData    <- df[-trainIDX, ]
  trainData   <- df[trainIDX, ]
  
  cvfit = cv.glmnet(trainData$X, trainData$Y,family = "binomial",nfolds=5) # LASSO regression: lamba is determined by 5 CV
  LASSO_coef[,i] = as.matrix(coef(cvfit, s = LambdaType)) #write the coeffcients
  result = (predict(cvfit, newx = testData$X, s = LambdaType,type="response")>threshold)*1
  cm <- confusionMatrix(as.factor(result), as.factor(testData$Y), positive="1") 
  Sens_perm[i,4] = cm$byClass['Sensitivity']
  Spec_perm[i,4] = cm$byClass['Specificity']
  Accu_perm[i,4] = cm$byClass['Balanced Accuracy']
}

# SUMMARY & STATISTICS --------------------------------------------------------------------

# LASSO items
items = rowSums(abs(LASSO_coef)>0)
labels = rownames(coef(cvfit, s = LambdaType))
LAS <- data.frame(items,labels)
LAS <- LAS[order(-LAS$items),]
print(paste0("LASSO variables in order ",LAS$labels[2:11]))
print(paste0("LASSO variables consistency ",LAS$items[2:11]/max(perms)))

# Statistics
library(reshape2)
source(paste(Path,"pairwisettestwithtanddf.R",sep=""))
lab <- c("BSSS clinical","BSSS","GEM+K6","LASSO")

# repeated-measures setup
SensModel <- data.frame(Sens_perm) # Sensitivity
names(SensModel) = lab
SensModel$ID = 1:nrow(SensModel)
SensModel <- melt(SensModel,id.vars = "ID")
model <- aov(value ~ variable, data = SensModel)
summary(model) #print results
SensStats = pairwisettestwithtanddf(SensModel$value,SensModel$variable, p.adjust.method = "bonferroni") #t-tests

SpecModel <- data.frame(Spec_perm) #Specificity
names(SpecModel) = lab
SpecModel$ID = 1:nrow(SpecModel)
SpecModel <- melt(SpecModel,id.vars = "ID")
model <- aov(value ~ variable, data = SpecModel)
summary(model)
SpecStats = pairwisettestwithtanddf(SpecModel$value,SpecModel$variable, p.adjust.method = "bonferroni")

# PLOTS ---------------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)

lims <- c(0, 1.05) #yaxis limits

dat <- stack(as.data.frame(Sens_perm))
means <- colMeans(Sens_perm)
p1 <-ggplot(dat) + 
  geom_jitter(aes(x = ind, y = values),alpha = 0.25,size=0.25) +
  geom_boxplot(aes(x = ind, y = values),width = 0.3,fill="lightsteelblue1",lwd=0.5,outlier.shape = NA) +
  labs(x = "",y = "Sensitivity") +
  scale_x_discrete(labels= lab) +
  scale_y_continuous(limits = lims) +
  geom_text(x=1,y=1.05,label = round(means[1],2),size = 3.5) +
  geom_text(x=2,y=1.05,label = round(means[2],2),size = 3.5) +
  geom_text(x=3,y=1.05,label = round(means[3],2),size = 3.5) +
  geom_text(x=4,y=1.05,label = round(means[4],2),size = 3.5) +
  theme_classic()

dat <- stack(as.data.frame(Spec_perm))
means <- colMeans(Spec_perm)
p2 <-ggplot(dat) + 
  geom_jitter(aes(x = ind, y = values),alpha = 0.25,size=0.25) +
  geom_boxplot(aes(x = ind, y = values),width = 0.3,fill="lightsteelblue1",lwd=0.5,outlier.shape = NA) +
  labs(x = "",y = "Specificity") +
  scale_x_discrete(labels= lab) +
  scale_y_continuous(limits = lims) +
  geom_text(x=1,y=1.05,label = round(means[1],2),size = 3.5) +
  geom_text(x=2,y=1.05,label = round(means[2],2),size = 3.5) +
  geom_text(x=3,y=1.05,label = round(means[3],2),size = 3.5) +
  geom_text(x=4,y=1.05,label = round(means[4],2),size = 3.5) +
  theme_classic()

dat <- stack(as.data.frame(Accu_perm))
p3 <-ggplot(dat) + 
  geom_boxplot(aes(x = ind, y = values),width = 0.4,fill="dodgerblue1",lwd=1) +
  labs(x = "",y = "Balanced Accuracy") +
  scale_x_discrete(labels= lab) +
  scale_y_continuous(limits = lims) +
  theme_classic()

grid.arrange(p1,p2,ncol=2)