# Analysis script for FBR identification
# reads in data then performs a few logistic regressions and
# a lasso regression using GLMnet.

library(caret)
library(e1071)
library(glmnet)
set.seed(17072018)
rm(list=ls())

# read in data and give appropriate variable names
Data <- read.csv("C:\\Users\\HearneL\\AnacondaProjects\\FBR_Screen\\Raw_data.csv", stringsAsFactors = FALSE)
K6 <- Data[,grepl("K6", names(Data))]
GEM <- Data[,grepl("GEM_", names(Data))]
TSCL <- Data[,grepl("TSCL", names(Data))]
PTCI <- Data[,grepl("PTCI", names(Data))]
BSSS <- Data[,grepl("ptsd_screen", names(Data))]
BSSStot <- Data[,grepl("PTSD_SCREEN", names(Data))]
PTSD <- Data[,grepl("ICD_PTSD12m", names(Data))]

# parameters & preallocation
threshold=0.5
perms <- 1:500
N <- 104
Sens_perm <- matrix(data=NA,nrow=max(perms),ncol=4)
Spec_perm <- matrix(data=NA,nrow=max(perms),ncol=4)
Accu_perm <- matrix(data=NA,nrow=max(perms),ncol=4)
LASSO_coef <- matrix(data=NA,nrow=83,ncol=max(perms))
S <- round(104*.90,0) # percentage to split train/test

#-[1]--- 
df = data.frame(Data$CODE)
df$Y <- PTSD

for (i in perms){
  print(i)
  trainIDX = sample(1:N,S) #same training IDX is used across all tests.
  
  # ensure a certain level of PTSD balance
  t <-sum(df$Y[-trainIDX])

  while((t< 3 || t > 7)){
    trainIDX = sample(1:N,S)
    t <-sum(df$Y[-trainIDX])
  }
  
  #--- BSSS clincal
  df$X <- as.matrix(BSSStot) #the predictor(s)
  
  testData    <- df[-trainIDX, ] #test predictors(x) and outcomes (y)
  trainData   <- df[trainIDX, ]  #training "" ""
  model       <- glm(Y~X,family="binomial",data=trainData) # do logistic regression on trainData
  result <- (predict(model, testData,type = "response")>threshold)*1 #apply to new data
  cm <- confusionMatrix(as.factor(result), as.factor(testData$Y), positive="1") #create Confusion matrix
  
  #save results
  Sens_perm[i,1] <-cm$byClass['Sensitivity'] 
  Spec_perm[i,1] <-cm$byClass['Specificity']
  Accu_perm[i,1] <-cm$byClass['Balanced Accuracy']
  
  #--- BSSS regression
  df$X <- as.matrix(BSSS)
  
  testData    <- df[-trainIDX, ]
  trainData   <- df[trainIDX, ]
  model       <- glm(Y~X,family="binomial",data=trainData)
  result <- (predict(model, testData,type = "response")>threshold)*1
  cm <- confusionMatrix(as.factor(result), as.factor(testData$Y), positive="1") 
  
  Sens_perm[i,2] <-cm$byClass['Sensitivity']
  Spec_perm[i,2] <-cm$byClass['Specificity']
  Accu_perm[i,2] <-cm$byClass['Balanced Accuracy']
  
  #--- K6/GEM
  df$X <- as.matrix(cbind(K6,GEM))
  
  testData    <- df[-trainIDX, ]
  trainData   <- df[trainIDX, ]
  model       <- glm(Y~X,family="binomial",data=trainData)
  result <- (predict(model, testData,type = "response")>threshold)*1
  cm <- confusionMatrix(as.factor(result), as.factor(testData$Y), positive="1") 
  
  Sens_perm[i,3] <-cm$byClass['Sensitivity']
  Spec_perm[i,3] <-cm$byClass['Specificity']
  Accu_perm[i,3] <-cm$byClass['Balanced Accuracy']
  
  #--- LASSO w/all
  df$X <- as.matrix(cbind(BSSS,TSCL,PTCI,K6,GEM))
  
  testData    <- df[-trainIDX, ]
  trainData   <- df[trainIDX, ]
  
  #model <- glmnet(trainData$X,trainData$Y, family = "binomial", alpha = 1)
  cvfit = cv.glmnet(trainData$X, trainData$Y,family = "binomial",nfolds=5) # the lamba is determined by cross validation
  
  #lambda.1se, which gives the most regularized model such that error is within one standard error of the minimum.
  LASSO_coef[,i] <- as.matrix(coef(cvfit, s = "lambda.min"))
  result <-(predict(cvfit, newx = testData$X, s = "lambda.min",type="response")>threshold)*1
  cm <- confusionMatrix(as.factor(result), as.factor(testData$Y), positive="1") 
  
  Sens_perm[i,4] <-cm$byClass['Sensitivity']
  Spec_perm[i,4] <-cm$byClass['Specificity']
  Accu_perm[i,4] <-cm$byClass['Balanced Accuracy']
}

## LASSO - what items are included?
items = rowSums(abs(LASSO_coef)>0)
labels = rownames(coef(cvfit, s = "lambda.min"))
LAS <- data.frame(items,labels)
LAS <- LAS[order(-LAS$items),]

print(paste0("LASSO variables in order ",LAS$labels[1:16]))
print(paste0("LASSO variables consistency ",LAS$items[1:16]/max(perms)))

## Comparison Statistics
# repeated measures ANOVA.
write.csv(Sens_perm, file = "SENS.csv") #save the weights
write.csv(Spec_perm, file = "SPEC.csv") #save the weights

## PLOT
library(ggplot2)
library(gridExtra)


lims <- c(0, 1)
lab <- c("BSSS clinical","BSSS","GEM+K6","LASSO")
dat <- stack(as.data.frame(Sens_perm))
p1 <-ggplot(dat) + 
  geom_boxplot(aes(x = ind, y = values),width = 0.4,fill="gray90",lwd=0.5) +
  labs(x = "",y = "Sensitivity") +
  scale_x_discrete(labels= lab) +
  scale_y_continuous(limits = lims) +
  theme_classic()

dat <- stack(as.data.frame(Spec_perm))
p2 <-ggplot(dat) + 
  geom_boxplot(aes(x = ind, y = values),width = 0.4,fill="gray90",lwd=0.5) +
  labs(x = "",y = "Specificity") +
  scale_x_discrete(labels= lab) +
  scale_y_continuous(limits = lims) +
  theme_classic()

dat <- stack(as.data.frame(Accu_perm))
p3 <-ggplot(dat) + 
  geom_boxplot(aes(x = ind, y = values),width = 0.4,fill="gray90",lwd=0.5) +
  labs(x = "",y = "Balanced Accuracy") +
  scale_x_discrete(labels= lab) +
  scale_y_continuous(limits = lims) +
  theme_classic()

grid.arrange(p1,p2,ncol=2)

