#
#
#

library(caret)
library(e1071)
Data <- read.csv("C:\\Users\\HearneL\\AnacondaProjects\\FBR_Screen\\Raw_data.csv", stringsAsFactors = FALSE)
K6 <- Data[,grepl("K6", names(Data))]
GEM <- Data[,grepl("GEM_", names(Data))]
TSCL <- Data[,grepl("TSCL", names(Data))]
PTCI <- Data[,grepl("PTCI", names(Data))]
BSSS <- Data[,grepl("ptsd_screen", names(Data))]

# 1) Logistic regression BSSS
#######################################
y <- Data$ICD_PTSD12m
x <- Data$PTSD_SCREEN

mod <- glm(y~x, family="binomial")
nullmod <- glm(y~1, family="binomial")
R2 <- 1-logLik(mod)/logLik(nullmod) #McFadden R2
summary(mod)
anova (mod,nullmod,test="Chisq") #Chi Squared stats

threshold=0.5
predicted_values<-ifelse(predict(mod,type="response")>threshold,1,0)
actual_values<-mod$y
confusionMatrix(as.factor(predicted_values), as.factor(y), positive="1") 

# 2) Logistic regression GEM and K6
#######################################
x <- cbind(K6,GEM)
x <- as.matrix(x)

mod2 <- glm(y~x, family="binomial")
R2 <- 1-logLik(mod2)/logLik(mod) #McFadden R2 with the BSSS as the null model
R2
summary(mod2)
anova (mod2,mod,test="Chisq") #Chi Squared stats

threshold=0.5
predicted_values<-ifelse(predict(mod2,type="response")>threshold,1,0)
actual_values<-mod2$y
confusionMatrix(as.factor(predicted_values), as.factor(y), positive="1") 

# 3) LASSO regression with all variables.
#######################################
library(glmnet)

# Lasso penalized regression with all scales.
x <- cbind(K6,GEM,TSCL,PTCI,BSSS)
x <- as.matrix(x)
set.seed(1)
lass <- glmnet(x, y, family = "binomial", alpha = 1)
plot(lass)

#lam <-0.075 # lam is set manually to derive 10 items (~ number of useful items)
lam <-0.13 # 5 items
lasso_coef = coef(lass, lam)

predicted_values <- predict(lass, x, s = lam, type = "class")
actual_values<-mod2$y
confusionMatrix(as.factor(predicted_values), as.factor(y), positive="1") 

tLL <- lass$nulldev - deviance(lass)
k <- lass$df
n <- lass$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
plot(AICc)

