####################### SHAP values code #################

### Libraries ###
library(caret)
library(shapviz)
library(ggplot2)
library(kernelshap)
### Load the data ###
load("./Models/Models_CADP_whole.Rda")
model<-Models_all$Fold1$rus$SVM

data<-read.csv("./Train_test_sets/Whole_set_CADP_rus_1.csv")
target<-"SYNTAX.BINARY"

# IMPORTANT: return numeric vector (probability for ONE class)
pfun <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, 1]
}
set.seed(123)

# Compute Shapley values
X<-data[,-14]
shap<-kernelshap(object=model,
                 X=X,
                 pred_fun = pfun,
                 bg_X=X[sample(nrow(X),50),])



# ----------------------------
# 5. Create shapviz object
# ----------------------------
sv <- shapviz(shap)


###############################
new_names <- c(
  "CHEST.PAIN" = "CHEST PAIN",
  "ATRIAL.FIBRILLATION" = "ATRIAL FIBRILLATION",
  "ATYPICAL.ANGINA.SYMPTOMS" = "ATYPICAL ANGINA SYMPTOMS",
  "EOS."="EOS%"
)



p<-sv_importance(sv, kind = "bar") +
  scale_y_discrete(labels=new_names)+
  theme(
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 14),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 16)
  )



library(tidyverse)
tiff('./Desktop/Shap_CADP.jpeg',units="in",width=10,height=4,res=300)
p
dev.off()







### Load the data ###
load("./Models/Models_CADS_whole.Rda")
model<-Models_all$Fold1$rus$SVM

data<-read.csv("./Train_test_sets/Whole_set_CADS_rus_1.csv")
target<-"SYNTAX.SEVERITY"

# IMPORTANT: return numeric vector (probability for ONE class)
pfun <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, 1]
}
set.seed(123)

# Compute Shapley values
X<-data[,-9]
shap<-kernelshap(object=model,
                 X=X,
                 pred_fun = pfun,
                 bg_X=X[sample(nrow(X),50),])



# ----------------------------
# 5. Create shapviz object
# ----------------------------
sv <- shapviz(shap)


###############################
new_names <- c(
  "DIABETES.MELLITUS" = "DIABETES MELLITUS",
  "CRUSADE.SCORE" = "CRUSADE SCORE",
  "CHEST.PAIN"="CHEST PAIN"
)



p<-sv_importance(sv, kind = "bar") +
  scale_y_discrete(labels=new_names)+
  theme(
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 14),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 16)
  )



library(tidyverse)
tiff('./Desktop/Shap_CADS.jpeg',units="in",width=10,height=4,res=300)
p
dev.off()





### Load the data ###
load("./Models/Models_CADHS_whole.Rda")
model<-Models_all$Fold1$rus$SVM

data<-read.csv("./Train_test_sets/Whole_set_CADHS_rus_1.csv")
target<-"SYNTAX.HIGH.SEVERITY"

# IMPORTANT: return numeric vector (probability for ONE class)
pfun <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, 1]
}
set.seed(123)

# Compute Shapley values
X<-data[,-7]
shap<-kernelshap(object=model,
                 X=X,
                 pred_fun = pfun,
                 bg_X=X[sample(nrow(X),50),])



# ----------------------------
# 5. Create shapviz object
# ----------------------------
sv <- shapviz(shap)


###############################
new_names <- c(
  "DIABETES.MELLITUS" = "DIABETES MELLITUS",
  "CRUSADE.SCORE" = "CRUSADE SCORE"
  #"CHEST.PAIN"="CHEST PAIN"
)



p<-sv_importance(sv, kind = "bar") +
  scale_y_discrete(labels=new_names)+
  theme(
    axis.text.x  = element_text(color = "black", size = 12),
    axis.text.y  = element_text(color = "black", size = 14),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 16)
  )



library(tidyverse)
tiff('./Desktop/Shap_CADHS.jpeg',units="in",width=10,height=4,res=300)
p
dev.off()
