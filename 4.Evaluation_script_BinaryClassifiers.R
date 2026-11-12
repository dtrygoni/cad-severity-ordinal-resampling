####################### Point estimates: evaluation of resampling in 1. CAD presence, 2. CAD severity, 3. CAD High severity ###############

#Libraries
library(caret)
library(pROC)
#Functions
Evaluation_binary_fnc<-function(model,test,target,train){
  #Metrics: Balanced Accuracy, Sensitivity, Specificity, ROC-AUC
  X_test<-test[predictors(model)]
  y_test<-test[,target]
  
  X_train<-train[predictors(model)]
  y_train<-train[,target]
  
  
  # Targets as factors
  y_train <- factor(y_train, levels = c(0,1),labels=c("neg","pos"))
  y_test  <- factor(y_test,  levels = c(0,1),labels=c("neg", "pos"))

  # Predicted probabilities on train and test
  train_prob <- predict(model, X_train, type = "prob")[, "pos"]
  test_prob  <- predict(model, X_test,  type = "prob")[, "pos"]
  
  # ROC on training set
  roc_train <- pROC::roc(
    response = y_train,
    predictor = train_prob,
    levels = c("neg", "pos"),
    direction = "<"
  )
  
  thr <- pROC::coords(
    roc_train,
    x = "best",
    best.method = "youden",
   ret = "threshold"
  )
  if(length(as.numeric(thr$threshold))==1){
    thr <- as.numeric(thr)
  }else{
    thr <- as.numeric(thr$threshold[1])
    
  }
  # Apply threshold to test probabilities
  y_pred <- ifelse(test_prob >= thr, "pos", "neg")
  y_pred <- factor(y_pred, levels = c("neg", "pos"))
  
  cm <- confusionMatrix(
    data = y_pred,
    reference = y_test,
    positive = "pos"
  )
  
  # Metrics
  sensitivity <- round(as.numeric(cm$byClass["Sensitivity"]),3)
  specificity <- round(as.numeric(cm$byClass["Specificity"]),3)
  bal_acc     <- round(as.numeric(cm$byClass["Balanced Accuracy"]),3)
  
  roc_test <- pROC::roc(
    response = y_test,
    predictor = test_prob,
    levels = c("neg", "pos"),
    direction = "<"
  )
  auc_val <- round(as.numeric(pROC::auc(roc_test)),3)
  
  
  res<-data.frame(Balanced.Accuracy=bal_acc,Sensitivity=sensitivity,Specificity=specificity,ROC.AUC=auc_val)
  return(res)
  }
#1. CAD presence
load('./Models/Models_CADP.Rda')
res_all<-NULL
target_col<-'SYNTAX.BINARY'
classifiers<-c("RF","SVM","kNN","LR")
resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
for(jj in 1:10){
  test<-read.csv(paste0('./Train_test_sets/Test_set_CADP_orig_',jj,'.csv'))
  for(resamp in resamplings){
    train<-read.csv(paste0('./Train_test_sets/Train_set_CADP_',resamp,'_',jj,'.csv'))
    for(classifier in classifiers){
      model_wanted<-Models_all[[jj]][[resamp]][[classifier]]
      res<-Evaluation_binary_fnc(model_wanted,test,target_col,train)
      res$Resampling<-resamp
      res$Fold<-jj
      res$Classifier<-classifier
      res_all<-rbind(res_all,res)
    }
  }
}
res_presence<-res_all
save(res_presence,file='./Results/Binary_CAD_presence.Rda')




#2. CAD severity
load('./Models/Models_CADS.Rda')
res_all<-NULL
target_col<-'SYNTAX.SEVERITY'
classifiers<-c("RF","SVM","kNN","LR")
resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
for(jj in 1:10){
  test<-read.csv(paste0('./Train_test_sets/Test_set_CADS_orig_',jj,'.csv'))
  for(resamp in resamplings){
    train<-read.csv(paste0('./Train_test_sets/Train_set_CADS_',resamp,'_',jj,'.csv'))
    for(classifier in classifiers){
      model_wanted<-Models_all[[jj]][[resamp]][[classifier]]
      res<-Evaluation_binary_fnc(model_wanted,test,target_col,train)
      res$Resampling<-resamp
      res$Fold<-jj
      res$Classifier<-classifier
      res_all<-rbind(res_all,res)
    }
  }
}
res_severity<-res_all
save(res_severity,file='./Results/Binary_CAD_severity.Rda')



#3. CAD severity
load('./Models/Models_CADHS.Rda')
res_all<-NULL
target_col<-'SYNTAX.HIGH.SEVERITY'
classifiers<-c("RF","SVM","kNN","LR")
resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
for(jj in 1:10){
  test<-read.csv(paste0('./Train_test_sets/Test_set_CADHS_orig_',jj,'.csv'))
  for(resamp in resamplings){
    train<-read.csv(paste0('./Train_test_sets/Train_set_CADHS_',resamp,'_',jj,'.csv'))
    for(classifier in classifiers){
      model_wanted<-Models_all[[jj]][[resamp]][[classifier]]
      res<-Evaluation_binary_fnc(model_wanted,test,target_col,train)
      res$Resampling<-resamp
      res$Fold<-jj
      res$Classifier<-classifier
      res_all<-rbind(res_all,res)
    }
  }
}
res_high_severity<-res_all
save(res_high_severity,file='./Results/Binary_CAD_high_severity.Rda')















############################################## For repeats #################################
for(rep in 1:10){
  print(rep)
  #1. CAD presence
  load(paste0('./Repeats/Rep',rep,'/Models/Models_CADP.Rda'))
  res_all<-NULL
  target_col<-'SYNTAX.BINARY'
  classifiers<-c("RF","SVM","kNN","LR")
  resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
  for(jj in 1:10){
    test<-read.csv(paste0('./Repeats/Rep',rep,'/Train_test_sets/Test_set_CADP_orig_',jj,'.csv'))
    for(resamp in resamplings){
      train<-read.csv(paste0('./Repeats/Rep',rep,'/Train_test_sets/Train_set_CADP_',resamp,'_',jj,'.csv'))
      for(classifier in classifiers){
        model_wanted<-Models_all[[jj]][[resamp]][[classifier]]
        res<-Evaluation_binary_fnc(model_wanted,test,target_col,train)
        res$Resampling<-resamp
        res$Fold<-jj
        res$Classifier<-classifier
        res_all<-rbind(res_all,res)
      }
    }
  }
  res_presence<-res_all
  save(res_presence,file=paste0('./Repeats/Rep',rep,'/Results/Binary_CAD_presence.Rda'))
  
  
  
  
  #2. CAD severity
  load(paste0('./Repeats/Rep',rep,'/Models/Models_CADS.Rda'))
  res_all<-NULL
  target_col<-'SYNTAX.SEVERITY'
  classifiers<-c("RF","SVM","kNN","LR")
  resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
  for(jj in 1:10){
    test<-read.csv(paste0('./Repeats/Rep',rep,'/Train_test_sets/Test_set_CADS_orig_',jj,'.csv'))
    for(resamp in resamplings){
      train<-read.csv(paste0('./Repeats/Rep',rep,'/Train_test_sets/Train_set_CADS_',resamp,'_',jj,'.csv'))
      for(classifier in classifiers){
        model_wanted<-Models_all[[jj]][[resamp]][[classifier]]
        res<-Evaluation_binary_fnc(model_wanted,test,target_col,train)
        res$Resampling<-resamp
        res$Fold<-jj
        res$Classifier<-classifier
        res_all<-rbind(res_all,res)
      }
    }
  }
  res_severity<-res_all
  save(res_severity,file=paste0('./Repeats/Rep',rep,'/Results/Binary_CAD_severity.Rda'))
  
  
  
  #3. CAD severity
  load(paste0('./Repeats/Rep',rep,'/Models/Models_CADHS.Rda'))
  res_all<-NULL
  target_col<-'SYNTAX.HIGH.SEVERITY'
  classifiers<-c("RF","SVM","kNN","LR")
  resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
  for(jj in 1:10){
    test<-read.csv(paste0('./Repeats/Rep',rep,'/Train_test_sets/Test_set_CADHS_orig_',jj,'.csv'))
    for(resamp in resamplings){
      train<-read.csv(paste0('./Repeats/Rep',rep,'/Train_test_sets/Train_set_CADHS_',resamp,'_',jj,'.csv'))
      for(classifier in classifiers){
        model_wanted<-Models_all[[jj]][[resamp]][[classifier]]
        res<-Evaluation_binary_fnc(model_wanted,test,target_col,train)
        res$Resampling<-resamp
        res$Fold<-jj
        res$Classifier<-classifier
        res_all<-rbind(res_all,res)
      }
    }
  }
  res_high_severity<-res_all
  save(res_high_severity,file=paste0('./Repeats/Rep',rep,'/Results/Binary_CAD_high_severity.Rda'))
  
}
