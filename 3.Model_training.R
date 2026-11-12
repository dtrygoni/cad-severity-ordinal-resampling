####################### Code for model training for each .csv #####################
#### Hyperparameter random optimization 50 4-fold cross validation with metric ROC-AUC

#Libraries
library(caret)
library(MLmetrics)
library(randomForest)
library(pROC)
library(mltools)
library(dplyr)
set.seed(1)
#Functions
summary_fnc <- function(data, lev = NULL, model = NULL) {
  # Ensure prediction and true labels are factors
  obs <- factor(data$obs, levels = lev)
  pred <- factor(data$pred, levels = lev)
  
  cm <- caret::confusionMatrix(pred, obs)
  sens <- cm$byClass["Sensitivity"]
  spec <- cm$byClass["Specificity"]
  
  bal_acc <- mean(c(sens, spec), na.rm = TRUE)
  acc <- MLmetrics::Accuracy(pred, obs)
  f1 <- FBeta_Score(obs, pred, positive = "pos", beta = 1)
  f2 <- FBeta_Score(obs, pred, positive = "pos", beta = 2)
  f05 <- FBeta_Score(obs, pred, positive = "pos", beta = 0.5)
  
  # AUC Calculation (assumes the column name for positive class is 'pos')
  if ("pos" %in% colnames(data)) {
    auc <- pROC::roc(response = obs, predictor = data[["pos"]], levels = rev(lev), direction = ">")$auc
  } else {
    auc <- NA
  }
  
  out <- c(
    balAcc = bal_acc,
    F1Score = f1,
    Spec = spec,
    Sens = sens,
    Accuracy = acc,
    F2Score = f2,
    F05Score = f05,
    AUC = auc
  )
  return(out)
}
Model_training<-function(clinical_train,target_col,model_selection='RF',number_iterations=50,optimization_metric='AUC'){
  set.seed(123)
  
  input_x_clinical <- as.matrix(dplyr::select(clinical_train, -c(target_col)))
  input_y_clinical <- clinical_train[,target_col]
  input_y_clinical <- factor(input_y_clinical, levels = c(0, 1), labels = c("neg", "pos"))
  
  ctrl <- trainControl(
    method = "cv",
    number = 4,
    search = "random",
    summaryFunction = summary_fnc,
    classProbs = TRUE,
    verboseIter = TRUE,
    allowParallel = TRUE
  )
  
  # Control for final model training
  tune_control <- trainControl(
    method = "cv",
    number = 4,
    classProbs = TRUE,
    summaryFunction = summary_fnc,
    verboseIter = FALSE,
    allowParallel = TRUE
  )
  if(model_selection=='RF'){
    rf_model_clinical <- train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),          
      method = "rf",
      trControl = ctrl,
      tuneLength =number_iterations,
      metric=optimization_metric
    )
    
    (final_grid<- expand.grid(
      mtry=rf_model_clinical$bestTune$mtry
    ))
    
    (model_Clinical <- caret::train(
      x = input_x_clinical,
      y = as.factor(input_y_clinical),
      trControl = tune_control,
      tuneGrid = final_grid,
      method = "rf",
      verbose = TRUE
    ))
    
    
    }else if(model_selection=='kNN'){
      
      
      model_Clinical <- tryCatch({
        final_grid <- expand.grid(
          k = seq(3, 15, by = 2)
        )
        
        knn_model_clinical <- caret::train(
          x = input_x_clinical,
          y = as.factor(input_y_clinical),
          method = "knn",
          trControl = ctrl,
          tuneGrid = final_grid,
          metric = optimization_metric
        )
        
        (final_grid<- expand.grid(
          k=knn_model_clinical$bestTune$k
        ))
        
        caret::train(
          x = input_x_clinical,
          y = as.factor(input_y_clinical),
          trControl = tune_control,
          tuneGrid = final_grid,
          method = "knn")
        
      }, error = function(e) {
        message("SVM training failed: ", conditionMessage(e))
        NULL
      })
      
      
    }else if(model_selection=='SVM'){
      
      
      model_Clinical <- tryCatch({
        svm_grid <- expand.grid(
          sigma = seq(0.001, 0.02, length.out = 5),
          C = seq(0.5, 5, length.out = 5)
        ) 
        
        svm_model_clinical <- caret::train(
          x = input_x_clinical,
          y = as.factor(input_y_clinical),
          method = "svmRadial",
          trControl = ctrl,
          tuneGrid = svm_grid,
          metric = optimization_metric
        )
        
        final_grid <- expand.grid(
          sigma = svm_model_clinical$bestTune$sigma,
          C = svm_model_clinical$bestTune$C
        )
        
        caret::train(
          x = input_x_clinical,
          y = as.factor(input_y_clinical),
          trControl = tune_control,
          tuneGrid = final_grid,
          method = "svmRadial",
          probs = TRUE
        )
        
      }, error = function(e) {
        message("SVM training failed: ", conditionMessage(e))
        NULL
      })
    
      
      
      
      
      
      
      
      
      
    }else if(model_selection=="LR"){
      log_model_clinical <- train(
        x = input_x_clinical,
        y = as.factor(input_y_clinical),
        method = "glm",
        family = binomial(),
        trControl = ctrl,
        metric = optimization_metric
      )
      
      model_Clinical <- caret::train(
        x = input_x_clinical,
        y = as.factor(input_y_clinical),
        trControl = tune_control,
        method = "glm",
        family = binomial()
      )
      
    }
  
  
  return(list(`ClinicalModel`=model_Clinical))
  
}




# Model training for CAD presence 
target<-"SYNTAX.BINARY"
classifiers<-c("RF","SVM","kNN","LR")
resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
Models_all<-list()




for(jj in 1:10){
  resamp_models<-list()
  for(resamp in resamplings){
  data<-read.csv(paste0('./Train_test_sets/Train_set_CADP_',resamp,'_',jj,'.csv'))
  data[,target]<-as.factor(data[,target])
  data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY")]
  
  classifiers_list<-list()
  for(classifier in classifiers){
    m<-Model_training(data,target,model_selection = classifier)
    model<-m$ClinicalModel
    classifiers_list[[classifier]]<-model
  }
  
  resamp_models[[resamp]]<-classifiers_list
  }
  Models_all[[paste0("Fold",jj)]]<-resamp_models
}
save(Models_all,file='./Models/Models_CADP.Rda')




# Model training for CAD severity 
target<-"SYNTAX.SEVERITY"
classifiers<-c("RF","SVM","kNN","LR")
resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
Models_all<-list()
for(jj in 1:10){
  resamp_models<-list()
  for(resamp in resamplings){
    data<-read.csv(paste0('./Train_test_sets/Train_set_CADS_',resamp,'_',jj,'.csv'))
    data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.BINARY","SYNTAX.HIGH.SEVERITY")]
    data[,target]<-as.factor(data[,target])
    classifiers_list<-list()
    for(classifier in classifiers){
      m<-Model_training(data,target,model_selection = classifier)
      model<-m$ClinicalModel
      classifiers_list[[classifier]]<-model
    }
    
    resamp_models[[resamp]]<-classifiers_list
  }
  Models_all[[paste0("Fold",jj)]]<-resamp_models
}
save(Models_all,file='./Models/Models_CADS.Rda')




# Model training for CAD high severity 
target<-"SYNTAX.HIGH.SEVERITY"
classifiers<-c("RF","SVM","kNN","LR")
resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
Models_all<-list()
for(jj in 1:10){
  resamp_models<-list()
  for(resamp in resamplings){
    data<-read.csv(paste0('./Train_test_sets/Train_set_CADHS_',resamp,'_',jj,'.csv'))
    data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.BINARY","SYNTAX.SEVERITY")]
    data[,target]<-as.factor(data[,target])
    classifiers_list<-list()
    for(classifier in classifiers){
      m<-Model_training(data,target,model_selection = classifier)
      model<-m$ClinicalModel
      classifiers_list[[classifier]]<-model
    }
    
    resamp_models[[resamp]]<-classifiers_list
  }
  Models_all[[paste0("Fold",jj)]]<-resamp_models
}
save(Models_all,file='./Models/Models_CADHS.Rda')



############## Complete #####################
# Model training for CAD presence 
target<-"SYNTAX.BINARY"
classifiers<-c("SVM")
resamplings<-c("rus")
Models_all<-list()
for(jj in 1){
  resamp_models<-list()
  for(resamp in resamplings){
    data<-read.csv(paste0('./Train_test_sets/Whole_set_CADP_',resamp,'_',jj,'.csv'))
    data[,target]<-as.factor(data[,target])
    data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY")]
    
    classifiers_list<-list()
    for(classifier in classifiers){
      m<-Model_training(data,target,model_selection = classifier)
      model<-m$ClinicalModel
      classifiers_list[[classifier]]<-model
    }
    
    resamp_models[[resamp]]<-classifiers_list
  }
  Models_all[[paste0("Fold",jj)]]<-resamp_models
}
save(Models_all,file='./Models/Models_CADP_whole.Rda')




# Model training for CAD severity 
target<-"SYNTAX.SEVERITY"
classifiers<-c("SVM")
resamplings<-c("rus")
Models_all<-list()
for(jj in 1){
  resamp_models<-list()
  for(resamp in resamplings){
    data<-read.csv(paste0('./Train_test_sets/Whole_set_CADS_',resamp,'_',jj,'.csv'))
    data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.BINARY","SYNTAX.HIGH.SEVERITY")]
    data[,target]<-as.factor(data[,target])
    classifiers_list<-list()
    for(classifier in classifiers){
      m<-Model_training(data,target,model_selection = classifier)
      model<-m$ClinicalModel
      classifiers_list[[classifier]]<-model
    }
    
    resamp_models[[resamp]]<-classifiers_list
  }
  Models_all[[paste0("Fold",jj)]]<-resamp_models
}
save(Models_all,file='./Models/Models_CADS_whole.Rda')




# Model training for CAD high severity 
target<-"SYNTAX.HIGH.SEVERITY"
classifiers<-c("SVM")
resamplings<-c("rus")
Models_all<-list()
for(jj in 1){
  resamp_models<-list()
  for(resamp in resamplings){
    data<-read.csv(paste0('./Train_test_sets/Whole_set_CADHS_',resamp,'_',jj,'.csv'))
    data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.BINARY","SYNTAX.SEVERITY")]
    data[,target]<-as.factor(data[,target])
    classifiers_list<-list()
    for(classifier in classifiers){
      m<-Model_training(data,target,model_selection = classifier)
      model<-m$ClinicalModel
      classifiers_list[[classifier]]<-model
    }
    
    resamp_models[[resamp]]<-classifiers_list
  }
  Models_all[[paste0("Fold",jj)]]<-resamp_models
}
save(Models_all,file='./Models/Models_CADHS_whole.Rda')


















###################################### For the repeats ############################
for(rep in 1:10){
  print(paste0("Repeat :",rep))
  setwd(paste0("./Repeats/Rep",rep))
  
  
  
  target<-"SYNTAX.BINARY"
  classifiers<-c("RF","SVM","kNN","LR")
  resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
  Models_all<-list()
  for(jj in 1:10){
    resamp_models<-list()
    for(resamp in resamplings){
      data<-read.csv(paste0('./Train_test_sets/Train_set_CADP_',resamp,'_',jj,'.csv'))
      data[,target]<-as.factor(data[,target])
      data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY")]
      
      classifiers_list<-list()
      for(classifier in classifiers){
        m<-Model_training(data,target,model_selection = classifier)
        model<-m$ClinicalModel
        classifiers_list[[classifier]]<-model
      }
      
      resamp_models[[resamp]]<-classifiers_list
    }
    Models_all[[paste0("Fold",jj)]]<-resamp_models
  }
  save(Models_all,file='./Models/Models_CADP.Rda')
  
  
  # Model training for CAD severity 
  target<-"SYNTAX.SEVERITY"
  classifiers<-c("RF","SVM","kNN","LR")
  resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
  Models_all<-list()
  for(jj in 1:10){
    resamp_models<-list()
    for(resamp in resamplings){
      data<-read.csv(paste0('./Train_test_sets/Train_set_CADS_',resamp,'_',jj,'.csv'))
      data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.BINARY","SYNTAX.HIGH.SEVERITY")]
      data[,target]<-as.factor(data[,target])
      classifiers_list<-list()
      for(classifier in classifiers){
        m<-Model_training(data,target,model_selection = classifier)
        model<-m$ClinicalModel
        classifiers_list[[classifier]]<-model
      }
      
      resamp_models[[resamp]]<-classifiers_list
    }
    Models_all[[paste0("Fold",jj)]]<-resamp_models
  }
  save(Models_all,file='./Models/Models_CADS.Rda')
  
  # Model training for CAD High Severity
  target<-"SYNTAX.HIGH.SEVERITY"
  classifiers<-c("RF","SVM","kNN","LR")
  resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
  Models_all<-list()
  for(jj in 1:10){
    resamp_models<-list()
    for(resamp in resamplings){
      data<-read.csv(paste0('./Train_test_sets/Train_set_CADHS_',resamp,'_',jj,'.csv'))
      data<-data[!colnames(data) %in% c("SYNTAX.SCORE","SYNTAX.BINARY","SYNTAX.SEVERITY")]
      data[,target]<-as.factor(data[,target])
      classifiers_list<-list()
      for(classifier in classifiers){
        m<-Model_training(data,target,model_selection = classifier)
        model<-m$ClinicalModel
        classifiers_list[[classifier]]<-model
      }
      
      resamp_models[[resamp]]<-classifiers_list
    }
    Models_all[[paste0("Fold",jj)]]<-resamp_models
  }
  save(Models_all,file='./Models/Models_CADHS.Rda')
  
  
}
