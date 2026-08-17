####################### Point estimates: evaluation of resampling in Multiclassification: Three classes, Four classes###############

#Libraries
library(caret)
library(pROC)
library(irr)
library(Iso)
library(tidyr)
library(dplyr)
#Functions
Probabilities_correction_with_Iso<-function(predictions_list,n_class=3){
  if(n_class==3){
    predictions.0<-predictions_list[[1]]
    predictions.1<-predictions_list[[2]]
    n<-nrow(predictions.0)
    q0<-c()
    q1<-c()
    for(i in 1:n){
      qs_new<-pava(c(predictions.0[i,2],predictions.1[i,2]),decreasing=TRUE)
      q0<-c(q0,qs_new[1])
      q1<-c(q1,qs_new[2])
    }
    
    p0<-1-q0
    p1<-q0-q1
    p2<-q1
    
    
    probs<-data.frame(Class.0=p0,Class.1=p1,Class.2=p2)
    pred_class <- colnames(probs)[max.col(probs, ties.method = "first")]
    
    preds<-ifelse(pred_class=="Class.0",0,ifelse(pred_class=="Class.1",1,2))

    return(preds)
  }else if(n_class==4){
    predictions.0<-predictions_list[[1]]
    predictions.1<-predictions_list[[2]]
    predictions.2<-predictions_list[[3]]
    n<-nrow(predictions.0)
    q0<-c()
    q1<-c()
    q2<-c()
    for(i in 1:n){
      qs_new<-pava(c(predictions.0[i,2],predictions.1[i,2],predictions.2[i,2]),decreasing=TRUE)
      q0<-c(q0,qs_new[1])
      q1<-c(q1,qs_new[2])
      q2<-c(q2,qs_new[3])
    }
    
    p0<-1-q0
    p1<-q0-q1
    p2<-q1-q2
    p3<-q2
    
    probs<-data.frame(Class.0=p0,Class.1=p1,Class.2=p2,Class.3=p3)
    pred_class <- colnames(probs)[max.col(probs, ties.method = "first")]
    
    preds<-ifelse(pred_class=="Class.0",0,ifelse(pred_class=="Class.1",1,
                                                 ifelse(pred_class=="Class.2",2,3)))
    
    return(preds)
  }
  
}
Multiclass_evalutation<-function(test,models_list,n_class=3){
  if(n_class==3){
    model.0<-models_list[[1]]
    test.0<-test[predictors(model.0)]
    predictions.0<-predict(model.0,newdata=test.0,type="prob")
    model.1<-models_list[[2]]
    test.1<-test[predictors(model.1)]
    predictions.1<-predict(model.1,newdata=test.1,type="prob")
    
    preds<-Probabilities_correction_with_Iso(list(predictions.0,predictions.1),3)
    
    
    
    target<-Target_class_construction(test,n_class)
    
    cm <- caret::confusionMatrix(
      factor(preds),
      factor(target)
    )

    byClass <- cm$byClass
    macro_recall <- round(mean(byClass[, "Sensitivity"], na.rm = TRUE),3)
    macro_specificity <- round(mean(byClass[, "Specificity"], na.rm = TRUE),3)
    macro_ballAcc<-round(mean(byClass[, "Balanced Accuracy"], na.rm = TRUE),3)
    
    mae <- round(mean(abs(preds - target)),3)
    mse <- round(mean((preds - target)^2),3)
    kappa<-round(kappa2(data.frame(preds, target), "squared")$value,3)
    
    df <- as.data.frame(round(cm$byClass,3))
    
    df$class <- rownames(df)
    
    library(tidyr)
    library(dplyr)
    
    df_long <- df %>%
      tidyr::pivot_longer(-class, names_to = "metric", values_to = "value")
    
    df_wide <- df_long %>%
      tidyr::unite(metric_class, metric, class) %>%
      tidyr::pivot_wider(names_from = metric_class, values_from = value)
    
    
    df_wide$Macro.Recall<-macro_recall
    df_wide$Macro.Specificity<-macro_specificity
    df_wide$Macro.BalancedAccuracy<-macro_ballAcc
    df_wide$MAE<-mae
    df_wide$MSE<-mse
    df_wide$Kappa2<-kappa
    
    return(data.frame(df_wide))
  }else if(n_class==4){
    model.0<-models_list[[1]]
    test.0<-test[predictors(model.0)]
    predictions.0<-predict(model.0,newdata=test.0,type="prob")
    model.1<-models_list[[2]]
    test.1<-test[predictors(model.1)]
    predictions.1<-predict(model.1,newdata=test.1,type="prob")
    model.2<-models_list[[3]]
    test.2<-test[predictors(model.2)]
    predictions.2<-predict(model.2,newdata=test.2,type="prob")
    
    
    preds<-Probabilities_correction_with_Iso(list(predictions.0,predictions.1,predictions.2),n_class)
    
    
    
    target<-Target_class_construction(test,n_class)
    
    cm <- caret::confusionMatrix(
      factor(preds),
      factor(target)
    )
    
    byClass <- cm$byClass
    macro_recall <- round(mean(byClass[, "Sensitivity"], na.rm = TRUE),3)
    macro_specificity <- round(mean(byClass[, "Specificity"], na.rm = TRUE),3)
    macro_ballAcc<-round(mean(byClass[, "Balanced Accuracy"], na.rm = TRUE),3)
    
    mae <- round(mean(abs(preds - target)),3)
    mse <- round(mean((preds - target)^2),3)
    kappa<-round(kappa2(data.frame(preds, target), "squared")$value,3)
    
    df <- as.data.frame(round(cm$byClass,3))
    
    df$class <- rownames(df)
    

    
    df_long <- df %>%
      tidyr::pivot_longer(-class, names_to = "metric", values_to = "value")
    
    df_wide <- df_long %>%
      tidyr::unite(metric_class, metric, class) %>%
      tidyr::pivot_wider(names_from = metric_class, values_from = value)
    
    
    df_wide$Macro.Recall<-macro_recall
    df_wide$Macro.Specificity<-macro_specificity
    df_wide$Macro.BalancedAccuracy<-macro_ballAcc
    df_wide$MAE<-mae
    df_wide$MSE<-mse
    df_wide$Kappa2<-kappa
    
    return(data.frame(df_wide))
  }
}
Target_class_construction<-function(test,nclass=3){
  if(nclass==3){
    wanted0<-'SYNTAX.BINARY'
    wanted1<-'SYNTAX.SEVERITY'
    test_target<-test[,c(wanted0,wanted1)]
    target_col<-ifelse(test_target$SYNTAX.SEVERITY==1,2,ifelse(test_target$SYNTAX.BINARY==1,1,0))
    return(target_col)
  }else{
    wanted0<-'SYNTAX.BINARY'
    wanted1<-'SYNTAX.SEVERITY'
    wanted2<-'SYNTAX.HIGH.SEVERITY'
    test_target<-test[,c(wanted0,wanted1,wanted2)]
    target_col<-ifelse(test_target$SYNTAX.HIGH.SEVERITY==1,3,ifelse(test_target$SYNTAX.SEVERITY==1,2,
                                                                    ifelse(test_target$SYNTAX.BINARY==1,1,0)))
  }
}
Test_set_preparation<-function(list_of_test_sets){
  test_first<-list_of_test_sets[[1]]
  for(i in 2:length(list_of_test_sets)){
    test_i<-list_of_test_sets[[i]]
    diff_cols<-setdiff(colnames(test_i),test_first)
    for(col in diff_cols){
      test_first[,col]<-test_i[,col]
    }
  }
  return(test_first)
}
Extract_test_sets<-function(test_path){
  test_all<-list()
  for(jj in 1:10){
    test<-read.csv(paste0(test_path,jj,'.csv'))
    test_all[[jj]]<-test
  }
  return(test_all)
}
#1. CAD presence
load('./Models/Models_CADP.Rda')
CADP_models<-Models_all
#2. CAD severity
load('./Models/Models_CADS.Rda')
CADS_models<-Models_all
#3. CAD high severity
load('./Models/Models_CADHS.Rda')
CADHS_models<-Models_all


#Extract the test sets for all
test_set_CADP<-Extract_test_sets('./Train_test_sets/Test_set_CADP_orig_')
test_set_CADS<-Extract_test_sets('./Train_test_sets/Test_set_CADS_orig_')
test_set_CADHS<-Extract_test_sets('./Train_test_sets/Test_set_CADHS_orig_')


f_test<-list()
for(jj in 1:10){
  test_p<-test_set_CADP[[jj]]
  test_s<-test_set_CADS[[jj]]
  test_hs<-test_set_CADHS[[jj]]
  test_whole<-Test_set_preparation(list(test_p,test_s,test_hs))
  
  f_test[[jj]]<-test_whole
}





####################### Four-class classification ##################
res_all<-NULL
classifiers<-c("RF","SVM","kNN","LR")
resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
for(jj in 1:10){
  test_jj<-f_test[[jj]]
  for(resamp in resamplings){
    for(classifier in classifiers){
      model_p<-CADP_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
      model_s<-CADS_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
      model_hs<-CADHS_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
      res<-Multiclass_evalutation(test_jj,list(model_p,model_s,model_hs),4)
      res$Sampling<-resamp
      res$Classifier<-classifier
      res$Fold<-jj
      res_all<-rbind(res_all,res)
    }
  }
  
}
Four.class.multiclass<-res_all
save(Four.class.multiclass,file='./Results/Multiclass_four_classes.Rda')



















