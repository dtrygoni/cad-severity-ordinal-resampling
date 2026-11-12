####################### Point estimates: evaluation of resampling in Multiclassification: Three classes, Four classes###############
#Libraries
library(caret)
library(pROC)
library(irr)
library(Iso)
library(tidyr)
library(dplyr)
library(kernlab)
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
Macro_MAE_fnc<-function(target,predicted){
  n_groups<-length(unique(target))
  maes<-c()
  for(val in unique(target)){
    mae_val<-mean(abs(target[target==val] - predicted[target==val]))
    maes<-c(maes,mae_val)
  }
  macro_mae<-round(sum(maes)/length(maes),3)
  return(macro_mae)
}
MAE_per_class<-function(target,predicted){
  n_groups<-length(unique(target))
  maes<-list()
  for(val in unique(target)){
    mae_val<-mean(abs(target[target==val] - predicted[target==val]))
    maes[[as.character(val)]]<-round(mae_val,3)
  }
  return(maes)
}
Multiclass_evalutation<-function(test,models_list,n_class=4){
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
    kappa<-round(kappa2(data.frame(preds, target), "equal")$value,3)
    
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
    
    
    macro_mae<-Macro_MAE_fnc(target,preds)
    df_wide$Macro.MAE<-macro_mae
    
    ############### MAE per class ############
    maes_list<-MAE_per_class(target,preds)
    for(class in names(maes_list)){
      nm<-paste0("MAE.",class)
      val<-maes_list[[class]]
      df_wide[[nm]]<-val
    }
    
    return(list(df=data.frame(df_wide),Preds=preds,Target=target))
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
df_whole<-NULL
classifiers<-c("RF","SVM","kNN","LR")
resamplings<-c("rus","ros","smote","smote_enn","orig","enn")



for(jj in 1:10){
  test_jj<-f_test[[jj]]
  for(resamp in resamplings){
    for(classifier in classifiers){
      model_p<-CADP_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
      model_s<-CADS_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
      model_hs<-CADHS_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
      
      
      res1<-Multiclass_evalutation(test_jj,list(model_p,model_s,model_hs),4)
      res<-res1$df
      res$Sampling<-resamp
      res$Classifier<-classifier
      res$Fold<-jj
      res_all<-rbind(res_all,res)
      
      
      target_jj<-res1$Target
      pred_jj<-res1$Preds
      fold_jj<-rep(jj,length(pred_jj))
      df_temp<-data.frame(Fold=fold_jj,Algorithm=rep(classifier,length(pred_jj)),Sampling=rep(resamp,length(pred_jj)),
                          Predicted=pred_jj,Target=target_jj)
      df_whole<-rbind(df_whole,df_temp)
    }
  }
  
}

Four.class.multiclass<-res_all
save(Four.class.multiclass,file='./Results/Multiclass_four_classes_f.Rda')



############## Evaluation of predictions #############
df_rf_rus<-df_whole[df_whole$Algorithm=="SVM" & df_whole$Sampling=="enn",]
df_temp<-df_rf_rus
maes<-c()
for(fold in 1:10){
  df_temp_fold<-df_temp[df_temp$Fold==fold,]
  macro.mae<-Macro_MAE_fnc(df_temp_fold$Target,df_temp_fold$Predicted)
  maes<-c(maes,macro.mae)
}
Macro_MAE_fnc(df_temp$Target,df_temp$Predicted)
mae<-c()
for(jj in 1:10){
  df_jj<-df_temp[df_temp$Fold==jj,]
  mae_jj<-mean(abs(df_jj$Target-df_jj$Predicted))
  mae<-c(mae,mae_jj)
}
confusionMatrix(factor(df_rf_rus$Predicted),factor(df_rf_rus$Target))






########## Evaluation of correct calculation ########
df_rf_rus<-df_whole[df_whole$Algorithm=="SVM" & df_whole$Sampling=="orig",]
target<-df_rf_rus$Target
predicted<-df_rf_rus$Predicted

n_groups<-length(unique(target))
maes<-c()
for(val in unique(target)){
  mae_val<-mean(abs(target[target==val] - predicted[target==val]))
  maes<-c(maes,mae_val)
}

df_rf_rus<-df_whole[df_whole$Algorithm=="SVM" & df_whole$Sampling=="rus",]
target<-df_rf_rus$Target
predicted<-df_rf_rus$Predicted

n_groups<-length(unique(target))
maes<-c()
for(val in unique(target)){
  mae_val<-mean(abs(target[target==val] - predicted[target==val]))
  maes<-c(maes,mae_val)
}
macro_mae<-round(sum(maes)/length(maes),3)






df_all<-NULL
for(fold in 1:10){
  for(classifier in classifiers){
    for(resamp in resamplings){
      df_temp<-df_whole[df_whole$Sampling==resamp & df_whole$Algorithm==classifier & df_whole$Fold==fold,]
      target<-df_temp$Target
      predicted<-df_temp$Predicted
      maes<-c()
      for(val in c(0,1,2,3)){
        mae_val<-mean(abs(target[target==val] - predicted[target==val]))
        maes<-c(maes,mae_val)
      }
      MAE.class0<-maes[1]
      MAE.class1<-maes[2]
      MAE.class2<-maes[3]
      MAE.class3<-maes[4]
      
      dfjj<-data.frame(Fold=fold,Algorithm=classifier,Sampling=resamp,MAE.Class0=round(MAE.class0,3),
                       MAE.Class1=round(MAE.class1,3),MAE.Class2=round(MAE.class2,3),MAE.Class3=round(MAE.class3,3),
                       Macro.MAE=round(mean(maes),3))
      
      df_all<-rbind(df_all,dfjj)
    }
  }
}
save(df_all,file='./Results/MAE_per_class_per_fold.Rda')





############################### For repeats #############################
for(rep in 1:10){
  print(rep)
  #1. CAD presence
  load(paste0('./Repeats/Rep',rep,'/Models/Models_CADP.Rda'))
  CADP_models<-Models_all
  #2. CAD severity
  load(paste0('./Repeats/Rep',rep,'/Models/Models_CADS.Rda'))
  CADS_models<-Models_all
  #3. CAD high severity
  load(paste0('./Repeats/Rep',rep,'/Models/Models_CADHS.Rda'))
  CADHS_models<-Models_all
  
  
  #Extract the test sets for all
  test_set_CADP<-Extract_test_sets(paste0('./Resampling/Code/Repeats/Rep',rep,'/Train_test_sets/Test_set_CADP_orig_'))
  test_set_CADS<-Extract_test_sets(paste0('./Repeats/Rep',rep,'/Train_test_sets/Test_set_CADS_orig_'))
  test_set_CADHS<-Extract_test_sets(paste0('./Repeats/Rep',rep,'/Train_test_sets/Test_set_CADHS_orig_'))
  
  
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
  df_whole<-NULL
  classifiers<-c("RF","SVM","kNN","LR")
  resamplings<-c("rus","ros","smote","smote_enn","orig","enn")
  
  
  
  for(jj in 1:10){
    test_jj<-f_test[[jj]]
    for(resamp in resamplings){
      for(classifier in classifiers){
        model_p<-CADP_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
        model_s<-CADS_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
        model_hs<-CADHS_models[[paste0("Fold",jj)]][[resamp]][[classifier]]
        
        
        res1<-Multiclass_evalutation(test_jj,list(model_p,model_s,model_hs),4)
        res<-res1$df
        res$Sampling<-resamp
        res$Classifier<-classifier
        res$Fold<-jj
        res_all<-rbind(res_all,res)
        
        
        target_jj<-res1$Target
        pred_jj<-res1$Preds
        fold_jj<-rep(jj,length(pred_jj))
        df_temp<-data.frame(Fold=fold_jj,Algorithm=rep(classifier,length(pred_jj)),Sampling=rep(resamp,length(pred_jj)),
                            Predicted=pred_jj,Target=target_jj)
        df_whole<-rbind(df_whole,df_temp)
      }
    }
    
  }
  
  Four.class.multiclass<-res_all
  save(Four.class.multiclass,file=paste0('./Repeats/Rep',rep,'/Results/Multiclass_four_classes_f.Rda'))
  
}
