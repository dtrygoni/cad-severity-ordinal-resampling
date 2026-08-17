######### Code for train/set construction: 1. Train/Test splitting, 2. Univariate Filtering, 3. MissForest Imputation
#### 4. Min-Max normalization, 5. Boruta Feature selection
#### Two scenarios: One for CAD presence and One for CAD severity
set.seed(1)
#Library
library(caret)
library(missForestPredict)
library(Boruta)
library(SNPassoc)
#Functions
remove_NA_values<-function(data,thresh=0.5){
  cols<-c()
  n<-round(nrow(data)*thresh)
  for(col in colnames(data)){
    if(sum(is.na(data[,col]))>n){
      cols<-c(cols,col)
    }
  }
  data<-data[!colnames(data) %in% cols]
  return(data)
  
}
Univariate_LR_filtering<-function(train,test,target_col,cols_to_exclude,LR_threshold=0.10){
  clinical_predictors<-c()
  predictors_clinical<-setdiff(colnames(train),c(cols_to_exclude,target_col))
  for(predictor in predictors_clinical){
    form <- as.formula(paste(target_col, '~', paste0('`', predictor, '`')))
    res<-glm(form,family='binomial',data=train)
    train_temp<-train[!is.na(train[[predictor]]),]
    res_null<-glm(as.formula(paste0(target_col,' ~ 1')),family='binomial',data=train_temp)
    res_anova<-anova(res_null,res,test="Chisq")
    pval <- res_anova$`Pr(>Chi)`[2]
    if(!is.na(pval)){
    if(pval<0.10){
      clinical_predictors<-c(clinical_predictors,predictor)
    }
    }
  }
  
  
  
  
  #### Clinical sets
  train_clinical<-train[colnames(train) %in% c(clinical_predictors,target_col,cols_to_exclude)]
  #train_clinical<-train_clinical[complete.cases(train_clinical),]
  test_clinical<-test[colnames(test) %in% c(clinical_predictors,target_col,cols_to_exclude)]
  #test_clinical<-test_clinical[complete.cases(test_clinical),]
 
  
  return(list(`Train_clinical`=train_clinical,`Test_clinical`=test_clinical,`Clinical_predictors`=clinical_predictors
              ))
  
}
missForestImputation_fnc<-function(train,test){
  mm<-missForest(train)
  train_imputed<-mm$ximp
  test_imputed<-missForestPredict(mm,newdata=test)
  
  return(list(Train=train_imputed,Test=test_imputed))
}
Scaling_fnc<-function(train,test,cols_to_exclude){
  
  X_train<-train[!colnames(train) %in% cols_to_exclude]
  X_test<-test[!colnames(test) %in% cols_to_exclude]
  
  pp <- preProcess(X_train, method = "range")
  
  x_train_scaled <- predict(pp, X_train)
  x_test_scaled  <- predict(pp, X_test)
 
  
  x_train_scaled[,cols_to_exclude]<-train[,cols_to_exclude]
  x_test_scaled[,cols_to_exclude]<-test[,cols_to_exclude]
  
  return(list(Train=x_train_scaled,Test=x_test_scaled))
}
Boruta_fnc<-function(train,test,target_col,cols_to_exclude,option=2){
 x_train<-train[setdiff(colnames(train),c(cols_to_exclude,target_col))]
 
 
 r<-Boruta(x_train,train[,target_col],pValue=0.05)
 feats<-names(r$finalDecision[r$finalDecision %in% c("Confirmed")])
 if(option==1){
   x_train<-train[setdiff(colnames(x_train),colnames(train)[grepl("rs",colnames(train))])]
   x_snp<-train[colnames(train)[grepl("rs",colnames(train))]]
   r<-Boruta(x_train,train[,target_col],pValue=0.05)
   feats<-names(r$finalDecision[r$finalDecision %in% c("Confirmed")])
   
   r<-Boruta(x_snp,train[,target_col],pValue=0.05)
   feats_snp<-names(r$finalDecision[r$finalDecision %in% c("Confirmed")])
   train_bor<-train[c(feats,feats_snp,target_col,cols_to_exclude)]
   test_bor<-test[c(feats,feats_snp,target_col,cols_to_exclude)]
 }else if(option==2){
   train_bor<-train[c(feats,target_col,cols_to_exclude)]
   test_bor<-test[c(feats,target_col,cols_to_exclude)]
 }
 return(list(Train=train_bor,Test=test_bor))
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
snp_numeric_convert <-function(dataset,snp_columns,levels){
  dataset_no_snp<-dataset[! colnames(dataset) %in% snp_columns]
  dataset_snp<-dataset[colnames(dataset) %in% snp_columns]
  for(snp in snp_columns){
    levels_snp<-levels[[snp]]
    if(length(levels_snp)==2){
      level_0<-strsplit(levels_snp[1],split="-")[[1]]
      level_1<-strsplit(levels_snp[2],split="-")[[1]]
      dataset_snp[[snp]]<-as.factor(ifelse(dataset_snp[[snp]] %in% level_0,0,1))
    }else{
      level_0<-levels_snp[1]
      level_1<-levels_snp[2]
      dataset_snp[[snp]]<-as.numeric(ifelse(dataset_snp[[snp]] %in% level_0,0,ifelse(dataset_snp[[snp]] %in% level_1,1,2)))
    }
    
  }
  
  df<-cbind(dataset_no_snp,dataset_snp)
  return(df)
}
quality_controlSNP<-function(train,test,target,genetic_modelling='Optimal',HWE_threshold=0.05,MAF_threshold=5){
  cols<-which(grepl("rs",colnames(train)))
  clinical_cols<-setdiff(colnames(train),colnames(train)[cols])
  snps<-setupSNP(data=train,colSNPs=cols,sep="")
  snps_test<-setupSNP(data=test,colSNPs = cols,sep="")
  ########### 1. HWE 0.05 p value
  hwe<-tableHWE(snps,snps[,target])
  #Multiple testing correction: FDR
  p_vals <- hwe[,2]
  fdr_pvals <- p.adjust(p_vals, method = "fdr")
  threshold_p<-HWE_threshold
  snp_filtered<-rownames(hwe)[fdr_pvals < threshold_p | is.na(fdr_pvals)]
  snp_list<-setdiff(rownames(hwe),snp_filtered)
  snp_df<-snps[colnames(snps) %in% c(clinical_cols,snp_list)]
  snp_test_df<-snps_test[colnames(snps) %in% c(clinical_cols,snp_list)]
  ## 2. MAF < 1% (https://doi.org/10.1038/s41598-024-53310-x)
  summary_ss<-summary(snp_df,print=FALSE)
  MAF_thresh<-100-MAF_threshold
  snps_after_MAF<-rownames(summary_ss)[summary_ss[,2] < MAF_thresh]
  #no snp filtered with MAF filter (they were filtered before because those with p-values in HWE test, were eliminated)
  snp_df<-snps[colnames(snps) %in% c(clinical_cols,snps_after_MAF)]
  snp_test_df<-snps_test[colnames(snps) %in% c(clinical_cols,snps_after_MAF)]
  
  ## 3. LD R^2>0.5 (10.1002/gepi.20309)
  res<-SNPassoc::LD(snp_df[cols[1]:ncol(snp_df)],which=c("r"))
  r_matrix <- data.frame(res[["r"]]^2)  # Square to get r²
  snp_ld<-list()
  for(snp in colnames(r_matrix)){
    snp_matrix_temp<-r_matrix[snp]
    if(sum(snp_matrix_temp>0.5,na.rm=TRUE)>1){
      
      snp_names<-row.names(snp_matrix_temp)
      
      snps_in_ld_all<-row.names(snp_matrix_temp)[snp_matrix_temp>0.5]
      
      snps_in_ld<-snps_in_ld_all[!is.na(snps_in_ld_all)]
      snp_ld[[snp]]<-snps_in_ld
    } 
  }
  
  groups_snp <- list()
  i <- 0
  for (snp in names(snp_ld)) {
    group_temp <- c(snp, snp_ld[[snp]])
    matched_groups <- c()
    
    for (group in names(groups_snp)) {
      if (any(group_temp %in% groups_snp[[group]])) {
        matched_groups <- c(matched_groups, group)
      }
    }
    
    if (length(matched_groups) == 0) {
      i <- i + 1
      groups_snp[[as.character(i)]] <- group_temp
    } else {
      # Merge all matching groups into one
      merged_snps <- unique(c(group_temp, unlist(groups_snp[matched_groups])))
      
      # Remove old groups
      groups_snp[matched_groups] <- NULL
      
      # Add new merged group
      i <- i + 1
      groups_snp[[as.character(i)]] <- merged_snps
    }
  }
  clean_snps<-setdiff(colnames(snp_df)[cols[1]:ncol(snp_df)],unlist(groups_snp))
  snps_to_check<-groups_snp
  snp_without_ld<-snps[colnames(snps) %in% c(clinical_cols,clean_snps)]
  snp_ld_df<-snps[colnames(snps) %in% c(clinical_cols,unique(unlist(snps_to_check)))]
  ### selected SNPs with high LD from the groups
  selected_snps<-c()
  if(genetic_modelling=='Codominant'){
    ind<-2
  }else if(genetic_modelling=='Dominant'){
    ind<-6
  }else if(genetic_modelling=='Recessive'){
    ind<-9
  }else if(genetic_modelling=='Optimal'){
    ind<-NA
  }
  
  
  for(group in names(snps_to_check)){
    group_snp<-snps_to_check[[group]]
    snp_temp<-snp_ld_df[colnames(snp_ld_df) %in% c(target,group_snp)]
    aic_vals<-c()
    for(snp in group_snp){
      r<-SNPassoc::association(as.formula(paste0(target,' ~ ',snp)),data=snp_ld_df)
      if(!is.na(ind)){
        aic_vals<-c(aic_vals,r[ind,"AIC"])
      }else{
        val<-min(r[,9],na.rm=TRUE)
        aic_vals<-c(aic_vals,val)
      }
    }
    snp_select<-group_snp[which.min(aic_vals)]
    selected_snps<-c(selected_snps,snp_select)
  }
  all_snps<-c(selected_snps,clean_snps)
  ###### Train/test snps set: 693 train / 171 test  (features depend on the seed)
  train_snps<-snps[colnames(snps) %in% c(clinical_cols,all_snps)]
  test_snps<-snps_test[colnames(snps) %in% c(clinical_cols,all_snps)]
  #train_unconverted<-snp_data_train[colnames(snp_data_train) %in% c(clinical_cols,all_snps)]
  
  #data_filtered<-train_unconverted
  
  
  # ######## Modelling - levels to numerical
  cols<-which(grepl("rs",colnames(train_snps)))
  snp_cols<-colnames(train_snps)[cols[1]:ncol(train_snps)]
  levels_modelling<-list()
  snps<-c()
  for(snp in snp_cols){
    form<-as.formula(paste(target,' ~',snp))
    res<-association(form,data=train_snps)
    if(length(res[,1])==15){
      snps<-c(snps,snp)
      levels_co<-names(res[,1])[2:4]
      levels_do<-names(res[,1])[6:7]
      levels_re<-names(res[,1])[9:10]
      levels_over<-names(res[,1])[12:13]

      if(genetic_modelling=='Codominant'){
        levels_modelling[[snp]]<-levels_co
      }else if(genetic_modelling=='Dominant'){
        levels_modelling[[snp]]<-levels_do
      }else if(genetic_modelling=='Recessive'){
        levels_modelling[[Snp]]<-levels_re
      }else if(genetic_modelling=='Optimal'){
        index<-which.min(res[1:12,9])
        if(index==2){
          levels_opt<-levels_co
        }else if(index==6){
          levels_opt<-levels_do
        }else if(index==9){
          levels_opt<-levels_re
        }else if(index==12){
          levels_opt<-levels_over
        }
        levels_modelling[[snp]]<-levels_opt
      }


    }
  }

  train_converted<-snp_numeric_convert(train_snps,snp_cols,levels_modelling)
  test_converted<-snp_numeric_convert(test_snps,snp_cols,levels_modelling)
  
  
  return(list(Train=train_converted,Test=test_converted))
}

#Load the dataset
load("./Data/Whole_dataset.Rda")


clinical_data<-full_dataset
clinical_data<-clinical_data[!is.na(clinical_data$SYNTAX.SCORE),]
clinical_data$SYNTAX.BINARY<-as.factor(ifelse(clinical_data$SYNTAX.SCORE>0,1,0))
clinical_data$SYNTAX.SEVERITY<-as.factor(ifelse(clinical_data$SYNTAX.SCORE>23,1,0))
clinical_data$SYNTAX.HIGH.SEVERITY<-as.factor(ifelse(clinical_data$SYNTAX.SCORE>32,1,0))

#Remove columns with more then 50% missing data
data<-remove_NA_values(clinical_data)
data<-data[!colnames(data) %in% c('ID','SYNTAX.RISK','TOTAL.SS')]
#Test indices
ids<-createFolds(data$SYNTAX.BINARY,k=10,list=TRUE)
data$MPV<-as.numeric(data$MPV)
data$PDW<-as.numeric(data$PDW)


target<-"SYNTAX.BINARY"
#Cross validation sets
for(jj in 1:10){
  print(jj)
  train_set<-data[-ids[[jj]],]
  test_set<-data[ids[[jj]],]
  
  ### 0. SNP QC ####
  res<-quality_controlSNP(train_set,test_set,target)
  train_set<-res$Train
  test_set<-res$Test
  
  
  ### 1. Univariate filtering
  res<-Univariate_LR_filtering(train_set,test_set,target,c(target,setdiff(c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY","SYNTAX.BINARY"),target)))
  train_un<-res$Train_clinical
  test_un<-res$Test_clinical
  
  target1<-"SYNTAX.SEVERITY"
  res<-Univariate_LR_filtering(train_set,test_set,target1,c(target1,setdiff(c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY","SYNTAX.BINARY"),target1)))
  train_un1<-res$Train_clinical
  test_un1<-res$Test_clinical

  target2<-"SYNTAX.HIGH.SEVERITY"
  res<-Univariate_LR_filtering(train_set,test_set,target2,c(target2,setdiff(c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY","SYNTAX.BINARY"),target2)))
  train_un2<-res$Train_clinical
  test_un2<-res$Test_clinical

  train_un<-Test_set_preparation(list(train_un,train_un1,train_un2))
  test_un<-Test_set_preparation(list(test_un,test_un1,test_un2))

  
  ### 2. MissForest Predict (fit to train set and then transform test)
  res<-missForestImputation_fnc(train_un,test_un)
  train_imp<-res$Train
  test_imp<-res$Test
  ### 3. Min-Max Scaling
  res<-Scaling_fnc(train_imp,test_imp,c(target,setdiff(c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY","SYNTAX.BINARY"),target)))
  train_scaled<-res$Train
  test_scaled<-res$Test
  
  ### 4. Boruta FS
  res<-Boruta_fnc(train_scaled,test_scaled,target,setdiff(c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY","SYNTAX.BINARY"),target),option=2)
  train_set<-res$Train
  test_set<-res$Test
  #train_set<-train_set[!colnames(train_set) %in% paste0(target,".1")]
  #test_set<-test_set[!colnames(test_set) %in% paste0(target,".1")]
  
  write.csv(train_set,row.names=FALSE,file=paste0("C:/Users/30697/Dropbox/PhD/Working/Resampling/Code/Train_test_sets/Train_set_CADP_orig_",jj,".csv"))
  write.csv(test_set,row.names=FALSE,file=paste0("C:/Users/30697/Dropbox/PhD/Working/Resampling/Code/Train_test_sets/Test_set_CADP_orig_",jj,".csv"))
  
  
  target1<-"SYNTAX.SEVERITY"
  res1<-Boruta_fnc(train_scaled,test_scaled,target1,setdiff(c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY","SYNTAX.BINARY"),target1),option=2)
  train_set1<-res1$Train
  test_set1<-res1$Test
  #train_set1<-train_set1[!colnames(train_set1) %in% paste0(target1,".1")]
  #test_set1<-test_set1[!colnames(test_set1) %in% paste0(target1,".1")]

  write.csv(train_set1,row.names=FALSE,file=paste0("./Train_test_sets/Train_set_CADS_orig_",jj,".csv"))
  write.csv(test_set1,row.names=FALSE,file=paste0("./Train_test_sets/Test_set_CADS_orig_",jj,".csv"))

  res2<-Boruta_fnc(train_scaled,test_scaled,target2,setdiff(c("SYNTAX.SCORE","SYNTAX.SEVERITY","SYNTAX.HIGH.SEVERITY","SYNTAX.BINARY"),target2),option=2)
  train_set2<-res2$Train
  test_set2<-res2$Test
  #train_set2<-train_set2[!colnames(train_set2) %in% paste0(target2,".1")]
  #test_set2<-test_set2[!colnames(test_set2) %in% paste0(target2,".1")]

  write.csv(train_set2,row.names=FALSE,file=paste0("./Train_test_sets/Train_set_CADHS_orig_",jj,".csv"))
  write.csv(test_set2,row.names=FALSE,file=paste0("./Train_test_sets/Test_set_CADHS_orig_",jj,".csv"))

  
}


