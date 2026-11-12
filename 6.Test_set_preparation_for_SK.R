################ Prepare the binary classifiers for ScottKnott ##########
load('./Results/Binary_CAD_presence.Rda')
library(dplyr)
Correct_dataset_format<-function(dataset){
  wanted<-c("Classifier","Resampling","Fold")
  data<-dataset[,c(wanted,setdiff(colnames(dataset),wanted))]
  data2<-data %>%
    rename(Algorithm=Classifier,
           Technique=Resampling)
  
  data2$Algorithm<-factor(data2$Algorithm)
  data2$Technique<-factor(data2$Technique)
  data2$Fold<-factor(data2$Fold)
  levels(data2$Technique)<-c("ENN","No-sampling","ROS","RUS","SMOTE","SMOTE ENN")
  return(data2)
}


dataset<-Correct_dataset_format(res_presence)
save(dataset,file='./Results/SK_ready_CADP_res.Rda')


load('./Results/Binary_CAD_severity.Rda')
dataset<-Correct_dataset_format(res_severity)
save(dataset,file='./Results/SK_ready_CADS_res.Rda')

load('./Results/Binary_CAD_high_severity.Rda')
dataset<-Correct_dataset_format(res_high_severity)
save(dataset,file='./Results/SK_ready_CADHS_res.Rda')






load('./Results/Multiclass_four_classes_f.Rda')
data<-Four.class.multiclass%>%
  rename(Resampling=Sampling)
dataset<-Correct_dataset_format(data)
save(dataset,file='./Results/SK_ready_fourClass_res_f.Rda')



############################### For repeats ################################
for(rep in 1:10){
  print(rep)
  load(paste0('./Repeats/Rep',rep,'/Results/Binary_CAD_presence.Rda'))
  dataset<-Correct_dataset_format(res_presence)
  save(dataset,file=paste0('./Repeats/Rep',rep,'/Results/SK_ready_CADP_res.Rda'))
  
  
  load(paste0('./Repeats/Rep',rep,'/Results/Binary_CAD_severity.Rda'))
  dataset<-Correct_dataset_format(res_severity)
  save(dataset,file=paste0('./Repeats/Rep',rep,'/Results/SK_ready_CADS_res.Rda'))
  
  load(paste0('./Repeats/Rep',rep,'/Results/Binary_CAD_high_severity.Rda'))
  dataset<-Correct_dataset_format(res_high_severity)
  save(dataset,file=paste0('./Repeats/Rep',rep,'/Results/SK_ready_CADHS_res.Rda'))
  
  
  
  
  
  
  load(paste0('./Repeats/Rep',rep,'/Results/Multiclass_four_classes_f.Rda'))
  data<-Four.class.multiclass%>%
    rename(Resampling=Sampling)
  dataset<-Correct_dataset_format(data)
  save(dataset,file=paste0('./Repeats/Rep',rep,'/Results/SK_ready_fourClass_res_f.Rda'))
  
  
}







