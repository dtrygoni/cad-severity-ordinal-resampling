########### Scott-Knot results############
### Libraries ####
library(reshape2)
library(psych)
library(ggplot2)
#library(ggExtra)#
library(ggthemes)
library(lmerTest)
library(lsmeans)
library(emmeans)
library("multcompView")
library(dplyr)
library("ScottKnott")
library(rlang)

############################

### Function #####
StatisticalComparison <- function (Dataset, Metric = "Balanced.Accuracy"){
  ### Exploratory Analytics for AUC
  ExploratoryTable <- describeBy(Dataset[,paste0(Metric)], list(Dataset$Algorithm,
                                                                Dataset$Technique),
                                 mat=T)
  
  ExploratoryTable <- ExploratoryTable[,c("group1","group2","n","mean","sd","median","min","max")]
  
  ExploratoryTable <- ExploratoryTable[order(ExploratoryTable$group1,decreasing = F),]
  
  
  ExploratoryTable[,c("mean","sd","median","min","max")]<- lapply(ExploratoryTable[,c("mean","sd","median","min","max")],
                                                                  round, 4)
  
  Boxplots <- ggplot(data = Dataset ,aes(y = Technique , x = Dataset[,Metric])) +
    geom_boxplot()+
    facet_wrap(.~Algorithm)+
    #scale_fill_viridis_d( option = "D", direction = -1)+  
    theme_bw()+
    #removeGrid(x=TRUE, y=TRUE)+
    guides(fill = guide_legend(override.aes = list(alpha = 1,color="black")))+
    theme(legend.position="none",
          strip.text.x = element_text(size = 16, hjust = 0.5, vjust = 0.5, face = 'bold'),
          strip.text.y = element_text(size = 16, hjust = 0.5, vjust = 0.5, face = 'bold'),
          
          legend.title=element_blank(),
          legend.text = element_text(size=30),
          #axis.text.x = element_blank(),
          axis.text.x = element_text(face="bold", color="black", size=15),
          axis.text.y = element_text(face="bold", color="black", size=15),
          #axis.title.x=element_blank(),
          axis.title.x = element_text(face="bold", color="black", size=20),
          axis.title.y = element_text(face="bold", color="black", size=20),
          plot.title = element_text(face="bold", color = "black", size=20))+
    labs(x = "AUC-ROC", y = "Technique")
  
  
  
  # LME model with interaction
  # Formulas
  form_interaction <- as.formula(paste(Metric, "~ Algorithm * Technique + (1 | Fold)"))
  form_main <- as.formula(paste(Metric, "~ Algorithm + Technique + (1 | Fold)"))
  
  # Model comparison
  ModelInteraction <- lmerTest::lmer(form_interaction, data = Dataset, REML = FALSE)
  ModelMain <- lmerTest::lmer(form_main, data = Dataset, REML = FALSE)
  ModelComparison <- anova(ModelInteraction, ModelMain)
  
  p_interaction <- ModelComparison$`Pr(>Chisq)`[2]
  
  
  
  
  
  if (p_interaction>0.05){
    print("Insignificant interaction term")
    
    ModelFinal <- lmer(Dataset[,Metric] ~ Algorithm+Technique+
                         # Preprocess*Dataset+
                         (1 | Fold), 
                       data = Dataset[c("Algorithm","Technique","Fold",Metric)],
                       REML=TRUE)
    
    AnovaTable <- anova(ModelFinal,ddf="Kenward-Roger")
    
    AnovaTable <- AnovaTable [,c(3:6)]
    
    plot(ModelFinal)
    qqnorm(resid(ModelFinal))
    
    # #### Scott Knott
    formula.SK <- paste(Metric,"~","Algorithm+Technique+(1|Fold)")
    
    anv <- lmer(as.formula(formula.SK), data=Dataset[c("Algorithm","Technique","Fold",Metric)])
    
    ### SK results for Algorithm #### 
    SKResults.Algorithm <- ScottKnott::SK(x=anv,
                                          y=Metric,
                                          which='Algorithm')
    
    
    
    SKResultsStatistics.Algorithm <- data.frame(SKResults.Algorithm$info)
    
    Algorithm <- SKResultsStatistics.Algorithm [,1]
    Measure <- SKResultsStatistics.Algorithm [,2]
    Lower <- SKResultsStatistics.Algorithm[,4]
    Upper <- SKResultsStatistics.Algorithm[,5]
    
    ClusterMatrix <- SKResults.Algorithm$out$Result
    ClusterMatrix <- ClusterMatrix[-1]
    ClusterMatrix <- cbind(ClusterMatrix, Group = apply(ClusterMatrix, 1, 
                                                        function(x) paste0(na.omit(x), collapse = "")))
    
    Groups.SK.Algorithm <- data.frame(Algorithm,Measure,Lower,Upper,ClusterMatrix$Group)
    
    colnames(Groups.SK.Algorithm)[5] <- c("Cluster")
    
    ### SK results for Technique #### 
    SKResults.Technique <- ScottKnott::SK(x=anv,
                                          y=Metric,
                                          which='Technique')
    
    SKResultsStatistics.Technique <- data.frame(SKResults.Technique$info)
    
    Technique <- SKResultsStatistics.Technique [,1]
    Measure <- SKResultsStatistics.Technique [,2]
    Lower <- SKResultsStatistics.Technique[,4]
    Upper <- SKResultsStatistics.Technique[,5]
    
    ClusterMatrix <- SKResults.Technique$out$Result
    ClusterMatrix <- ClusterMatrix[-1]
    ClusterMatrix <- cbind(ClusterMatrix, Group = apply(ClusterMatrix, 1, 
                                                        function(x) paste0(na.omit(x), collapse = "")))
    
    Groups.SK.Technique <- data.frame(Technique,Measure,Lower,Upper,ClusterMatrix$Group)
    
    colnames(Groups.SK.Technique)[5] <- c("Cluster")
    
    list(
      ExploratoryTable = ExploratoryTable,
      Boxplots = Boxplots, 
      Groups.SK.Algorithm = Groups.SK.Algorithm,
      Groups.SK.Technique = Groups.SK.Technique)
    
    
  } else {
    print("Significant interaction term was identified")
    
    ModelFinal <- lmer(Dataset[,Metric] ~ Algorithm*Technique+
                         # Preprocess*Dataset+
                         (1 | Fold), 
                       data = Dataset[c("Algorithm","Technique","Fold",Metric)],
                       REML=TRUE)
    
    AnovaTable <- anova(ModelFinal,ddf="Kenward-Roger")
    
    AnovaTable <- AnovaTable [,c(3:6)]
    
    plot(ModelFinal)
    qqnorm(resid(ModelFinal))
    
    # #### Scott Knott
    formula.SK <- paste(Metric,"~","Algorithm*Technique+(1 | Fold)")
    
    anv <- lme4::lmer(as.formula(formula.SK), data=Dataset[c("Algorithm","Technique","Fold",Metric)])
    ### SK results for Algorithm #### 
    Groups.SK.Algorithm <- data.frame()
    
    
    ############################################ Error is here #################################
    # for(i in 1:nlevels(Dataset$Algorithm)){
    #   SKResults.Algorithm <- ScottKnott::SK(x=anv,
    #                                         y=Metric,
    #                                         which='Algorithm:Technique',
    #                                         fl1=i)
    # 
    # 
    #   SKResultsStatistics.Algorithm <- data.frame(SKResults.Algorithm$info)
    #   
    #   Algorithm <- SKResultsStatistics.Algorithm [,1]
    #   Technique <- SKResultsStatistics.Algorithm [,2]
    #   Measure <- SKResultsStatistics.Algorithm [,3]
    #   Lower <- SKResultsStatistics.Algorithm[,6]
    #   Upper <- SKResultsStatistics.Algorithm[,7]
    #   
    #   ClusterMatrix <- SKResults.Algorithm$out$Result
    #   ClusterMatrix <- ClusterMatrix[-1]
    #   ClusterMatrix <- cbind(ClusterMatrix, Group = apply(ClusterMatrix, 1, 
    #                                                       function(x) paste0(na.omit(x), collapse = "")))
    #   
    #   Groups.SK.Algorithm0 <- data.frame(Algorithm,Technique,Measure,Lower,Upper,ClusterMatrix$Group)
    #   
    #   colnames(Groups.SK.Algorithm0)[6] <- c("Cluster")
    #   Groups.SK.Algorithm <- rbind(Groups.SK.Algorithm, Groups.SK.Algorithm0)
    # }
    #############################################
    for(i in 1:nlevels(Dataset$Algorithm)){
      lv<-levels(Dataset$Algorithm)[i]
      Dataset_i<-Dataset[Dataset$Algorithm==lv,]
      formula.SK <- paste(Metric,"~","Technique+(1 | Fold)")
      
      anv <- lmer(as.formula(formula.SK), data=Dataset_i[c("Algorithm","Technique","Fold",Metric)])
      
      SKResults.Algorithm <- ScottKnott::SK(x=anv,
                                            y=Metric,
                                            which='Technique')
      
      
      SKResultsStatistics.Algorithm <- data.frame(SKResults.Algorithm$info)
      
      Algorithm <- lv
      Technique <- SKResultsStatistics.Algorithm [,1]
      Measure <- SKResultsStatistics.Algorithm [,2]
      Lower <- SKResultsStatistics.Algorithm[,4]
      Upper <- SKResultsStatistics.Algorithm[,5]
      
      ClusterMatrix <- SKResults.Algorithm$out$Result
      ClusterMatrix <- ClusterMatrix[-1]
      ClusterMatrix <- cbind(ClusterMatrix, Group = apply(ClusterMatrix, 1,
                                                          function(x) paste0(na.omit(x), collapse = "")))
      
      Groups.SK.Algorithm0 <- data.frame(Algorithm,Technique,Measure,Lower,Upper,ClusterMatrix$Group)
      
      colnames(Groups.SK.Algorithm0)[6] <- c("Cluster")
      Groups.SK.Algorithm <- rbind(Groups.SK.Algorithm, Groups.SK.Algorithm0)
    }
    
    
    
    ########################################## Rest of the code ###############################
    
    
    Dataset$Interaction <- interaction(Dataset$Algorithm, 
                                       Dataset$Technique, 
                                       sep = "/",
                                       drop=T)
    
    ModelFinalInteraction <- lmer(Dataset[,Metric] ~ Interaction+
                                    # Preprocess*Dataset+
                                    (1 | Fold), 
                                  data = Dataset[c("Interaction","Fold",Metric)],
                                  REML=TRUE)
    
    AnovaTableInteraction <- anova(ModelFinalInteraction,ddf="Kenward-Roger")
    
    AnovaTableInteraction <- AnovaTableInteraction [,c(3:6)]
    
    # #### Scott Knott
    formula.SK.Interaction <- paste(Metric,"~","Interaction+(1 | Fold)")
    
    anvInteraction <- lmer(as.formula(formula.SK.Interaction), 
                           data=Dataset[c("Interaction","Fold",Metric)])
    
    ### SK results for Algorithm #### 
    SKResults.Interaction <- ScottKnott::SK(x=anvInteraction,
                                            y=Metric,
                                            which='Interaction')
    
    
    SKResultsStatistics.Interaction <- data.frame(SKResults.Interaction$info)
    
    AlgorithmxTechnique.Interaction <- SKResultsStatistics.Interaction [,1]
    Measure.Interaction <- SKResultsStatistics.Interaction [,2]
    Lower.Interaction <- SKResultsStatistics.Interaction[,4]
    Upper.Interaction <- SKResultsStatistics.Interaction[,5]
    
    ClusterMatrix.Interaction <- SKResults.Interaction$out$Result
    ClusterMatrix.Interaction <- ClusterMatrix.Interaction[-1]
    ClusterMatrix.Interaction <- cbind(ClusterMatrix.Interaction, Group = apply(ClusterMatrix.Interaction, 1, 
                                                                                function(x) paste0(na.omit(x), collapse = "")))
    
    Groups.SK.Interaction <- data.frame(AlgorithmxTechnique.Interaction,
                                        Measure.Interaction,
                                        Lower.Interaction,
                                        Upper.Interaction,
                                        ClusterMatrix.Interaction$Group)
    
    colnames(Groups.SK.Interaction)[5] <- c("Cluster")
    
    
    list(
      ExploratoryTable = ExploratoryTable,
      Boxplots = Boxplots,
      Groups.SK.Algorithm =Groups.SK.Algorithm,
      Groups.SK.Interaction = Groups.SK.Interaction) 
  }
  
}



#rm(list=ls())
rep<-8
load(paste0('./Repeats/Rep',rep,'/Result/SK_ready_fourClass_res_f.Rda'))
res_balAcc<-StatisticalComparison(dataset,Metric="Macro.BalancedAccuracy")

res_Recall<-StatisticalComparison(dataset,Metric="Macro.Recall")




res<-StatisticalComparison(dataset,Metric="MAE")
res$Groups.SK.Interaction

format_fnc<-function(dataset,Metric="Macro.BalancedAccuracy"){
  vals<-c()
  combs<-c()
  for(algorithm in unique(dataset$Algorithm)){
    for(technique in unique(dataset$Technique)){
      mean_val<-round(mean(dataset[dataset$Algorithm==algorithm & dataset$Technique==technique,Metric]),3)
      sd_val<-round(sd(dataset[dataset$Algorithm==algorithm & dataset$Technique==technique,Metric]),3)
      comb_jj<-paste0(algorithm,"/",technique)
      
      #means<-c(means,mean_val)
      vals<-c(vals,paste0(mean_val," (",sd_val,")"))
      combs<-c(combs,comb_jj)
    }
  }
  
  df<-data.frame(Combination=combs,Values=vals)
  return(df)
}

res<-format_fnc(dataset,Metric="MAE.3")
res



load('./Result/SK_ready_CADS_res.Rda')
format_fnc_non_important<-function(dataset,Metric="Balanced.Accuracy",type="Algorithm"){
  if(type=="Algorithm"){
    vals<-c()
    algs<-c()
    for(algorithm in unique(dataset$Algorithm)){
      mean_val<-round(mean(dataset[dataset$Algorithm==algorithm,Metric]),3)
      sd_val<-round(sd(dataset[dataset$Algorithm==algorithm,Metric]),3)
      comb_jj<-algorithm
      
      #means<-c(means,mean_val)
      vals<-c(vals,paste0(mean_val," (",sd_val,")"))
      algs<-c(algs,comb_jj)
    }
    
    res<-data.frame(Algorithm=algs,Value=vals)
  }else{
    vals<-c()
    algs<-c()
    
    for(algorithm in unique(dataset$Technique)){
      mean_val<-round(mean(dataset[dataset$Technique==algorithm,Metric]),3)
      sd_val<-round(sd(dataset[dataset$Technique==algorithm,Metric]),3)
      comb_jj<-algorithm
      
      #means<-c(means,mean_val)
      vals<-c(vals,paste0(mean_val," (",sd_val,")"))
      algs<-c(algs,comb_jj)
    }
    
    res<-data.frame(Technique=algs,Value=vals)
  }
  
  return(res)
}
res_balAcc<-StatisticalComparison(dataset,Metric="Sensitivity")
res_balAcc$Groups.SK.Algorithm
res_balAcc$Groups.SK.Technique
res_balAcc$Groups.SK.Interaction
#r<-format_fnc_non_important(dataset,Metric="Sensitivity",type="Algorithm")
#r

r<-format_fnc(dataset,Metric="Sensitivity")
r







###################### Confusion matrix ##############
predictions<-df_whole #From 5.Multiclass_evaluation
preds_rf_enn<-predictions[predictions$Sampling=="enn" & predictions$Algorithm=="RF",]

target<-preds_rf_enn$Target
predicted<-preds_rf_enn$Predicted
class_levels <- 0:3
class_labels <- c("No CAD", "LS CAD", "IS CAD", "HS CAD")

target <- factor(
  preds_rf_enn$Target,
  levels = class_levels,
  labels = class_labels
)

predicted <- factor(
  preds_rf_enn$Predicted,
  levels = class_levels,
  labels = class_labels
)

cm <- caret::confusionMatrix(
  data = predicted,
  reference = target
)
cm_df <- as.data.frame(cm$table)

p<-ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "#F0F5FA", high = "#3182BD") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(limits = rev(levels(predicted)), drop = FALSE) +
  coord_equal() +
  labs(
    x = "Observed CAD severity",
    y = "Predicted CAD severity",
    fill = "Patients"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
library(tidyverse)
tiff('./confusion_matrix.jpeg',units="in",width=10,height=5,res=300)
p
dev.off()
