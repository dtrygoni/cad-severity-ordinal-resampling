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


##### Load CADP  results ######
load('./Results/SK_ready_CADP_res.Rda')

res_balAcc<-StatisticalComparison(dataset,Metric="Balanced.Accuracy")
res_ROCAUC<-StatisticalComparison(dataset,Metric="ROC.AUC")
res_Recall<-StatisticalComparison(dataset,Metric="Sensitivity")
save(res_balAcc,file='./Results/ScottKnott/Scott_Knott_res_CADP_balAcc.Rda')
save(res_ROCAUC,file='./Results/ScottKnott/Scott_Knott_res_CADP_ROCAUC.Rda')

##### Load CADS  results ######
load('./Results/SK_ready_CADS_res.Rda')
res_balAcc<-StatisticalComparison(dataset,Metric="Balanced.Accuracy")
res_ROCAUC<-StatisticalComparison(dataset,Metric="ROC.AUC")
res_Recall<-StatisticalComparison(dataset,Metric="Sensitivity")
save(res_balAcc,file='./Results/ScottKnott/Scott_Knott_res_CADS_balAcc.Rda')
save(res_ROCAUC,file='./Results/ScottKnott/Scott_Knott_res_CADS_ROCAUC.Rda')
#save(res_Recall,file='./Results/ScottKnott/Scott_Knott_res_CADS_Recall.Rda')

##### Load CADHS validation results ######
load('./Results/SK_ready_CADHS_res.Rda')
res_balAcc<-StatisticalComparison(dataset,Metric="Balanced.Accuracy")
res_ROCAUC<-StatisticalComparison(dataset,Metric="ROC.AUC")
res_Recall<-StatisticalComparison(dataset,Metric="Sensitivity")
#res_Recall<-StatisticalComparison(dataset,Metric="Sensitivity")
save(res_balAcc,file='./Results/ScottKnott/Scott_Knott_res_CADHS_balAcc.Rda')
save(res_ROCAUC,file='./Results/ScottKnott/Scott_Knott_res_CADHS_ROCAUC.Rda')
#save(res_Recall,file='./Results/ScottKnott/Scott_Knott_res_CADHS_Recall.Rda')


####### Load Multiclass results ###########



#load("C:/Users/dtryg/Dropbox/Resampling_results/SK_ready_fourClass_res.Rda")
load('./Results/SK_ready_fourClass_res_f.Rda')
res_balAcc<-StatisticalComparison(dataset,Metric="Macro.BalancedAccuracy")
res_Recall<-StatisticalComparison(dataset,Metric="Macro.Recall")
res_kappa<-StatisticalComparison(dataset,Metric="Kappa2")

res_rmse<-StatisticalComparison(dataset,Metric="MSE")
save(res_balAcc,file='./Results/ScottKnott/Scott_Knott_Four_class_res_macrobalAcc.Rda')
save(res_kappa,file='./Results/ScottKnott/Scott_Knott_Four_class_res_kappa.Rda')
save(res_Recall,file='./Results/ScottKnott/Scott_Knott_Four_class_res_macroRecall.Rda')


### MAE metrics ##
res_mae<-StatisticalComparison(dataset,Metric="MAE")

r<-res_mae$Groups.SK.Interaction
r$Final<-paste0(round(r$Measure.Interaction,2),'[',
                round(r$Lower.Interaction,2),'-',
                round(r$Upper.Interaction,2),']')

r[c("AlgorithmxTechnique.Interaction","Final","Cluster")]


res_mae0<-StatisticalComparison(dataset,Metric="MAE.0")
r<-res_mae0$Groups.SK.Interaction
r$Final<-paste0(round(r$Measure.Interaction,2),'[',
                round(r$Lower.Interaction,2),'-',
                round(r$Upper.Interaction,2),']')

r[c("AlgorithmxTechnique.Interaction","Final","Cluster")]

res_mae1<-StatisticalComparison(dataset,Metric="MAE.1")


r<-res_mae1$Groups.SK.Interaction
r$Final<-paste0(round(r$Measure.Interaction,2),'[',
                round(r$Lower.Interaction,2),'-',
                round(r$Upper.Interaction,2),']')

r[c("AlgorithmxTechnique.Interaction","Final","Cluster")]


res_mae2<-StatisticalComparison(dataset,Metric="MAE.2")


r<-res_mae2$Groups.SK.Interaction
r$Final<-paste0(round(r$Measure.Interaction,2),'[',
                round(r$Lower.Interaction,2),'-',
                round(r$Upper.Interaction,2),']')

r[c("AlgorithmxTechnique.Interaction","Final","Cluster")]


res_mae3<-StatisticalComparison(dataset,Metric="MAE.3")


r<-res_mae3$Groups.SK.Interaction
r$Final<-paste0(round(r$Measure.Interaction,2),'[',
                round(r$Lower.Interaction,2),'-',
                round(r$Upper.Interaction,2),']')

r[c("AlgorithmxTechnique.Interaction","Final","Cluster")]

