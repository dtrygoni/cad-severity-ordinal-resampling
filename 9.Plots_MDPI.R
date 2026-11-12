rm(list=ls())
load("./Results/SK_ready_CADS_res.Rda")
library(dplyr)




########### Scott-Knot results############
### Libraries ####
library(reshape2)
library(psych)
library(ggplot2)
library(ggExtra)
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
    removeGrid(x=TRUE, y=TRUE)+
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


load('./Results/SK_ready_CADHS_res.Rda')
Dataset<-dataset
levels(Dataset$Technique)<-c("ENN","wo-resampling","ROS","RUS","SMOTE","SMOTE ENN" )

Metric<-'Balanced.Accuracy'



technique_order <- c("wo-resampling","RUS","ENN","ROS","SMOTE","SMOTE ENN")
Dataset$Technique <- factor(
  Dataset$Technique,
  levels = technique_order
)
plot_data <- Dataset %>%
  mutate(Technique = factor(Technique, levels = technique_order)) %>%
  group_by(Technique, Algorithm) %>%
  summarise(
    mean_score = mean(.data[[Metric]], na.rm = TRUE),
    min_score  = min(.data[[Metric]], na.rm = TRUE),
    max_score  = max(.data[[Metric]], na.rm = TRUE),
    .groups = "drop"
  )

dodge <- position_dodge(width = 0.6)

p1<-ggplot(plot_data, aes(x = Technique, y = mean_score,
                      shape = Algorithm, group = Algorithm)) +
  geom_errorbar(aes(ymin = min_score, ymax = max_score),
                width = 0.15, position = dodge) +
  geom_point(size = 3.5, position = dodge) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(x = "Technique", y = "Balanced Accuracy" , shape = "Algorithm") +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15),
    axis.text.y = element_text(face = "bold", color = "black", size = 15),
    axis.title  = element_text(face = "bold", color = "black", size = 20),
    legend.position = "right"
  )


library(tidyverse)
tiff('C:/Users/dtryg/Desktop/CADHS_BA.jpeg',units="in",width=10,height=5,res=300)
p1
dev.off()





load('./Results/SK_ready_CADHS_res.Rda')
Dataset<-dataset
levels(Dataset$Technique)<-c("ENN","wo-resampling","ROS","RUS","SMOTE","SMOTE ENN" )
Metric<-'ROC.AUC'

technique_order <- c("wo-resampling","RUS","ENN","ROS","SMOTE","SMOTE ENN")
Dataset$Technique <- factor(
  Dataset$Technique,
  levels = technique_order
)
plot_data <- Dataset %>%
  mutate(Technique = factor(Technique, levels = technique_order)) %>%
  group_by(Technique, Algorithm) %>%
  summarise(
    mean_score = mean(.data[[Metric]], na.rm = TRUE),
    min_score  = min(.data[[Metric]], na.rm = TRUE),
    max_score  = max(.data[[Metric]], na.rm = TRUE),
    .groups = "drop"
  )

dodge <- position_dodge(width = 0.6)

p1<-ggplot(plot_data, aes(x = Technique, y = mean_score,
                          shape = Algorithm, group = Algorithm)) +
  geom_errorbar(aes(ymin = min_score, ymax = max_score),
                width = 0.15, position = dodge) +
  geom_point(size = 3.5, position = dodge) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(x = "Technique", y = "AUC-ROC" , shape = "Algorithm") +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15),
    axis.text.y = element_text(face = "bold", color = "black", size = 15),
    axis.title  = element_text(face = "bold", color = "black", size = 20),
    legend.position = "right"
  )


library(tidyverse)
tiff('C:/Users/dtryg/Desktop/CADHS_AUCROC.jpeg',units="in",width=10,height=5,res=300)
p1
dev.off()
















###################### Sensitivity ##################
load('./Results/SK_ready_CADHS_res.Rda')
Dataset<-dataset
levels(Dataset$Technique)<-c("ENN","wo-resampling","ROS","RUS","SMOTE","SMOTE ENN" )
Metric<-'Sensitivity'

technique_order <- c("wo-resampling","RUS","ENN","ROS","SMOTE","SMOTE ENN")
Dataset$Technique <- factor(
  Dataset$Technique,
  levels = technique_order
)
plot_data <- Dataset %>%
  mutate(Technique = factor(Technique, levels = technique_order)) %>%
  group_by(Technique, Algorithm) %>%
  summarise(
    mean_score = mean(.data[[Metric]], na.rm = TRUE),
    min_score  = min(.data[[Metric]], na.rm = TRUE),
    max_score  = max(.data[[Metric]], na.rm = TRUE),
    .groups = "drop"
  )

dodge <- position_dodge(width = 0.6)

p1<-ggplot(plot_data, aes(x = Technique, y = mean_score,
                          shape = Algorithm, group = Algorithm)) +
  geom_errorbar(aes(ymin = min_score, ymax = max_score),
                width = 0.15, position = dodge) +
  geom_point(size = 3.5, position = dodge) +
  scale_shape_manual(values = c(16, 17, 15, 18)) +
  labs(x = "Technique", y = "Recall" , shape = "Algorithm") +
  theme_bw() +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 15),
    axis.text.y = element_text(face = "bold", color = "black", size = 15),
    axis.title  = element_text(face = "bold", color = "black", size = 20),
    legend.position = "right"
  )


library(tidyverse)
tiff('C:/Users/dtryg/Desktop/CADHS_Recall.jpeg',units="in",width=10,height=5,res=300)
p1
dev.off()

