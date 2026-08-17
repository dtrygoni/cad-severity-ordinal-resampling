############# Descriptive statistics ########
library(arsenal)
rm(list=ls())
load("./Data/Whole_dataset.Rda")

dataset<-dataset[!is.na(dataset$SYNTAX.SCORE),]
dataset$CAD.Risk.Group<-as.factor(ifelse(dataset$SYNTAX.SCORE==0,"No CAD",ifelse(
  dataset$SYNTAX.SCORE<23,"LS CAD",ifelse(dataset$SYNTAX.SCORE<32,"IS CAD","HS CAD")
)))

dataset$CAD.Risk.Group <- factor(
  dataset$CAD.Risk.Group,
  levels = c("No CAD","LS CAD","IS CAD", "HS CAD")
)
props <- prop.table(table(dataset$CAD.Risk.Group))



rs_cols<-colnames(dataset)[which(grepl("rs",colnames(dataset)))]
dataset_clinical<-dataset[!colnames(dataset) %in% rs_cols]
dataset_clinical<-dataset_clinical[!colnames(dataset_clinical) %in% c("SYNTAX.SCORE","SYNTAX.Tertiles")]

# Convert character variables to factors
dataset_clinical[] <- lapply(dataset_clinical, function(x) {
  if (is.character(x)) factor(x) else x
})

dataset_clinical<-dataset_clinical[!colnames(dataset_clinical) %in% c('ID')]
tb_ctrl <- tableby.control(
  test = TRUE,                   # compute p-values
  numeric.stats = c("medianq1q3"),
  numeric.test = "kwt",        # nonparametric test
  cat.stats = c("countpct"),
  cat.test = "chisq"
)
table1 <- tableby(
  CAD.Risk.Group ~ .,
  data = dataset_clinical,
  control = tb_ctrl
)

summary(table1,text=TRUE)
summary(
  table1,
  text = TRUE,
  digits = 2
)
tab_sum <- summary(
  table1,
  digits = 2,
  digits.nonzero = 2,
  digits.integer = 0,
  text = FALSE
)
write2word(tab_sum, file = "./SYNTAX_resampling_suppl.docx")


