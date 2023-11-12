library(tidyverse)
library(tableone)
library(DataExplorer)
rm(list = ls())
df=read.csv("ChronicKidneyDisease_EHRs_from_AbuDhabi.csv",header = T) 

df1=df %>% select(-TimeToEventMonths,-TIME_YEAR) %>% 
  mutate(Sex=factor(Sex),
         ACEIARB=factor(ACEIARB),
         DLDmeds=factor(DLDmeds),
         DMmeds=factor(DMmeds),
         #EventCKD35=factor(EventCKD35),
         HistoryDiabetes=factor(HistoryDiabetes),
         HistoryCHD=factor(HistoryCHD),
         HistoryVascular=factor(HistoryVascular),
         HistorySmoking=factor(HistorySmoking),
         HistoryHTN=factor(HistoryHTN),
         HistoryDLD=factor(HistoryDLD),
         HistoryObesity=factor(HistoryObesity)
         ) %>% 
  mutate(EventCKD35=factor(ifelse(EventCKD35==1,"CKD-Yes","CKD-No")))




## Table one
T1=CreateTableOne(data=df1,strata = c("EventCKD35"))
T3=print(T1,showAllLevels = TRUE)
DT::datatable(T3)


## 可视化
dfa=df1 %>% select(2,12,14:20) %>% 
  gather(x, y,-EventCKD35) %>% 
  na.omit()


dfa %>% 
  ggplot(aes(x = y, color = EventCKD35, fill = EventCKD35)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 3)

ggsave("EDA.pdf",height = 12,width = 12,dpi=500)


## ML

###### Machine learning
## 变量筛选
library(caret)
names(getModelInfo())

set.seed(355)
trainIndex = createDataPartition(df1$EventCKD35, p = 0.8, list = FALSE)
trainingSet = df1[trainIndex,]
testSet = df1[-trainIndex,]

####
set.seed(42)
fitControl <- trainControl(method = "repeatedcv", 
                           number = 5, 
                           repeats = 10, 
                           verboseIter = FALSE)
## model
model_rf = caret::train(EventCKD35 ~ .,
                        data = trainingSet,
                        method = "rf",
                        preProcess = c("scale", "center"),
                        trControl = fitControl)

# 1.LR-linear regression
set.seed(12)
model_lmx =train(EventCKD35 ~ .,
                data = trainingSet,
                method = "glm",
                family = "binomial",
                preProcess = c("scale", "center"),
                trControl = fitControl)


# 1.LR-linear regression
set.seed(12)
model_lm =train(EventCKD35 ~ .,
                data = trainingSet,
                method = "xgbTree",
                family = "binomial",
                trControl = fitControl)

# 7.nnet
set.seed(12)
model_nnet = train(EventCKD35 ~ .,
                   data = trainingSet,
                   method='nnet',
                   trControl = fitControl)  
model_nnet



##### list all
all <- list(LM = model_lmx,
            RF = model_rf,
            NNET=model_nnet,
            XGBOOST = model_lm
)

resampling <- resamples(all)
bwplot(resampling)


library(pROC)
library(ROCit)
testingSet=testSet
dfm=tibble()
for (i in c(1:4)) {
  fit=all[[i]]
  result.predicted.prob = predict(fit, newdata=testingSet, type="prob") # Prediction
  ROCit_obj1 <- rocit(score=result.predicted.prob$`CKD-Yes`,class=testingSet$EventCKD35)
  x=ciAUC(ROCit_obj1)
  auc <- round(auc(testingSet$EventCKD35, result.predicted.prob$`CKD-Yes`),3)

  m1=tibble(
    name=names(all)[i],
    TPR=ROCit_obj1$TPR,
    FPR=ROCit_obj1$FPR,
    AUC=paste0(name," AUC=",auc,"(95%CI:",round(x$lower,3),"-",round(x$upper,3),")")
  )
  dfm=m1 %>% bind_rows(dfm)
  print(head(m1))
}

# fit2
## ggplot
ggplot(dfm,aes(x = FPR, y = TPR,color=AUC)) +
  geom_path() +
  labs(
    x = "False Positive Rate (1-Specificity)", 
    y = "True Positive Rate (Sensitivity)")+
  geom_abline(lty = 3) +
  #coord_equal() +
  theme_classic()+
  theme(
    legend.position = c(0.7, 0.2),
    legend.title = element_blank(),
    strip.background = element_blank(),
    text = element_text(size = 12),
  )
ggsave("ROC-logistic.pdf",height = 8,width = 10)


###############################################################
##### 3.Calibration Curves
###############################################################
## Generate the test set results
lift_results <- data.frame(Class = testingSet$EventCKD35)
lift_results$FDA <- predict(model_lmx, testingSet, type = "prob")[,"CKD-Yes"]
lift_results$LDA <- predict(model_rf, testingSet, type = "prob")[,"CKD-Yes"]
lift_results$C5.0 <- predict(model_lm, testingSet,type = "prob")[,"CKD-Yes"]
head(lift_results)



trellis.par.set(caretTheme())
cal_obj <- calibration(Class ~  FDA + LDA + C5.0,
                       data = lift_results,
                       cuts = 13)

ggplot(cal_obj)


###############################################################
##### 1. Variable Importance
###############################################################
df_imp=c()
for (i in c(1:4)) {
  fit=all[[i]]
  # bind data 
  a=varImp(fit)
  print(names(all)[i])
  b=as.data.frame(a$importance)
  df= tibble(variables=rownames(b),
             Score=b$Overall,
             model=names(all)[i]) %>% 
    arrange(desc(Score)) %>% mutate(id=1:n())
  df_imp=rbind(df_imp,df)
}


figimp=function(i=1){
  fit=all[[i]]
  # bind data 
  a=varImp(fit)
  print(names(all)[i])
  b=as.data.frame(a$importance)
  df= tibble(variables=rownames(b),
             Score=b$Overall,
             model=names(all)[i]) %>% 
    arrange(desc(Score)) %>% mutate(id=1:n())
  #df_imp=rbind(df_imp,df)
  p=ggplot(df)+
    geom_bar(aes(x=reorder(id,Score),Score),
             stat = "identity")+
    scale_x_discrete(labels =rev(df$variables))+
    coord_flip()+
    theme_bw() +
    theme(
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      axis.text.y = element_text( face="bold", colour = "black"),  
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"))+
    labs(y="Variable importance",x="",title =unique(df$model) )
  
  return(p)
}

library(patchwork)
x1=figimp(i=1)
x2=figimp(i=2)
x3=figimp(i=3)
x4=figimp(i=4)

(x1+x2)/(x3+x4)
ggsave("VIP-logistic.pdf",height = 12,width = 16)


#########################
##### ALE ##### 
#########################
library(pdp)
library(iml)


mod <- Predictor$new(model_lm, data = trainingSet,type = "prob")
# Compute the accumulated local effects for the first feature
eff <- FeatureEffect$new(mod, feature = "AgeBaseline")
eff$plot()

#########################
##### PdP ##### 
#########################
# Again, but this time with a partial dependence plot and ice curves
eff <- FeatureEffect$new(mod,feature = "eGFRBaseline", method = "pdp+ice",grid.size = 30)
p1=plot(eff)+ ggtitle("")+xlab("eGFRBaseline")+
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    axis.text.y = element_text( face="bold", colour = "black"),  
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))


eff <- FeatureEffect$new(mod,feature = "AgeBaseline", method = "pdp+ice",grid.size = 30)
p2=plot(eff)+ ggtitle("")+xlab("AgeBaseline")+
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    axis.text.y = element_text( face="bold", colour = "black"),  
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))



eff <- FeatureEffect$new(mod,feature = "CholesterolBaseline", method = "pdp+ice",grid.size = 30)
p3=plot(eff)+ggtitle("")+xlab("CholesterolBaseline")+
  theme_bw() +
  theme(
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    axis.text.y = element_text( face="bold", colour = "black"),  
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))

x=eff$results %>% rename("class"=2) %>% filter(class=="CKD-Yes") %>% filter(.type=="pdp")

p1/p2/p3+ plot_annotation(tag_levels = 'A')
 
ggsave("PDP-logistic.pdf",height = 14,width = 12)

#################### #################### 
## RCS
#################### #################### 

# RCS curves for a Cox regression model
x=rcsplot(data = trainingSet,
          outcome = "EventCKD35",
          exposure = "AgeBaseline",
          covariates = c("DMmeds", "eGFRBaseline", "CreatinineBaseline", "AgeBaseline"))




### Explainer for survival

library("DALEX")
explain_titanic_rf <- explain(model_lm, 
                              data = trainingSet[,-20],
                              y = trainingSet$EventCKD35== "CKD-Yes", 
                              label = "XGBOOST")

library("iBreakDown")
rf_la <- local_attributions(explain_titanic_rf, trainingSet[19,])
rf_la
plot(rf_la)

rf_la <- local_attributions(explain_titanic_rf, trainingSet[71,],path = "average")
plot(rf_la)



xa=rbind(trainingSet[19,],trainingSet[71,])



### HR
library(compareGroups)
df=read.csv("ChronicKidneyDisease_EHRs_from_AbuDhabi.csv",header = T) 

df1=df %>% select(-TIME_YEAR) %>% 
  mutate(Sex=factor(Sex),
         ACEIARB=factor(ACEIARB),
         DLDmeds=factor(DLDmeds),
         DMmeds=factor(DMmeds),
         #EventCKD35=factor(EventCKD35),
         HistoryDiabetes=factor(HistoryDiabetes),
         HistoryCHD=factor(HistoryCHD),
         HistoryVascular=factor(HistoryVascular),
         HistorySmoking=factor(HistorySmoking),
         HistoryHTN=factor(HistoryHTN),
         HistoryDLD=factor(HistoryDLD),
         HistoryObesity=factor(HistoryObesity)
  ) %>% 
  mutate(EventCKD35=factor(ifelse(EventCKD35==1,"CKD-Yes","CKD-No")))
