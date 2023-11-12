library(tidyverse)
library(tableone)
library(compareGroups)
library(survminer)
rm(list = ls())
df=read.csv("ChronicKidneyDisease_EHRs_from_AbuDhabi.csv",header = T) 
df1=df %>% select(-TIME_YEAR) %>% 
  mutate(Sex=factor(ifelse(Sex==1,"Male","Female")),
         ACEIARB=factor(ifelse(ACEIARB==1,"yes","no")),
         DLDmeds=factor(ifelse(DLDmeds==1,"yes","no")),
         DMmeds=factor(ifelse(DMmeds==1,"yes","no")),
         #EventCKD35=factor(EventCKD35),
         HistoryDiabetes=factor(ifelse(HistoryDiabetes==1,"yes","no")),
         HistoryCHD=factor(ifelse(HistoryCHD==1,"yes","no")),
         HistoryVascular=factor(ifelse(HistoryVascular==1,"yes","no")),
         HistorySmoking=factor(ifelse(HistorySmoking==1,"yes","no")),
         HistoryHTN=factor(ifelse(HistoryHTN==1,"yes","no")),
         HistoryDLD=factor(ifelse(HistoryDLD==1,"yes","no")),
         HistoryObesity=factor(ifelse(HistoryObesity==1,"yes","no"))
  ) %>% 
  mutate(EventCKD35=factor(ifelse(EventCKD35==1,"yes","no")))


######################################################################
##### Survival analysis
######################################################################
library(survival)
df1$tmain <- with(df1, Surv(TimeToEventMonths, EventCKD35 == "yes"))

x=createTable(compareGroups(tmain ~ ., data = df1), 
              show.ratio = TRUE)
export2word(x, file='table1.docx')



x=createTable(compareGroups(tmain ~ AgeBaseline+HistoryDiabetes+
                              eGFRBaseline+HistoryHTN, data = df1), 
              show.ratio = TRUE)

df2=df1 %>% select(-BMIBaseline,-dBPBaseline,-HistoryObesity)
##
x=coxph(Surv(TimeToEventMonths, EventCKD35 == "yes") ~ ., data = df2)
summary(x)


library(ggforce)
ggforest(x)

## for all data
fit = survfit(Surv(TimeToEventMonths, EventCKD35 == "yes") ~ 1, data = df1)
print(fit)

## plot
ggsurvplot(fit)

cox_fit <- survfit(x)
plot(cox_fit)

ggsurvplot(fit,data = df1,
           conf.int = TRUE,  # 增加置信区间
           risk.table = TRUE)


## log-rank test
fit=survfit(Surv(TimeToEventMonths, EventCKD35 == "yes") ~ sex, data = df1)

ggsurvplot(fit,data = df1,
           conf.int = TRUE,  # 增加置信区间
           risk.table = TRUE)



## add p-value in fig

dfx=df1 %>% select(1,3:11,13,20,21)

## for gender
for (i in 1:11) {
  dfx2=dfx %>% select(i,TimeToEventMonths, EventCKD35) 
  dfx1=dfx %>% select(i,TimeToEventMonths, EventCKD35) %>% 
    rename("v"=1)
  ## plot
  fit=survfit(Surv(TimeToEventMonths, EventCKD35 == "yes") ~ v, data = dfx1)
  ggsurvplot(fit,data=dfx1,legend.title = colnames(dfx2)[1],pval = TRUE,
             conf.int = TRUE,  # 增加置信区间
             pval.method=T,palette=c("red","blue"))+
    labs(x="Time (Month)")
  print(colnames(dfx2)[1])
  print(fit)
  ## plot save
  ggsave(paste0(colnames(dfx2)[1],".pdf"),height = 4,width = 6,dpi = 500)
}

######################################################################
##### XGBOOST
######################################################################











