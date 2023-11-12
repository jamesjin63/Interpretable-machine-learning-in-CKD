library(tidyverse)
library(survival)
library(survminer)
library(rms)
library(survminer)
#install.packages("plotRCS")
library(plotRCS)
# View data
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
  ) 

# RCS curves for a Cox regression model
x1=rcsplot(data = df1, conf.int = F,
          outcome = "EventCKD35",
          time = "TimeToEventMonths",
          exposure = "AgeBaseline",
          covariates = c("DMmeds", "eGFRBaseline", "CreatinineBaseline", "AgeBaseline"))

# RCS curves for a Cox regression model
x2=rcsplot(data = df1, conf.int = F,
          outcome = "EventCKD35",
          time = "TimeToEventMonths",
          exposure = "CreatinineBaseline",
          covariates = c("DMmeds", "eGFRBaseline", "CreatinineBaseline", "AgeBaseline"))


# RCS curves for a Cox regression model
x3=rcsplot(data = df1, conf.int = F,
          outcome = "EventCKD35",
          time = "TimeToEventMonths",
          exposure = "eGFRBaseline",
          covariates = c("DMmeds", "eGFRBaseline", "CreatinineBaseline", "AgeBaseline"))

library(patchwork)
x1/x2/x3+ plot_annotation(tag_levels = 'A')
ggsave("rcs.pdf",width = 12,height = 10,dpi=300)

x1+ylim(0,30)

############################################################
#################### Compare table
####################
############################################################
library(compareGroups)

## compare
res <- compareGroups(EventCKD35 ~ ., data =df1)
createTable(res,show.all = T)
export2csv(createTable(res,show.all = T), file='table1.csv')
export2word(createTable(res,show.all = T), file='table1.docx')

xa=read.csv('table1.csv',header=T) %>% 
  set_names("Variables","ALL","SP_Mets","SP_No_Mets","p-value")
DT::datatable(xa)



############################################################
#################### HR
####################
############################################################












