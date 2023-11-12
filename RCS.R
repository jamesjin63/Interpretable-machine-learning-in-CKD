# Load required packages
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


# Split data into training (70%) and test (30%) datasets
set.seed(123)
train_idx <- sample(1:nrow(df1), 0.7 * nrow(df1))
train_data <- df1[train_idx, ]
test_data <- df1[-train_idx, ]

# Fit the Cox model to the training data
#关键变量根据一定原则自己灵活定义
cox_model <- coxph(Surv(TimeToEventMonths,EventCKD35) ~ ., data=train_data,x=TRUE)
summary(cox_model)

# Calculate the C-index for train data
train_pred <- predict(cox_model, train_data)
train_cindex <- rcorr.cens(-train_pred, Surv(train_data$TimeToEventMonths, train_data$EventCKD35))
train_cindex 

# Calculate the C-index for test data
test_pred <- predict(cox_model, test_data)
test_cindex <- rcorr.cens(-test_pred, Surv(test_data$TimeToEventMonths, test_data$EventCKD35))
test_cindex 


cox_model <- coxph(Surv(TimeToEventMonths,EventCKD35) ~ DMmeds+eGFRBaseline+AgeBaseline,
                   data=df1,x=TRUE)
summary(cox_model)
ggforest(cox_model,df1)




# RCS curves for a Cox regression model
rcsplot(data = df1,
        outcome = "EventCKD35",
        time = "TimeToEventMonths",
        exposure = "eGFRBaseline",
        covariates = c("DMmeds", "eGFRBaseline", "CreatinineBaseline", "AgeBaseline"))

# RCS curves for a Cox regression model
x=rcsplot(data = df1,
        outcome = "EventCKD35",
        time = "TimeToEventMonths",
        exposure = "AgeBaseline",
        covariates = c("DMmeds", "eGFRBaseline", "CreatinineBaseline", "AgeBaseline"))



mv_fit <- coxph(Surv(TimeToEventMonths,EventCKD35) ~ ., data=df1,x=TRUE)
ccox <- cox.zph(mv_fit)
print(ccox)
#                    chisq df     p
# age                 2.260  1 0.133
# serum_creatinine    3.646  1 0.056
# serum_sodium        1.531  1 0.216
# high_blood_pressure 0.313  1 0.576
# GLOBAL              7.324  4 0.120
options(repr.plot.width=10, repr.plot.height=40)
ggcoxzph(ccox)



mo <- coxph(Surv(TimeToEventMonths,EventCKD35) ~ ., data=df1,x=TRUE)
ggforest(mo,df1)

mo <- coxph(Surv(TimeToEventMonths,EventCKD35) ~ DMmeds+eGFRBaseline+CreatinineBaseline+AgeBaseline, data=df1,x=TRUE)
ggforest(mo,df1)


xa=tibble()
for (i in 1:19) {
  dfa=df1 %>% select(i,20,21) %>% rename("x"=1)
  mv_fit <- coxph(Surv(TimeToEventMonths,EventCKD35) ~ x, data=dfa,x=TRUE)
  x=summary(mv_fit )
  x$coefficients[5]
  xb=tibble(name=colnames(df1)[i],p=  x$coefficients[5])
  xa=rbind(xa,xb)
}
xa %>% filter(p>0.05)

# Load necessary libraries
library(glm)

# Create a logistic regression model
model <- glm(EventCKD35 ~ ., data = df1, family = binomial(link = "logit"))
ggforest(model,df1)


df2=df1 %>% select(- HistoryObesity,-dBPBaseline,-BMIBaseline)
mo <- coxph(Surv(TimeToEventMonths,EventCKD35) ~ ., data=df2,x=TRUE)
ggforest(mo,df1)







