# Load required packages
library(tidyverse)
library(survival)
library(survminer)
library(rms)
library(survminer)
library(randomForestSRC)
rm(list = ls())
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

df2=df1 %>% select(-HistoryObesity,-dBPBaseline,-BMIBaseline)
threedata=df2
library(pec)
set.seed(19910220)
index1 <- sample(1:nrow(threedata), round(0.8*nrow(threedata)))
train_set_1 <- threedata[index1,]
test_set_1 <- threedata[-index1,]
cox1 <- coxph(Surv(TimeToEventMonths, EventCKD35) ~., x = T, y = T, 
              data = train_set_1)

rsf1 <- rfsrc(Surv(TimeToEventMonths, EventCKD35)~.,
              data=train_set_1,
              ntree=1500,forest=TRUE,
              tree.err=T, importance=T,
              na.action = "na.impute")
cf1=calPlot(list(cox1,   rsf1),
            col=c("black",'red'),
            time=36,
            type="survival",
            legend=F, 
            data = train_set_1,
            splitMethod = "cv",
            B=10)
legend("topleft", legend=c("Cox","randomForestSRC"),
       col=c("black",'red'), lty=1, cex=0.8)



# overtime cindex using pec package ###
# overtime C index in complete cases for three year cohort (training set)
ApparrentCindex1 <- pec::cindex(list(cox1,rsf1),
                                formula=Surv(TimeToEventMonths, EventCKD35) ~ .,
                                data=train_set_1,
                                splitMethod="cv",
                                B=10,
                                eval.times=seq(0,100,1))
plot(ApparrentCindex1,legend=c(2,1),xlim=c(0,100))

xa=ApparrentCindex1$AppCindex %>% as.data.frame() %>% set_names("Cox","RSF") %>% 
  mutate(time=1:n()) %>% gather("cindex","V",-time)
p11=ggplot(xa, aes(x = time, y = V,fill=cindex,color=cindex)) +
  geom_line() +ylim(0.4,1)+
  geom_hline(yintercept=c(0.5), linetype='dashed')+
  theme_minimal() +
  labs(title = "Train set", x = "Survival time (month)", y = "Concordance Index") +
  theme(
    axis.text.y = element_text( face="bold", colour = "black"),  
    axis.line = element_line(colour = "black"))



ApparrentCindex1 <- pec::cindex(list(cox1,                    
                                     rsf1),
                                formula=Surv(TimeToEventMonths, EventCKD35) ~ .,
                                data=test_set_1,
                                splitMethod="cv",
                                B=10,
                                eval.times=seq(0,100,1))


plot(ApparrentCindex1,legend=c(2,1),xlim=c(0,100))

xa=ApparrentCindex1$AppCindex %>% as.data.frame() %>% 
  set_names("Cox","RSF") %>% 
  mutate(time=1:n()) %>% gather("cindex","V",-time)
p22=ggplot(xa, aes(x = time, y = V,fill=cindex,color=cindex)) +
  geom_line() +ylim(0.4,1)+
  geom_hline(yintercept=c(0.5), linetype='dashed')+
  theme_minimal() +
  labs(title = "Test set", x = "Survival time (month)", y = "Concordance Index") +
  theme(
    axis.text.y = element_text( face="bold", colour = "black"),  
    axis.line = element_line(colour = "black"))

library(patchwork)
p11+p22+ plot_annotation(tag_levels = 'A')
ggsave("Cindex-x.pdf",width = 12,height = 8,dpi=300)

#(p11+p22)/(p1+p2)+ plot_annotation(tag_levels = 'A')

ggsave("Cindex-x.pdf",width = 12,height = 8,dpi=300)

