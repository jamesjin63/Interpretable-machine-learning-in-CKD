# Load required packages
library(survival)
library(survminer)
library(rms)
library(survminer)
# Split data into training (70%) and test (30%) datasets
set.seed(123)

df2=df1 %>% select(- HistoryObesity,-dBPBaseline,-BMIBaseline)
train_idx <- sample(1:nrow(df2), 0.7 * nrow(df2))
train_data <- df2[train_idx, ]
test_data <- df2[-train_idx, ]

# Fit the Cox model to the training data
#关键变量根据一定原则自己灵活定义
cox_model <- coxph(Surv(TimeToEventMonths,EventCKD35) ~ ., data=train_data,x=TRUE)
summary(cox_model)

# Calculate the C-index for train data
train_pred <- predict(cox_model, train_data)
train_cindex <- rcorr.cens(-train_pred, Surv(train_data$TimeToEventMonths, train_data$EventCKD35))
train_cindex[1] 

# Calculate the C-index for test data
test_pred <- predict(cox_model, test_data)
test_cindex <- rcorr.cens(-test_pred, Surv(test_data$TimeToEventMonths, test_data$EventCKD35))
test_cindex[1] 

############################################################
# plot Brier scores for the train and test data over time
############################################################
library(survival)
library(pec)
library(prodlim)
# Estimate Brier score for train data
train_brier <- pec(cox_model, newdata = train_data, 
                   formula = Surv(TimeToEventMonths, EventCKD35) ~ 1, data = train_data)

# Estimate Brier score for test data
test_brier <- pec(cox_model, newdata = test_data, 
                  formula = Surv(TimeToEventMonths, EventCKD35) ~ 1, data = test_data)

plot(train_brier, col = "blue", main = "Brier Score for Train and Test Data", xlab = "Time", ylab = "Brier Score")
lines(test_brier, col = "red")
legend("topright", legend = c("Train", "Test"), col = c("blue", "red"), lty = 1, cex = 0.8)



train_brier_df <- as.data.frame(train_brier$AppErr) %>%
  mutate(Dataset = "Train") %>% 
  mutate(Time=train_brier$time)

test_brier_df <- as.data.frame(test_brier$AppErr) %>%
  mutate(Dataset = "Test",Time=test_brier$time)

brier_df <- rbind(train_brier_df, test_brier_df)
  

p1=ggplot(brier_df, aes(x = Time, y = coxph, color = Dataset)) +
  geom_line() +
  #geom_line(data=brier_df,aes(Time, Reference ,color="black")) +
  labs( x = "Time", y = "Brier Score",title = "Cox") +
  theme_minimal() +ylim(0,0.15)+
  theme(
    axis.text.y = element_text( face="bold", colour = "black"),  
    axis.line = element_line(colour = "black"))




ggsave("brier_scorecox.pdf",width = 8,height = 6,dpi=300)


############################################################
# plot Brier scores for the train and test data over time
############################################################
# https://www.randomforestsrc.org/articles/survival.html#c-index-calculation-1
set.seed(362)
rf_model <- rfsrc(Surv(TimeToEventMonths, EventCKD35) ~ ., 
                  #nodesize = 20,mtry = 8,
                  data = train_data, ntree = 900, importance = TRUE)
rf_train_prediction <- predict.rfsrc(rf_model, newdata = train_data,outcome="train")
yvar_train <- rf_train_prediction$yvar

rf_test_prediction <- predict.rfsrc(rf_model, newdata = test_data,outcome="test")
yvar_test <- rf_test_prediction$yvar

train_c_index <- get.cindex(yvar_train$TimeToEventMonths, yvar_train$EventCKD35,rf_model$predicted.oob)
1-train_c_index
test_c_index <-  get.cindex(yvar_test$TimeToEventMonths, yvar_test$EventCKD35, rf_test_prediction$predicted.oob)
1-test_c_index
# 

bs.rsf <- get.brier.survival(rf_model, cens.mode = "rfsrc")$brier.score
dim(bs.rsf)
bs.rsf1 <- get.brier.survival(rf_test_prediction, cens.mode = "rfsrc")$brier.score
dim(bs.rsf1)

brier_rf= rbind(bs.rsf %>% mutate(Dataset = "Train"),
                bs.rsf1 %>% mutate(Dataset = "Test"))

## plot the brier score

p2=ggplot(brier_rf, aes(x = time, y = brier.score, color = Dataset)) +
  geom_line() +
  #geom_line(data=brier_df, aes(x = Time, y = coxph, color = Dataset))+
  labs(x = "Time", y = "Brier Score",title = "RSF") +
  theme_minimal() +ylim(0,0.15)+
  theme(
    axis.text.y = element_text( face="bold", colour = "black"),  
    axis.line = element_line(colour = "black"))


p1+p2+ plot_annotation(tag_levels = 'A')
ggsave("brier_scorecox.pdf",width = 8,height = 6,dpi=300)




ggplot(brier_df, aes(x = Time, y = coxph, color = Dataset)) +
  geom_line() +
  #geom_line(data=brier_df,aes(Time, Reference ,color="black")) +
  labs(title = "Brier Score for Train and Test Data", x = "Time", y = "Brier Score") +
  theme_minimal() +
  theme(legend.title = element_blank())


jk.obj <- subsample(rf_model)
plot(jk.obj , xlab = "Variable Importance (x 100)", cex = 1.2)


pdf("VIMPsur.pdf", width = 15, height = 20)
par(oma = c(0.5, 10, 0.5, 0.5))
par(cex.axis = 2.0, cex.lab = 2.0, cex.main = 2.0, mar = c(6.0,17,1,1), mgp = c(4, 1, 0))
plot(jk.obj, xlab = "Variable Importance", cex = 1.2)
dev.off()




