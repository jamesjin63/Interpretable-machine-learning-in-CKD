library(tidyverse)
library(DataExplorer)
df=read.csv("./kidney_disease.csv",header = T) %>% 
  select(-1) %>% 
  set_names('age', 'blood_pressure', 'specific_gravity', 'albumin', 'sugar', 'red_blood_cells', 'pus_cell',
             'pus_cell_clumps', 'bacteria', 'blood_glucose_random', 'blood_urea', 'serum_creatinine', 'sodium',
             'potassium', 'haemoglobin', 'packed_cell_volume', 'white_blood_cell_count', 'red_blood_cell_count',
             'hypertension', 'diabetes_mellitus', 'coronary_artery_disease', 'appetite', 'peda_edema',
             'aanemia', 'class') %>% 
  mutate(class=ifelse(class=="ckd\t","ckd",class)) %>% 
  mutate(coronary_artery_disease=ifelse(coronary_artery_disease=="\tno","no",coronary_artery_disease)) %>% 
  mutate(diabetes_mellitus=ifelse(diabetes_mellitus=="\tyes","yes",diabetes_mellitus)) %>% 
  mutate(diabetes_mellitus=ifelse(diabetes_mellitus==" yes","yes",diabetes_mellitus)) %>% 
  mutate(diabetes_mellitus=ifelse(diabetes_mellitus=="\tno","no",diabetes_mellitus)) %>% 
  mutate(red_blood_cells=factor(red_blood_cells),
         pus_cell=factor(pus_cell),
         pus_cell_clumps=factor(pus_cell_clumps),
         bacteria=factor(bacteria),
         hypertension=factor(hypertension),
         diabetes_mellitus=factor(diabetes_mellitus),
         coronary_artery_disease=factor(coronary_artery_disease),
         appetite=factor(appetite),
         peda_edema =factor(peda_edema ),
         aanemia=factor(aanemia),
         class=factor(class),
         packed_cell_volume=as.numeric(packed_cell_volume),
         white_blood_cell_count=as.numeric(white_blood_cell_count),
         red_blood_cell_count=as.numeric(red_blood_cell_count))

## 数据概览
plot_missing(df)

## Filter missing over 15%
x= profile_missing(df) %>% 
  filter(pct_missing<0.15)
x$feature

## 更改分类变量
df_feature=df %>% select(x$feature) %>% 
  fill(age,blood_pressure,specific_gravity,    
       albumin,sugar,blood_glucose_random,
       blood_urea,serum_creatinine,haemoglobin,.direction = "up",)

## 补全缺失值-uppper
library(tableone)

T1=CreateTableOne(data=df_feature,strata = c("class"))
T3=print(T1,showAllLevels = TRUE)
DT::datatable(T3)



## 可视化
dfa=df_feature %>% select(1:5,10:13,20) %>% 
  gather(x, y,-class) %>% 
  na.omit()


dfa %>% 
  ggplot(aes(x = y, color = class, fill = class)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 3)

ggsave("EDA.pdf",height = 12,width = 12,dpi=500)




## 单因素分析
library(ggpubr)
xa=df_feature %>% select(blood_pressure,class) %>% 
  set_names("value","group")
# Global test
ggboxplot(xa, x = "group", y = "value",
          color = "group", palette = "jco")+
  stat_compare_means(method = "t.test")+
  labs(y="blood_pressure")


## 模型参数选择
subsets <- c(1:5, 10, 15, 20,25,35,55)
set.seed(355)
rfeCtrl <- rfeControl(functions = rfFuncs,
                      method = "cv",
                      verbose = FALSE)

rfProfile <- rfe(x = trainingSet[,-20], 
                 y = trainingSet$class, 
                 sizes = subsets,
                 rfeControl = rfeCtrl)
plot(rfProfile)





###### Machine learning
## 变量筛选
library(caret)
names(getModelInfo())

set.seed(355)
trainIndex = createDataPartition(df_feature$class, p = 0.9, list = FALSE)
trainingSet = df_feature[trainIndex,]
testSet = df_feature[-trainIndex,]

####
set.seed(42)
fitControl <- trainControl(method = "repeatedcv", 
                           number = 5, 
                           repeats = 10, 
                           verboseIter = FALSE)
## model
model_rf = caret::train(class ~ .,
                        data = trainingSet,
                        method = "rf",
                        preProcess = c("scale", "center"),
                        trControl = fitControl)

# 1.LR-linear regression
set.seed(12)
model_lm =train(class ~ .,
                data = trainingSet,
                method = "glm",
                family = "binomial",
                trControl = fitControl)

## probability
final_smote <- data.frame(actual = testSet$class,predict(model_rf, newdata = testSet, type = "prob"))
x=predict(model_rf, newdata = testSet)
confusionMatrix(table(testSet$class, x))


# Visualize the curve using ggplot2 manually
library(yardstick)
library(dplyr)
roc_curve(two_class_example, truth, Class1) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw()

## probability
final_smote <- tibble(truth = testSet$class,predict(model_rf, newdata = testSet, type = "prob")) %>% 
  mutate(predicted=x)

roc_curve(final_smote, truth, ckd) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  geom_abline(
    lty = 2, alpha = 0.5,
    color = "gray50",
    size = 1.2
  ) + 
  geom_path(size = 1.2, color = "cornflowerblue") +
  theme_minimal() +
  labs(x = "False Positive Rate",
       y = "True Positive Rate",
       title = "ROC for RF")




#########################
##### ALE ##### 
#########################
library(pdp)
# Compute the accumulated local effects for the first feature
eff <- FeatureEffect$new(mod, feature = "Vascular")
eff$plot()

#########################
##### PdP ##### 
#########################
# Again, but this time with a partial dependence plot and ice curves
eff <- FeatureEffect$new(mod,feature = "Vascular", method = "pdp+ice",grid.size = 30)
plot(eff)


##IML
library(iml)
mod <- Predictor$new(model_rf, data = df_feature,type = "prob")

# Compute the accumulated local effects for the first feature
eff1 <-eff <- FeatureEffect$new(mod,feature = "albumin", method = "pdp+ice",grid.size = 30)
plot(eff1)

