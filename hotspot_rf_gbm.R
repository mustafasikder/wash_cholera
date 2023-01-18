# This script runs random forest and GBM hotspot model using country fold cross validation 
# and using the full data. 
#
# Paper output: 




rm(list= ls())
setwd("X:/Spatial Stat/WASH Cholera/clean_repo") 
output_dir<- 'X:/Spatial Stat/WASH Cholera/clean_repo/results/'

library(randomForest)
library(permimp)
library(tidyverse)
library(plotROC)
library(pROC)
library(ROSE)
library(data.table)
library(viridis)
library(ggthemes)
library(boot)
library(gridExtra)
library(vip)
library(patchwork)

theme_set(theme_minimal())
theme_update(legend.position = "bottom")
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")
options(digits=2)


dist_df<- read.csv(paste0(output_dir, 'district_6v.csv'))
dist_df<- dist_df[complete.cases(dist_df), ]

########################################################################
### -------------------------- Hotspot ------------------------------###
########################################################################

###############################################################################
##### ------------------- Country Fold Cross-Validation -----------------######
###############################################################################
hotspot_df<- dist_df[,c(4:9, 11)]
dist_dt<- as.data.table(dist_df)
adm0<- dist_dt[,.N,by= NAME_0.x]
folds.h<- rep(1:nrow(adm0), adm0$N)
hotspot_df$hotspot<- as.factor(hotspot_df$hotspot)

############################# Random Forest CV Model #############################
test.list.h<- list()
train.list.h<- list()
rf.list.h<- list()
for(i in 1:nrow(adm0)){
  testIndexes <- which(folds.h==i,arr.ind=TRUE)
  testData <- hotspot_df[testIndexes, ]
  trainData <- hotspot_df[-testIndexes, ]
  trainData<- ovun.sample(hotspot~., data=trainData,
                          N=nrow(trainData),  p= 0.5,
                          seed=1, method="both")$data
  test.list.h[[i]]<- testData
  train.list.h[[i]]<- trainData
  rf<- randomForest(
    hotspot~.,
    data = trainData,
    ntree=500,
    localImp= TRUE, 
    keep.forest= TRUE, 
    keep.inbag= TRUE
  )
  rf.list.h[[i]]<- rf
  rf.list.h[[i]]$pred.newdata<- predict(rf, newdata = testData)
  rf.list.h[[i]]$pred.newdata.vote<- predict(rf, newdata = testData, type = "vote",norm.votes=TRUE)
  rf.list.h[[i]]$testdata<- testData
  rf.list.h[[i]]$conf.matrix<- table(testData$hotspot ,rf.list.h[[i]]$pred.newdata)
}

saveRDS(rf.list.h, "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf.list.h_6vars.RDS")

# data preparation for the plot RoC plot 
hotspot_predict<- data.table()
full.hotspot.data<- data.table()
for (i in 1:length(rf.list.h)){
  hotspot_predict<- data.table() # since the lengths are different
  temp<- hotspot_predict[ , c("cv.number", "Predicted", "Observed", "vote.0", "vote.1") := 
                              list(rep(paste0("cv.", i), length(rf.list.h[[i]]$pred.newdata)), 
                                   rf.list.h[[i]]$pred.newdata, 
                                   rf.list.h[[i]]$testdata$hotspot, 
                                   rf.list.h[[i]]$pred.newdata.vote[,1], 
                                   rf.list.h[[i]]$pred.newdata.vote[,2]),]
  full.hotspot.data<- rbind(full.hotspot.data, temp)
}

full.hotspot.data$Predicted<- as.numeric(as.character(full.hotspot.data$Predicted))
full.hotspot.data$Observed<- as.numeric(as.character(full.hotspot.data$Observed))

saveRDS(full.hotspot.data, "X:/Spatial Stat/WASH Cholera/clean_repo/results/full.hotspot.data_6vars.RDS")
full.hotspot.data<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/full.hotspot.data_6vars.RDS")

### Bootstrap for cvACU and 95% CI
actual<- data.table()
predicted<- data.table()
for (i in 1:nrow(adm0)){
  actual<- rbind(actual, rf.list.h[[i]]$testdata$hotspot)
  predicted<- rbind(predicted, rf.list.h[[i]]$pred.newdata.vote[,2])
}
auc(as.factor(actual$x), predicted$x)

# Bootstrap for cvAUC for 95% CI
auc_data<- data.frame(actual= as.factor(actual$x), predicted= predicted$x)
fc<- function(d, i){
  d2<- d[i,]
  return(auc(d2$actual, d2$predicted))
}
cvAUC<- boot(auc_data, statistic=fc, R=5000)
rf_cvAUC_cl<- boot.ci(cvAUC, conf=0.95, type="bca")

#################################### GBM CV Model ###################################
hotspot_df_gbm<- hotspot_df
hotspot_df_gbm$hotspot<- as.numeric(hotspot_df_gbm$hotspot)-1 # as GBM wants non factor for two categories 

test.list.h<- list()
train.list.h<- list()
gbm.list.h<- list()
for(i in 1:nrow(adm0)){
  testIndexes <- which(folds.h==i,arr.ind=TRUE)
  testData <- hotspot_df_gbm[testIndexes, ]
  trainData <- hotspot_df_gbm[-testIndexes, ]
  trainData<- ovun.sample(hotspot~., data=trainData,
                          N=nrow(trainData),  p= 0.5,
                          seed=1, method="both")$data
  test.list.h[[i]]<- testData
  train.list.h[[i]]<- trainData
  gbm.hp<- gbm(
    hotspot~.,
    data = trainData,
    distribution="bernoulli",  
    n.trees=500, 
    verbose=T
  )
  gbm.list.h[[i]]<- gbm.hp
  gbm.list.h[[i]]$pred.newdata<- predict(gbm.hp, newdata = testData)
  gbm.list.h[[i]]$pred.newdata.vote<- predict(gbm.hp, newdata = testData, type= "response")
  gbm.list.h[[i]]$testdata<- testData
  gbm.list.h[[i]]$conf.matrix<- table(testData$hotspot ,gbm.list.h[[i]]$pred.newdata)
}

############### ADM0 as folds for cvRMSE - hotspot model ####################
actual<- data.table()
predicted<- data.table()
for (i in 1:nrow(adm0)){
  actual<- rbind(actual, gbm.list.h[[i]]$testdata$hotspot)
  predicted<- rbind(predicted, gbm.list.h[[i]]$pred.newdata.vote)
}
auc(as.factor(actual$x), predicted$x)

# Bootstrap for cvAUC for 95% CI
auc_data<- data.frame(actual= as.factor(actual$x), predicted= predicted$x)
fc<- function(d, i){
  d2<- d[i,]
  return(auc(d2$actual, d2$predicted))
}
cvAUC_gbm<- boot(auc_data, statistic=fc, R=5000)
gbm_cvAUC_cl<- boot.ci(cvAUC_gbm, conf=0.95, type="bca")

### --------- data preparation for the RoC plot (Fig 3B) -------------------- ###
## prepare data for ROC plot
## NOTE: the plot also requires full model data 
cv.hotspot.data.gbm<- data.table()
for (i in 1:length(gbm.list.h)){
  hotspot_predict.gbm<- data.table() # since the lengths are different
  temp<- hotspot_predict.gbm[ , c("cv.number", "Observed", "vote.1") := 
                                list(rep(paste0("cv.", i), length(gbm.list.h[[i]]$pred.newdata)), 
                                     gbm.list.h[[i]]$testdata$hotspot, 
                                     gbm.list.h[[i]]$pred.newdata.vote),]
  cv.hotspot.data.gbm<- rbind(cv.hotspot.data.gbm, temp)
}

cv.hotspot.data.gbm$Observed<- as.numeric(as.character(cv.hotspot.data.gbm$Observed))
saveRDS(cv.hotspot.data.gbm, "X:/Spatial Stat/WASH Cholera/clean_repo/results/cv.hotspot.data.gbm_6vars.RDS")
cv.hotspot.data.gbm<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/cv.hotspot.data.gbm_6vars.RDS")

###############################################################################
##### ---------------------------- Full data model  ---------------------######
###############################################################################

########################## GBM full data Model #############################
gbm.hp<- gbm(hotspot~.,  
             data = hotspot_df_gbm, 
             distribution="bernoulli",  
             n.trees=500, 
             verbose=T)
gbm.hp.full.pred<- predict.gbm(gbm.hp, hotspot_df_gbm, n.trees = 500, distribution="bernoulli", type= "response")

gbm_hp_full.df<- data.frame(Observed.f= hotspot_df_gbm$hotspot, vote.1.f= gbm.hp.full.pred)
saveRDS(gbm_hp_full.df, file = "X:/Spatial Stat/WASH Cholera/clean_repo/results/gbm_hp_full.df_6vars.RDS")
gbm_hp_full.df<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/gbm_hp_full.df_6vars.RDS")



##### combine the full model and cv model results for ROC plot
## not reported since RF cvAUC is greatter 
combined_df_hp.gbm<- data.frame(cv.hotspot.data.gbm, gbm_hp_full.df)
combined_df_hp.gbm$Observed<- as.numeric(as.character(combined_df_hp.gbm$Observed))
combined_df_hp.gbm$Observed.f<- as.numeric(as.character(combined_df_hp.gbm$Observed.f))
saveRDS(combined_df_hp.gbm, "X:/Spatial Stat/WASH Cholera/clean_repo/results/combined_df_hp.gbm_6vars.RDS")
combined_df_hp.gbm<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/combined_df_hp.gbm_6vars.RDS")

auc_plot<- ggplot(data=combined_df_hp.gbm)+ 
  geom_roc(aes(m= vote.1, d= Observed, color= cv.number), n.cuts=0, linealpha = .3) + 
  geom_roc(aes(m= vote.1.f, d= Observed), n.cuts=0, linealpha = 1) + 
  coord_equal()+
  theme(legend.position = 'none')+
  scale_color_grey(end= 0) +
  labs(x= "False positive rate", y= "True positive rate")

saveRDS(auc_plot, paste0(output_dir, 'gbm_hotspot_auc_6vr_plot.RDS'))
png('X:/Spatial Stat/WASH Cholera/clean_repo/results/hotspot_auc_6vr.png', 
    width = 6, height = 6, res = 400, units = 'in')
auc_plot
dev.off()

################## GBM full model importance plot hotspot (Fig: SI) ########################
vip_full_hp<- vip(gbm.hp)
v_imp_hp_gbm<- data.table(vip_full_hp$data)
setorder(v_imp_hp_gbm, Variable)
v_imp_hp_gbm$var<- c("Piped or other improved sanitation", 
                     "Open defecation", 
                     "Piped sanitation", 
                     "Piped or other improved water", 
                     "Piped water", 
                     "Surface water")
saveRDS(v_imp_hp_gbm, "X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_hp_gbm_6vars.RDS")
v_imp_hp_gbm<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_hp_gbm_6vars.RDS")

mean_impurity_decrease.hp.gbm <-
  ggplot(data = v_imp_hp_gbm, aes(Importance, reorder(var, Importance), fill = var)) +
  geom_bar(stat = "identity") +
  scale_fill_manual( values = c('#2d588a', '#58954c', '#e9a044', '#c12f32', '#723e77', '#7d807f'))+#("#264653","#2a9d8f","#F16745","#FFC65D","#f4a261","#4CC3D9","#93648D","#457b9d")) + 
  labs(y = NULL, x = "Variable Importance (GBM)") + theme(
    legend.title = element_blank(),
    legend.position = 'none',
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size= 7, face = "plain"),
    axis.text.x = element_text(size = 6),
    #plot.background = element_rect(color = "White"),
    panel.grid.major = element_blank(),
    panel.grid.minor =  element_blank(), 
    panel.grid.minor.x = element_blank())+ 
  scale_x_continuous(position = "top") #, limits = c(min(v_imp_hp$mean_accuracy_decrease), max(v_imp_hp$mean_accuracy_decrease)), breaks = c(as.vector(summary(v_imp_hp$mean_accuracy_decrease)[c(1,3,6)])))

png('X:/Spatial Stat/WASH Cholera/clean_repo/results/hotspot_vimp_gbm_6vars.png', 
    width = 6, height = 6, res = 400, units = 'in')
mean_impurity_decrease.hp.gbm
dev.off()


########################## Random forest full data Model #############################

rf_hp_full<- randomForest(
  hotspot~.,
  data = hotspot_df,#_ovun,
  ntree=500,
  localImp= TRUE, 
  keep.forest= TRUE, 
  keep.inbag= TRUE
)

rf_hp_full.df<- data.frame(Predicted.f= rf_hp_full$predicted, Observed.f= rf_hp_full$y, vote.0.f= rf_hp_full$votes[,1], vote.1.f= rf_hp_full$votes[,2])
saveRDS(rf_hp_full.df, file = "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_hp_full.df_6vars.RDS")
rf_hp_full.df<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_hp_full.df_6vars.RDS")

imp.h_full<- permimp::permimp(rf_hp_full, conditional = T, do_check = FALSE)
saveRDS(imp.h_full, "X:/Spatial Stat/WASH Cholera/clean_repo/results/permutationImpFullData_hotspot_6vars.rds")
imp.h_full<- readRDS("results/permutationImpFullData_hotspot_6vars.rds")

## AUC  for full data model
auc(rf_hp_full$y , rf_hp_full$votes[,1], quiet= T)
#95% CI of AUC
ci.auc(auc(rf_hp_full$y , rf_hp_full$votes[,1], quiet= T))
# Note: we will report RF model results in the main text since the cvAUC is larger than GBM cvAUC 

#### combining the full model and cv model  data frames for the ROC plot
combined_df_hp<- data.frame(full.hotspot.data, rf_hp_full.df)
combined_df_hp$Observed<- as.numeric(as.character(combined_df_hp$Observed))
combined_df_hp$Observed.f<- as.numeric(as.character(combined_df_hp$Observed.f))
combined_df_hp$Predicted<- as.numeric(as.character(combined_df_hp$Predicted))
combined_df_hp$Predicted.f<- as.numeric(as.character(combined_df_hp$Predicted.f))

saveRDS(combined_df_hp, "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_combined_df_hp_6vars.RDS")
combined_df_hp<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_combined_df_hp_6vars.RDS")
### to make sure combined_df_hp$Observed!=combined_df_hp$Observed.f to use in the main2 plot 
sum(combined_df_hp$Observed!=combined_df_hp$Observed.f)

################################## Hotspot ROC AUC Plot (Fig: 3B) ##########################
#### Using RF AUC plot in the main text as RF cvAUC> GBM cvAUC
rf_auc_plot<- ggplot(data=combined_df_hp)+ 
  geom_roc(aes(m= vote.1, d= Observed, color= cv.number), n.cuts=0, linealpha = .3, size= .2) + 
  geom_roc(aes(m= vote.1, d= Observed), n.cuts=0, linealpha = 1) + # using CV model data 
  coord_equal()+
  theme(legend.position = 'none')+
  scale_color_grey(end= 0) +
  #scale_color_tableau()+ 
  # scale_color_viridis_d()+
  geom_point(aes(x= 0.447, y= 0.796), color= 'red', size= 2)+ # x, y from roc_cutpoint.R script; max Youden's index
  labs(x= "False positive rate", y= "True positive rate")

# saveRDS(rf_auc_plot, paste0(output_dir, 'rf_hotspot_auc_6vr_plot.RDS'))
saveRDS(rf_auc_plot, paste0(output_dir, 'rf_hotspot_auc_6vr_plot_light.RDS'))
# saveRDS(rf_auc_plot, paste0(output_dir, 'rf_hotspot_auc_6vr_plot.RDS'))
png('X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_hotspot_roc_6vars.png', 
    width = 6, height = 6, res = 400, units = 'in')
rf_auc_plot
dev.off()

################## RF full model importance plot hotspot (Fig: SI) ########################
v_imp_hp<- cbind(imp.h_full$values)
rownames(v_imp_hp)<- c("Piped water", 
                       "Piped or other improved water", 
                       "Surface water",
                       "Piped sanitation", 
                       "Piped or other improved sanitation", 
                       "Open defecation")
v_imp_hp<- cbind(rownames(v_imp_hp), v_imp_hp)
v_imp_hp<- as.data.frame(v_imp_hp)
names(v_imp_hp)<- c("var", 'mean_accuracy_decrease')
v_imp_hp$mean_accuracy_decrease<- as.numeric(v_imp_hp$mean_accuracy_decrease)
saveRDS(v_imp_hp, "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_v_imp_hp_6vr.RDS")
v_imp_hp<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_v_imp_hp_6vr.RDS")

mean_impurity_decrease.hp <-
  ggplot(data = v_imp_hp,
         aes(
           mean_accuracy_decrease,
           reorder(var, mean_accuracy_decrease),
           fill = var
         )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c('#2d588a', '#58954c', '#e9a044', '#c12f32', '#723e77', '#7d807f')
  ) + labs(y = NULL, x = "Conditional variable importance (RF)") + theme(
    legend.title = element_blank(),
    legend.position = 'none',
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size= 8),
    axis.text.x = element_text(size = 6),
    #plot.background = element_rect(color = "White"),
    panel.grid.major = element_blank(),
    panel.grid.minor =  element_blank()
  )+ 
  scale_x_continuous(position = "top")

png('X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_hotspot_vimp_6vars.png', 
    width = 5, height = 3, res = 300, units = 'in')
mean_impurity_decrease.hp
dev.off()


################### FIGURE SI Hotspot Variable importance plot: RF + GBM ###################
hp_imp_plot<- mean_impurity_decrease.hp + mean_impurity_decrease.hp.gbm

hp_var_imp_plot<- hp_imp_plot+ plot_annotation(tag_levels = 'A') & #, tag_prefix = '3'
  theme(plot.tag.position = c(.4, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

png(paste0(output_dir, 'plots/fig_si_hotspot_var_imp_6vs.png'), 
    width = 8, height = 3, res = 400, units = 'in')
hp_var_imp_plot
dev.off()







plot4<- main2+annotation_custom(ggplotGrob(mean_impurity_decrease.hp), xmin = .3 , xmax = 1.5, ymin = -.2, ymax = .7)

plot3.v1<- ggpubr::ggarrange(plot3, plot4, ncol = 1, labels = c("A", "B"),  heights = c(1, 1))

save.image(file = "revised.model.results.RData")


# making ROC manually -----------------------------
# discrimination threshold
d<- seq(0, 1, by= 0.05)

# for detecting  hotspots 
sensetivity<- function(dis_thr){ # TP/TP+FN
  length(full.hotspot.data[cv.number=="cv.9" & Observed == 1 & vote.1>=dis_thr, ]$Predicted)/ # TP= when the vote was >= dis_thr how many were predicted as 1
    (length(full.hotspot.data[cv.number=="cv.9" & Observed == 1 & vote.1>=dis_thr, ]$Predicted) + length(full.hotspot.data[cv.number=="cv.9" & Observed == 1 & vote.1<dis_thr, ]$Predicted))
}

specificity<- function(dis_thr){ # TN/TN+ FP
  length(full.hotspot.data[cv.number=="cv.9" & vote.1<=dis_thr & Observed ==0,]$Predicted)/
    (length(full.hotspot.data[cv.number=="cv.9" & vote.1<=dis_thr & Observed ==0, ]$Predicted) + length(full.hotspot.data[cv.number=="cv.9" & vote.1>dis_thr & Observed==0,]$Predicted))
}
temp_roc<- NULL
for(i in d){
  roc1<- c("FPR"= 1-specificity(i), "TPR"= sensetivity(i) )
  temp_roc<- rbind(temp_roc, roc1)
}
ggplot(data = data.frame(temp_roc), aes(x= FPR, y= TPR))+ geom_point()

## ROC for cv 9 
ggplot(data=full.hotspot.data[cv.number=="cv.9",], aes(m= vote.1, d= Observed, color= cv.number))+
  geom_roc(n.cuts=0, linealpha = .5) +
  coord_equal()+
  theme(legend.position = 'none')+
  scale_color_tableau()+
  # scale_color_viridis_d()+
  labs(x= "False positive rate", y= "True positive rate")

# full plot with cv9 highlighted
ggplot(data=full.hotspot.data, aes(m= vote.1, d= Observed, color= cv.number))+ 
  geom_roc(n.cuts=0, linealpha = .5) + 
  coord_equal()+
  theme(legend.position = 'none', panel.grid = element_line(linetype = 3, size = .5), panel.grid.minor = element_blank())+
  scale_color_tableau(direction = -1)+ 
  # scale_color_viridis_d()+
  labs(x= "False positive rate", y= "True positive rate")+ 
  geom_roc(data = full.hotspot.data[cv.number=="cv.9",], 
           aes(m= vote.1, d= Observed), size= 1.2, 
           cutoffs.at = c(.145), # this is the discrimination threshold 
           cutoff.labels = "FPR= 0.40,     \nTPR= 0.85      ", # from sensetivity(.145) and 1- specificity(.145)
           pointsize = 1, labelsize = 2)


save.image("Jan312022.RData")



################################################################################
## hotspot models for countries w least 1 admin_2 with incidence > 1 per 1000 ##
################################################################################
# get the country vector 
country_names<- unique(dist_dt[incidence_in_thousan>1, NAME_0])
hotspot_df2<- dist_dt[NAME_0 %in% country_names, ]

##### -------------------------- cross-validation ------------------------ #####
adm0<- hotspot_df2[,.N,by= NAME_0]
folds.h<- rep(1:nrow(adm0), adm0$N)

hotspot_df2<- hotspot_df2[,c(10:12, 14:16, 20)]
hotspot_df2$hotspot<- as.factor(hotspot_df2$hotspot)

test.list.h2<- list()
train.list.h2<- list()
rf.list.h2<- list()


for(i in 1:nrow(adm0)){
  testIndexes <- which(folds.h==i,arr.ind=TRUE)
  testData <- hotspot_df2[testIndexes, ]
  trainData <- hotspot_df2[-testIndexes, ]
  trainData<- ovun.sample(hotspot~., data=trainData,
                          N=nrow(trainData),  p= 0.5,
                          seed=1, method="both")$data
  test.list.h2[[i]]<- testData
  train.list.h2[[i]]<- trainData
  rf<- randomForest(
    hotspot~.,
    data = trainData,
    ntree=500,
    localImp= TRUE, 
    keep.forest= TRUE, 
    keep.inbag= TRUE
  )
  rf.list.h2[[i]]<- rf
  rf.list.h2[[i]]$pred.newdata<- predict(rf, newdata = testData)
  rf.list.h2[[i]]$pred.newdata.vote<- predict(rf, newdata = testData, type = "vote",norm.votes=TRUE)
  rf.list.h2[[i]]$testdata<- testData
  rf.list.h2[[i]]$conf.matrix<- table(testData$hotspot ,rf.list.h2[[i]]$pred.newdata)
}

saveRDS(rf.list.h2, "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf.list.h2.RDS")

### Predict using CV models
hotspot_predict2<- data.table()
full.hotspot.data2<- data.table()
for (i in 1:length(rf.list.h2)){
  hotspot_predict2<- data.table() # since the lengths are different
  temp<- hotspot_predict2[ , c("cv.number", "Predicted", "Observed", "vote.0", "vote.1") := 
                            list(rep(paste0("cv.", i), length(rf.list.h2[[i]]$pred.newdata)), 
                                 rf.list.h2[[i]]$pred.newdata, 
                                 rf.list.h2[[i]]$testdata$hotspot, 
                                 rf.list.h2[[i]]$pred.newdata.vote[,1], 
                                 rf.list.h2[[i]]$pred.newdata.vote[,2]),]
  full.hotspot.data2<- rbind(full.hotspot.data2, temp)
}

full.hotspot.data2$Predicted<- as.numeric(as.character(full.hotspot.data2$Predicted))
full.hotspot.data2$Observed<- as.numeric(as.character(full.hotspot.data2$Observed))

saveRDS(full.hotspot.data2, "X:/Spatial Stat/WASH Cholera/clean_repo/results/full.hotspot.data2.RDS")
full.hotspot.data2<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/full.hotspot.data2.RDS")


############ --------------- hotspot full model ------------------##############
### Note: these results are used to produce ROC and variable importance plots


rf_hp_full2<- randomForest(
  hotspot~.,
  data = hotspot_df2,#_ovun,
  ntree=500,
  localImp= TRUE, 
  keep.forest= TRUE, 
  keep.inbag= TRUE
)


rf_hp_full.df2<- data.frame(Predicted.f= rf_hp_full2$predicted, Observed.f= rf_hp_full2$y, vote.0.f= rf_hp_full2$votes[,1], vote.1.f= rf_hp_full2$votes[,2])
saveRDS(rf_hp_full.df2, file = "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_hp_full.df2.RDS")
rf_hp_full.df2<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_hp_full.df2.RDS")



#### combining the full model and cv model  data frames for hte ROC plot
combined_df_hp2<- data.frame(full.hotspot.data2, rf_hp_full.df2)
combined_df_hp2$Observed<- as.numeric(as.character(combined_df_hp2$Observed))
combined_df_hp2$Observed.f<- as.numeric(as.character(combined_df_hp2$Observed.f))
combined_df_hp2$Predicted<- as.numeric(as.character(combined_df_hp2$Predicted))
combined_df_hp2$Predicted.f<- as.numeric(as.character(combined_df_hp2$Predicted.f))

saveRDS(combined_df_hp2, "X:/Spatial Stat/WASH Cholera/clean_repo/results/combined_df_hp2.RDS")
combined_df_hp2<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/combined_df_hp2.RDS")
### to make sure combined_df_hp$Observed!=combined_df_hp$Observed.f to use in the main2 plot 
sum(combined_df_hp2$Observed!=combined_df_hp2$Observed.f)


# variable importance and 95% CI

## Will check later  
imp.h_full2<- permimp::permimp(rf_hp_full2, conditional = T, do_check = FALSE)
saveRDS(imp.h_full2, "X:/Spatial Stat/WASH Cholera/clean_repo/results/imp.h_full2.RDS")
# imp.h_full2<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/imp.h_full2.RDS")
# 
# ## AUC  for full data model
# auc(rf_hp_full$y , rf_hp_full$votes[,1], quiet= T)
# 
# #95% CI of AUC
# ci.auc(auc(rf_hp_full$y , rf_hp_full$votes[,1], quiet= T))


# v_imp_hp<- cbind(imp.h$values) commented since we are using full data model now 12/28/2021
v_imp_hp2<- cbind(imp.h_full2$values)
#rownames(v_imp_hp)<- c("Water Improved", "Water Piped", "Water Surface", "Water Unimproved", "Sanitation Improved", "Open Defecation", "Sanitation Piped", "Sanitation Unimproved")
rownames(v_imp_hp2)<- c("Water Piped & Improved", "Water Piped", "Water Surface", "Sanitation Piped & Improved", "Open Defecation", "Sanitation Piped")
v_imp_hp2<- cbind(rownames(v_imp_hp2), v_imp_hp2)
v_imp_hp2<- as.data.frame(v_imp_hp2)
names(v_imp_hp2)<- c("var", 'mean_accuracy_decrease')
v_imp_hp2$mean_accuracy_decrease<- as.numeric(v_imp_hp2$mean_accuracy_decrease)
saveRDS(v_imp_hp2, "X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_hp2.RDS")
v_imp_hp2<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_hp2.RDS")




# ---------------------

main22<- ggplot(data=combined_df_hp2)+ 
  geom_roc(aes(m= vote.1, d= Observed, color= cv.number), n.cuts=0, linealpha = .3) + 
  geom_roc(aes(m= vote.1.f, d= Observed), n.cuts=0, linealpha = 1) + 
  coord_equal()+
  theme(legend.position = 'none')+
  scale_color_grey(end= 0) +
  #scale_color_tableau()+ 
  # scale_color_viridis_d()+
  labs(x= "False positive rate", y= "True positive rate")


mean_impurity_decrease.hp2 <-
  ggplot(data = v_imp_hp2,
         aes(
           mean_accuracy_decrease,
           reorder(var, mean_accuracy_decrease),
           fill = var
         )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c(
      "#264653",
      "#2a9d8f",
      "#F16745",
      "#FFC65D",
      "#f4a261",
      "#4CC3D9",
      "#93648D",
      "#457b9d"
    )
  ) + labs(y = NULL, x = "Conditional permutation importance") + theme(
    legend.title = element_blank(),
    legend.position = 'none',
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size= 8),
    axis.text.x = element_text(size = 6),
    #plot.background = element_rect(color = "White"),
    panel.grid.major = element_blank(),
    panel.grid.minor =  element_blank()
  )+ 
  scale_x_continuous(position = "top")



###############################################################################
##### ------------------- District Fold Cross-Validation -----------------######
###############################################################################
hotspot_df<- dist_df[,c(4:9, 11)]
dist_dt<- as.data.table(dist_df)
adm0<- dist_dt[,.N,by= NAME_0.x]
folds.h<- rep(1:nrow(adm0), adm0$N)
hotspot_df$hotspot<- as.factor(hotspot_df$hotspot)
saveRDS(hotspot_df, 'results/hotspot_df.RDS')

############################# Random Forest CV Model #############################
test.list.h<- list()
train.list.h<- list()
rf.list.h.dist<- list()
for(i in 1:nrow(hotspot_df)){
  print(i)
  testIndexes <- i
  testData <- hotspot_df[testIndexes, ]
  trainData <- hotspot_df[-testIndexes, ]
  trainData<- ovun.sample(hotspot~., data=trainData,
                          N=nrow(trainData),  p= 0.5,
                          seed=1, method="both")$data
  test.list.h[[i]]<- testData
  train.list.h[[i]]<- trainData
  rf<- randomForest(
    hotspot~.,
    data = trainData,
    ntree=500,
    localImp= TRUE, 
    keep.forest= TRUE, 
    keep.inbag= TRUE
  )
  rf.list.h.dist[[i]]<- i
  # rf.list.h.dist[[i]]$pred.newdata<- predict(rf, newdata = testData)
  rf.list.h.dist[[i]]$pred.newdata.vote<- predict(rf, newdata = testData, type = "vote",norm.votes=TRUE)[, 2]
  rf.list.h.dist[[i]]$testdata<- testData$hotspot
  # rf.list.h.dist[[i]]$conf.matrix<- table(testData$hotspot ,rf.list.h.dist[[i]]$pred.newdata)
}

saveRDS(rf.list.h.dist, "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf.list.h_dist_6vars.RDS")
rf.list.h.dist<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/rf.list.h_dist_6vars.RDS")

### Bootstrap for cvACU and 95% CI
actual<- data.table()
predicted<- data.table()
for (i in 1:length(rf.list.h.dist)){
  actual<- rbind(actual, rf.list.h.dist[[i]]$testdata)
  predicted<- rbind(predicted, rf.list.h.dist[[i]]$pred.newdata.vote)
}
auc(as.factor(actual$x), predicted$x)

# Bootstrap for cvAUC for 95% CI
auc_data<- data.frame(actual= as.factor(actual$x), predicted= predicted$x)
fc<- function(d, i){
  d2<- d[i,]
  return(auc(d2$actual, d2$predicted))
}
cvAUC<- boot(auc_data, statistic=fc, R=5000)
rf_cvAUC_cl<- boot.ci(cvAUC, conf=0.95, type="bca")

### Predict using CV models
hotspot_predict2_dis<- data.table()
full.hotspot.data2_dis<- data.table()
for (i in 1:length(rf.list.h.dist)){
  hotspot_predict2_dis<- data.table() # since the lengths are different
  temp<- hotspot_predict2_dis[ , c("cv.number", "Observed", "vote.1") := 
                             list(paste0("cv.", i), 
                                  # rf.list.h.dist[[i]]$pred.newdata, 
                                  rf.list.h.dist[[i]]$testdata, 
                                  rf.list.h.dist[[i]]$pred.newdata.vote), ]
  full.hotspot.data2_dis<- rbind(full.hotspot.data2_dis, temp)
}

# full.hotspot.data2_dis$Predicted<- as.numeric(as.character(full.hotspot.data2_dis$Predicted))
full.hotspot.data2_dis$Observed<- as.numeric(as.character(full.hotspot.data2_dis$Observed))

saveRDS(full.hotspot.data2_dis, "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf.hotspot.cv_leave_one_dist.RDS")
full.hotspot.data2_dis<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/rf.hotspot.cv_leave_one_dist.RDS")

### combine full data model results with CV for ROC plot 
combined_df_hp_dist<- data.frame(full.hotspot.data2_dis, rf_hp_full.df)
combined_df_hp_dist$Observed<- as.numeric(as.character(combined_df_hp_dist$Observed))
combined_df_hp_dist$Observed.f<- as.numeric(as.character(combined_df_hp_dist$Observed.f))
# combined_df_hp_dist$Predicted<- as.numeric(as.character(combined_df_hp_dist$Predicted))
combined_df_hp_dist$Predicted.f<- as.numeric(as.character(combined_df_hp_dist$Predicted.f))

saveRDS(combined_df_hp_dist, "X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_combined_df_hp_dist_6vars.RDS")
combined_df_hp_dist<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_combined_df_hp_dist_6vars.RDS")


################################## Hotspot ROC AUC Plot with leave one district out (Fig: 3B) ##########################

rf_auc_plot_dist<- ggplot(data=combined_df_hp_dist)+ 
  # geom_roc(aes(m= vote.1, d= Observed, color= cv.number), n.cuts=0, linealpha = .3, size= .2) + 
  geom_roc(aes(m= vote.1, d= Observed), n.cuts=0, linealpha = 1) + # using CV model data 
  coord_equal()+
  theme(legend.position = 'none')+
  scale_color_grey(end= 0) +
  #scale_color_tableau()+ 
  # scale_color_viridis_d()+
  geom_point(aes(x= 0.26036, y= 0.71521 ), color= 'red', size= 2)+ # x, y from roc_cutpoint.R script; max Youden's index
  labs(x= "False positive rate", y= "True positive rate") +
  annotate("text", x= .28, y= .72, 
           label= sprintf("(FPR = 0.26, TPR = 0.72)"), hjust= 0, size = 3)

saveRDS(rf_auc_plot_dist, paste0(output_dir, 'rf_hotspot_auc_6vr_plot_distCV.RDS'))
# saveRDS(rf_auc_plot, paste0(output_dir, 'rf_hotspot_auc_6vr_plot_light.RDS'))
# saveRDS(rf_auc_plot, paste0(output_dir, 'rf_hotspot_auc_6vr_plot.RDS'))
png('X:/Spatial Stat/WASH Cholera/clean_repo/results/rf_hotspot_roc_6vars.png', 
    width = 6, height = 6, res = 400, units = 'in')
rf_auc_plot
dev.off()






#################################### GBM CV Model leave one District out ###################################
hotspot_df_gbm<- hotspot_df
hotspot_df_gbm$hotspot<- as.numeric(hotspot_df_gbm$hotspot)-1 # as GBM wants non factor for two categories 

test.list.h<- list()
train.list.h<- list()
gbm.list.h.dist<- list()
for(i in 1:nrow(hotspot_df_gbm)){
  testIndexes <- i
  testData <- hotspot_df_gbm[testIndexes, ]
  trainData <- hotspot_df_gbm[-testIndexes, ]
  trainData<- ovun.sample(hotspot~., data=trainData,
                          N=nrow(trainData),  p= 0.5,
                          seed=1, method="both")$data
  test.list.h[[i]]<- testData
  train.list.h[[i]]<- trainData
  gbm.hp<- gbm(
    hotspot~.,
    data = trainData,
    distribution="bernoulli",  
    n.trees=500, 
    verbose=T
  )
  gbm.list.h.dist[[i]]<- i
  # gbm.list.h.dist[[i]]$pred.newdata<- predict(gbm.hp, newdata = testData)
  gbm.list.h.dist[[i]]$pred.newdata.vote<- predict(gbm.hp, newdata = testData, type= "response")
  gbm.list.h.dist[[i]]$testdata<- testData$hotspot
  # gbm.list.h.dist[[i]]$conf.matrix<- table(testData$hotspot ,gbm.list.h[[i]]$pred.newdata)
}

############### district as folds for cvRMSE - hotspot model ####################
actual<- data.table()
predicted<- data.table()
for (i in 1:length(gbm.list.h.dist)){
  actual<- rbind(actual, gbm.list.h.dist[[i]]$testdata)
  predicted<- rbind(predicted, gbm.list.h.dist[[i]]$pred.newdata.vote)
}
auc(as.factor(actual$x), predicted$x)

# Bootstrap for cvAUC for 95% CI
auc_data<- data.frame(actual= as.factor(actual$x), predicted= predicted$x)
fc<- function(d, i){
  d2<- d[i,]
  return(auc(d2$actual, d2$predicted))
}
cvAUC_gbm<- boot(auc_data, statistic=fc, R=5000)
gbm_cvAUC_cl<- boot.ci(cvAUC_gbm, conf=0.95, type="bca")


