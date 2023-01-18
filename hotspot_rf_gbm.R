# This script runs random forest and GBM hotspot model using district fold cross validation 
# and full data. 
#


rm(list= ls())
output_dir<- 'DIRECTORY_PATH'

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

hotspot_df<- dist_df[,c(4:9, 11)]
hotspot_df$hotspot<- as.factor(hotspot_df$hotspot)

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
    panel.grid.major = element_blank(),
    panel.grid.minor =  element_blank(), 
    panel.grid.minor.x = element_blank())+ 
  scale_x_continuous(position = "top") #, limits = c(min(v_imp_hp$mean_accuracy_decrease), max(v_imp_hp$mean_accuracy_decrease)), breaks = c(as.vector(summary(v_imp_hp$mean_accuracy_decrease)[c(1,3,6)])))



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

imp.h_full<- permimp::permimp(rf_hp_full, conditional = T, do_check = FALSE)

## AUC  for full data model
auc(rf_hp_full$y , rf_hp_full$votes[,1], quiet= T)
#95% CI of AUC
ci.auc(auc(rf_hp_full$y , rf_hp_full$votes[,1], quiet= T))

################################## Hotspot ROC AUC Plot (Fig: 3B) ##########################
#### Using RF AUC plot in the main text as RF cvAUC> GBM cvAUC
rf_auc_plot<- ggplot(data=combined_df_hp)+ 
  geom_roc(aes(m= vote.1, d= Observed), n.cuts=0, linealpha = 1) + # using CV model data 
  coord_equal()+
  theme(legend.position = 'none')+
  scale_color_grey(end= 0) +
  geom_point(aes(x= 0.447, y= 0.796), color= 'red', size= 2)+ # x, y from roc_cutpoint.R; max Youden's index
  labs(x= "False positive rate", y= "True positive rate")

saveRDS(rf_auc_plot, paste0(output_dir, 'rf_hotspot_auc_6vr_plot_light.RDS'))

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
    panel.grid.major = element_blank(),
    panel.grid.minor =  element_blank()
  )+ 
  scale_x_continuous(position = "top")


################### FIGURE SI Hotspot Variable importance plot: RF + GBM ###################
hp_imp_plot<- mean_impurity_decrease.hp + mean_impurity_decrease.hp.gbm

hp_var_imp_plot<- hp_imp_plot+ plot_annotation(tag_levels = 'A') & #, tag_prefix = '3'
  theme(plot.tag.position = c(.4, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

png(paste0(output_dir, 'plots/fig_si_hotspot_var_imp_6vs.png'), 
    width = 8, height = 3, res = 400, units = 'in')
hp_var_imp_plot
dev.off()


###############################################################################
##### ------------------- District Fold Cross-Validation -----------------######
###############################################################################
hotspot_df<- dist_df[,c(4:9, 11)]
hotspot_df$hotspot<- as.factor(hotspot_df$hotspot)

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
  rf.list.h.dist[[i]]$pred.newdata.vote<- predict(rf, newdata = testData, type = "vote",norm.votes=TRUE)[, 2]
  rf.list.h.dist[[i]]$testdata<- testData$hotspot
}

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

full.hotspot.data2_dis$Observed<- as.numeric(as.character(full.hotspot.data2_dis$Observed))

### combine full data model results with CV for ROC plot 
combined_df_hp_dist<- data.frame(full.hotspot.data2_dis, rf_hp_full.df)
combined_df_hp_dist$Observed<- as.numeric(as.character(combined_df_hp_dist$Observed))
combined_df_hp_dist$Observed.f<- as.numeric(as.character(combined_df_hp_dist$Observed.f))
combined_df_hp_dist$Predicted.f<- as.numeric(as.character(combined_df_hp_dist$Predicted.f))

################################## Hotspot ROC AUC Plot  ##########################

rf_auc_plot_dist<- ggplot(data=combined_df_hp_dist)+ 
  geom_roc(aes(m= vote.1, d= Observed), n.cuts=0, linealpha = 1) + # using CV model data 
  coord_equal()+
  theme(legend.position = 'none')+
  scale_color_grey(end= 0) +
  geom_point(aes(x= 0.26036, y= 0.71521 ), color= 'red', size= 2)+ # x, y from roc_cutpoint.R script; max Youden's index
  labs(x= "False positive rate", y= "True positive rate") +
  annotate("text", x= .28, y= .72, 
           label= sprintf("(FPR = 0.26, TPR = 0.72)"), hjust= 0, size = 3)

saveRDS(rf_auc_plot_dist, paste0(output_dir, 'rf_hotspot_auc_6vr_plot_distCV.RDS'))
png('rf_hotspot_roc_6vars.png', 
    width = 6, height = 6, res = 400, units = 'in')
rf_auc_plot
dev.off()

#################################### GBM CV Model leave one District out ###################################
hotspot_df_gbm<- hotspot_df
hotspot_df_gbm$hotspot<- as.numeric(hotspot_df_gbm$hotspot)-1 # GBM expects non-factor for two categories 

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
  gbm.list.h.dist[[i]]$pred.newdata.vote<- predict(gbm.hp, newdata = testData, type= "response")
  gbm.list.h.dist[[i]]$testdata<- testData$hotspot
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

