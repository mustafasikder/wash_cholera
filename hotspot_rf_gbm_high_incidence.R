rm(list= ls())

library(randomForest)
library(gbm)
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
library(patchwork)

theme_set(theme_minimal())
theme_update(legend.position = "bottom")
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")
options(digits=2)

# Read and subset the data 
dist_df<- readRDS("results/dist_df.rds")
dist_df<- dist_df[complete.cases(dist_df), ]
dist_dt<- as.data.table(dist_df)
dist_dt[hotspot==1, .N, by= NAME_0]
# subset countries with at least 1 adm2 with incidence > 1 per 1000
subset_countries<- unique(dist_dt$NAME_0[dist_dt$hotspot==1])
dist_dt_sub<- dist_dt[NAME_0 %in% c(subset_countries), ]
hotspot_df<- dist_dt_sub[,c(10:12, 14:16, 20)]

###############################################################################
##### ------------------- District Fold Cross-Validation -----------------######
###############################################################################

############################### Random Forest district fold #################################
hotspot_df$hotspot<- as.factor(hotspot_df$hotspot)

test.list.h<- list()
train.list.h<- list()
rf.list.h<- list()

for(i in 1:nrow(hotspot_df)){
  testIndexes <- i
  print(i)
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
  
  rf.list.h[[i]]<- i
  rf.list.h[[i]]$pred.newdata.vote1<- predict(rf, newdata = testData, type = "vote",norm.votes=TRUE)[, 2]
  rf.list.h[[i]]$testdata<- testData$hotspot
}

saveRDS(rf.list.h, "rf.list.h_6vars_1caseAdm2_distCV.RDS")

# Data prep for ROC plot
hotspot_predict<- data.table()
full.hotspot.data<- data.table()
for (i in 1:length(rf.list.h)){
  temp<- hotspot_predict[ , c("cv.number", "vote.1", "Observed") := 
                            .(rf.list.h[[i]][[1]], 
                                 rf.list.h[[i]]$pred.newdata.vote1, 
                                 rf.list.h[[i]]$testdata)]
  full.hotspot.data<- rbind(full.hotspot.data, temp)
}

full.hotspot.data$Observed<- as.numeric(as.character(full.hotspot.data$Observed))

saveRDS(full.hotspot.data, "full.hotspot.data_6vars_1caseAdm2_distCV.RDS")


### Bootstrap for cvACU and 95% CI
actual<- data.table()
predicted<- data.table()
for (i in 1:nrow(hotspot_df)){
  actual<- rbind(actual, rf.list.h[[i]]$testdata)
  predicted<- rbind(predicted, rf.list.h[[i]]$pred.newdata.vote1)
}
auc(as.factor(actual$x), predicted$x)

# Bootstrap for cvAUC for 95% CI
auc_data_rf<- data.frame(actual= as.factor(actual$x), predicted= predicted$x)
fc<- function(d, i){
  d2<- d[i,]
  return(auc(d2$actual, d2$predicted))
}
cvAUC_sub_rf<- boot(auc_data_rf, statistic=fc, R=5000)
rf_cvAUC_cl_sub<- boot.ci(cvAUC_sub_rf, conf=0.95, type="bca")




#################### GBM cv model with subset data district fold ##################################
hotspot_df_gbm<- hotspot_df
hotspot_df_gbm$hotspot<- as.numeric(hotspot_df_gbm$hotspot)-1 # as GBM wants non factor for two categories 

test.list.h<- list()
train.list.h<- list()
gbm.list.h<- list()
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
  gbm.list.h[[i]]<- i
  gbm.list.h[[i]]$pred.newdata<- predict(gbm.hp, newdata = testData)
  gbm.list.h[[i]]$pred.newdata.vote<- predict(gbm.hp, newdata = testData, type= "response")
  gbm.list.h[[i]]$testdata<- testData$hotspot
}


saveRDS(gbm.list.h, "gbm.hotspot.data_6vars_1caseAdm2_distCV.RDS")


############### cvRMSE - hotspot model ####################
actual_gbm<- data.table()
predicted_gbm<- data.table()
for (i in 1:length(gbm.list.h)){
  actual_gbm<- rbind(actual, gbm.list.h[[i]]$testdata)
  predicted_gbm<- rbind(predicted, gbm.list.h[[i]]$pred.newdata.vote)
}
auc(as.factor(actual_gbm$x), predicted_gbm$x)

# Bootstrap for cvAUC for 95% CI
auc_data_gbm<- data.frame(actual= as.factor(actual_gbm$x), predicted= predicted_gbm$x)
fc<- function(d, i){
  d2<- d[i,]
  return(auc(d2$actual, d2$predicted))
}
cvAUC_sub_gbm<- boot(auc_data_gbm, statistic=fc, R=5000)
gbm_cvAUC_sub_cl<- boot.ci(cvAUC_sub_gbm, conf=0.95, type="bca")



# --------------- hotspot full model ------------------

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

#### combining the full model and cv model  data frames for hte ROC plot
combined_df_hp<- data.frame(full.hotspot.data, rf_hp_full.df)
combined_df_hp$Observed<- as.numeric(as.character(combined_df_hp$Observed))
combined_df_hp$Predicted.f<- as.numeric(as.character(combined_df_hp$Predicted.f))

# variable importance 
v_imp_hp<- cbind(imp.h_full$values)
rownames(v_imp_hp)<- c("Piped or other improved water", 
                       "Piped water",
                       "Surface water",
                       "Piped or other improved sanitation", 
                       "Open defecation", 
                       "Piped sanitation")
v_imp_hp<- cbind(rownames(v_imp_hp), v_imp_hp)
v_imp_hp<- as.data.frame(v_imp_hp)
names(v_imp_hp)<- c("var", 'mean_accuracy_decrease')
v_imp_hp$mean_accuracy_decrease<- as.numeric(v_imp_hp$mean_accuracy_decrease)
saveRDS(v_imp_hp, "X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_hp_6vars_1caseAdm2.RDS")
v_imp_hp<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_hp_6vars_1caseAdm2.RDS")

main2_sub<- ggplot(data=combined_df_hp)+ 
  # geom_roc(aes(m= vote.1, d= Observed, color= cv.number), n.cuts=0, linealpha = .3) + 
  geom_roc(aes(m= vote.1, d= Observed), n.cuts=0, linealpha = 1) + 
  coord_equal()+
  theme(legend.position = 'none')+
  scale_color_grey(end= 0) +
  labs(x= "False positive rate", y= "True positive rate")

png('hotspot_roc_6vars_1caseAdm2.png', 
    width = 6, height = 6, res = 400, units = 'in')
main2_sub
dev.off()

mean_impurity_decrease.hp_sub <-
  ggplot(data = v_imp_hp,
         aes(
           mean_accuracy_decrease,
           reorder(var, mean_accuracy_decrease),
           fill = var
         )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values = c('#2d588a', '#58954c', '#e9a044', '#c12f32', '#723e77', '#7d807f')
  ) + labs(y = NULL, x = "Conditional permutation importance") + theme(
    legend.title = element_blank(),
    legend.position = 'none',
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size= 8),
    axis.text.x = element_text(size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor =  element_blank()
  )+ 
  scale_x_continuous(position = "top")

png('hotspot_vimp_6vars_1caseAdm2.png', 
    width = 5, height = 3, res = 300, units = 'in')
mean_impurity_decrease.hp_sub
dev.off()


# Figure SI 

sub_plot<- main2_sub+ mean_impurity_decrease.hp_sub
fig_SI3<- sub_plot+ plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(.3, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

png(paste0(output_dir, 'plots/fig_SI_3_sub_auc_Vimp_6vs.png'), 
    width = 6, height = 3.5, res = 400, units = 'in')
fig_SI3
dev.off()




