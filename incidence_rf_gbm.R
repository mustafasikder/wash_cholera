# Incidence random forest (RF) and GBM model using (1) district as cross-validation fold
# and (2) the full data 
# 
#
# Paper input from cross-validation model:  1. cvRMSE (95% CI RMSE)
#                                           2. Figure 3A scatter plot (each fold predicts the hold out data)
#                                           3. Spearman's rho (compare observed vs hold out prediction)
#                
# Paper input from full data model:         3. SI Figure variable importance plot
#                                            


rm(list= ls())
output_dir<- 'DIRECTORY_PATH'


library(randomForest)
library(gbm)
library(ggplot2)
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
library(scales)
library(doParallel)
library(patchwork)

theme_set(theme_minimal())
theme_update(legend.position = "bottom")
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")
options(digits=5)

#### read district data 
dist_df<- read.csv(paste0(output_dir, 'district_6v.csv'))
dist_df$log_incidence<- log10(dist_df$incidence_in_thousan)
dist_df<- dist_df[complete.cases(dist_df), ]
cv_df<- dist_df[,c(4:9, 12)]
###############################################################################
##### ------------------- District Fold Cross-Validation -----------------######
###############################################################################

############################# Random Forest CV Model #############################
# set up and register cluster
cl = makeCluster(8)
registerDoParallel(cl)

test.list<- list()
train.list<- list()
rf.list.dist<- data.frame()

system.time(
  rf.list.dist<- foreach(i = 1:nrow(cv_df), .packages='randomForest', .combine= rbind) %dopar%{
    print(i)
    testIndexes <- i
    testData <- cv_df[testIndexes, ]
    trainData <- cv_df[-testIndexes, ]
    rf<- randomForest(
      log_incidence~.,
      data = trainData,
      ntree=500
    )
    rf.list.dist[i, c('cv', 'pred.incidence', 'actual.incidence')]<- c(i, predict(rf, newdata = testData), testData$log_incidence)
  }
)

stopCluster(cl)

saveRDS(rf.list.dist, paste0(output_dir, 'incidence_rf_model_w_district_cv_6vr.RDS'))
# rf.list.dist <- readRDS(paste0(output_dir, 'incidence_rf_model_w_district_cv_6vr.RDS'))

# convert to dataframe
rf.cv.dist<- data.frame(rf.list.dist)
colnames(rf.cv.dist)<- c('fold', 'pred.incidence', 'actual.incidence')
rownames(rf.cv.dist)<- NULL

# Spearman's rho 
rho<- round(cor.test(rf.cv.dist$actual.incidence, rf.cv.dist$pred.incidence, method = "spearman")$estimate, 2)
# RMSE
rmse<- round((sum((rf.cv.dist$actual.incidence- rf.cv.dist$pred.incidence)^2)/nrow(rf.cv.dist))^.5, 2)

################################## Incidence Scatter Plot ##########################
# Scatter plot (Figure 3A) of district CV results observed vs predicted
scatter_plot_cv.dist<- ggplot(data= rf.cv.dist, aes(pred.incidence, actual.incidence))+ 
  geom_point(size= .5, alpha= .5, color= 'grey40')+ 
  labs(x= 'Predicted', y= 'Observed')+ 
  theme(legend.position = "none", panel.grid = element_line(linetype = 1, size = .5))+
  scale_x_continuous(limits = c(-5, 3), breaks = c(-4, -2, 0, 2),labels = c("0.0001", "0.01", "0.0", "100"), position = "bottom")+
  scale_y_continuous(limits = c(-5, 2), breaks = c(-4, -2, 0, 2),labels = c("0.0001", "0.01", "0.0", "100"), position = "left")+
  annotate("text", x= 1, y= -1, 
           label= sprintf("atop (italic(rho) == %1.2f, 
                          RMSE == %f)", 
                          rho, rmse), parse = TRUE, size = 3)+ 
  coord_equal()
saveRDS(scatter_plot_cv.dist, paste0(output_dir, 'incidence_rf_scatter_6vr_distCV.RDS'))

#### --------- k fold cross validation RMSE for RF model --------------#####
# Bootstrap for cvRMSE with 95% CI
cvRMSE<- data.frame()
for (i in 1:length(rf.list.dist)) {
  cvRMSE<- rbind(cvRMSE, cbind(m_id= paste0("model", i), rmse= rf.list.dist[[i]]$rmse))
}
cvRMSE$rmse<- as.numeric(cvRMSE$rmse)
fc<- function(d, i){
  d2<- d[i, ]
  return((mean((d2$rmse)^2))^.5)
}
#95% CI of cvRMSE
rf_cvRMSE_CI<- boot::boot.ci(boot(cvRMSE, statistic=fc, R=5000), conf=0.95, type="bca")

#################################### GBM CV Model ###################################

test.list<- list()
train.list<- list()
gbm.list.dist<- list()
for(i in 1:nrow(cv_df)){
  print(i)
  testIndexes <- i
  testData <- cv_df[testIndexes, ]
  trainData <- cv_df[-testIndexes, ]
  gbm.i<- gbm(
    log_incidence~.,
    data = trainData,
    distribution="gaussian",
    n.trees=500
  )
  gbm.list.dist[[i]]<- i
  gbm.list.dist[[i]]$rmse<- mean((predict(gbm.i, newdata = testData)-testData$log_incidence)^2)^0.5
}
saveRDS(gbm.list.dist, paste0(output_dir, 'incidence_gbm_model_w_district_cv_6vr.RDS'))
gbm.list.dist <- readRDS(paste0(output_dir, 'incidence_gbm_model_w_district_cv_6vr.RDS'))

#### --------- k fold cross validation RMSE for GBM model --------------#####
(sum(unlist(lapply(gbm.list.dist, function(x) (x$rmse)^2)))/length(gbm.list.dist))^.5
# Bootstrap for cvRMSE with 95% CI
cvRMSE_GBM<- data.frame()
for (i in 1:length(gbm.list.dist)) {
  cvRMSE_GBM<- rbind(cvRMSE_GBM, cbind(m_id= paste0("model", i), rmse= gbm.list.dist[[i]]$rmse))
}
cvRMSE_GBM$rmse<- as.numeric(cvRMSE_GBM$rmse)
gbm_cvRMSE_CI<- boot::boot.ci(boot(cvRMSE_GBM, statistic=fc, R=5000), conf=0.95, type="bca")


###############################################################################
##### ---------------------------- Full data model  ---------------------######
###############################################################################

########################## Random Forest full data Model #############################
incidence_df<- dist_df[,c(4:9, 12)]
set.seed(123)

rf_full<- randomForest(
  log_incidence~.,
  data = incidence_df,
  ntree=500,
  localImp= TRUE, 
  keep.forest= TRUE, 
  keep.inbag= TRUE
)


imp_full<- permimp::permimp(rf_full, conditional = T, do_check = FALSE)
saveRDS(imp_full, paste0(output_dir, 'permutationImpFullData.RDS'))
#imp_full<- readRDS("results/permutationImpFullData.RDS")

v_imp_inc<- data.frame( var= c("Piped water", 
                          "Piped or other improved water", 
                          "Surface water",
                          "Piped sanitation", 
                          "Piped or other improved sanitation", 
                          "Open defecation"), 
                   mean_accuracy_decrease= imp_full$values)
v_imp_inc$mean_accuracy_decrease<- as.numeric(v_imp_inc$mean_accuracy_decrease)

saveRDS(v_imp_inc, "X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_inc_6vars.RDS")
v_imp_inc<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_inc_6vars.RDS")

#### RF variable importance plot for SI 
mean_impurity_decrease <-
  ggplot(data = v_imp_inc,
         aes(
           mean_accuracy_decrease,
           reorder(var, mean_accuracy_decrease),
           fill = var
         )) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('#2d588a', '#58954c', '#e9a044', '#c12f32', '#723e77', '#7d807f')) +
  labs(y = NULL, x = "Conditional variable importance (RF)") + theme(
    legend.title = element_blank(),
    legend.position = 'none',
    axis.text.y = element_text(size = 8), 
    axis.title.x = element_text(size= 8),
    axis.text.x = element_text(size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor =  element_blank()
  )+ 
  scale_x_continuous(position = "top")

saveRDS(mean_impurity_decrease, paste0(output_dir, 'incidence_rf_vimp_6vars.RDS'))
png('X:/Spatial Stat/WASH Cholera/clean_repo/results/incidence_rf_vimp_6vars.png', 
    width = 5, height = 3, res = 300, units = 'in')
mean_impurity_decrease
dev.off()


########################## GBM full data Model #############################
gbm.inc<- gbm(log_incidence~.,  
              data = incidence_df, 
              distribution="gaussian",  
              n.trees=500, 
              verbose=T)

######################## full model importance plot  ########################## 
vip_full<- vip(gbm.inc)
v_imp_inc_gbm<- data.table(vip_full$data)
setorder(v_imp_inc_gbm, Variable)
v_imp_inc_gbm$var<- c("Piped or other improved sanitation",
                      "Open defecation", 
                      "Piped sanitation",
                      "Piped or other improved water", 
                      "Piped water",
                      "Surface water"
                      )

saveRDS(v_imp_inc_gbm, file = "X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_inc_gbm_6vars.RDS")
v_imp_inc_gbm<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/v_imp_inc_gbm_6vars.RDS")

#### GBM variable importance plot for SI 
mean_impurity_decrease_gbm <-
  ggplot(data = v_imp_inc_gbm,aes( Importance, reorder(var, Importance), fill = var)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('#2d588a', '#58954c', '#e9a044', '#c12f32', '#723e77', '#7d807f')) +
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
  scale_x_continuous(position = "top")#, limits = c(min(v_imp_inc$mean_accuracy_decrease), max(v_imp_inc$mean_accuracy_decrease)), breaks = c(as.vector(summary(v_imp_inc$mean_accuracy_decrease)[c(1,3,6)])))



################### FIGURE SI Incidence Variable importance plot: RF + GBM ###################
incidence_imp_plot<- mean_impurity_decrease + mean_impurity_decrease_gbm

incidence_var_imp_plot<- incidence_imp_plot+ plot_annotation(tag_levels = 'A') & #, tag_prefix = '3'
  theme(plot.tag.position = c(.4, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

png(paste0(output_dir, 'plots/fig_si_incidence_var_imp_6vs.png'), 
    width = 8, height = 3, res = 400, units = 'in')
incidence_var_imp_plot
dev.off()

