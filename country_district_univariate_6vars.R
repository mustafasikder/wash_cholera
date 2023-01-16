# This script runs country level Quasi-Poisson regression and district level 
# Poisson generalized estimating equations (GEE) using the six variable and 
# saves results and produces Figure 2.  
#
# If only editing Figure 2, go to the end
#
# Paper output: (1) Figure 2
#               (2) Table S2
#               (3) Risk ratios for main text 



rm(list= ls())
output_dir<- 'X:/Spatial Stat/WASH Cholera/clean_repo/results/'

# library(raster)
# library(exactextractr)
# library(sf)
#library(lme4)
library(gee)
library(data.table)
library(tidyrules)
library(dplyr)
library(ggplot2)
#library(MASS)
#library(AER)

theme_set(theme_minimal())
theme_update(legend.position = "bottom")
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

#####################################################################################
###### ------------------------------ Country Level --------------------------- #####
#####################################################################################
#legacy
# ctry_p<- readRDS("X:/Spatial Stat/WASH Cholera/new_data/country.data.rds")
# country_df<- ctry_p%>%select(mean_pop_sum, W_Imp, W_Pip, W_Sur, S_Imp, S_OD, S_Pip, total_cases.int)
# fwrite(country_df, paste0(output_dir, 'country_6v.csv'))
country_df<- read.csv(paste0(output_dir, 'country_6v.csv')) # created from ctry_p

# function to get quasipoisson cooef and CI
convert_to_CI<- function(data, var){
  res<- c(
    round((exp(coef(glm(data[,"total_cases.int"]~data[,var], data = data, family = "quasipoisson", offset = log(mean_pop_sum)))[2])), 2), 
    suppressMessages(round((exp(confint(glm(data[,"total_cases.int"]~data[,var], data = data, family = "quasipoisson", offset = log(mean_pop_sum)))[2,])), 2)))
  print(res)
}
# get country results
country_result<- NULL
for (i in c(2:7)) { # updated to include new + % better chnage 
  tmp <- c(print(names(country_df)[i]), convert_to_CI(country_df, i))
  country_result <- rbind(country_result, tmp)
}

country_result<- as.data.frame(country_result)
colnames(country_result)<- c("var", "est", "lower", "upper")
country_result<- transform(country_result, est= as.numeric(est), upper= as.numeric(upper), lower= as.numeric(lower))
fwrite(country_result, paste0(output_dir, 'country_univariate_results.csv'))
# country_result<- read.csv(paste0(output_dir, 'country_univariate_results.csv'))
# ggplot(data = filter(country_result, var!='incdn__'), aes(x= as.numeric(est), y=var))+geom_pointrange(aes(xmin = as.numeric(upper), xmax = as.numeric(lower)))

# country multivariate regression
summary(glm(total_cases.int~W_Imp+W_Sur+W_Pip+S_Imp+S_Pip+S_OD, data = country_df, family = "quasipoisson", offset = log(mean_pop_sum)))


#####################################################################################
###### ----------------------------- District Level --------------------------- #####
#####################################################################################
# legacy
# dist_p<- readRDS("X:/Spatial Stat/WASH Cholera/clean_repo/results/dist_p.RDS")
# district_df<- dist_p%>%select( total_cases.int, NAME_0.x, mean_pop_sum, W_Pip, W_Imp, W_Sur, S_Pip, S_Imp, S_OD, incidence_in_thousan, hotspot)
# fwrite(district_df, paste0(output_dir, 'district_6v.csv'))

# read data 
district_df<- read.csv(paste0(output_dir, 'district_6v.csv'))

# GEE
# function to get cooef and CI
convert_to_CI_dist<- function(est, stde){
  e<- round((exp(est)), 4)
  ci<- cbind(upper= round((exp(est - 1.96*stde)), 4), lower= round((exp(est + 1.96*stde)), 4))
  cbind.data.frame(e, ci)
}
# get results 
dist_result<- NULL
for (i in 4:9) { 
  temp2<- cbind(names(district_df)[i], convert_to_CI_dist(summary(gee(district_df[,"total_cases.int"]~ district_df[,i]+ offset(log(district_df [, "mean_pop_sum"])), data = district_df , id= as.factor(NAME_0.x),  family = poisson, silent = T))$coefficients[2,1][1], 
                                                          summary(gee(district_df[,"total_cases.int"]~ district_df[,i]+ offset(log(district_df [, "mean_pop_sum"])), data = district_df , id= as.factor(NAME_0.x),  family = poisson, silent = T))$coefficients[2,4][1]))
  dist_result <- rbind(dist_result, temp2)
}

dist_result<- as.data.frame(dist_result)
colnames(dist_result)<- c("var", "est", "lower", "upper")

#ggplot(data = dist_result, aes(x= as.numeric(est), y=var))+geom_pointrange(aes(xmin = as.numeric(upper), xmax = as.numeric(lower)))


# ----------- multivariate model - District leve
summary(gee(total_cases.int~ 
              W_Pip+ 
              W_Imp+ 
              W_Sur+ 
              S_Pip+ 
              S_Imp+ 
              S_OD+ 
              offset(log(mean_pop_sum)), 
            data = district_df , 
            id= as.factor(NAME_0.x),  
            family = poisson, 
            silent = T))$coefficients


######  ----------- Combine country and district regression data --------------###### 
country_result<- cbind(country_result, source=rep("Country-level", 6))
dist_result<- cbind(dist_result, source= rep("Second-level", 6))
rownames(country_result)<- NULL
com_df<- full_join(country_result, dist_result)
com_df$var<- factor(com_df$var, 
                    levels = c("W_Pip", "W_Imp", "W_Sur", "S_Pip", "S_Imp", "S_OD"), 
                    labels = c("Piped water", 
                               "Piped or other improved water", 
                               "Surface water",
                               "Piped sanitation", 
                               "Piped or other improved sanitation", 
                               "Open defecation"))

saveRDS(object = com_df, file = 'X:/Spatial Stat/WASH Cholera/clean_repo/results/country_district_univariate_results.RDS')
fwrite(com_df, paste0(output_dir, 'country_district_univariate_results.csv'))


#####################################################################################
##### -------------------------- Figure 2 --------------------------------------##### 
#####################################################################################
com_df<- read.csv(paste0(output_dir, 'country_district_univariate_results.csv'))

levels(com_df$var) <- c("Piped water", 
                          "Piped or other improved water", 
                          "Surface water",
                          "Piped sanitation", 
                          "Piped or other improved sanitation", 
                          "Open defecation")
# change Second-level to District-level
com_df$source<- case_when(com_df$source=='Second-level' ~ 'District-level', com_df$source=='Country-level' ~ 'Country-level')

fig2<- ggplot(data = com_df, aes(x= est, y=var, color= source))+
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 3.5, color= "grey", linetype= 2)+ 
  geom_errorbar(aes(xmin = upper, xmax = lower), size= 1, width= .5, show.legend = T, position=position_dodge(width=0.5))+
  geom_point(aes(x= as.numeric(est), y=var), size= 2, position=position_dodge(width=0.5))+ 
  scale_color_manual(values = c("#2a9d8f", '#f9c74f'))+
  theme(panel.grid.major.y =  element_blank(),legend.title = element_blank())+
  xlab("Risk Ratio")+ylab("") + 
  #annotate("text", x= c(0.95, 1.05), y= c(8.5, 8.5), label= c("Increase of incidence", "Decrease of incidence"), size= 3)+
  scale_y_discrete(limits= rev(levels(com_df$var))) #+
  #annotate("text", x = c(.88), y=c(0,3.5), label = c("Sanitation \n predictor", "Water \npredictor"), angle= 90, size= 3, hjust= -.4)

# ggsave(paste0(output_dir, "plots/Poisson results_6v.jpg"), width = 6, height = 4, units = "in", dpi= 400)

png(paste0(output_dir, 'plots/Poisson results_6v.png'), 
    width = 6, height = 4, res = 400, units = 'in')
fig2
dev.off()

pdf(paste0(output_dir, 'plots/Poisson results_6v.pdf'), 
    width = 6, height = 4, bg= 'white')
fig2
dev.off()




#--------------------- Table S2 and risk ratios for main text ----------------------
View(com_df%>%mutate_if(is.numeric, round, 2)) # results saved in new_data folder Excel
# Table S2
com_df%>%
  mutate(out= paste0(round(est, 2),' (', round(lower, 2),'-',round(upper, 2),')'))%>%
  select(var, out, source)
# 1- est risk ratios for risk factor categories; doing 1- est to convert increase in cholera to decrease in cholera   
com_df%>% filter(var=='Open defecation' | var== 'Surface water')%>%
  mutate(out= paste0(scales::percent(1-est, accuracy = 0.01),' (', scales::percent(1-lower, accuracy = 0.01),'-',scales::percent(1-upper, accuracy = 0.01),')'))%>%
  select(var, out, source)





