#   Figure 1
library(tidyr)
library(ggplot2)
library(ggnewscale)
library(ggExtra)
library(ggpubr)

theme_set(theme_minimal())
theme_update(legend.position = "bottom")
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")
options(digits=2)

output_dir<- 'DIRECTORY_PATH'



################################### FIGURE 1 TILE PLOT ############################## 
ctry_0<- sf::st_read("country_wash_incidence.shp")
country_data<- sf::st_drop_geometry(ctry_0)
# remove "piped or other" categories 
country_data<- base::subset(country_data, select= -c(W_Imp, S_Imp))
country_data<- country_data%>%pivot_longer(cols = !c(GID_0, NAME_0, incdn__))#incdn__ is incidence in 1,000 population 
country_data$NAME_0<- ifelse(country_data$NAME_0=="Democratic Republic of the Congo", "DR Congo", country_data$NAME_0)
country_data$NAME_0<- ifelse(country_data$NAME_0=="Central African Republic", "Central African Rep.", country_data$NAME_0)
country_data$NAME_0<- factor(country_data$NAME_0)
country_data$Incidence<- log10(country_data$incdn__)


country_data$name<- factor(country_data$name, 
                           levels = c("S_ImpOt",
                                      "S_Pip", 
                                      "S_OD", 
                                      "S_Uni", 
                                      "W_ImpOt", 
                                      "W_Pip", 
                                      "W_Sur", 
                                      "W_Uni"), 
                           labels = c("Other improved sanitation",
                                      "Piped sanitation", 
                                      "Open defecation", 
                                      "Unimproved sanitation", 
                                      "Other improved water", 
                                      "Piped water", 
                                      "Surface water", 
                                      "Unimproved water"))

levels(country_data$name) <- c("Piped water",
                               "Other improved water", 
                               "Piped sanitation",
                               "Other improved sanitation",
                               "Unimproved water",
                               "Surface water",
                               "Unimproved sanitation",
                               "Open defecation")



rasterplot <-
  ggplot() +
  geom_tile(data = subset(country_data, name %in% c("Other improved sanitation",
                                                    "Piped sanitation", 
                                                    "Other improved water", 
                                                    "Piped water" 
                                                    )),
            aes(y = NAME_0, x = name, fill = value)) +
  scale_fill_gradient(
    low = '#006d66',
    high = "#cff3eb",
    guide = guide_colorbar(
      title.position = "top",
      title.theme = element_text(angle = 0, size = 10), 
      order = 1
    )
  ) +
  labs(fill = "Access (%)") +
  ggnewscale::new_scale("fill") +
  geom_tile(data = subset(country_data, name %in% c("Open defecation", 
                                                    "Unimproved sanitation", 
                                                    "Surface water", 
                                                    "Unimproved water"
  )),
            aes(y = NAME_0, x = name, fill = value)) +
  scale_fill_gradient(
    low = '#fff75e',
    high = "#9d0208",
    guide = guide_colorbar(
      title.position = "top",
      title.theme = element_text(angle = 0, size = 10), 
      order = 2
    )
  ) +
  labs(fill = "Reliance (%)") +
  
  geom_vline(xintercept = 4.5, color = "#edae49") +
  scale_y_discrete(limits = rev(levels(country_data$NAME_0))) +
  scale_x_discrete(limits = levels(country_data$name)) +
  
  ggnewscale::new_scale("fill") +
  geom_tile(data = country_data , aes(x = 9.2, y = NAME_0, fill = Incidence)) +
  scale_fill_gradient(
    low = 'white',
    high = 'red',
    guide = guide_colorbar(
      title.position = "top",
      title.theme = element_text(angle = 0, size = 10)
    ),
    labels = c("1e-3", "0.01", "0.1", "1")
  ) +
  labs(fill = " Incidence ") +
  
  annotate(
    "text",
    label = "Protective",
    x = 2.5,
    y = 41.5,
    size = 3,
    hjust = .5,
    vjust = 1
  ) +
  annotate(
    "text",
    label = "Risk",
    x = 6.5,
    y = 41.5,
    size = 3,
    hjust = .5,
    vjust = 1
  ) +
  
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1.1
    ),
    panel.grid =   element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.height = unit(.3, "cm"),
    legend.key.width = unit(.5, 'cm'),
    legend.title.align = .5
  )



#######################FIGURE 1 BIVARIATE PLOT WITH DENSITY #########################
# bivariate plot with density 

district_df<- read.csv(paste0(output_dir, 'district_6v.csv'))
# viridis points
imp_wat_san<- ggplot(data = district_df, aes(x= W_Imp, y= S_Imp, color= (log10(incidence_in_thousan))))+
  geom_point(alpha= .3)+ 
  labs(y= "Piped or other improved sanitation", x= "Piped or other improved water", color= "Incidence")+
  theme(legend.title = element_text(vjust = -.1, size = 9, hjust = .5), 
        aspect.ratio = 1/1, 
        legend.key.height =  unit(.3, "cm"), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 9))+ 
  scale_color_continuous (guide= guide_colorbar(title.position = "top"), labels= c("1e-05", "", "", "0.01", "", "", "10")) + 
  scale_x_continuous(labels = function(x) paste0(x, "%"))+ 
  scale_y_continuous(labels = function(x) paste0(x, "%"))
imp_wat_san<- ggMarginal(imp_wat_san, type = "histogram" ,xparams = list(binwidth= 2, color= "white", fill= "#00b4d8"), yparams = list(binwidth= 2, color= "white", fill= "#f3722c"))

# grey points
imp_wat_san<- ggplot(data = district_df, aes(x= W_Imp, y= S_Imp, color= (log10(incidence_in_thousan))))+
  geom_point(alpha= .3)+ 
  labs(y= "Piped or other improved sanitation", x= "Piped or other improved water", color= "Incidence")+
  theme(legend.title = element_text(vjust = -.1, size = 9, hjust = .5), 
        aspect.ratio = 1/1, 
        legend.key.height =  unit(.3, "cm"), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 9))+ 
  scale_color_gradient (low= 'white', high= 'black', guide= guide_colorbar(title.position = "top"), labels= c("1e-05", "", "", "0.01", "", "", "10")) + 
  scale_x_continuous(labels = function(x) paste0(x, "%"))+ 
  scale_y_continuous(labels = function(x) paste0(x, "%"))
imp_wat_san<- ggMarginal(imp_wat_san, type = "histogram" ,xparams = list(binwidth= 2, color= "white", fill= "#00b4d8"), yparams = list(binwidth= 2, color= "white", fill= "#f3722c"))

# viridis points
surf_odef<- ggplot(data = district_df, aes(W_Sur, S_OD, color= log10(incidence_in_thousan)))+
  geom_point(alpha= .3)+
  labs(x= "Surface Water", y= "Open Defecation", color= "Log of incidence")+
  theme(legend.position = "none", 
        aspect.ratio = 1/1, 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 9))+ 
  scale_color_continuous (guide= guide_colorbar(title.position = "top"))+ 
  scale_x_continuous(labels = function(x) paste0(x, "%"))+ 
  scale_y_continuous(labels = function(x) paste0(x, "%"))
surf_odef<- ggMarginal(surf_odef, type = "histogram" ,xparams = list(binwidth= 2, color= "white", fill= "#00b4d8"), yparams = list(binwidth= 2, color= "white", fill= "#f3722c"))

# grey points
surf_odef<- ggplot(data = district_df, aes(W_Sur, S_OD, color= log10(incidence_in_thousan)))+
  geom_point(alpha= .3)+
  labs(x= "Surface Water", y= "Open Defecation", color= "Log of incidence")+
  theme(legend.position = "none", 
        aspect.ratio = 1/1, 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 9))+ 
  scale_color_gradient (low= 'white', high= 'black', guide= guide_colorbar(title.position = "top"))+ 
  scale_x_continuous(labels = function(x) paste0(x, "%"))+ 
  scale_y_continuous(labels = function(x) paste0(x, "%"))
surf_odef<- ggMarginal(surf_odef, type = "histogram" ,xparams = list(binwidth= 2, color= "white", fill= "#00b4d8"), yparams = list(binwidth= 2, color= "white", fill= "#f3722c"))


plot1<- ggarrange(rasterplot, ncol = 2, labels = c("A"), ggarrange(surf_odef, imp_wat_san, ncol = 1, labels = c("B", "C"), heights = c(.8, 1)))
ggsave(paste0(output_dir, "plots/country_district_paner.jpeg"), plot= plot1, width = 7, height = 7, units = "in", dpi= 600)




