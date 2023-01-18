library(gridExtra)
library(ggplot2)
library(plotROC)
library(pROC)
library(patchwork)

theme_set(theme_minimal())
theme_update(legend.position = "bottom")
options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

output_dir<- 'DIRECTORY_PATH'

scatter_plot<- readRDS(paste0(output_dir, 'incidence_rf_scatter_6vr_distCV.RDS'))
auc_plot<- readRDS(paste0(output_dir, 'rf_hotspot_auc_6vr_plot_distCV.RDS')) # district CV 



patchwork <- (scatter_plot + auc_plot)

fig3<- patchwork+ plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(.5, 1),
        plot.tag = element_text(size = 12, hjust = 0, vjust = 0))

png(paste0(output_dir, 'plots/fig3_scatter_auc_6vs_distCV.png'), 
    width = 6, height = 4, res = 400, units = 'in')
fig3
dev.off()

