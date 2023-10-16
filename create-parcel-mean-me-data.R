library(data.table)
source('get_roi_data.R')
d <- get_roi_data(roi = 1, cope = 'fearGTcalm', ncpus = 16, tsv_dir = 'tsv_nosmooth')

library(ggplot2)
ggplot(d, aes(x = y)) + 
  geom_histogram(fill = 'black', position = position_identity(), alpha = 1, bins = 100) + 
  facet_grid(sess~id) + 
  coord_cartesian(xlim = c(-50, 50), y = c(0, 75)) + 
  theme_minimal() + 
  theme(strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_blank(), panel.grid = element_blank(),
        axis.title = element_blank())

library(sandwich)
library(lmtest)
