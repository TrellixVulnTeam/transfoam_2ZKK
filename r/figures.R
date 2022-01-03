library(dplyr)
library(ggplot2)

c <- read.csv('robustness/controlFlux.csv', header = FALSE)
colnames(c) <- 'controlFlux'
o <- read.csv('robustness/objFlux.csv', header = FALSE)
colnames(o) <- 'objFlux'

df <- cbind(c, o)

p <- ggplot(data = df, aes(x = controlFlux, y = objFlux))

p + 
    geom_line(color = 'aquamarine3') +
    theme_bw() +
    xlab('Styrene Uptake (mmol/(gDw*hr))') +
    ylab('PHB Output (mmol/(gDw*hr))') +
    coord_cartesian(xlim = c(-13, 0), ylim = c(0, 10)) +
    scale_x_continuous(expand = c(0,0), breaks = c(-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"), panel.grid.minor = element_blank()) 

ggsave(filename = "robustness.png", path = "figures/", width = 6, height = 4, device='png', dpi=1400)
