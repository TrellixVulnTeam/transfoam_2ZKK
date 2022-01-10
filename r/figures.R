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

# ggsave(filename = "robustness.png", path = "figures/", width = 6, height = 4, device='png', dpi=1400)

## glucose production envelope
biomass <- read.csv('p_e/biomass_wt.csv', header = FALSE)
phb <- read.csv('p_e/phb_wt.csv', header = FALSE)

wt <- cbind(biomass, phb)
colnames(wt) <- c('biomass', 'phb_lb', 'phb_ub')
wt['color'] <- 'wt'

biomass <- read.csv('p_e/biomass_ac.csv', header = FALSE)
phb <- read.csv('p_e/phb_ac.csv', header = FALSE)

ac <- cbind(biomass, phb)
colnames(ac) <- c('biomass', 'phb_lb', 'phb_ub')
ac['color'] <- 'ac'

biomass <- read.csv('p_e/biomass_tca.csv', header = FALSE)
phb <- read.csv('p_e/phb_tca.csv', header = FALSE)

tca <- cbind(biomass, phb)
colnames(tca) <- c('biomass', 'phb_lb', 'phb_ub')
tca['color'] <- 'tca'

pe_wt <- ggplot()

pe_wt +
    geom_line(data = wt[1:2,], aes(x = biomass, y = phb_lb, color = color), size = 1, alpha = 0.9) +
    geom_line(data = wt[2:20,], aes(x = biomass, y = phb_lb, color = color), size = 1, alpha = 0.9) +
    geom_line(data = wt, aes(x = biomass, y = phb_ub, color = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(wt), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'lightgoldenrod', alpha = 0.2) +
    geom_line(data = ac[1:2,], aes(x = biomass, y = phb_lb, color = color), size = 1, alpha = 0.9) +
    geom_line(data = ac[2:20,], aes(x = biomass, y = phb_lb, color = color), size = 1, alpha = 0.9) +
    geom_line(data = ac, aes(x = biomass, y = phb_ub, color = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(ac), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'aquamarine2', alpha = 0.2) +
    geom_line(data = tca[1:2,], aes(x = biomass, y = phb_lb, color = color), size = 1, alpha = 0.9) +
    geom_line(data = tca[2:20,], aes(x = biomass, y = phb_lb, color = color), size = 1, alpha = 0.9) +
    geom_line(data = tca, aes(x = biomass, y = phb_ub, color = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(tca), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'salmon3', alpha = 0.2) +
    scale_color_manual(name = 'Phenotype', values = c('wt' = 'lightgoldenrod', 'ac' = 'aquamarine2', 'tca' = 'salmon3'), labels = c('Wild Type', 'Acetate Modifications', 'TCA Modifications')) +
    theme_minimal() +
    ggtitle('Styrene Minimal') +
    xlab(expression ('Growth Rate '~(h^-1))) +
    ylab('Citrate Synthase Flux (mmol/gDWh)') +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = c(0.8, 0.815)) +
    theme(legend.background = element_rect(fill='white',
                                           size=0.2, linetype='solid', 
                                           colour ='black')) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue")) 

ggsave(filename = "cs_sty_envelope.png", path = "figures/", width = 6, height = 5, device='png', dpi=1400)

## styrene production envelope
biomass <- read.csv('p_e_sty/biomass_wt.csv', header = FALSE)
phb <- read.csv('p_e_sty/phb_wt.csv', header = FALSE)

wt <- cbind(biomass, phb)
colnames(wt) <- c('biomass', 'phb_lb', 'phb_ub')
wt['color'] <- 'wt'

biomass <- read.csv('p_e_sty/biomass_ac.csv', header = FALSE)
phb <- read.csv('p_e_sty/phb_ac.csv', header = FALSE)

ac <- cbind(biomass, phb)
colnames(ac) <- c('biomass', 'phb_lb', 'phb_ub')
ac['color'] <- 'ac'

biomass <- read.csv('p_e_sty/biomass_tca.csv', header = FALSE)
phb <- read.csv('p_e_sty/phb_tca.csv', header = FALSE)

tca <- cbind(biomass, phb)
colnames(tca) <- c('biomass', 'phb_lb', 'phb_ub')
tca['color'] <- 'tca'

pe_wt <- ggplot()

pe_wt +
    geom_line(data = wt, aes(x = biomass, y = phb_lb, color = color), size = 2, alpha = 0.9) +
    geom_line(data = wt, aes(x = biomass, y = phb_ub, color = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(wt), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'lightgoldenrod', alpha = 0.2) +
    geom_line(data = ac, aes(x = biomass, y = phb_lb, color = color), size = 2, alpha = 0.9) +
    geom_line(data = ac, aes(x = biomass, y = phb_ub, color = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(ac), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'aquamarine2', alpha = 0.2) +
    geom_line(data = tca, aes(x = biomass, y = phb_lb, color = color), size = 2, alpha = 0.9) +
    geom_line(data = tca, aes(x = biomass, y = phb_ub, color = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(tca), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'salmon3', alpha = 0.2) +
    scale_color_manual(name = 'Phenotype', values = c('wt' = 'lightgoldenrod', 'ac' = 'aquamarine2', 'tca' = 'salmon3'), labels = c('Wild Type', 'Acetate Modifications', 'TCA Modifications')) +
    theme_minimal() +
    ggtitle('Styrene M9 Medium') +
    xlab(expression ('Growth Rate '~(h^-1))) +
    ylab('PHB Production (mmol/gDWh)') +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = c(0.8, 0.815)) +
    theme(legend.background = element_rect(fill='white',
                                           size=0.2, linetype='solid', 
                                           colour ='black')) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))

# ggsave(filename = "sty_envelope.png", path = "figures/", width = 6, height = 5, device='png', dpi=1400)

