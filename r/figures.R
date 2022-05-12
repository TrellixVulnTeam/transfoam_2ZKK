library(dplyr)
library(ggplot2)


## sty robustness
c <- read.csv('robustness/controlFlux.csv', header = FALSE)
colnames(c) <- 'controlFlux'
o <- read.csv('robustness/objFlux.csv', header = FALSE)
colnames(o) <- 'objFlux'

df <- cbind(c, o)

p <- ggplot(data = df, aes(x = controlFlux, y = objFlux))

p + 
    geom_line(color = 'aquamarine3') +
    theme_bw() +
    xlab('Styrene Uptake (mmol/gDWh)') +
    ylab('PHB Output (mmol/gDWh)') +
    coord_cartesian(xlim = c(-13, 0), ylim = c(0, 10)) +
    scale_x_continuous(expand = c(0,0), breaks = c(-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"), panel.grid.minor = element_blank(), plot.background = element_rect(fill='transparent', color = NA)) 

ggsave(filename = "robustness.png", path = "figures/", width = 6, height = 4, device='png', dpi=1400)


## glc robustness
c <- read.csv('robustness/controlFlux_glc.csv', header = FALSE)
colnames(c) <- 'controlFlux'
o <- read.csv('robustness/objFlux_glc.csv', header = FALSE)
colnames(o) <- 'objFlux'

df <- cbind(c, o)

p <- ggplot(data = df, aes(x = controlFlux, y = objFlux))

p + 
    geom_line(color = 'aquamarine3') +
    theme_bw() +
    xlab('Glucose Uptake (mmol/gDWh)') +
    ylab('PHB Output (mmol/gDWh)') +
    coord_cartesian(xlim = c(-10, 0), ylim = c(0, 13)) +
    scale_x_continuous(expand = c(0,0), breaks = c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0)) +
    scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"), panel.grid.minor = element_blank(), plot.background = element_rect(fill='transparent', color = NA)) 

ggsave(filename = "robustness_glc.png", path = "figures/", width = 6, height = 4, device='png', dpi=1400)

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

# ggsave(filename = "cs_sty_envelope.png", path = "figures/", width = 6, height = 5, device='png', dpi=1400)

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
    scale_color_manual(name = 'Phenotype', values = c('wt' = 'lightgoldenrod', 'ac' = 'aquamarine2', 'tca' = 'salmon3'), labels = c('Control', 'Acetate Modifications', 'TCA Modifications')) +
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

## atpm production envelope
biomass <- read.csv('p_e_atpm/obj_wt.csv', header = FALSE)
phb <- read.csv('p_e_atpm/phb_wt.csv', header = FALSE)

wt <- cbind(biomass, phb)
colnames(wt) <- c('biomass', 'phb_lb', 'phb_ub')
wt['color'] <- 'wt'

biomass <- read.csv('p_e_atpm/obj_o.csv', header = FALSE)
phb <- read.csv('p_e_atpm/phb_o.csv', header = FALSE)

o <- cbind(biomass, phb)
colnames(o) <- c('biomass', 'phb_lb', 'phb_ub')
o['color'] <- 'o'

pe_wt <- ggplot()

pe_wt +
    geom_line(data = wt, aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 2, alpha = 0.9) +
    geom_line(data = wt, aes(x = biomass, y = phb_ub, color = color, linetype = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(wt), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'lightgoldenrod', alpha = 0.2) +
    geom_line(data = o[1:19,], aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 2, alpha = 0.9) +
    geom_line(data = o[19:20,], aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 1, alpha = 0.9) +
    geom_line(data = o, aes(x = biomass, y = phb_ub, color = color, linetype = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(o), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'aquamarine2', alpha = 0.2) +
    scale_color_manual(name = 'Phenotype', values = c('wt' = 'lightgoldenrod', 'o' = 'aquamarine2'), labels = c('Control', 'ICDHyr KO')) +
    scale_linetype_manual(name = 'Phenotype', values = c('wt' = 'solid', 'o' = '42'), labels = c('Control', 'ICDHyr KO')) +
    guides(color = guide_legend(override.aes = list(size = .5))) +
    theme_minimal() +
    ggtitle('Coupling PHB Production to ATPM (Glucose M9)') +
    xlab('ATP Synthesis (mmol/gDWh)') +
    ylab('PHB Production (mmol/gDWh)') +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = c(0.855, 0.865)) +
    theme(legend.background = element_rect(fill='white',
                                           size=0.2, linetype='solid', 
                                           colour ='black')) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))

ggsave(filename = "atpm_envelope.png", path = "figures/", width = 6, height = 5, device='png', dpi=1400)

## ac+tca new production envelope
biomass <- read.csv('res_bm/obj_wt_glc.csv', header = FALSE)
phb <- read.csv('res_bm/phb_wt_glc.csv', header = FALSE)

wt <- cbind(biomass, phb)
colnames(wt) <- c('biomass', 'phb_lb', 'phb_ub')
wt['color'] <- 'wt'

biomass <- read.csv('res_bm/obj_pta_glc.csv', header = FALSE)
phb <- read.csv('res_bm/phb_pta_glc.csv', header = FALSE)

pta <- cbind(biomass, phb)
colnames(pta) <- c('biomass', 'phb_lb', 'phb_ub')
pta['color'] <- 'pta'

biomass <- read.csv('res_bm/obj_ack_glc.csv', header = FALSE)
phb <- read.csv('res_bm/phb_ack_glc.csv', header = FALSE)

ack <- cbind(biomass, phb)
colnames(ack) <- c('biomass', 'phb_lb', 'phb_ub')
ack['color'] <- 'ack'

biomass <- read.csv('res_bm/obj_ex_glc.csv', header = FALSE)
phb <- read.csv('res_bm/phb_ex_glc.csv', header = FALSE)

ex <- cbind(biomass, phb)
colnames(ex) <- c('biomass', 'phb_lb', 'phb_ub')
ex['color'] <- 'ex'

pe_wt <- ggplot()

pe_wt +
    geom_line(data = wt, aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 2, alpha = 0.9) +
    geom_line(data = wt, aes(x = biomass, y = phb_ub, color = color, linetype = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(wt), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'lightgoldenrod', alpha = 0.2) +
    # geom_line(data = pta, aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 2, alpha = 0.9) +
    # geom_line(data = pta, aes(x = biomass, y = phb_ub, color = color, linetype = color), size = 1, alpha = 0.9) +
    # geom_ribbon(data = subset(pta), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'aquamarine2', alpha = 0.2) +
    # geom_line(data = ack[1:19,], aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 2, alpha = 0.9) +
    # geom_line(data = ack[19:20,], aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 1, alpha = 0.9) +
    # geom_line(data = ack, aes(x = biomass, y = phb_ub, color = color, linetype = color), size = 1, alpha = 0.9) +
    # geom_ribbon(data = subset(ack), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'salmon3', alpha = 0.2) +
    geom_line(data = ex[1:19,], aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 2, alpha = 0.9) +
    geom_line(data = ex[19:20,], aes(x = biomass, y = phb_lb, color = color, linetype = color), size = 1, alpha = 0.9) +
    geom_line(data = ex, aes(x = biomass, y = phb_ub, color = color, linetype = color), size = 1, alpha = 0.9) +
    geom_ribbon(data = subset(ex), aes(x = biomass, ymin = phb_lb, ymax = phb_ub), fill = 'aquamarine2', alpha = 0.1) +
    # scale_color_manual(name = 'Phenotype', values = c('wt' = 'lightgoldenrod', 'pta' = 'aquamarine2', 'ack' = 'salmon3', 'ex' = 'darkorchid'), labels = c('Control', 'ACKr + PTAr + SUCOAS KO', 'AC Production + SUCOAS KO', 'AC Exchange + SUCOAS KO')) +
    # scale_linetype_manual(name = 'Phenotype', values = c('wt' = 'solid', 'pta' = '42', 'ack' = '45', 'ex' = '13'), labels = c('Control', 'ACKr + PTAr + SUCOAS KO', 'AC Production + SUCOAS KO', 'AC Exchange + SUCOAS KO')) +
    
    # scale_color_manual(name = 'Phenotype', values = ('wt' = 'lightgoldenrod'), labels = 'Control') +
    # scale_linetype_manual(name = 'Phenotype', values = ('wt' = 'solid'), labels = 'Control') +
    
    scale_color_manual(name = 'Phenotype', values = c('wt' = 'lightgoldenrod', 'ex' = 'aquamarine2'), labels = c('Control', 'AC Production & SUCOAS KO')) +
    scale_linetype_manual(name = 'Phenotype', values = c('wt' = 'solid', 'ex' = '42'), labels = c('Control', 'AC Production & SUCOAS KO')) +
    # sty: 42, 45, 13
    # glc: 42, 11, 35
    guides(color = guide_legend(override.aes = list(size = .5))) +
    theme_minimal() +
    ggtitle('Coupling PHB Production to Growth') +
    xlab(expression ('Growth Rate '~(h^-1))) +
    # xlab('ATP Synthesis (mmol/gDWh)') +
    ylab('PHB Production (mmol/gDWh)') +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position = c(0.785, 0.81)) +
    # c(0.215, 0.195))
    # c(0.785, 0.81))
    theme(legend.background = element_rect(fill='white',
                                           size=0.2, linetype='solid', 
                                           colour ='black')) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"), plot.background = element_rect(fill='transparent', color = NA), panel.background = element_rect(fill='white'))

ggsave(filename = "new_both.png", path = "figures/", width = 6, height = 5, device='png', dpi=1400)
 
