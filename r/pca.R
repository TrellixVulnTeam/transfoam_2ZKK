library(dplyr)
library(ggplot2)

s <- read.csv('sampling/s.csv')
t <- read.csv('sampling/t.csv')

test <- bind_rows(s, t)

pca <- prcomp(test)

pca_out <- as.data.frame(pca$x)
pca_out$medium <- 'glc'
pca_out <- pca_out[,c(ncol(pca_out), 1:ncol(pca_out))]
pca_out[21:nrow(pca_out), 1] <- 'phleth'

p <- ggplot(pca_out, aes(x=PC1, y=PC2, color=medium))
p +
    geom_point(aes(color = medium), alpha = 0.5,subset(pca_out,medium != 'glc')) +
    geom_point(aes(color = medium), alpha = 0.5,subset(pca_out,medium != 'phleth')) +
    ggtitle('growth medium flux samples') +
    theme_minimal() +
    xlab('PC1') + 
    ylab('PC2') + 
    labs(color = 'growth medium')  +
    theme(legend.position = c(0.725, 0.225)) +
    theme(legend.background = element_rect(fill='white',
                                           size=0.2, linetype='solid', 
                                           colour ='black')) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), text = element_text(family = "Helvetica Neue"))
