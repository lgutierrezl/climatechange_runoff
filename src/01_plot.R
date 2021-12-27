theme_gg <- function(){ 
  font <- "URWTimes"  
  
  theme_bw() %+replace%    
    
    theme(
      panel.background = element_rect(fill = "#E9EDF0"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = '#FFFFFF', linetype = 'solid', size = 0.5),
      axis.text.y = element_text(size = 10,hjust = 0.5, vjust = 0.5, colour = 'black'),
      axis.text.x = element_text(size = 10,hjust = 0.5, vjust = 0.5, colour = 'black', angle = 20),
      axis.title.y = element_text(size = 9,hjust = 0.5, vjust = 0.5, angle = 90),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 15,hjust = 0.5, vjust = 0.5, colour = 'black', face = 'bold'),
      strip.background = element_rect(fill = 'white'),
      strip.text = element_text(face = "bold", size=10),
      strip.placement = "outside",
      legend.box.background = element_rect(color="white", size=0.5),
      legend.text = element_text(size = 11, colour = "black"),
      legend.title = element_text(face = "bold", size=11, hjust = 0.5, vjust = 0.5),
      legend.box.margin = margin(rep(0.5,4)),
      legend.margin = margin(rep(0.5,4)),
      legend.position='bottom'
    )
}

fig_chn <- function(datab = rnf_cc){
  barras2 <- ggplot(datab, aes(x = MOD, y = RNF, fill = MOD)) + 
    geom_hline(yintercept = 0, size=1.2, alpha=0.75, col='tomato')+
    geom_boxplot(aes(fill = MOD), outlier.alpha = 0, size = 0.2, alpha = 0.4)+
    labs(x = "Modelos", y = 'Cambio (%)')+
    guides(fill = guide_legend(title.position = "top",
                               reverse=T,title = 'Modelo',
                               title.vjust = 0.5, nrow = 2))+
    scale_fill_manual(values = c('Presente' = "#DA2E20",
                                 'ACCESS 1.0' = "#94ED3C",
                                 'HadGEM2-ES' = "#4192CF",
                                 'MPI-ESM-LR' = "#FFB036"))+
    stat_summary(aes(group=MOD),fun="mean", geom="point", size=1.4, position="identity", color="#AE2D10")+
    theme_gg()+
    theme(legend.position='none')
  return(barras2)
}

fig_abs <- function(datab = REG_DEPc){
  barras1 <- ggplot(datab, aes(x = MOD, y = RNF, fill = MOD)) + 
    geom_boxplot(aes(fill = MOD), outlier.alpha = 0, size = 0.2, alpha = 0.4)+
    labs(x = "Modelos", y = 'Escurrimiento (mm)')+
    guides(fill = guide_legend(title.position = "top",
                               reverse=T,title = 'Modelo',
                               title.vjust = 0.5, nrow = 2))+
    scale_fill_manual(values = c('Presente' = "#DA2E20",
                                 'ACCESS 1.0' = "#94ED3C",
                                 'HadGEM2-ES' = "#4192CF",
                                 'MPI-ESM-LR' = "#FFB036"))+
    stat_summary(aes(group=MOD),fun="mean", geom="point", size=1.4, position="identity", color="#AE2D10")+
    theme_gg()+
    theme(legend.position='top')
  return(barras1)
}

