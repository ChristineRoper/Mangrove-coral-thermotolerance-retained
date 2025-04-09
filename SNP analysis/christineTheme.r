library(ggplot2)
christineTheme = theme(
    panel.background = element_blank(),
    panel.border=element_rect(fill=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background=element_blank(),
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),
    axis.ticks=element_line(colour="black"),
    plot.margin=unit(c(0.1,0.1,0.1,0.1), "line")
)

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set3"))

saveGraph = function(name){
    ggsave(paste0(name,'-81.svg'), device='svg',  width=82, units='mm')
    ggsave(paste0(name,'-105.svg'), device='svg',  width=105, units='mm')
    ggsave(paste0(name,'-169.svg'), device='svg',  width=169, units='mm')
}