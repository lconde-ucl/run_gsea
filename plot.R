library(ggplot2) 
library(dplyr)


data<-read.table("final_table.txt", sep="\t", header=T)

#- Get range of NES for ylim
rangeNES=round(max(abs(max(data$NES)), abs(min(data$NES)))+0.5)

#- Replace FDR = 0 to FDR = 0.001 (as this was run on 1000 iterations)
data$FDR_p[data$FDR_p == 0] <- 0.001

#- and readjust pvalues
data<-data %>%  mutate(FDR=p.adjust(FDR_p, method="BH"))

#- get plot (based on code from https://www.biostars.org/p/168044/)

png("final_plot.png", width=800)
p <- ggplot(data, aes(NES, GENESET)) + 
    geom_point(aes(colour=FDR, size=RATIO)) +
    scale_color_gradient(limits=c(0, 0.10), low="red", high="white") +
    geom_vline(xintercept=0, size=0.5, colour="gray50") +
    theme(strip.text.x = element_text(size = 8), 
	  panel.background=element_rect(fill="gray95", colour="gray95"),
          panel.grid.major=element_line(size=0.25,linetype='solid', colour="gray90"), 
          panel.grid.minor=element_line(size=0.25,linetype='solid', colour="gray90"),
          axis.title.y=element_blank()) +
    expand_limits(x=c(-rangeNES,rangeNES)) +
    scale_x_continuous(breaks=seq(-rangeNES, rangeNES, 2)) +
    facet_grid(.~RANK+TYPE)
print(p)
dev.off()

