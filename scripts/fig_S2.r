# make figure for mapping efficiency

library(ggplot2)
library(plyr)

map_stats<-read.table("map_stats.txt", header=TRUE)

ggplot(map_stats, aes(x=State, y=Sym, fill=State, colour=State))+
  geom_boxplot(aes(alpha=0.2), outlier.size = NA)+
  geom_point(size=3)+
  geom_line(aes(group=Colony))+
  theme_minimal()+
  scale_fill_manual(values=c("#80CDC1", "#BF812D"))+
  scale_colour_manual(values=c("#003C30", "#543005"))+
  annotate("text", label="p=0.012", x=1.5, y=.23, size=4, fontface="italic")+
  theme(axis.text = element_text(face="bold", size=8, colour="black"), 
        axis.title = element_text(face="bold", size=10, colour="black"),
        panel.grid.minor=element_blank(),
        legend.position = "none")+
  guides(alpha=FALSE)+guides(colour=FALSE)+
  annotate("segment", x=1, xend=1, y=0.205, yend=0.22)+
  annotate("segment", x=1, xend=2, y=0.22, yend=0.22)+
  annotate("segment", x=2, xend=2, y=0.205, yend=0.22)+
  ylab("Symbiont Mapped Reads")+scale_y_continuous(breaks=c(0.05, 0.1,0.15,.2),limits = c(0,.23), labels=c("5%", "10%", "15%", "20%"))+
  scale_x_discrete(labels=c("Aposymbiotic", "Symbiotic"))

ggsave(filename="/Users/hannyrivera/Documents/BU/DaviesLab/_Data/Oculina/Mapped_counts.nosync/Coral_collapsed_transc_best_map/figures/mapping_stats/figS2.png", width = 3, height=4, units= "in", dpi=300)      

t.test(subset(map_stats, State=="Apo")$Sym, subset(map_stats, State=="Sym")$Sym, alternative="less", paired=TRUE)

