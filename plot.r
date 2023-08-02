#!/usr/bin/env Rscript
library("ggplot2")
library("ggrepel")
library("ggsci")

data <- read.table('data.tsv', header=T)
data$Pval<-dbinom(data$Overlaps,data$CatalogNumber,data$Fraction,log=F)*data$CatalogNumber
data[order(data[,6]),]->data
data$rank<-nrow(data):1

png(file='output.png',height = 8, width = 8, res=600, units = "in", family="Arial")
ggplot(data, aes(x=Odds, y=rank,label=Disease)) + 
geom_point(shape=19, alpha=0.5, aes(size=Overlaps,color=Disease)) + 
xlab("Odds Ratio") + ylab("P-value Rank") + 
ggtitle ("GWAS SNPs enrichment - binomial test") + 
scale_size(range = c(5, 20))+
theme(axis.text=element_text(size=18), 
axis.title=element_text(size=18))+
guides(color = guide_legend(override.aes = list(size = 5)))+
scale_color_npg()+
geom_text_repel(aes(label=Disease, color = Disease), size=6,show.legend=F)
dev.off()

write.table(data,'data.tsv',sep='\t',quote=F,row.names=F, col.names=T)

