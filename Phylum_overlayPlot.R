

# overlay plot of phylum data read out of MPA
# and phylum level OTU table from 16S rRNA data
setwd("~/Desktop/Non work bits/Co-ex_R")
library(vegan)
library(Rcmdr)
library(reshape2)

## ===== with just bacteria and archaea and 16S data======

DNA2 = read.table('Phylum_dnaRNA.txt', header=T, sep='\t')
DNAx = DNA2[,2:7]
DNAxPC = as.data.frame(colPercents(DNAx))
DNAPC = DNAxPC[1:15,]
DNAPC$tax = DNA2$tax



d2 = as.data.frame(colPercents(DNAx))
d2 = d2[1:15,]
d2$tax = DNA2$tax

prot2 = read.table('Soil2_protein.txt', header=T, sep='\t')
p2a = prot2[,2:4]
p2 = colPercents(p2a)
p2 = as.data.frame(p2)
p2 = p2[1:13,]
p2$tax = prot2$tax

prot4 = read.table('Soil4_protein.txt', header=T, sep='\t')
p4a = prot4[,2:4]
p4 = as.data.frame(colPercents(p4a))
p4 = p4[1:11,]
p4$tax = prot4$tax

prot7 = read.table('Soil8_protein.txt', header=T, sep='\t')
p7a = prot7[,2:4]
p7 = as.data.frame(colPercents(p7a))
p7 = p7[1:9,]
p7$tax = prot7$tax


## ===== chosing the top OTUs ========

#  chose the taxa from DNA, RNA and protein that represented 5% or above of
# total 'hits' at level of phylum 
# these were then 'merged' into one list which is what will be used. 
# DNA and RNA had same top 5 
# two of the 'top 5' were the same in nucleic and protein samples (proteobac and actinobac)

# Proteobacteria|Actinobacteria|Verrucomicrobia|Acidobacteria|Euryarchaeota|Cyanobacteria|Bacteroidetes|Firmicutes

nucTop5 = DNAPC[grep("Proteobacteria|Actinobacteria|Verrucomicrobia|Acidobacteria|Euryarchaeota|Cyanobacteria|Bacteroidetes|Firmicutes", DNAPC$tax),]

# read in files with zeros
eury = read.table('Eury-DNA.txt', header=T, sep='\t')
nucTop5full = rbind(nucTop5, eury)

DNATop5 = nucTop5full[,c(1,2,3,7)]
mDNA = melt(DNATop5)
#mDNA$tax <- factor(mDNA$tax, levels = unique(mDNA$tax)) # assign order

RNATop5 = nucTop5full[,c(4:7)]
mRNA = melt(RNATop5)
#mRNA$tax <- factor(mRNA$tax, levels = unique(mDNA$tax)) # assign order  according to DNA plot


# now for protein
p2Top5 = p2[grep("Proteobacteria|Actinobacteria|Verrucomicrobia|Acidobacteria|Euryarchaeota|Cyanobacteria|Bacteroidetes|Firmicutes", p2$tax),]
verr = read.table('verru-prot.txt', header=T, sep='\t')
protTop5 = rbind(p2Top5, verr)
mprot= melt(protTop5)
#mprot$tax <- factor(mprot$tax, levels = unique(mDNA$tax)) # assign order according to DNA plot


# run setFactorOrder.R and assign the same order to phyla in each dataset
mDNA[["tax"]] <- setFactorOrder(mDNA[["tax"]], c("Verrucomicrobia", "Proteobacteria", "Firmicutes", "Euryarchaeota",  "Cyanobacteria", "Bacteroidetes", "Actinobacteria","Acidobacteria" ))
mRNA[["tax"]] <- setFactorOrder(mRNA[["tax"]], c("Verrucomicrobia", "Proteobacteria", "Firmicutes", "Euryarchaeota",  "Cyanobacteria", "Bacteroidetes", "Actinobacteria","Acidobacteria" ))
mprot[["tax"]] <- setFactorOrder(mprot[["tax"]], c("Verrucomicrobia", "Proteobacteria", "Firmicutes", "Euryarchaeota",  "Cyanobacteria", "Bacteroidetes", "Actinobacteria","Acidobacteria" ))

##  ===== finally to actually plotting =======

library(ggplot2)

FacetLabs <- c(
  'DNA1' = 'DNA - Replicate 1', 
  'DNA2' = 'DNA - Replicate 2',
  'DNA3' = 'DNA - Replicate 3')

DNAplot <- ggplot(mDNA, aes(x=tax, y=value, fill=tax), show_guide=FALSE) +
  geom_bar(position="dodge",stat="identity") + coord_flip() +
  scale_fill_manual(values = rev(c("darkred", "orange", "seagreen3","darkgreen", "cornflowerblue",'mediumorchid1', 'purple4', "black"))) +
  theme_bw() + labs(y="") +ylim(0,75) +
  theme(legend.key.size = unit(2.5,"line"),
        axis.text.x=element_text(size=14, colour='black'),
        axis.text.y=element_text(size=13, colour='black'),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, colour='black'),
        legend.text = element_text(size=13, colour='black'),
        legend.title = element_text(size=15, colour='black'),
        strip.text.x = element_text(size = 13),
        strip.background =  element_rect(fill = "white")) + guides(fill=FALSE) +
        facet_wrap(~variable, labeller = as_labeller(FacetLabs))
DNAplot 

RNALabs <- c(
  'cDNA1' = 'RNA - Replicate 1', 
  'cDNA2' = 'RNA - Replicate 2',
  'cDNA3' = 'RNA - Replicate 3')

RNAplot <- ggplot(mRNA, aes(x=tax, y=value, fill=tax), show_guide=FALSE) +
  geom_bar(position="dodge",stat="identity") + coord_flip() +
  scale_fill_manual(values = rev(c("darkred", "orange", "seagreen3","darkgreen", "cornflowerblue",'mediumorchid1', 'purple4', "black"))) +
  theme_bw() + labs(y="") + ylim(0,75) +
  theme(legend.key.size = unit(2.5,"line"),
        axis.text.x=element_text(size=14, colour='black'),
        axis.text.y=element_text(size=13, colour='black'),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, colour='black'),
        legend.text = element_text(size=13, colour='black'),
        legend.title = element_text(size=15, colour='black'),
        strip.text.x = element_text(size = 13),
        strip.background =  element_rect(fill = "white")) + guides(fill=FALSE) +
  facet_wrap(~variable, labeller = as_labeller(RNALabs))
RNAplot 

protLabs <- c(
  'X2a_protein' = 'Protein - Replicate 1', 
  'X2b_protein' = 'Protein - Replicate 2',
  'X2c_protein' = 'Protein - Replicate 3')

protPlot <- ggplot(mprot, aes(x=tax, y=value, fill=tax), show_guide=FALSE) +
  geom_bar(position="dodge",stat="identity") + coord_flip() +
  scale_fill_manual(values = rev(c("darkred", "orange", "seagreen3","darkgreen", "cornflowerblue",'mediumorchid1', 'purple4', "black"))) +
  theme_bw() + labs(y="Relative Abundance (%)") + ylim(0,75) +
  theme(legend.key.size = unit(2.5,"line"),
        axis.text.x=element_text(size=14, colour='black'),
        axis.text.y=element_text(size=13, colour='black'),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=14, colour='black'),
        legend.text = element_text(size=13, colour='black'),
        legend.title = element_text(size=15, colour='black'),
        strip.text.x = element_text(size = 13),
        strip.background =  element_rect(fill = "white")) + guides(fill=FALSE) +
  facet_wrap(~variable, labeller = as_labeller(protLabs))
protPlot 

# now arrange these three plots above each other with ggpubr
library(ggpubr)

ggarrange(DNAplot, RNAplot, protPlot,
          nrow=3)

