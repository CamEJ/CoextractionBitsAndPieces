
# nmds of phylum data read out of MPA
setwd("~/Desktop/Non work bits/Co-ex_R")
library(vegan)


## ===== with just bacteria and archaea and 16S data======
library(Rcmdr)

DNA2 = read.table('Phylum_dnaRNA.txt', header=T, sep='\t')
DNAx = DNA2[,2:7]


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

all1 = merge(p4, p7, by = "tax", all = "TRUE")

all2 = merge(d2, p2, by = "tax", all = "TRUE")

all = merge(all1, all2, by = "tax", all = "TRUE")

all[is.na(all)] <- 0

write.csv(all, 'phylumNMDSdata.csv')

all = read.csv('phylumNMDSdata.csv', header=T, sep=',')

soil2 = all[,c(9:17)]

protsy = all[,3:17]
pt = t(protsy)
pt = as.data.frame(pt)

meTa <- read.table("nmds_tmts.txt", header = TRUE, row.names = 1, sep='\t')
little = meTa[7:15,]

# v.dist <- vegdist(pt, method="jaccard", binary=FALSE, diag=FALSE, upper=FALSE)

v.dist <- vegdist(pt, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE)


NMDS <- metaMDS(v.dist, distance = "euclidean", k = 2, trymax = 20, autotransform =TRUE,
                noshare = 0.1, wascores = TRUE, expand = TRUE, trace = 1,
                old.wa = FALSE)

plot(NMDS, type="t")
plot(NMDS$points)

ordiellipse(NMDS, meTa$Fraction, display = "sites", kind = "sd", label = T)


# then basic plot of nmds to check it worked
plot(NMDS, type="t")


library(ggplot2)


plot(NMDS$points, col = meTa$Soil)

data.scores <- as.data.frame(scores(NMDS)) # have a look at it and check ok

data.scores$Soil<- meTa$Soil
data.scores$Fraction<- meTa$Fraction



# and replot with prettier colours etc

library(ggplot2)
P <- ggplot(data.scores, aes(x= NMDS1, y= NMDS2, fill=Soil, shape=Fraction))+
  geom_point(size=8, stroke=1)  + # size and line thickness of plotted points
  scale_fill_manual(values = c("darkred", "olivedrab4", "cornflowerblue")) +
  scale_shape_manual(values=c(23, 25, 21))+
  theme_bw() + 
  theme(legend.key.size = unit(2.5,"line"),
        axis.text.x=element_text(size=14, colour='black'),
        axis.text.y=element_text(size=13, colour='black'),
        legend.text = element_text(size=13, colour='black'),
        legend.title = element_text(size=15, colour='black'))  +
guides(fill = guide_legend(override.aes=list(shape=21, size=10)),
       shape = guide_legend(override.aes = list(size=10))) 
               
P + annotate("text", x = -32, y = 30, label = "Stress = 0.02", size = 5) 


# and replot with prettier colours etc

library(ggplot2)
P <- ggplot(data.scores, aes(x= NMDS1, y= NMDS2, colour=Fraction))+
  geom_point(size=8, stroke=1, shape=21)  + # size and line thickness of plotted points
  scale_colour_manual(values = c("darkred", "olivedrab4", "cornflowerblue")) +
  theme_bw() + 
  theme(legend.key.size = unit(2.5,"line"),
        axis.text.x=element_text(size=14, colour='black'),
        axis.text.y=element_text(size=13, colour='black'),
        legend.text = element_text(size=13, colour='black'),
        legend.title = element_text(size=15, colour='black'))  +
  guides(fill = guide_legend(override.aes=list(shape=21, size=10)),
         shape = guide_legend(override.aes = list(size=10))) 

P + annotate("text", x = -32, y = 30, label = "Stress = 0.02", size = 5) 

## ======== rel abund plots ==========
library(dplyr)
library("phyloseq")



sharedsubFile = read.table('stability.opti_mcc.0.03.subsample.0.03.pick.shared')
sharedsubFile = t(sharedsubFile)
rownames(sharedsubFile) = sharedsubFile[,1]
colnames(sharedsubFile) = sharedsubFile[2,]
sharedsubFile = sharedsubFile[,2:7]
sharedsubFile = sharedsubFile[4:2083,]
class(sharedsubFile) <- "numeric"
head(sharedsubFile)

taxFile = read.table('stability.cons copy.taxonomy.txt', header=T, sep='\t')
rownames(taxFile) = taxFile[,1]
taxFile = taxFile[,2:8]
taxFile = as.matrix(taxFile)
head(taxFile)


OTU = otu_table(sharedsubFile, taxa_are_rows = TRUE)
TAX = tax_table(taxFile)
#META = sample_data(metaFile)
physeqSub = phyloseq(OTU, TAX)

OTURel = transform_sample_counts(physeqSub, function(x) x / sum(x) )

## readout phyloseq object, untrimmed
dat <- psmelt(OTURel)
write.csv(dat, file='16S_OTU_Coex.csv')

dat = read.csv("16S_OTU_Coex.csv")
dat[["Sample"]] <- setFactorOrder(dat[["Sample"]], c("DNA1", "DNA2", "DNA3", "cDNA1", "cDNA2", "cDNA3"))

m_df <- tbl_df(dat)

library(RColorBrewer)
# using colour brewer make a palette
mypal <- colorRampPalette( brewer.pal( 11 , "Set1" ) )
intercalate <- function(n){ #my crude attempt to shuffle the colors
  c(rbind(1:(n/2), n:(n/2+1))) #it will only work for even numbers
}

## now plot
# change as appropritae:
# x=
# y=
# fill = 
# facet_wrap(SamplingTime) # this is if you want subplots of tmt or day eg.

# family

ggplot(m_df, aes(x=Sample, y=Abundance, fill=Genus)) + 
  geom_bar(stat="identity") + 
 # facet_wrap(~SamplingTime, scales="free") + 
  scale_fill_manual( values = mypal(571)[intercalate(571)] )+
# change this 30 as appropriate for the number of phylum (or whatever)
# in bar plot. Put very high ~100 at first if you dont know. 
theme_bw() + ylab("Relative Abundance") + xlab("") +
  theme(axis.text.x=element_text(size=12, colour='black'),
        axis.text.y=element_text(size=11, colour='black'),
        axis.title.y=element_text(size=14, colour='black')) 

+ guides(fill=FALSE)



# ==== old bits
prots = read.table('Phylum_Table3Reps.txt', header=T, sep='\t')

protsy = prots[,2:10]
pt = t(protsy)
pt = as.data.frame(pt)
v.dist <- vegdist(pt, method="bray", binary=FALSE, diag=FALSE, upper=FALSE)

v.dist <- vegdist(pt, method="euclidean", binary=FALSE, diag=FALSE, upper=FALSE)


NMDS <- metaMDS(v.dist, distance = "euclidean", k = 2, trymax = 20, autotransform =TRUE,
                noshare = 0.1, wascores = TRUE, expand = TRUE, trace = 1,
                old.wa = FALSE)

# then basic plot of nmds to check it worked
plot(NMDS, type="t")



