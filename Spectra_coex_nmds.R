
spec = read.csv('SoilNMDS.csv', header=T, sep=',')
library(vegan)
dim(spec)
all = spec[,c(3:11)]

pt = t(all)
pt = as.data.frame(pt)

meTa <- read.table("nmds_tmts.txt", header = TRUE, row.names = 1, sep='\t')
little = meTa[7:15,]

v.dist <- vegdist(pt, method="raup", binary=TRUE, diag=FALSE, upper=FALSE)

v.dist <- vegdist(pt, method="chao", binary=TRUE, diag=FALSE, upper=FALSE)


NMDS <- metaMDS(v.dist, distance = "raup", k = 2, trymax = 20, autotransform =TRUE,
                noshare = 0.1, wascores = TRUE, expand = TRUE, trace = 1,
                old.wa = FALSE)

plot(NMDS, type="t")


meTa <- read.table("treatments_nmds2.txt", header = TRUE, row.names = 1, sep='\t')


dune.ano <- anosim(v.dist, meTa$soil)
summary(dune.ano)
plot(dune.ano)



data.scores <- as.data.frame(scores(NMDS)) # have a look at it and check ok
data.scores$Soil<- meTa$soil
data.scores$Soil<- as.factor(data.scores$Soil)

library(ggplot2)

P <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2, fill=Soil))+
  geom_point(size=9, stroke=1, shape=21)  + # size and line thickness of plotted points
  scale_fill_manual(values = c("darkred", "olivedrab4", "cornflowerblue")) +
  theme_bw() +
  theme(legend.key.size = unit(2.5,"line"),
        axis.text.x=element_text(size=14, colour='black'),
        axis.text.y=element_text(size=13, colour='black'),
        legend.text = element_text(size=13, colour='black'),
        legend.title = element_text(size=15, colour='black')) +
  guides(fill = guide_legend(override.aes=list(shape=21, size=12))) 
P


P + annotate("text", x = -0.34, y = 0.25, label = "2D Stress = 0.03", size = 5)

# SpectraBasednMDS_3soils
