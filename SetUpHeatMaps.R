mouseFeatureLab <- read.delim2("C:/Users/shivang/Desktop/mouse.matrix.out", row.names=1)
zfaFeatureLab <- read.csv("C:/Users/shivang/Desktop/xenopus_matrix.out", row.names=1)
xenopusFeaturesLab <-read.csv("C:/Users/shivang/Desktop/zfa.matrix.out", row.names=1)

mouseFeatureLab <- read.delim2("//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/ShivangSingh/AlignmentTool Project/mouse.matrix.out", row.names=1)
zfaFeatureLab<- read.csv("//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/ShivangSingh/AlignmentTool Project/zfa.matrix.out", row.names=1)
xenopusFeaturesLab <- read.csv("//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/ShivangSingh/AlignmentTool Project/XenopusExcel.txt", row.names=1)

mus.sample.table <- read.csv("Z:/Shared/Undergraduate_Students/ShivangSingh/Sample Data/Mus musculus/Mmus.sample.table.3.csv")
Mus.sp.OGG.exp<- read.csv("//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/ShivangSingh/AlignmentTool Project/OGG Filtered/MicroArray Filtered/Mus.sp.OGG.exp.ma", row.names=1, sep="")
drer.sample.table <- read.csv("//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/ShivangSingh/Sample Data/Danio rerio/drer.sample.table.csv")
Dan.sp.OGG.exp <- read.csv("//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/ShivangSingh/AlignmentTool Project/OGG Filtered/MicroArray Filtered/Dan.sp.OGG.exp.ma", row.names=1, sep="")
xtrop.sample.table <- read.csv("//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/ShivangSingh/Sample Data/Xenopus tropicalis/xtrop.sample.table.csv", row.names=1)
Xtrop.sp.OGG.exp <- read.csv("//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/ShivangSingh/AlignmentTool Project/OGG Filtered/MicroArray Filtered/Xtrop.sp.OGG.exp.ma", row.names=1, sep="")

source("https://bioconductor.org/biocLite.R")
if (!require("data.table", character.only=T, quietly=T)) biocLite("data.table")
library(data.table, lib.loc="~/R/R-3.2.1/library")
if (!require("ggplot2", character.only=T, quietly=T)) biocLite("ggplot2")
library(ggplot2, lib.loc="~/R/R-3.2.1/library")
source('//files.ccbb.utexas.edu/hofmannLab/Shared/Undergraduate_Students/Pranav_Bhamidipati/Rscripts/MA_lib-AvgParalogousGeneFIX.r')

install.packages("flexclust")
library("flexclust")

install.packages( pkgs= "gplots" )
library("gplots")

rownames(Dan.shared.OGG.exp) <- drer.sample.table$time.point
Dan.shared.OGG.exp <- Dan.filtered
dat.ma<- Dan.shared.OGG.exp
#clears the data up
dat <- dat.ma
time.points <- drer.sample.table$time.point
# Import necessary methods from data.table
rbindlist <- data.table::rbindlist

dat.list <- list()
for (time.point in unique(as.character(time.points))) {
  replicate_chips <- time.points == time.point
  print(time.point)
  dat.list[[time.point]] <- as.data.frame(t(colMeans(dat[replicate_chips,])))
}

# Concatenate data.frame rows from dat.list to make the data frame dat.stages.
dat.stages <- as.data.frame(rbindlist(dat.list ))
rownames(dat.stages) <- names(dat.list)
dat.ma <- dat.stages
Dan.shared.OGG.exp <- dat.ma

#Publishes the shared matrices
write.table(Xtrop.shared.OGG.exp, "C:/Users/ss74244/Desktop/Xtrop.filtered.Averaged" ,sep="\t") 
write.table(Mus.shared.OGG.exp, "C:/Users/ss74244/Desktop/Mus.filtered.Averaged" ,sep="\t") 
write.table(Dan.shared.OGG.exp, "C:/Users/ss74244/Desktop/Dan.filtered.Averaged" ,sep="\t") 

write.table(t(zfa.obo), "C:/Users/nidhi/Desktop/zfa.obo" ,sep="\t") 
write.table(t(mouse.obo), "C:/Users/nidhi/Desktop/mouse.obo" ,sep="\t") 
write.table(t(xenopus_anatomy), "C:/Users/nidhi/Desktop/xenopus.obo" ,sep="\t") 


#Inverts the matrices
Xtrop.sp.OGG.exp <- t(Xtrop.sp.OGG.exp)
Mus.sp.OGG.exp <- t(Mus.sp.OGG.exp)
Dan.sp.OGG.exp <- t(Dan.sp.OGG.exp)

#Sorts the stages 
Xtrop.shared.OGG.exp <- Xtrop.shared.OGG.exp[order(as.numeric(rownames(Xtrop.shared.OGG.exp))),,drop=FALSE]
Dan.shared.OGG.exp <- Dan.shared.OGG.exp[order(as.character(rownames(Dan.shared.OGG.exp))),,drop=FALSE]

rownames(Xtrop.shared.OGG.exp) <- c('S02', 'S08', 'S09' , 'S10' , 'S12', 'S13', 'S14', 'S16', 'S18', 'S20', 'S23', 'S25', 'S30', 'S33', 'S40')
rownames(Dan.shared.OGG.exp) <- c('01c', '09h', '10h', '11.7h', '16h', '24h' , '36h', '48h', 'd04', 'd05', '06h', '08h', 'd14', 'd30', 'd90')

#Generates Red heatmaps for the eucleidean distance matrices
XtropDan_heatmap <- heatmap.2(dist2(Xtrop.shared.OGG.exp, Dan.shared.OGG.exp , "euclidean"),main = "D. rer vs. X. trop Genetic based Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256) , scale="column", margins=c(5,10))
DanMus_heatmap <- heatmap.2(dist2(Dan.shared.OGG.exp , Mus.shared.OGG.exp , "euclidean"),main = "M. mus vs. D. rer Genetic based Stage Comparison", density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
XtropMus_heatmap <- heatmap.2(dist2(Xtrop.shared.OGG.exp, Mus.shared.OGG.exp , "euclidean"),main = "M. mus vs. X. trop Genetic based Stage Comparison", density.info="none", trace="none", Rowv=NA, Colv=NA,  col = heat.colors(256), scale="column", margins=c(5,10))



Dan_heatmap <- heatmap.2(dist2(Dan.shared.OGG.exp, Dan.shared.OGG.exp , "euclidean"),main = "D. rer vs. D.rer Genetic based Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256) , scale="column", margins=c(5,10))
Mus_heatmap <- heatmap.2(dist2(Mus.shared.OGG.exp , Mus.shared.OGG.exp , "euclidean"),main = "M. mus vs. M. mus Genetic based Stage Comparison", density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
Xtrop_heatmap <- heatmap.2(dist2(Xtrop.shared.OGG.exp, Xtrop.shared.OGG.exp , "euclidean"),main = "X. trop vs. X. trop Genetic based Stage Comparison", density.info="none", trace="none", Rowv=NA, Colv=NA,  col = heat.colors(256), scale="column", margins=c(5,10))

nrow(t(Dan.shared.OGG.exp))
nrow(t(Mus.shared.OGG.exp))

kmeans(Xtrop.shared.OGG.exp)

xenopus_anatomy.obo$X <- NULL
mouse.obo$X <- NULL
zfa.obo$X <- NULL

zfa.obo <- t(zfa.obo)
mouse.obo <- t(mouse.obo)
xenopus_anatomy.obo <- t(xenopus_anatomy.obo)

XtropDan_heatmap <- heatmap.2(dist2(xenopus_anatomy.obo , zfa.obo , "euclidean"),main = "D. rer vs. X. trop Genetic based Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256) , scale="column", margins=c(5,10))
DanMus_heatmap <- heatmap.2(dist2(zfa.obo , mouse.obo , "euclidean"),main = "M. mus vs. D. rer Genetic based Stage Comparison", density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
XtropMus_heatmap <- heatmap.2(dist2(xenopus_anatomy.obo, mouse.obo , "euclidean"),main = "M. mus vs. X. trop Genetic based Stage Comparison", density.info="none", trace="none", Rowv=NA, Colv=NA,  col = heat.colors(256), scale="column", margins=c(5,10))



xenopusFeaturesLab$X <- NULL
zfaFeatureLab$X <- NULL

#Alphabatizes by rownames
CommonClassMouse <- CommonClassMouse[ order(row.names(CommonClassMouse)), ]
CommonClassesFish  <- CommonClassesFish[ order(row.names(CommonClassesFish)), ]
CommonClassesXeno  <- CommonClassesXeno[ order(row.names(CommonClassesXeno)), ]


CommonClassesXeno <- t(CommonClassesXeno)
CommonClassesFish <- t(CommonClassesFish)
CommonClassMouse <- t(CommonClassMouse)


zfaFeatureLab <- t(zfaFeatureLab)
xenopusFeaturesLab <- t(xenopusFeaturesLab)
mouseFeatureLab <- t(mouseFeatureLab)

write.csv(xenopus_anatomy,file="E:/independent research/xenopus_anatomy.obo")
write.csv(t(CommonClassesXeno), file="F:/CommonClassesXeno")
write.text(vegdist(rbind(unlist(mouse.obo, use.names=F),unlist(xenopus_anatomy, use.names=F)), method = "jaccard"),file="E:/independent research/mouse.obo")

colnames(xenopus_anatomy)
ncol(CommonClassMouse)
ncol(CommonClassesFish)
ncol(CommonClassesXeno)

mouse.obo <- t(mouse.obo)
xenopus_anatomy.obo <- t(xenopus_anatomy.obo)
zfa.obo <- t(zfa.obo)

ncol(zfa.obo)
ncol(xenopus_anatomy)
ncol(mouse.obo)


install.packages("proxy")

library("proxy")

library("vegdist")
dist(mouse.obo, xenopus_anatomy, method = "Jaccard")

vegdist(rbind(unlist(mouse.obo, use.names=F),unlist(xenopus_anatomy, use.names=F)), method = "jaccard")

# Generates heatmaps

XtropDanFeatures_heatmap <- heatmap.2(dist(mouse.obo, xenopus_anatomy, method = "jaccard"),main = "M. mus vs. D. rer feature (SES-based) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
DanMusFeatures_heatmap <- heatmap.2(dist(xenopus_anatomy, zfa.obo , "Jaccard"),main = "M. mus vs. X. trop feature (SES-based) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
XtropMusFeatures_heatmap <- heatmap.2(dist(zfa.obo, xenopus_anatomy , "Jaccard"),main = "D. rer vs. X. trop feature (SES-based) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA,  col = heat.colors(256), scale="column", margins=c(5,10))

XtropDanFeatures_heatmap <- heatmap.2(dist(CommonClassesXeno, CommonClassesFish , "Jaccard"),main = "D. rer vs. X. trop feature (Common-Features) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
DanMusFeatures_heatmap <- heatmap.2(dist(CommonClassMouse , CommonClassesFish , "Jaccard"), density.info="none" ,main = "M. mus vs. D. rer feature (Common-Features) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
XtropMusFeatures_heatmap <- heatmap.2(dist(CommonClassMouse, CommonClassesXeno , "Jaccard"), main = "M. mus vs. X. trop feature (Common-Features) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA,  col = heat.colors(256), scale="column", margins=c(5,10))

DanFeatures_heatmap <- heatmap.2(dist(CommonClassMouse, CommonClassesFish , method = "Jaccard"),main = "D. rer vs. D. rer Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
MusFeatures_heatmap <- heatmap.2(dist(CommonClassMouse , CommonClassMouse , method = "Jaccard"), main = "M. mus vs. M. mus Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))
XtropFeatures_heatmap <- heatmap.2(dist(TempoMouse, Style1Xeno , method = "Jaccard"),main = "M. mus vs. X. trop feature (Common-Features) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA,  col = heat.colors(256), scale="column", margins=c(5,10))
XtropFeatures_heatmap <- heatmap.2(dist(TempoMouse, TempoZfa , method = "Jaccard"),main = "M. mus vs. D. rer feature (Common-Features) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA,  col = heat.colors(256), scale="column", margins=c(5,10))
XtropFeatures_heatmap <- heatmap.2(dist(Style1Xeno, TempoZfa , method = "Jaccard"),main = "X. trop vs. D. rer feature (Common-Features) Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA,  col = heat.colors(256), scale="column", margins=c(5,10))
MusFeatures_heatmap <- heatmap.2(dist(Style1Xeno , TempoMouse , method = "Jaccard"), main = "X.trop vs. M. mus Stage Comparison" , density.info="none", trace="none", Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))

#K - means clustering 


distance <- as.dist(dist(Xeno, Zfa.adjust.into.29, method = "Jaccard")
)

plot(hclust(distance))

X <- t(Dan.shared.OGG.exp)
Z <- t(Mus.shared.OGG.exp)
Y <- t(Xtrop.shared.OGG.exp)

cool <- rbind(t(X),t(Z), t(Y))

plot(hclust(dist(cool)))

#Generates the PC analysis

library(cluster)

clusplot(cool, fit$cluster, color=TRUE, shade=TRUE,
         labels=2, lines=0)


ncol(Style1Xeno)
ncol(TempoMouse)
ncol(TempoZfa)

Style1Xeno <- t(AppearDisappearXenopus.out)
TempoMouse <- t(Lezzz)
TempoZfa <- t(zfaNew)

Lezzz <- t(TempoMouse)
t <- zfaNew[order(rownames(zfaNew)),] 

# Heat maps


XtropFeatures_heatmap <- heatmap.2(dist(Style1Xeno, TempoMouse, method = "jaccard"), Rowv=FALSE, Colv=FALSE, 
          dendrogram="none", heat.colors(256), 
          key=T, keysize=1.5, density.info="none", 
          trace="none", labRow=NA)


hist(unlist(dist2(CommonClassesXeno, CommonClassesFish , method = "euclidean")), breaks = 100)
