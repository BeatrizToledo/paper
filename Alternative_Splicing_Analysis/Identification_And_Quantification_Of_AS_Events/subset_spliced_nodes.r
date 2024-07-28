# load libraries
library(dplyr)
library(stringr)

# read the Whippet PSI files
L222.whippet.psi <- read.delim("/../L222.whippet.psi.gz")
L223.whippet.psi <- read.delim("/../L223.whippet.psi.gz")
L224.whippet.psi <- read.delim("/../L224.whippet.psi.gz")
L355.whippet.psi <- read.delim("/../L355.whippet.psi.gz")
L356.whippet.psi <- read.delim("/../L356.whippet.psi.gz")
L357.whippet.psi <- read.delim("/../L357.whippet.psi.gz")
L393.whippet.psi <- read.delim("/../L393.whippet.psi.gz")
L394.whippet.psi <- read.delim("/../L394.whippet.psi.gz")
L395.whippet.psi <- read.delim("/../L395.whippet.psi.gz")

#get a unique identifier for each node (useful for later processing and fusing of nodes into a single AS event)
L222.whippet.psi$name<-gsub(" ", "",paste(L222.whippet.psi$Gene, L222.whippet.psi$Node, sep="-"))
L223.whippet.psi$name<-gsub(" ", "",paste(L223.whippet.psi$Gene, L223.whippet.psi$Node, sep="-"))
L224.whippet.psi$name<-gsub(" ", "",paste(L224.whippet.psi$Gene, L224.whippet.psi$Node, sep="-"))
L355.whippet.psi$name<-gsub(" ", "",paste(L355.whippet.psi$Gene, L355.whippet.psi$Node, sep="-"))
L356.whippet.psi$name<-gsub(" ", "",paste(L356.whippet.psi$Gene, L356.whippet.psi$Node, sep="-"))
L357.whippet.psi$name<-gsub(" ", "",paste(L357.whippet.psi$Gene, L357.whippet.psi$Node, sep="-"))
L393.whippet.psi$name<-gsub(" ", "",paste(L393.whippet.psi$Gene, L393.whippet.psi$Node, sep="-"))
L394.whippet.psi$name<-gsub(" ", "",paste(L394.whippet.psi$Gene, L394.whippet.psi$Node, sep="-"))
L395.whippet.psi$name<-gsub(" ", "",paste(L395.whippet.psi$Gene, L395.whippet.psi$Node, sep="-"))

#select a narrow confidence interval 
L222.c<-subset(L222.whippet.psi, CI_Width<=0.3)
L223.c<-subset(L223.whippet.psi, CI_Width<=0.3)
L224.c<-subset(L224.whippet.psi, CI_Width<=0.3)
L355.c<-subset(L355.whippet.psi, CI_Width<=0.3)
L356.c<-subset(L356.whippet.psi, CI_Width<=0.3)
L357.c<-subset(L357.whippet.psi, CI_Width<=0.3)
L393.c<-subset(L393.whippet.psi, CI_Width<=0.3)
L394.c<-subset(L394.whippet.psi, CI_Width<=0.3)
L395.c<-subset(L395.whippet.psi, CI_Width<=0.3)

PP.DP.DPSI.diff <- read.csv("/../PP.DP.DPSI.diff.gz", sep="")
PP.DP.DPSI.diff$name<-gsub(" ", "",paste(PP.DP.DPSI.diff$Gene, PP.DP.DPSI.diff$Node, sep="-"))

Whippet.expressed.genes.at.05TpM <- read.delim("~/Desktop/Whippet_corrected_analysis/Whippet.expressed.genes.at.05TpM.txt")
PP.DP.expressed<-subset(Whippet.expressed.genes.at.05TpM, TpM_L222>=0.5 & TpM_L355>=0.5 & TpM_L393>=0.5 & TpM_L223>=0.5 & TpM_L356>=0.5 & TpM_L394>=0.5)
expressed<-subset(PP.DP.DPSI.diff, PP.DP.DPSI.diff$Gene %in% PP.DP.expressed$Gene)
confidence<-subset(expressed, expressed$name %in% L222.c$name & expressed$name %in% L355.c$name & expressed$name %in% L393.c$name & expressed$name %in% L223.c$name & expressed$name %in% L356.c$name & expressed$name %in% L394.c$name)
whippet.spliced<-subset(confidence, Probability>=0.9 & DeltaPsi >=0.1 | Probability>=0.9 & DeltaPsi <=-0.1)
whip.act.spliced.PP.DP<-subset(whippet.spliced, Type=="CE" | Type=="AA" | Type=="AD" | Type=="RI")
write.table(whip.act.spliced.PP.DP, "Pb.stringtie.corrected.PP.DP.txt",sep="\t", row.names = FALSE, col.names=TRUE, quote=FALSE)


DP.N.DPSI.diff <- read.csv("/../DP.N.DPSI.diff.gz", sep="")
DP.N.DPSI.diff$name<-gsub(" ", "",paste(DP.N.DPSI.diff$Gene, DP.N.DPSI.diff$Node, sep="-"))
DP.N.expressed<-subset(Whippet.expressed.genes.at.05TpM, TpM_L224>=0.5 & TpM_L357>=0.5 & TpM_L395>=0.5 & TpM_L223>=0.5 & TpM_L356>=0.5 & TpM_L394>=0.5)
expressed2<-subset(DP.N.DPSI.diff, DP.N.DPSI.diff$Gene %in% DP.N.expressed$Gene)
confidence2<-subset(DP.N.DPSI.diff, DP.N.DPSI.diff$name %in% L223.c$name & DP.N.DPSI.diff$name %in% L356.c$name & DP.N.DPSI.diff$name %in% L394.c$name & DP.N.DPSI.diff$name %in% L224.c$name & DP.N.DPSI.diff$name %in% L357.c$name & DP.N.DPSI.diff$name %in% L395.c$name)
whippet.spliced2<-subset(confidence2, Probability>=0.9 & DeltaPsi >=0.1 | Probability>=0.9 & DeltaPsi <=-0.1)
whip.act.spliced.DP.N<-subset(whippet.spliced2, Type=="CE" | Type=="AA" | Type=="AD" | Type=="RI")
write.table(whip.act.spliced.DP.N, "Pb.stringtie.corrected.DP.N.txt",sep="\t", row.names = FALSE, col.names=TRUE, quote=FALSE)



