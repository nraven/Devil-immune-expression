##Code for manuscript titled 'Ecoimmunology of an endangered marsupial, the Tasmanian devil (Sarcophilus harrisii)

##Note: Not every detail of the script has been included in this document.
##Gene expression ratios were created and data was log transformed to improve normality fit.
##Gene expression matrix based on script from Prof Marcel Klaassen


library(MuMIn) ##model selection
library(gplots)  ##heatmap
library(plotrix)  ##heatmap legend
library(RColorBrewer) ##heatmap colours
library(ggeffects) ##marginal means
library(ggplot2)  ##marginal means plots
library(patchwork) ##plot layout



##Removing missing values from the explainatory variables
WPPall$OK <- 1
WPPall$OK[which(is.na(WPPall$Age)|is.na(WPPall$Sex)|is.na(WPPall$DFTD)|is.na(WPPall$Season)|is.na(WPPall$Total.ticks))]<- 0
WPPdf <- WPPall[which(WPPall$OK==1),]


###Remove any missing variables from individual genes.  
WPPdf$OK <- 1
WPPdf$OK[which(is.na(WPPdf$IgAL))]<- 0
WPPdf_IgA <- WPPdf[which(WPPdf$OK==1),]

##Create the full model. This script was written as I was learning R, so rather than looping, the script was
##copied and pasted to run each gene (I know, bad practice). 
##Will only paste one gene as an example of the scripts used.
ModIgA <- lm(IgAL~scale(AgeL) + Sex + DFTD + Season + scale(Total.ticks) + 
                DFTD:Sex + DFTD:scale(Total.ticks) + DFTD:scale(AgeL) + DFTD:Season +
                Sex:scale(Total.ticks) + Sex:scale(AgeL) + Sex:Season + scale(Total.ticks):Season + 
                scale(AgeL):scale(Total.ticks), data=WPPdf_IgA)

##Run the modelselection
#setting maximum number of expl variables at 8 (10 samples to one explainatory variable)
options(na.action = "na.fail")
results_IgA <- dredge(ModIgA, trace=2, beta=FALSE, evaluate=TRUE, m.lim= c(1,8)) 
#subset results to top 2 AICc
top_IgA <- subset(results_IgA, delta <= 2, recalc.weights=TRUE)
#Extract the importance values (model weights/variable)
Import_top_IgA <- importance(top_IgA)


#Extract and format the Importance values
Import_top_IgAd <- as.data.frame(t(as.data.frame(Import_top_IgA)))
Import_top_IgAd$gene <- "IgA"

#merge all the genes
Import_top <- merge(Import_top_IgAd, Import_top_IgEd, all=TRUE)

#formating the data for the heatmap
Import3 <- Import_top
Import3[is.na(Import3)] <- 0
row.names(Import3) = Import3$gene
Import3 <- Import3[,-3] #remove the 'gene' column, it is sometimes a different column
Import4 <- as.data.frame(t(Import3))

Import5 <- as.matrix(sapply(Import4, as.numeric)) 
row.names(Import5) =  rownames(Import4)


##Imputs for the heatmap
col <- colorRampPalette(brewer.pal(9, "YlOrRd"))(256)
fon <- round(Import5, 2)

##create the heatmap
pdf("Importance_heatmap.pdf", height = 7.5, width = 8.5, pointsize = 12)
heatmap.2(Import5, cellnote = fon, notecex=1.0, 
          notecol="black", scale = "none", main = "Importance values",
          dendrogram = "column", col = col, margins = c(5,14), key=F, lhei = c(2,7),
          keysize = 1.2, key.title = "Key",
          trace = "none", density.info = "none")
gradient.rect(0.115,0.25,0.135,0.65,border=F,gradient="y",col=smoothColors(col))
text(x=rep(0.108,7),y=seq(0.28,0.62,by=0.34),adj=1,cex=1,labels=c("0","1"))
dev.off()



###I assume there is a way to extract the importance values and import them into new models.
##Again, as I had no idea how to do this, I remade the new models manually.

##This is to make sure the seasons were in the right order in the marginal mean graphs and that healthy was the
##intercept value, it made the second heatmap more logical
WPPall$Season <- factor(WPPall$Season, levels= c("Summer","Autumn","Winter","Spring"))
WPPall$DFTD <- factor(WPPall$DFTD, levels= c("Healthy","Disease"))


##new models
ModIgG <- lm(IgGL~DFTD+scale(log(Age))+Sex+scale(log(Age)):Sex+scale(Total.ticks)+
               scale(Total.ticks):Sex+DFTD:Sex, data = WPPall)

ModIgM <- lm(IgML~DFTD+scale(log(Age))+Sex+Season+scale(log(Age)):Sex+
               DFTD:Sex+DFTD:Season, data = WPPall)

ModIgE <- lm(IgEL~DFTD+Season+Sex+DFTD:Sex+DFTD:Season, data = WPPall)

ModIgA <- lm(IgAL~DFTD+scale(log(Age))+Season+DFTD:scale(log(Age)), 
             data = WPPall)

ModCD4 <- lm(CD4L~Season, data = WPPall)

ModCD8 <- lm(CD8L~DFTD+scale(log(Age))+Season+Sex+scale(log(Age)):Sex+
               scale(Total.ticks), data = WPPall)

ModCD11 <- lm(CD11L~DFTD+scale(log(Age))+Season, data = WPPall)

ModCD16 <- lm(CD16L~DFTD+scale(log(Age))+Season+
                DFTD:scale(log(Age)), data = WPPall)

ModIgMG <- lm(IgM.IgGL~DFTD+scale(log(Age))+Season+scale(Total.ticks)+DFTD:Season, data = WPPall)

ModCD48 <- lm(CD4.CD8L~DFTD+scale(log(Age))+Season+Sex+scale(log(Age)):Sex+
                scale(Total.ticks), data = WPPall)

ModNKG2D <- lm(NKG2DL~DFTD+scale(log(Age))+Season+Sex+scale(log(Age)):Sex+
                 DFTD:scale(log(Age)), data = WPPall)

ModMHCC2 <- lm(MHCC2L~DFTD+Season, data = WPPall)


##To extract the coefficients from the model, this code extracts the p-values, but change the '4' in the third line
##to a '1' and you  get the estimates.
IgAP <- as.data.frame(summary(ModIgA)$coefficients)
IgAP2 <- IgAP
IgAP3 <- IgAP2[,4]
IgAP3 <- as.data.frame(IgAP3)
rownames(IgAP3) = rownames(IgAP)
IgAP3 <- t(IgAP3)
IgAP4 <- as.data.frame(IgAP3)
IgAP4$gene <- "IgA"

##The code chunk gets everything in the correct format that it can be merged into a single dataset, just as it did 
##in the previous heatmap. Make one dataset with the p-values and convert the values to the corresponding significant 
##sympbols (*, **, ***). Make a second dataset with the estimates.

##Again, both datasets need to be formatted for the heatmap 
Est2 <- Est2[,-1]
row.names(Est2) = Est2$gene
Est2 <- Est2[,-1]

Est_2 <- as.data.frame(t(Est2))
Est_3 <- as.matrix(sapply(Est_2, as.numeric)) 
row.names(Est_3) =  rownames(Est_2)

Pval <- Pval[,-1]
row.names(Pval) = Pval$gene
Pval <- Pval[,-1]
Pval2 <- as.data.frame(t(Pval))

##Both datasets need to be ordered the same way, otherwise there will be a mismatch
Pvalue <- Pval2
Est_3 <- Est_3[,order(colnames(Est_3))]
Pvalue <- Pvalue[,order(colnames(Pvalue))]

##Because there are negative and positive values, a new 'col' variable is created
col <- rev(colorRampPalette(brewer.pal(10, "RdYlBu"))(256))

##Create the heatmap
pdf("Estimate_heatmap.pdf", height = 7.5, width = 8.5)
heatmap.2(Est_3, Rowv = "False", cellnote = Pvalue, notecol = "black",  notecex=1.5,
          scale = "none", dendrogram = "column", main = "Model slope",
          col = col, margins = c(5,14), key=F, lhei = c(2,7), cexRow = 1,
          keysize = 1, key.title = "Key",
          trace = "none", density.info = "none")
gradient.rect(0.080,0.25,0.105,0.65,border=F,gradient="y",col=smoothColors(col))
text(x=rep(0.075,7),y=seq(0.28,0.62,by=0.17),adj=1,cex=0.8,labels=c("-2.4","0","2.4"))
dev.off()


###To create the marginal means plots, only significant values are used. The ggeffects package can be combined
##with ggplot2 to customise the graphs.


###This the code for season as a significant main effect.
z<-plot(ggpredict(ModMHCC2, terms = "Season"), colors = "reefs", limits = c(-3.7, 2), add.data = T)+
  labs(title = "MHC-II", x="Season", y="MHC-II expression (log)")
a<-plot(ggpredict(ModCD4, terms = "Season"), colors = "reefs",limits = c(-3.7, 2), add.data = T)+
  labs(title = "CD4", x="Season", y="CD4 expression (log)")
b<-plot(ggpredict(ModCD11, terms = "Season"), colors = "reefs", limits = c(-3.7, 2), add.data = T)+
  labs(title = "CD11", x="Season", y="CD11 expression (log)")
c<-plot(ggpredict(ModIgA, terms = "Season"), colors = "reefs", limits = c(-3.7, 2), add.data = T)+
  labs(title = "IgA", x="Season", y="IgA expression (log)")
d<-plot(ggpredict(ModIgM, terms = "Season"), colors = "reefs", limits = c(-3.7, 2), add.data = T)+
  labs(title = "IgM", x="Season", y="IgM expression (log)")
e<-plot(ggpredict(ModCD8, terms = "Season"), colors = "reefs", limits = c(-3.7, 2), add.data = T)+
  labs(title = "CD8", x="Season", y="CD8 expression (log)")
f<-plot(ggpredict(ModCD48, terms ="Season"), colors = "reefs", limits = c(-3.7, 2), add.data = T)+
  labs(title = "CD4:CD8", x="Season", y="CD4:CD8 expression (log)")
g<-plot(ggpredict(ModNKG2D, terms ="Season"), colors = "reefs", limits = c(-3.7, 2), add.data = T)+
  labs(title = "NKG2D", x="Season", y="NKG2D expression (log)")
h<-plot(ggpredict(ModCD16, terms ="Season"), colors = "reefs", limits = c(-3.7, 2), add.data = T)+
  labs(title = "CD16", x="Season", y="CD16 expression (log)")

pdf("Season_MMplots.pdf", width = 10, height = 10)
(a|f|z)/
  (e|b|h)/
  (c|g|d)
dev.off()






