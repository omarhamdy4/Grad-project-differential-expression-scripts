# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0, DESeq2 1.38.3
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
# load series and platform data from GEO

gset <- getGEO("GSE106817", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21263", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("00000000000000000000000000000000000000000000000000",
               "00000000000000000000000000000000000000000000000000",
               "000000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111111111111111",
               "11111111111111111111111111111111111111XXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]
# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }
# normalize data
exprs(gset) <- normalizeBetweenArrays(exprs(gset))

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("bc","no"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)
#check for NAs presence
sum(is.na(exprs(gset)))
# skip missing values
gset <- gset[complete.cases(exprs(gset)), ] 
# fit linear model
fit <- lmFit(gset, design)

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01, trend = T)

dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)
summary(dT)
tT <- topTable(fit2, adjust.method="fdr", sort.by="B", p.value = 0.05, lfc = 1, number = 1000)
tT$label <- ifelse(tT$logFC >= 1 & tT$adj.P.Val <= 0.05, "Up", 
                   ifelse(tT$logFC <= -1 & tT$adj.P.Val <= 0.05, "Down", "Not significant"))

#save tT
write.csv(tT, file = 'Top Table(1000)_GSE7trend.csv')

#filter UP and DOWN for 1000----
model = read.csv('Top Table(1000)_GSE7trend.csv')
str(model)
Uptrend <- subset(model, label == "Up")
write.csv(Uptrend, file = 'up-filtered from top table(1000)GSE7trend.csv')
downtrend <- subset(model, label == "Down")
write.csv(downtrend, file = 'down-filtered from top table(1000)GSE7trend.csv')
notsigtrend <- subset(model, label == "Not significant")
write.csv(notsigtrend, file = 'not significant-filtered from top table(1000)GSE7trend.csv')


#################################################################
# volcano plot (log P-value vs log fold change)
tT$log = -log10(tT$P.Value)
colnames(tT)[11] = "log10(P.value)"
tT3 = tT[tT$`log10(P.value)` != Inf,]
ggplot(data = tT3, aes(x = logFC, y = `log10(P.value)` , col = label)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("blue", "black", "red")) +
  theme(text = element_text(size = 20)) +
  labs(title = "Top table volcano Plot")
#################################################################
#ROC Curve
  # Generate the 1,0 from the expression set
x <- select(gset@phenoData@data, "description")
colnames(x) <- "disease_state"
x <- mutate(x,
            value=case_when(
              disease_state =="Breast Cancer" ~ "1",
              disease_state =="non-Cancer" ~ "0",
              TRUE ~ "X" #this replaces the non-mentioned other values to be x
            ))
x$value <- as.numeric(x$value)

  # convert mirbase accession to mirbase ids to generate the expression data for each id
    #BiocManager::install("miRBaseConverter",force = TRUE)
    # install with install.packages("pROC")
library(miRBaseConverter)
Nx2 <- t(exprs(gset))
mirbase_acc <- colnames(Nx2)
mirbase_ids <- miRNA_AccessionToName(mirbase_acc,targetVersion = "v22")
colnames(Nx2) <- mirbase_ids[,2]

CallMe <- function(q) {
  x$ex <- Nx2[,q]
  glm.fit=glm(value~ex,data=x, family=binomial)
  library(pROC)
  par(pty="s")
  roc(x$value, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE,
      xlab="100-Specificity",ylab="Sensitivity", col="#377eb8",lwd=3, print.auc=TRUE, print.auc.x=40)    
  legend("bottomright" ,legend=c(q), col=c("#377eb8"), lwd=4)
}
CallMe("hsa-miR-155-5p")
CallMe("hsa-miR-335-5p")
