library(GEOquery)
library(limma)
library(umap)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)
#set your directory
# read processed data downloaded from GEO database
data = read.delim('GSE211692_processed_data_matrix.txt', header = T)
#Adusting the data dimensions
  # get metadata
gse= getGEO(GEO = "GSE211692",GSEMatrix = TRUE)

metadata = pData(phenoData(gse[[1]]))
metadata.subset = select(metadata,c(1,10))

  # organize
metadata.subset = metadata %>%
  select(1,10) %>%
  rename(state = characteristics_ch1 ) %>%
  mutate(state = gsub("disease state: ","", state))%>%
  mutate(state = gsub("disease in the ","", state))


  # reshaping the data
dat.long = data %>%
  gather(key= "sample", value = "FPKM", -ID_REF)


# join 
dat.long = dat.long %>%
  left_join(.,metadata.subset,by = c("sample" = "title"))

#filter
bcvsno = dat.long %>%
  filter(state == "breast cancer" | state == "no cancer" )
#group_by(state)

#head(bcvsno)
bcvsno$sample <- paste(bcvsno$sample,bcvsno$state)
bcvsno = bcvsno[,-4]
bc= bcvsno %>% 
  spread(., key = sample, value = FPKM)
#head(bc)
row.names(bc) = bc$ID_REF
bc = bc[,-1]
# log2 transformation
qx <- as.numeric(quantile(bc, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { bc[which(bc <= 0)] <- NaN
exprs(bc) <- log2(bc) }
# force normalization
Nx <- normalizeBetweenArrays(bc) # normalize data
#check for NAs presence
sum(is.na(Nx))
pData.gse = pData(gse[[1]])
#head(pData.gse)
#pData.gse$`disease state:ch1`
GR = pData.gse$`disease state:ch1`
FG = factor(GR, levels = c('breast cancer','no cancer'))
designM = model.matrix(~0+FG)
colnames(designM) = c('breast','no')
contrast.matrix = makeContrasts(breast - no, levels = designM)
fit <- lmFit(Nx, designM)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, 0.01, trend=TRUE)
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=1)
summary(dT)
tT <- topTable(fit2, adjust.method="fdr", sort.by="B", p.value = 0.05, lfc = 1, number = 1000)
tT$label <- ifelse(tT$logFC >= 1 & tT$adj.P.Val <= 0.05, "Up", 
                   ifelse(tT$logFC <= -1 & tT$adj.P.Val <= 0.05, "Down", "Not significant"))

#save tT
write.csv(tT, file = 'Top Table(1000)_GSE6trend.csv')

#filter UP and DOWN for 1000----
model = read.csv('Top Table(1000)_GSE6trend.csv')
Uptrend <- subset(model, label == "Up")
write.csv(Uptrend, file = 'Up-filtered from top table(1000)GSE6trend.csv')
downtrend <- subset(model, label == "Down")
write.csv(downtrend, file = 'Down-filtered from top table(1000)GSE6trend.csv')
notsigtrend <- subset(model, label == "Not significant")
write.csv(notsigtrend, file = 'Not significant-filtered from top table(1000)GSE6trend.csv')

#################################################################
# volcano plot (log P-value vs log fold change)
tT$log = -log10(tT$P.Value)
colnames(tT)[8] = "log10(P.value)"
tT3 = tT[tT$`log10(P.value)` != Inf,]
ggplot(data = tT3, aes(x = logFC, y = `log10(P.value)` , col = label)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("blue", "black", "red")) +
  theme(text = element_text(size = 20)) +
  labs(title = "Top table volcano Plot")

#################################################################
#ROC Curve can be generated from script 2 since in this script we couldn't generate equal dimensions for the pheno data versus the normalized data to plot the ROC curve
gse1 <- gse[[1]]
x <- select(gse1@phenoData@data, "disease state:ch1")
colnames(x) <- "disease_state"
x <- mutate(x,
            value=case_when(
              disease_state =="breast cancer" ~ "1",
              disease_state =="no cancer" ~ "0",
              TRUE ~ "X" #this replaces the non-mentioned other values to be x
            ))
x$value <- as.numeric(x$value)
Nx2 <- t(Nx)
CallMe <- function(q) {
  x$ex <- Nx2[,q]
  glm.fit=glm(value~ex,data=x, family=binomial)
  library(pROC)
  par(pty="s")
  roc(x$value, glm.fit$fitted.values, plot=TRUE, legacy.axes=TRUE, percent=TRUE,
      xlab="100-Specificity",ylab="Sensitivity", col="#377eb8",lwd=4, print.auc=TRUE, print.auc.x=45)    
  legend("bottomright" ,legend=c(q), col=c("#377eb8"), lwd=4)
}
CallMe("hsa-miR-335-5p")
  