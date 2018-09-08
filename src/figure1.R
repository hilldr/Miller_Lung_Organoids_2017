## PCA analysis ------------------------------------------------------
## Load dataset from file
data.dir <- "../results/full-dataset"
df <- readr::read_csv(file = file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv"))

## rename for SeqRetriever compatibility

num.data <- df[,sapply(df,is.numeric)]

## read in table with sample metadata
samples <- readr::read_csv(file = file.path(data.dir,"AJMsample_key.csv"))
group <- samples$sample_name
## pca.data <- num.data[apply(num.data, 1, sd, na.rm=TRUE) != 0,]
## calculate variance by row (gene)
var <- apply(num.data, 1, sd, na.rm = TRUE)
## adjust cut off according to variance percentile
pca.data <- num.data[var > quantile(var, 0.9) & var != 0,]
pca <- prcomp(t(pca.data),scale = TRUE,center = TRUE)
scores <- data.frame(colnames(pca.data), pca$x[,1:ncol(pca$x)],group)


## PCA plot
## function to format decimals as precentage
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

library(ggplot2)
library(RColorBrewer)
source("../../../HIO_Ecoli_paper/DOC/eLife/src/ggplot2-themes.R")
fig1c <- qplot(x = PC1, y = PC2, data = scores) +  
    scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Paired"))(length(unique(scores$group)))) +
    scale_color_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Paired"))(length(unique(scores$group)))) +
    geom_point(shape = 21, aes(fill = factor(group), color = factor(group)), size = 12) +
    theme1 +
    theme(legend.position = "bottom",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white",fill = "white")) +
    geom_hline(yintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0,
               size = 1, linetype = "dashed", color = "grey70") +
    coord_fixed(ratio = 1) +
    xlab(paste("PC1 (",percent(round(summary(pca)$importance[2,1],4)),")",sep = "")) +
    ylab(paste("PC2 (",percent(round(summary(pca)$importance[2,2],4)),")",sep = ""))

print(fig1c)
print(length(rownames(pca.data)))

## TSNE analysis ----------------------------------------------------------------
##http://distill.pub/2016/misread-tsne/
## PCA analysis ------------------------------------------------------
## Load dataset from file
data.dir <- "../results/full-dataset"
df <- readr::read_csv(file = file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv"))

## rename for SeqRetriever compatibility

num.data <- df[,sapply(df,is.numeric)]
samples <- readr::read_csv(file = file.path(data.dir,"AJMsample_key.csv"))
group <- samples$sample_name
## pca.data <- num.data[apply(num.data, 1, sd, na.rm=TRUE) != 0,]
## calculate variance by row (gene)
var <- apply(num.data, 1, sd, na.rm=TRUE)
## adjust cut off according to variance percentile
pca.data <- num.data[var > quantile(var, 0.1) & var != 0,]
library(Rtsne)
set.seed(3210)
tsne <- Rtsne(as.matrix(unique(t(pca.data))),
              dims = 2,
	      perplexity = 10,
              theta = 0.5)
scores <- data.frame(tsne$Y,group)
#scores$fmt <- as.factor(sapply(strsplit(x = as.character(scores$group), split = "[_]"),"[",2))
#scores$day <- as.factor(sapply(strsplit(x = as.character(scores$group), split = "[_]"),"[",3))


library(RColorBrewer) 
library(ggplot2)
library(grid)

plot <- ggplot(data = scores, aes(x = X1, y = X2, fill = group)) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Set1"))(length(unique(scores$group)))) +
    geom_point(shape = 21, size = 12) +
  #  scale_fill_brewer(palette = "Paired") +
    theme1 + 
    theme(legend.position = "bottom",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white", fill = "white")) +
    coord_fixed(ratio = 1)

print(plot)

## TSNE analysis ----------------------------------------------------------------
##http://distill.pub/2016/misread-tsne/
## PCA analysis ------------------------------------------------------
## Load dataset from file
data.dir <- "../results/figure1-dataset"
df <- readr::read_csv(file = file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv"))

## rename for SeqRetriever compatibility

num.data <- df[,sapply(df,is.numeric)]
samples <- readr::read_csv(file = file.path(data.dir,"AJMsample_key.csv"))
samples <- samples[grep("fetal|enteroid", samples$sample_name),]
group <- samples$sample_name
## pca.data <- num.data[apply(num.data, 1, sd, na.rm=TRUE) != 0,]
## calculate variance by row (gene)
var <- apply(num.data, 1, sd, na.rm=TRUE)
## adjust cut off according to variance percentile
pca.data <- num.data[var > quantile(var, 0.5) & var != 0,]
library(Rtsne)
set.seed(3210)
tsne <- Rtsne(as.matrix(unique(t(pca.data))),
              dims = 2,
	      perplexity = 5,
              theta = 0.5)
scores <- data.frame(tsne$Y,group)
#scores$fmt <- as.factor(sapply(strsplit(x = as.character(scores$group), split = "[_]"),"[",2))
#scores$day <- as.factor(sapply(strsplit(x = as.character(scores$group), split = "[_]"),"[",3))


library(RColorBrewer) 
library(ggplot2)
library(grid)
source("../../../HIO_Ecoli_paper/DOC/eLife/src/ggplot2-themes.R")
plot <- ggplot(data = scores, aes(x = X1, y = X2, fill = group)) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Set1"))(length(unique(scores$group)))) +
    geom_point(shape = 21, size = 12) +
  #  scale_fill_brewer(palette = "Paired") +
    theme1 + 
    theme(legend.position = "bottom",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white", fill = "white")) +
    coord_fixed(ratio = 1)

print(plot)

## TSNE analysis ----------------------------------------------------------------
##http://distill.pub/2016/misread-tsne/
## PCA analysis ------------------------------------------------------
## Load dataset from file
data.dir <- "../results/figure1-hES-dataset"
df <- readr::read_csv(file = file.path(data.dir,"complete-dataset_DESeq2-normalized-counts.csv"))

## rename for SeqRetriever compatibility

num.data <- df[,sapply(df,is.numeric)]
samples <- readr::read_csv(file = file.path(data.dir,"AJMsample_key.csv"))
samples <- samples[grep("fetal|hPSC", samples$sample_name),]
samples[grep("5d|9d", samples$sample_name),]$sample_name <- "fetal lung"
group <- samples$sample_name
## pca.data <- num.data[apply(num.data, 1, sd, na.rm=TRUE) != 0,]
## calculate variance by row (gene)
var <- apply(num.data, 1, sd, na.rm=TRUE)
## adjust cut off according to variance percentile
pca.data <- num.data[var > quantile(var, 0.5) & var != 0,]
library(Rtsne)
set.seed(3210)
tsne <- Rtsne(as.matrix(unique(t(pca.data))),
              dims = 2,
	      perplexity = 5,
              theta = 0.5)
scores <- data.frame(tsne$Y,group)
#scores$fmt <- as.factor(sapply(strsplit(x = as.character(scores$group), split = "[_]"),"[",2))
#scores$day <- as.factor(sapply(strsplit(x = as.character(scores$group), split = "[_]"),"[",3))


library(RColorBrewer) 
library(ggplot2)
library(grid)
source("../../../HIO_Ecoli_paper/DOC/eLife/src/ggplot2-themes.R")
plot <- ggplot(data = scores, aes(x = X1, y = X2, fill = group)) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(length(unique(scores$group)), "Set1"))(length(unique(scores$group)))) +
    geom_point(shape = 21, size = 12) +
  #  scale_fill_brewer(palette = "Paired") +
    theme1 + 
    theme(legend.position = "bottom",
          legend.background = element_rect(colour = "white"),
          legend.key = element_rect(color = "white", fill = "white")) +
    coord_fixed(ratio = 1)

print(plot)
