## export complete dataset with DESeq2
## create directory to deposit results
## create directory to deposit results
data.dir <- "../results/full-dataset"
dir.create(path = data.dir, recursive = TRUE)

## read in table with sample metadata
samples <- readr::read_csv(file = file.path(data.dir,"AJMsample_key.csv"))

## setup access to kallisto read files
files <- file.path(samples$directory,
                   samples$file_name,
                   "abundance.h5") 

## set sample names as description_rep#_seq_rep#
names(files) <- samples$short_name
## check that all files are found
if (all(file.exists(files)) == FALSE) {
    print("kallisto files not found")
    stop()
}

## associate transcripts with gene IDs
## create biomart reference
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')
## create index of gene names
tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                         "external_gene_name"),
                          mart = mart)
## import kallisto data and generate count dataframe (dds)
## http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
library(readr)
txi <- tximport::tximport(files, type = "kallisto", tx2gene = tx2gene)

## export abundance counts
write.csv(txi$abundance, file = file.path(data.dir, "complete_dataset_txi.csv"))

library(DESeq2)
## https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
dds <- DESeq2::DESeqDataSetFromTximport(txi,
                                        colData = samples,
                                        design = ~ code_name)
## pre-filter out counts < 1
dds <- dds[rowSums(counts(dds)) > 1, ]

## write out normalized expression counts
dds <- DESeq2::estimateSizeFactors(dds)
ddscounts <- DESeq2::counts(dds, normalized = TRUE)

## write expression matrix to file
write.csv(ddscounts, file =  file.path(data.dir, "complete-dataset_DESeq2-normalized-counts.csv"))
save(dds, file = file.path(data.dir, "dds.Rdata"))
rm(list = ls())

## enable parallel processes
library("BiocParallel")

register(MulticoreParam(4))

## setup multifactor design
colData(dds)$code_name <- as.factor(colData(dds)$code_name)
ddsMF <- dds
DESeq2::design(ddsMF) <- ~ code_name

## Likelihood ratio test (ANOVA-like)
ddsMF <- DESeq2::DESeq(ddsMF, test = "LRT", reduced = ~1, parallel = TRUE)
res <- DESeq2::results(ddsMF)
write.csv(res, file = file.path(data.dir, "LRT.csv"))

## DESeq2 differential expression output
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "adult_lung","lung_buds_d0"))
write.csv(res, file = file.path(data.dir, "adult-lung_over_fresh-lung-buds.csv"))
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "adult_lung","lung_buds_d14"))
write.csv(res, file = file.path(data.dir, "adult-lung_over_lung-bud-isolated-organoids.csv"))
res <- DESeq2::results(ddsMF, test = "Wald",
                       contrast = c("code_name", "adult_lung","HLO"))
write.csv(res, file = file.path(data.dir, "adult-lung_over_PSC-derived-HLO.csv"))
