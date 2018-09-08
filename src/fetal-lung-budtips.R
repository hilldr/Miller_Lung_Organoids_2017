## export complete dataset with DESeq2
## create directory to deposit results
## create directory to deposit results
data.dir <- "../results/fetal_lung_uncultured/"
src.dir <- "./"
dir.create(path = data.dir, recursive = TRUE)

## read in table with sample metadata
samples <- readr::read_csv(file = file.path(src.dir,"AJMsample_key.csv"))
samples <- samples[grep("adult_lung|uncultured", samples$sample_name),]

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
write.csv(ddscounts, file =  file.path(data.dir, "complete-dataset_DESeq2-normalized.csv"))
save(dds, file = file.path(data.dir, "dds.Rdata"))
rm(list = ls())
