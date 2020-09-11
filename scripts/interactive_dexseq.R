library(DEXSeq)
library(tidyverse)
library(magrittr)
library(BiocParallel)
library(tximport)


files <- list.files("data/salmon", "quant.sf", full.names=TRUE, recursive=TRUE) %>%
    { .[str_detect(.,"(ery_d|lsk_tip60_wt_[12])")] } %>%
    set_names(str_extract(., "(lsk|ery)[^/]+"))

txi <- tximport(files, type="salmon", txOut=TRUE, 
                countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]

df_samples <- colnames(cts) %>%
    { data.frame(sample_id=.,
                 cell_type=ifelse(str_detect(., "ery"), "EryD", "HSC")) }

ENSEMBL_FILE="/data/mayrc/db/mm10/ensembl_identifiers_v101.txt"

df_ensembl <- read_tsv(ENSEMBL_FILE, skip=1,
                       col_names=c("Gene", "tid", "Gene_Type", "Transcript_Type", "Gene_Name"))

################################################################################
## Methods from QAPA
## License: GNU GPLv3
################################################################################

#### Split Ensembl field ####
format_multi_ensembl_ids <- function(ids) {
    ## Format Ensembl Transcript and Ensembl Gene IDs if there are multiple
    ## e.g.
    ## ENSMUST00000111043_ENSMUSG00000048482,ENSMUST00000111044_ENSMUSG00000048482_mm9_chr1
    ## becomes
    ## ENSMUST00000111043,ENSMUST00000111044_ENSMUSG00000048482_mm9_chr1
    ## Test regex: https://regex101.com/r/zuDsy1/1
    split_ids <- str_match(ids, "^(([^_]+_[^_,]+)(,[^_]+_[^_,]+)*)_(([^_]+).+)")

    if (is.na(split_ids[1])) {
        stop("Unable to format Ensembl ID by regex")
    }

    ## Separate multiple Transcript_Gene name
    ens <- strsplit(split_ids[,2], ",")
    ## Split transcript and gene names, then re-arrange to combine transcripts and genes
    ens <- lapply(ens, function(e) {
        strsplit(e, "_") %>% do.call("rbind", .) %>%
            apply(., 2, function(y) paste(unique(y), collapse=",")) %>%
            paste(., collapse="_")
    })
    stopifnot(length(ens) == length(ids))
    apply(cbind(ens, split_ids[,5]), 1, paste, collapse="_")
}

separate_ensembl_field <- function(utr_ids) {
    tx_pattern <- "^([^_]+_[^_,]+)(,[^_]+_[^_,]+)*_[^_]+_(chr)*\\w+_\\d+_\\d+_[-+]_utr_\\d+_\\d+"
    if (grepl(tx_pattern, utr_ids[1], perl = TRUE)) {
        ## Format Ensembl Transcript and Ensembl Gene IDs if there are multiple
        ## Remove utr tag
        str_extract(utr_ids, tx_pattern) %>%
            format_multi_ensembl_ids() %>%
            str_replace("_utr", "") %>%
            { tibble(UTR_ID=utr_ids, Transcript=.) } %>%
            separate(Transcript, sep="_",
                     into=c("Transcript", "Gene", "Species", "Chr", "LastExon.Start",
                            "LastExon.End", "Strand", "UTR3.Start", "UTR3.End")) %>%
            dplyr::select(-c("Species", "Gene")) %>%
            mutate(across(matches("(Start|End)$"), as.numeric),
                   Length=abs(UTR3.End - UTR3.Start),
                   tid=str_extract(Transcript, "^[^,]+"))
    } else {
        warning("Unable to find Ensembl IDs by regex")
        utr_ids
    }
}

################################################################################

utr_type <- function (utr_rank) {
    if (length(utr_rank) == 1) { "S" }
    else { ifelse(utr_rank > 1, "D", "P") }
}

df_utrs <- rownames(cts) %>%
    separate_ensembl_field %>%
    inner_join(df_ensembl, by='tid') %>%
    filter(Gene_Type == 'protein_coding') %>%
    group_by(Gene) %>%
    mutate(APA_ID=rank(Length) %>% { str_c(Gene, ., utr_type(.), sep="_") }) %>%
    ungroup() %>%
    dplyr::select(APA_ID, Gene_Name, Chr, UTR3.Start, UTR3.End, Strand, Length, everything())

utr_to_apa_id <- deframe(df_utrs[, c("UTR_ID", "APA_ID")])

cts %<>% `[`(rownames(.) %in% names(utr_to_apa_id),)

rownames(cts) %<>% { utr_to_apa_id[.] }

dxd <- DEXSeqDataSet(countData=round(cts), sampleData=df_samples,
                     design=~sample + exon + cell_type:exon,
                     featureID=rownames(cts),
                     groupID=str_extract(rownames(cts), "^[^_]+"))

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd, reducedModel=~sample + exon)

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)

columns <- c("featureID","groupID","pvalue")
dxr.t <- as.data.frame(dxr[,columns])
head(dxr.t)

saveRDS(dxd, "data/dexseq/dxd_hsc_eryd.Rds")

write_tsv(dxr.g, "data/dexseq/dxr_hsc_eryd_gene.tsv.gz")
write_tsv(dxr.t, "data/dexseq/dxr_hsc_eryd_tx.tsv.gz")

library(stageR)

dxr <- dxr[!is.na(dxr$pvalue),]
qval <- perGeneQValue(dxr)

pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(dxr$featureID,"transcript")
pScreen <- qval
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
    dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                   onlySignificantGenes=TRUE)
    })

write_tsv(dex.padj, "data/dexseq/staged_qvals_hsc_eryd_fdr5.tsv.gz")

