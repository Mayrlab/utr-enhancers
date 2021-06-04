library(DEXSeq)
library(tximport)
library(tidyverse)
library(magrittr)


files <- list.files("data/salmon", "quant.sf", full.names=TRUE, recursive=TRUE) %>%
    { .[str_detect(.,"/(pe/ery_d_gata1_wt|pe/lsk_sf3b1_wt_srsf2_wt)")] } %>%
    set_names(str_extract(., "(lsk|ery)[^/]+"))

txi <- tximport(files, type="salmon", txOut=TRUE, 
                countsFromAbundance="scaledTPM")
cts <- txi$counts

celltype_map <- c(hspc="HSC", ery_d="EryD")
df_samples <- read_tsv("metadata/sample-sheet.tsv") %>%
    filter(sample_id %in% colnames(cts)) %>%
    mutate(cell_type=celltype_map[cell_type]) %>%
    set_rownames(.$sample_id) %>%
    { .[colnames(cts),] } %>%
    set_rownames(.$sample_id)

df_qapa <- read_tsv("data/qapa/ery_hsc_select_results.txt",
                    col_types="cccccddcdddd____________________________")

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

## Use annots from QAPA output to get APA_ID
df_utrs <- rownames(cts) %>%
    separate_ensembl_field %>%
    select(UTR_ID, Chr, Strand, UTR3.Start, UTR3.End) %>%
    inner_join(df_qapa, by=c("Chr", "Strand", "UTR3.Start", "UTR3.End")) %>%
    select(APA_ID, Gene_Name, Chr, UTR3.Start, UTR3.End, Strand, Length, everything())

utr_to_apa_id <- deframe(df_utrs[, c("UTR_ID", "APA_ID")])

cts <- cts[rowSums(cts) > 0,]
cts %<>% `[`(rownames(.) %in% names(utr_to_apa_id),)

rownames(cts) %<>% { utr_to_apa_id[.] }

dxd <- DEXSeqDataSet(countData=round(cts), sampleData=df_samples,
                     design=~sample + exon + cell_type:exon,
                     featureID=rownames(cts),
                     groupID=str_extract(rownames(cts), "^[^_]+"))

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd)
dxd <- testForDEU(dxd, reducedModel=~sample + exon)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar='cell_type', denominator='HSC')
dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval), qval)
columns <- c("featureID","groupID","pvalue")
dxr.t <- as.data.frame(dxr[,columns])
head(dxr.t)
dxr.full <- dxr[,-c(11,12)] %>%
    as_tibble() %>%
    rename(Gene=groupID, APA_ID=featureID)
saveRDS(dxd, "data/dexseq/dxd_hsc_eryd_ling19_lee18.Rds")
write_tsv(dxr.g, "data/dexseq/dxr_hsc_eryd_ling19_lee18_gene.tsv.gz")
write_tsv(dxr.t, "data/dexseq/dxr_hsc_eryd_ling19_lee18_tx.tsv.gz")
write_tsv(dxr.full, "data/dexseq/dxr_hsc_eryd_ling19_lee18_full.tsv.gz")
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
write_tsv(dex.padj, "data/dexseq/staged_qvals_hsc_eryd_ling19_lee18_cdr5.tsv.gz")
