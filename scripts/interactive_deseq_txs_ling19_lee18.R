library(DESeq2)
library(Matrix)
library(tximport)
library(tidyverse)
library(magrittr)


files <- list.files("data/salmon", "quant.sf", full.names=TRUE, recursive=TRUE) %>%
    { .[str_detect(.,"/(pe/ery_d_gata1_wt|pe/lsk_sf3b1_wt_srsf2_wt)")] } %>%
    set_names(str_extract(., "(lsk|ery)[^/]+"))

txi <- tximport(files, type="salmon", txOut=TRUE, 
                countsFromAbundance="scaledTPM")
cts <- txi$counts
cts <- cts[rowSums(cts) > 0,]

celltype_map <- c(hspc="HSC", ery_d="EryD")
df_samples <- read_tsv("metadata/sample-sheet.tsv") %>%
    filter(sample_id %in% colnames(cts)) %>%
    mutate(cell_type=celltype_map[cell_type]) %>%
    set_rownames(.$sample_id) %>%
    DataFrame() %>%
    { .[colnames(cts),] }
    

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
apa_id_to_gene <- deframe(df_utrs[, c("APA_ID", "Gene")])

cts %<>% `[`(rownames(.) %in% names(utr_to_apa_id),)

rownames(cts) %<>% { utr_to_apa_id[.] } %>% unname

mode(cts) <- 'integer'

#M_gene_X_tx <- rownames(cts) %>% { apa_id_to_gene[.] } %>% fac2sparse

#cts_gene <- drop0(round(M_gene_X_tx %*% cts))

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=df_samples,
                              design=~cell_type)

dds <- DESeq(dds)

dres <- DESeq2::results(dds, contrast=c("cell_type", "EryD", "HSC"),
                        independentFiltering=FALSE, cooksCutoff=FALSE)

saveRDS(dds, "data/deseq/dds_hsc_eryd_ling19_lee18_txs.Rds")

saveRDS(dres, "data/deseq/dres_hsc_eryd_ling19_lee18_txs.Rds")
