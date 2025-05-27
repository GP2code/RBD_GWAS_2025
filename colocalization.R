#### Colocalization using GTEx v7 eQTL data 
# Load libraries
library(data.table)
library(dplyr)
library(coloc)

## Load and format GWAS
gwas <- read.table("gwas_RCOR1_hg19.txt", header = T)

gwas_formatted <- gwas %>%
  rename(snp = SNP, chr = CHR, pos = BP) %>%
  mutate(
    varbeta = se^2,
    type = "cc",
    N = 143957,             # total sample size
    s = 3564/143957,         # proportion of cases
    snp = paste0(chr, ":", pos, ":", A2, ":", A1)
  ) %>%
  select(snp, chr, pos, A1, A2, beta, varbeta, MAF, N, type, s, P)

## Check if any MAFs are > 0.5, if so use code below
# gwas_formatted <- gwas_formatted %>%
  # mutate(beta = ifelse(MAF > 0.5, -1 * beta, beta))

## Load and format eQTL data
eqtl <- fread("tissue.txt")

eqtl_formatted <- eqtl %>%
  filter(!is.na(slope), slope_se > 0) %>%
  mutate(
        varbeta = slope_se^2,
        type = "quant",
        N = ma_samples
    ) %>%
  select(snp = variant_id, chr, pos, A1, A2, gene_id, beta = slope, varbeta, MAF = maf, N, type, P = pval_nominal)

## Filter both datasets for +/- 500kb of lead RCOR1 snp
# lead var is rs35785423 
eqtl_formatted <- eqtl_formatted %>% filter(pos >= 102574798 & pos <= 103574798)
gwas_formatted <- gwas_formatted %>% filter(pos >= 102574798 & pos <= 103574798)

## Match SNPs and prep input
common_snps <- intersect(gwas_formatted$snp, eqtl_formatted$snp)
length(common_snps)

# Flip snp id allele coding to capture more matches
gwas_unmatched <- gwas_formatted %>% filter(!(snp %in% common_snps))
eqtl_unmatched <- eqtl_formatted %>% filter(!(snp %in% common_snps))

gwas_unmatched_flipped <- gwas_unmatched %>%
  mutate(
    snp_flipped = paste0(chr, ":", pos, ":", A1, ":", A2)
  )

flipped_matches <- intersect(gwas_unmatched_flipped$snp_flipped, eqtl_formatted$snp)

# Replace snp names with flipped versions where needed
gwas_unmatched_flipped <- gwas_unmatched_flipped %>%
  filter(snp_flipped %in% flipped_matches) %>%
  mutate(snp = snp_flipped) %>%
  select(-snp_flipped)

# Combine correct GWAS data
gwas_final <- bind_rows(
  gwas_formatted %>% filter(snp %in% common_snps),
  gwas_unmatched_flipped
)

# Pull only matches from eQTL data
eqtl_final <- eqtl_formatted %>% filter(snp %in% gwas_final$snp)

# Log snp counts
cat("SNPs before flipping:", length(common_snps), "\n")
cat("SNPs rescued by flipping:", length(flipped_matches), "\n")
cat("Total overlapping SNPs after flipping:", nrow(gwas_final), "\n")

## Gene-specific colocalization within 1 Mb region
# Remove duplicated SNPs
gwas_final <- gwas_final %>% distinct(snp, .keep_all = TRUE)

# combine datasets
df <- dplyr::inner_join(gwas_final, eqtl_final, by = "snp", suffix = c(".gwas", ".eqtl"))

# Get list of unique genes in eQTL data
genes <- unique(df$gene_id)

# apply colocalization function based on gene ID
results <- lapply(genes, function(gene) {
  df2 <- df[df$gene_id == gene, ]

    # remove dup snps
    df2 <- df2[!duplicated(df2$snp), ]

    # Skip genes with too few SNPs
    if (nrow(df2) < 10) return(NULL)

    # Prep GWAS input
    dataset1 <- list(
        beta = df2$beta.gwas,
        varbeta = df2$varbeta.gwas,
        MAF = df2$MAF.gwas,
        N = df2$N.gwas[1],
        s = df2$s[1],
        type = "cc",
        snp = df2$snp
    )

    # Prep eQTL input
    dataset2 <- list(
        beta = df2$beta.eqtl,
        varbeta = df2$varbeta.eqtl,
        MAF = df2$MAF.eqtl,
        N = df2$N.eqtl[1],        # use eQTL sample size
        type = "quant",
        snp = df2$snp
    )

  res <- coloc.abf(dataset1 = dataset1, dataset2 = dataset2)
  summary_df <- as.data.frame(t(res$summary))
  return(data.frame(gene = gene, summary_df))
})

# Combine results
results_table <- do.call(rbind, results)
results_table <- results_table %>%
  mutate(across(where(is.numeric), ~ round(., 4)))
write.table(results_table,"tissue_results.txt",col.names = T, row.names = F, quote = F)




