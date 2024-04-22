pacman::p_load(tidyverse, vroom, qqman, data.table, TwoSampleMR, fastman, patchwork)

d <- vroom("/Users/fh6520/data/non-linear/combined_pre_sub.tsv.gz")
names <- vroom("/Users/fh6520/data/non-linear/01.regenie.Ydict", col_names = c("column", "pheno"))

names
bmi_plot <- d %>%
  filter(A1FREQ >0.01, INFO >0.8) %>%
  mutate(p = 10^-LOG10P.Y1) %>%
  fastman::fastman_gg(chr = "CHROM", bp = "GENPOS", p = "p", maxP = 100) +
  ggtitle("Raw BMI")



ranked_plot <- d %>%
  filter(A1FREQ >0.01, INFO >0.8) %>%
  mutate(p = 10^-LOG10P.Y3) %>%
  fastman::fastman_gg(chr = "CHROM", bp = "GENPOS", p = "p", maxP = 50) +
  ggtitle("Ranked strata")

resid_plot <- d %>%
  filter(A1FREQ >0.01, INFO >0.8) %>%
  mutate(p = 10^-LOG10P.Y4) %>%
  fastman::fastman_gg(chr = "CHROM", bp = "GENPOS", p = "p", maxP = 50)+
  ggtitle("Residual strata")

extreme_ranked <- d %>%
  filter(A1FREQ >0.01, INFO >0.8) %>%
  mutate(p = 10^-LOG10P.Y5) %>%
  fastman::fastman_gg(chr = "CHROM", bp = "GENPOS", p = "p", maxP = 50)+
  ggtitle("Extreme strata (ranked)")

extreme_resid <- d %>%
  filter(A1FREQ >0.01, INFO >0.8) %>%
  mutate(p = 10^-LOG10P.Y6) %>%
  fastman::fastman_gg(chr = "CHROM", bp = "GENPOS", p = "p", maxP = 50) +
  ggtitle("Extreme strata (residual)")

overall <- bmi_plot / resid_plot / ranked_plot / extreme_resid / extreme_ranked +
  plot_annotation(tag_levels = "A")

ggsave("GWAS/overall_gwas.tiff", compression = "lzw", height = 14, width = 10)


# bring in BMI exposure from GIANT
bmi_e <- fread("https://portals.broadinstitute.org/collaboration/giant/images/1/15/SNP_gwas_mc_merge_nogc.tbl.uniq.gz") 


bmi_e %>%
  filter(SNP == "rs17234097")
clumped_bmi_exposure <- bmi_e %>%
  rename(pval.exposure = p) %>%
  filter(pval.exposure < 5e-8) %>%
  clump_data() %>%
 as_tibble() %>%
  mutate(effect_allele.exposure = A1,
         other_allele.exposure = A2,
         eaf.exposure = Freq1.Hapmap,
         beta.exposure = b,
         se.exposure = se,
         pval.exposure,
         exposure = "BMI")
  
  
bmi_outcome <- bmi_e %>%
  rename(pval.outcome = p) %>%
  as_tibble() %>%
  mutate(effect_allele.outcome = A1,
         other_allele.outcome = A2,
         eaf.outcome = Freq1.Hapmap,
         beta.outcome = b,
         se.outcome = se,
         pval.outcome,
         outcome = "BMI",
         id.outcome = "BMI")


# First, run MR from BMI onto these outcomes




run_mr_bmi_onto_outcome <- function( col = "Y1", name = "bmi", data = d) {
  beta_col <- paste0("BETA.", col)
  se_col <- paste0("SE.", col)
  pval_col <- paste0("LOG10P.", col)
  
  outcome <- data %>%
    transmute(
      SNP = ID,
      !!paste0("beta.outcome") := .data[[beta_col]],
      !!paste0("se.outcome") := .data[[se_col]],
      !!paste0("pval.outcome") := 10^(-.data[[pval_col]]),
      effect_allele.outcome = ALLELE1,
      other_allele.outcome = ALLELE0,
      eaf.outcome = A1FREQ,
      id.outcome = name,
      outcome = name,
      CHROM,
      GENPOS)
      
      
    dat <- harmonise_data(clumped_bmi_exposure, outcome, action = 1)
    mr(dat, method = "mr_ivw")
    
}



mr_bmi_onto_outcome_res <- map2_dfr(.x = names$column, .y = names$pheno, run_mr_bmi_onto_outcome)


# now run the reverse


run_mr_outcome_onto_bmi <- function( col = "Y1", name = "bmi", data = d) {
  beta_col <- paste0("BETA.", col)
  se_col <- paste0("SE.", col)
  pval_col <- paste0("LOG10P.", col)
  
  exposure <- data %>%
    transmute(
      SNP = ID,
      !!paste0("beta.outcome") := .data[[beta_col]],
      !!paste0("se.outcome") := .data[[se_col]],
      !!paste0("pval.outcome") := 10^(-.data[[pval_col]]),
      effect_allele.outcome = ALLELE1,
      other_allele.outcome = ALLELE0,
      eaf.outcome = A1FREQ,
      id.outcome = name,
      outcome = name,
      CHROM,
      GENPOS) %>%
    convert_outcome_to_exposure() %>%
    filter(pval.exposure < 5e-8) %>%
    clump_data()
  
  
  dat <- harmonise_data(exposure, bmi_outcome, action = 1)
  mr(dat, method = "mr_ivw")
  
}

run_mr_outcome_onto_bmi()



# now run interaction analyses

interaction <- vroom("/Users/fh6520/data/non-linear/prs_interaction_bmi.tsv.gz")
p_interaction <- interaction %>%
  filter(A1FREQ >0.01 & INFO >0.8 & TEST == "ADD-INT_SNPxbmi_prs") %>%
  mutate(p = 10^-LOG10P) %>%
  fastman::fastman_gg(chr = "CHROM", bp = "GENPOS", p = "p") +
  ggtitle("68 SNP BMI PRS x SNP interaction")

interaction %>%
  head(1e6) %>%
  filter(A1FREQ >0.01 & INFO >0.8 & TEST == "ADD-INT_SNP") %>%
  mutate(p = 10^-LOG10P) %>%
  fastman::fastman_gg(chr = "CHROM", bp = "GENPOS", p = "p", maxP = 50) +
  ggtitle("68 SNP BMI PRS x SNP interaction")


interaction %>%
  filter(A1FREQ >0.01 & INFO >0.8 & TEST == "ADD-INT_SNPxbmi_prs") %>%
  mutate(p = 10^-LOG10P) %>%
  arrange(p)

int_qq <- interaction %>%
  filter(A1FREQ >0.01 & INFO >0.8 & TEST == "ADD-INT_SNPxbmi_prs") %>%
  mutate(p = 10^-LOG10P) 




interaction %>%
  filter(ID == "rs77867957") %>%
  mutate(p = 10^-LOG10P) %>%
  arrange(p)

snps_to_remove = c("rs2228145", "rs555500075")

interaction_only <- interaction %>%
  filter(A1FREQ >0.01 & INFO >0.8 & TEST == "ADD-INT_SNPxbmi_prs")
# Define a window size
window_size <- 500000



# Create a list of intervals (chromosome and position range) to exclude
exclusion_intervals <- interaction %>%
  filter(ID %in% snps_to_remove) %>%
  transmute(CHROM, start = GENPOS - window_size, end = GENPOS + window_size) %>%
  distinct(CHROM, start, end)



# Initialize a vector to collect SNPs to potentially remove
snps_within_window <- vector("character")

# Loop through each SNP to find those within 500kb either side
for (snp in snps_to_remove) {
  # Get the chromosome and position of the current SNP
  snp_info <- interaction_only %>%
    filter(ID == snp) %>%
    select(CHROM, GENPOS)
  
  # If snp_info is empty (SNP not found), skip to next iteration
  if (nrow(snp_info) == 0) next
  
  # Find all SNPs within the window on the same chromosome
  current_snps <- interaction_only %>%
    filter(CHROM == snp_info$CHROM & 
             GENPOS >= snp_info$GENPOS - window_size & 
             GENPOS <= snp_info$GENPOS + window_size) %>%
    pull(ID)
  
  # Append to the list of SNPs to potentially remove
  snps_within_window <- c(snps_within_window, current_snps)
}

# Make the vector unique to remove duplicates
snps_within_window <- unique(snps_within_window)

  