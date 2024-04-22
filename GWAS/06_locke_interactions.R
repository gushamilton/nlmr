pacman::p_load(tidyverse, vroom, qqman, data.table, TwoSampleMR, broom, data.table)

files <- list.files("/Users/fh6520/data/non-linear/locke_interactions", full.names = T)
overall_prs <- fread(files[1])
all_snps <- fread("/Users/fh6520/data/non-linear/locke_interactions/raw_locke_snps_genotype.tsv.gz")

dat <- vroom("~/data/non-linear/pheno/analysis_data_R1.tsv.gz")

#ensure participants have measured BMI and BMI prs

exclusions <-vroom("~/data/non-linear/geno/exclusions/data.minimal_relateds.plink.txt.gz", col_names = c("ieu", "IID"))
exclusions_highly <-vroom("~/data/non-linear/geno/exclusions/data.highly_relateds.plink.txt.gz", col_names = c("ieu", "IID"))

dat <- dat %>% 
  anti_join(exclusions) %>%
  anti_join(exclusions_highly)

# highest hit SNPS

list_snps <- vroom("/Users/fh6520/R/non_linear_testing/inputs/BMI_risk_score_pub.tsv")
list_snps %>%
  arrange(-abs(b))

# FTO SNP is rs1558902
just_bmi <- dat %>%
  select(body_mass_index_bmi, bmi_prs, IID = ieu, contains("PC"), sex, age_at_recruitment) 

run_all_interactions <- function(name){
to_bring_in <- files %>%
  str_subset(name)

snp <- fread(to_bring_in)

d <- all_snps %>%
  select(IID, contains(name)) %>%
  select(!contains("HET")) %>%
  rename("SNP" = 2) %>%
  filter(SNP == 1 | SNP == 2 | SNP == 0) %>%
  right_join(just_bmi) %>%
  left_join(snp)


tidy(lm(body_mass_index_bmi ~ SNP*SCORE1_AVG + PC1 + PC2 + PC3 + PC4 + PC5 + age_at_recruitment + sex, data = d))  %>%
  mutate(SNP = name)

}

res <- map_dfr(list_snps$SNP, run_all_interactions)

res %>%
  filter(str_detect(term, ":")) %>%
  arrange(p.value)
