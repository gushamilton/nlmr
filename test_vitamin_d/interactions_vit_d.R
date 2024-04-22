pacman::p_load(tidyverse, vroom, janitor, data.table, nlmr, SUMnlmr, broom, ggforestplot, patchwork, OneSampleMR)
dat <- vroom("~/data/non-linear/pheno/analysis_data_R1.tsv.gz")

#ensure participants have measured BMI and BMI prs

exclusions <-vroom("~/data/non-linear/geno/exclusions/data.minimal_relateds.plink.txt.gz", col_names = c("ieu", "IID"))
exclusions_highly <-vroom("~/data/non-linear/geno/exclusions/data.highly_relateds.plink.txt.gz", col_names = c("ieu", "IID"))

dat <- dat %>% 
  anti_join(exclusions) %>%
  anti_join(exclusions_highly)

relevel(factor(dat$uk_biobank_assessment_centre), ref = "Leeds")

dat <- dat %>%
  mutate(ukbac = relevel(factor(dat$uk_biobank_assessment_centre), ref = "Leeds"))
d <- tidy(lm(log(vitd) ~ scale(prs_vitd)*townsend_deprivation_index_at_recruitment, data = dat))
