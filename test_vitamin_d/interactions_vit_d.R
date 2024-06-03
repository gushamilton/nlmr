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


dat %>%
  ggplot(aes(x = body_mass_index_bmi, y = bmi_prs)) +
  geom_point() +
  geom_smooth() +
  geom_smooth(method = "lm")

  tidy(lm(body_mass_index_bmi ~ bmi_prs, data = dat))
  
library(DRMR)
  
d_nom <- dat %>%
  drop_na(body_mass_index_bmi, bmi_prs, age_at_recruitment)

s <- SUMnlmr::create_nlmr_summary(d_nom$age_at_entry,x = d_nom$body_mass_index_bmi, d_nom$bmi_prs, q = 100, strata_method = "ranked") 
s$summary %>%
  ggplot(aes(x = bxse)) +
  geom_histogram() +
  scale_x_log10()
