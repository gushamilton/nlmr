pacman::p_load(tidyverse, vroom, janitor, data.table, nlmr, SUMnlmr, broom, ggforestplot, patchwork, OneSampleMR)
dat <- vroom("~/data/non-linear/pheno/analysis_data_R1.tsv.gz")

#ensure participants have measured BMI and BMI prs

exclusions <-vroom("~/data/non-linear/geno/exclusions/data.minimal_relateds.plink.txt.gz", col_names = c("ieu", "IID"))
exclusions_highly <-vroom("~/data/non-linear/geno/exclusions/data.highly_relateds.plink.txt.gz", col_names = c("ieu", "IID"))

dat <- dat %>% 
  anti_join(exclusions) %>%
  anti_join(exclusions_highly)

relevel(factor(dat$uk_biobank_assessment_centre), ref = "Leeds")

dat$pr
dat <- dat %>%
  mutate(ukbac = relevel(factor(dat$uk_biobank_assessment_centre), ref = "Leeds")) %>%
  drop_na(body_mass_index_bmi, prs_vitd)
d <- (lm(ldl_direct ~ scale(prs_LDL), data = dat))
bp_test_r <- bptest(d)
bp_test_r$statistic
pchisq(bp_test_r$statistic, df = 1, lower.tail = F, log.p = F)

skedastic::glejser(d)

# Print the R-squared value
print(paste("R-squared value for variance explained:", r_squared_residuals))

dat %>%
  mutate(strat = ntile(prs_LDL, 10)) %>%
  group_by(strat) %>%
  summarise(mean_vitd = mean(ldl_direct, na.rm = T), sd_vit_d = sd(ldl_direct, na.rm = T)^2) %>%
  gt::gt()


dat %>%
  ggplot(aes(x = body_mass_index_bmi, y = bmi_prs)) +
  geom_point() +
  geom_smooth() +
  geom_smooth(method = "lm")

  tidy(lm(body_mass_index_bmi ~ bmi_prs, data = dat))
  
library(DRMR)
  
d_nom <- dat %>%
  drop_na(body_mass_index_bmi, bmi_prs, age_at_recruitment)


