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


m <- lm(body_mass_index_bmi ~ SNP, data = d)
tidy(m) %>%
  bind_cols(skedastic::glejser(m) %>%
              transmute(gles = p.value)) %>%
  filter(term == "SNP") %>%
  mutate(SNP = name)

}

run_all_interactions("rs1000940")
res <- map_dfr(list_snps$SNP, run_all_interactions)

res %>%
  arrange(-gles) %>%
  left_join(list_snps)


internal_f <- function(d){
t <- tidy(lm(body_mass_index_bmi ~ bmi_prs, data = d)) %>%
  filter(term == "bmi_prs")}


dat %>%
  mutate(q_prs = ntile(bmi_prs,3)) %>%
  group_by(q_prs) %>%
  summarise(mean = mean(body_mass_index_bmi, na.rm = T))


tibble(bmi_prs = rnorm(mean = 5, n = 1e5),
       lbmi = 0.2* bmi_prs + rnorm(1e5) ,
       bmi = exp(lbmi)) %>%
  ggplot(aes(y = bmi, x = bmi_prs)) +
  geom_point() +
  geom_smooth() 
  
  
  
  # Create the dataset
  set.seed(124) # for reproducibility
d <- tibble(bmi_prs = rnorm(mean = 5, n = 3e5),
            lbmi = 0.4 * bmi_prs + rnorm(3e5),
            bmi = exp(lbmi))

# Fit Ordinary Least Squares (OLS) Regression
ols_model <- lm(bmi ~ bmi_prs, data = d)

# Fit Fractional Polynomial Regression
fracpoly_model <- mfp(bmi ~ fp(bmi_prs), data = d)

comparison <- data.frame(
  Model = c("OLS", "Fractional Polynomial"),
  AIC = c(AIC(ols_model), AIC(fracpoly_model)),
  BIC = c(BIC(ols_model), BIC(fracpoly_model)),
  Adjusted_R2 = c(summary(ols_model)$adj.r.squared, summary(fracpoly_model)$adj.r.squared)
)
comparison


# Fit Ordinary Least Squares (OLS) Regression
ols_model <- lm((vitd) ~ (prs_vitd), data = dat)


# Fit Fractional Polynomial Regression
fracpoly_model <- mfp(vitd ~ fp(prs_vitd), data = dat)

comparison <- data.frame(
  Model = c("OLS", "Fractional Polynomial"),
  AIC = c(AIC(ols_model), AIC(fracpoly_model)),
  BIC = c(BIC(ols_model), BIC(fracpoly_model)),
  Adjusted_R2 = c(summary(ols_model)$adj.r.squared, summary(fracpoly_model)$adj.r.squared)
)
comparison

to_bring_in <- files %>%
  str_subset("rs11165643")

snp <- fread(to_bring_in)

d <- all_snps %>%
  select(IID, contains("rs11165643")) %>%
  select(!contains("HET")) %>%
  rename("SNP" = 2) %>%
  filter(SNP == 1 | SNP == 2 | SNP == 0) %>%
  right_join(just_bmi) %>%
  left_join(snp)

tidy(lm(body_mass_index_bmi ~ SNP, data = d))
test <- d %>%
  drop_na(SNP, body_mass_index_bmi)

z <- test %>%
  mutate(g = scale(SNP),
         x = body_mass_index_bmi,
         y =age_at_recruitment) %>%
  as.data.frame() %>%
  generate_all_sumstats(k = 10)

z <- test %>%
  mutate(g = scale(SNP),
         u = rnorm(nrow(test)),
         x = scale(body_mass_index_bmi) + 2* u ,
         y = 2 *u +rnorm(nrow(test))) %>%
  as.data.frame() %>%
  generate_all_sumstats(k = 10)

z$ss %>%
  ggforestplot::forestplot(estimate = beta_mr, se = se_mr, name = strata) +
  facet_wrap(~ strata_method) +
  coord_flip() +
  theme_bw()


test %>%
       transmute(u = rnorm(377713),
                 x = -0.2* u + body_mass_index_bmi,
                 y = 2*u + rnorm(383793),
                 g = prs_vitd) %>%
       drop_na()
