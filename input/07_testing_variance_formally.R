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

just_bmi <- dat %>%
  select(body_mass_index_bmi, bmi_prs, IID = ieu, contains("PC"), sex, age_at_recruitment) 

d <- all_snps %>%
  left_join(just_bmi)

# function to extract the betas and glesjer test

get_betas<- function(name){

  
  d <- all_snps %>%
    select(IID, contains(name)) %>%
    select(!contains("HET")) %>%
    rename("SNP" = 2) %>%
    right_join(just_bmi) %>%
    left_join(snp)
  
  
  m <- lm(body_mass_index_bmi ~ SNP, data = d)
  tidy(m) %>%
    bind_cols(skedastic::glejser(m) %>%
                transmute(gles = p.value)) %>%
    filter(term == "SNP") %>%
    mutate(SNP = name)
  
}

all_snp_names <- colnames(all_snps) %>%
  str_subset("rs") %>%
  str_subset("HET", negate = T)

betas <- map_dfr(all_snp_names, get_betas)

#  now we have all the betas here:

betas

#  generate the master PRS

# Ensure SNP columns in all_snps have the same order as in betas
snps_data <- all_snps %>% select(IID, one_of(betas$SNP))

# Convert SNP data and betas to matrices
snp_matrix <- as.matrix(snps_data %>% select(-IID))
beta_vector <- as.numeric(betas$estimate)

# Calculate the PRS for each participant using matrix multiplication
prs_values <- snp_matrix %*% beta_vector
working <- tibble(prs = prs_values[,1],
       IID = all_snps$IID) %>%
  left_join(just_bmi)

full_prs <- working %>% select(prs_full = prs, IID)

#  yep, explains 2% of the BMI variance, and clearly fails glejeset

m <- (lm(body_mass_index_bmi ~ prs, data = working))
glance(m)
glejser(m)

# now, let's create some IV's which don't fail glesjer

invariant_thresh_0.05 <- betas %>%
  filter(gles >0.05) %>%
  filter(p.value < 0.05) 


snps_data <- all_snps %>% select(IID, one_of(invariant_thresh_0.05$SNP))

# Convert SNP data and betas to matrices
snp_matrix <- as.matrix(snps_data %>% select(-IID))
beta_vector <- as.numeric(invariant_thresh_0.05$estimate)


# Calculate the PRS for each participant using matrix multiplication
prs_values <- snp_matrix %*% beta_vector
working <- tibble(prs = prs_values[,1],
                  IID = all_snps$IID) %>%
  left_join(just_bmi)

prs_0.05 <- working %>% select(prs_0.05 = prs, IID)

#  yep, much weaker IV, p = 10^-35, and still fails glejser but weakly

m <- (lm(body_mass_index_bmi ~ prs, data = working))
glance(m)
glejser(m)
tidy(m)


# now, let's create some IV's which are a bit looser

invariant_thresh_0.01 <- betas %>%
  filter(gles >0.01) %>%
  filter(p.value < 0.05) 


snps_data <- all_snps %>% select(IID, one_of(invariant_thresh_0.01$SNP))

# Convert SNP data and betas to matrices
snp_matrix <- as.matrix(snps_data %>% select(-IID))
beta_vector <- as.numeric(invariant_thresh_0.01$estimate)


# Calculate the PRS for each participant using matrix multiplication
prs_values <- snp_matrix %*% beta_vector
working <- tibble(prs = prs_values[,1],
                  IID = all_snps$IID) %>%
  left_join(just_bmi)

prs_0.01 <- working %>% select(prs_0.01 = prs, IID)

#  yep, much weaker IV, p = 10^-35, and still fails glejser but weakly

m <- (lm(body_mass_index_bmi ~ prs, data = working))
glance(m)
glejser(m)

betas %>%
  arrange(-gles) %>%s
  mutate(f_stat = estimate^2/std.error^2)

rs11165643_T <-  betas %>%
  filter(SNP == "rs11165643_T") 



snps_data <- all_snps %>% select(IID, one_of(rs11165643_T$SNP))

# Convert SNP data and betas to matrices
snp_matrix <- as.matrix(snps_data %>% select(-IID))
beta_vector <- as.numeric(invariant_thresh_0.01$estimate)


# Calculate the PRS for each participant using matrix multiplication
prs_values <- snp_matrix %*% beta_vector
working <- tibble(prs = prs_values[,1],
                  IID = all_snps$IID) %>%
  left_join(just_bmi)

prs_rs11165643 <- working %>% select(prs_rs11165643 = prs, IID)

#  yep, much weaker IV, p = 10^-35, and still fails glejser but weakly

m <- (lm(body_mass_index_bmi ~ prs, data = working))
glance(m)
glejser(m)

final_data_set <- just_bmi %>%
  distinct(IID, .keep_all = T) %>%
  left_join(full_prs) %>%
  distinct(IID, .keep_all = T) %>%
  left_join(prs_0.01) %>%
  distinct(IID, .keep_all = T) %>%
  left_join(prs_0.05) %>%
  distinct(IID, .keep_all = T) %>%
  left_join(prs_rs11165643) 

just_bmi$sex
z <- final_data_set %>%
  mutate(y = rnorm(nrow(final_data_set))) %>%
  mutate(sex = as.numeric(sex== "Male")) %>%
  select(g = prs_0.01, x = body_mass_index_bmi, y = sex) %>%
  drop_na() %>%
  as.data.frame() %>%
  generate_all_sumstats()



z$ss %>%
  ggforestplot::forestplot(estimate = beta_mr, se = se_mr, name = strata) +
  facet_wrap(~ strata_method) +
  coord_flip() +
  theme_bw()



z <- final_daprs_0.01z <- final_data_set %>%
  mutate(u = rnorm(nrow(final_data_set))) %>%
  mutate(y = rnorm(nrow(final_data_set)) + 0.5*u) %>%
  transmute(g = prs_full, x = body_mass_index_bmi + 2 *u, y) %>%
  drop_na() %>%
  as.data.frame() %>%
  generate_all_sumstats()


z$ss %>%
  ggforestplot::forestplot(estimate = beta_mr, se = se_mr, name = strata) +
  facet_wrap(~ strata_method) +
  coord_flip() +
  theme_bw()



