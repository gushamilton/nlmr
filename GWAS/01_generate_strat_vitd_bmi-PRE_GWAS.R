pacman::p_load(tidyverse, vroom, janitor, lubridate, data.table, survival, nlmr, SUMnlmr, broom, faux, ggforestplot, patchwork, OneSampleMR)
dat <- vroom("~/data/non-linear/pheno/analysis_data.tsv.gz")


generate_ranked_strata <- function(g, x, q) {
  z = rank(g, ties.method = "random")
  
  strata1 = floor((z-1)/q)+1
  
  # check GR statistic
  
  set.seed(125)
  id= seq(x)
  temp<- data.frame(x=x,strata1=strata1,id=id, g=g)
  temp<- arrange(.data=temp, x)
  temp<-group_by(.data=temp, strata1)
  temp<-mutate(.data=temp, x0q= rank(x, ties.method = "random"))
  temp<-arrange(.data=temp, id)
  x0q <- temp$x0q
  return(x0q)
}

#ensure participants have measured BMI and BMI prs

exclusions <-vroom("~/data/non-linear/geno/exclusions/data.minimal_relateds.plink.txt.gz", col_names = c("ieu", "IID"))
exclusions_highly <-vroom("~/data/non-linear/geno/exclusions/data.highly_relateds.plink.txt.gz", col_names = c("ieu", "IID"))

dat <- dat %>% 
  anti_join(exclusions) %>%
  anti_join(exclusions_highly)

initial <- dat %>%
  select(FID = ieu, IID = ieu, body_mass_index_bmi,bmi_prs)




initial %>%
  drop_na(body_mass_index_bmi, bmi_prs) %>%
  mutate(ranked_bmi = generate_ranked_strata(bmi_prs, body_mass_index_bmi, 10),
         resid_bmi = ntile(lm(body_mass_index_bmi ~ bmi_prs, data = .,)$resid, 10)) %>%
  mutate(extreme_ranked_bmi = if_else(ranked_bmi ==10 | ranked_bmi == 1, 1,0),
          extreme_resid_bmi = if_else(resid_bmi ==10 | resid_bmi == 1, 1,0)) %>%
  distinct(IID, .keep_all = T) %>%
  write_tsv("~/data/non-linear/pheno/stratified_pre_submission_tsv.gz")


