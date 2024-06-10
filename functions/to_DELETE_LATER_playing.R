


run_glesjer_confound <- function(bgux, seed = 123, bux =0.5, buy = 0.5, bv = 0) {

s <- sim_mydata_outcomes_confounder(bgux = bgux,bux = bux, buy = buy,bv = bv, seed = runif(1,1,1e5) + bgux + seed)
d <-generate_all_sumstats(s)[2]$ss

get_pvals = function(df){
  p_quadratic <- metafor::rma(beta_mr ~ mean, (se_mr)^2, method="FE", data = df, control=list(maxiter=1000))$pval[2]
  p_Q <- 1 - pchisq(metafor::rma(beta_mr, vi=(se_mr)^2, data = df, control=list(maxiter=1000))$QE, df=(9))
  tibble(quadratic_p = p_quadratic, Q_p = p_Q)
}


p_res <- d %>%
    filter(strata_method != "full") %>%
    group_by(strata_method) %>%
    nest() %>%
    mutate(res = map(data, get_pvals)) %>%
    unnest(res) %>%
    ungroup()


ranked_ps <- p_res %>%
   filter(strata_method == "ranked") %>%
   mutate(quadratic_p_rank = if_else(quadratic_p ==0, 1e-300, quadratic_p)) %>%
   transmute(Q_p_rank = if_else(Q_p ==0, 0, Q_p))

resid_ps <- p_res %>%
    filter(strata_method == "residual") %>%
    mutate(quadratic_p_res = if_else(quadratic_p ==0, 1e-300, quadratic_p)) %>%
    transmute(Q_p_res= if_else(Q_p == 0, 0, Q_p))

bind_cols(resid_ps, ranked_ps, tibble(glesjer = skedastic::glejser(lm(x ~ g, data = s))$p.value), bgux = bgux, buy = buy, bux = bux) }


run_glesjer_confound_seperate_interact <- function(bgux, seed = 123, bux = 0.3, buy = 0, bv = 1) {
  
  s <- sim_mydata_outcomes_confounder(bgux = bgux,bux = bux, buy = buy,bv = bv, seed = runif(1,1,1e5) + bgux*1000 + seed)
  d <-generate_all_sumstats(s)[2]$ss
  
  get_pvals = function(df){
    p_quadratic <- metafor::rma(beta_mr ~ mean, (se_mr)^2, method="FE", data = df, control=list(maxiter=1000))$pval[2]
    p_Q <- 1 - pchisq(metafor::rma(beta_mr, vi=(se_mr)^2, data = df, control=list(maxiter=1000))$QE, df=(9))
    tibble(quadratic_p = p_quadratic, Q_p = p_Q)
  }
  
  
  p_res <- d %>%
    filter(strata_method != "full") %>%
    group_by(strata_method) %>%
    nest() %>%
    mutate(res = map(data, get_pvals)) %>%
    unnest(res) %>%
    ungroup()
  
  
  ranked_ps <- p_res %>%
    filter(strata_method == "ranked") %>%
    mutate(quadratic_p_rank = if_else(quadratic_p ==0, 1e-300, quadratic_p)) %>%
    transmute(Q_p_rank = if_else(Q_p ==0, 0, Q_p))
  
  resid_ps <- p_res %>%
    filter(strata_method == "residual") %>%
    mutate(quadratic_p_res = if_else(quadratic_p ==0, 1e-300, quadratic_p)) %>%
    transmute(Q_p_res= if_else(Q_p == 0, 0, Q_p))
  
  bind_cols(resid_ps, ranked_ps, tibble(glesjer = skedastic::glejser(lm(x ~ g, data = s))$p.value), bgux = bgux, buy = buy, bux = bux) }

run_glesjer_variance <- function(g_sd, seed = 123, buy = 1, bux = 1, bx = 0.3) {
  
  s <- sim_mydata_outcomes_vQTL(g_sd = g_sd,bux = bux, buy = buy,bx = bx, seed = runif(1,1,1e5) + g_sd*1000 + seed)
  d <-generate_all_sumstats(s)[2]$ss
  
  get_pvals = function(df){
    p_quadratic <- metafor::rma(beta_mr ~ mean, (se_mr)^2, method="FE", data = df, control=list(maxiter=1000))$pval[2]
    p_Q <- 1 - pchisq(metafor::rma(beta_mr, vi=(se_mr)^2, data = df, control=list(maxiter=1000))$QE, df=(9))
    tibble(quadratic_p = p_quadratic, Q_p = p_Q)
  }
  
  
  p_res <- d %>%
    filter(strata_method != "full") %>%
    group_by(strata_method) %>%
    nest() %>%
    mutate(res = map(data, get_pvals)) %>%
    unnest(res) %>%
    ungroup()
  
  
  ranked_ps <- p_res %>%
    filter(strata_method == "ranked") %>%
    mutate(quadratic_p_rank = if_else(quadratic_p ==0, 1e-300, quadratic_p)) %>%
    transmute(Q_p_rank = if_else(Q_p ==0, 0, Q_p))
  
  resid_ps <- p_res %>%
    filter(strata_method == "residual") %>%
    mutate(quadratic_p_res = if_else(quadratic_p ==0, 1e-300, quadratic_p)) %>%
    transmute(Q_p_res= if_else(Q_p == 0, 0, Q_p))
  
  bind_cols(resid_ps, ranked_ps, tibble(glesjer = skedastic::glejser(lm(x ~ g, data = s))$p.value), g_sd = g_sd, buy = buy, bux = bux, bx = bx) }



glesjer_confound_res <- map_dfr(seq(-1,1, 0.01), run_glesjer_confound_seperate_interact)


glesjer_confound_res %>%
  ggplot(aes(x = Q_p_rank, y =glesjer)) +
  geom_point(aes(colour = abs(bgux))) +
  scale_x_log10() +
  scale_y_log10()

glesjer_variance <- map_dfr(seq(0,0.5, 0.003), run_glesjer_variance)



glesjer_variance %>%
  ggplot(aes(x = Q_p_rank, y =glesjer)) +
  geom_point(aes(colour = abs(g_sd))) +
  scale_x_log10() +
  scale_y_log10()
