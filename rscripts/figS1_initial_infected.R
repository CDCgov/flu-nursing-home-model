## Set-up --------------------------------

# load packages
pacman::p_load(odin, dde, magrittr, tidyverse, here, patchwork)

# load eqns etc
source(here("rscripts", "helper_funs.R"))
source(here("rscripts", "eqns.R"))
source(here("rscripts", "parameters.R"))


## Baseline (i): Short-stay resident ---------------------

# 1. current guidelines
sims_df <- simulate_model(sir, n_sims = n_sims)

# 2. no prophy
sir$set_user(uptake_R = 0, uptake_HC = 0, uptake_R_prop = 0, uptake_HC_prop = 0)
sims_df_noP <- simulate_model(sir, n_sims = n_sims)

ssr <- sims_df_noP %>% mutate(scenario = "None") %>%
    rbind(., { sims_df %>% mutate(scenario = "Until 7d no cases") }) %>%
    mutate(seed = "Short-stay resident")


## Sensitivity (ii): Long-stay resident ---------------------

source(here("rscripts", "parameters.R"))

sir$set_user(S_ini_SR = 78  * N_factor,     A_ini_SR = 0,
             S_ini_LR = 22  * N_factor - 1, A_ini_LR = 1,
             S_ini_HC = 100 * N_factor,     A_ini_HC = 0)

# 1. current guidelines
sims_df <- simulate_model(sir, n_sims = n_sims)

# 2. no prophy
sir$set_user(uptake_R = 0, uptake_HC = 0, uptake_R_prop = 0, uptake_HC_prop = 0)
sims_df_noP <- simulate_model(sir, n_sims = n_sims)

lsr <- sims_df_noP %>% mutate(scenario = "None") %>%
    rbind(., { sims_df %>% mutate(scenario = "Until 7d no cases") }) %>%
    mutate(seed = "Long-stay resident")


## Sensitivity (iii): Healthcare personnel ---------------------

source(here("rscripts", "parameters.R"))

sir$set_user(S_ini_SR = 78  * N_factor,     A_ini_SR = 0,
             S_ini_LR = 22  * N_factor,     A_ini_LR = 0,
             S_ini_HC = 100 * N_factor - 1, A_ini_HC = 1)

# 1. current guidelines
sims_df <- simulate_model(sir, n_sims = n_sims)

# 2. no prophy
sir$set_user(uptake_R = 0, uptake_HC = 0, uptake_R_prop = 0, uptake_HC_prop = 0)
sims_df_noP <- simulate_model(sir, n_sims = n_sims)

hcp <- sims_df_noP %>% mutate(scenario = "None") %>%
    rbind(., { sims_df %>% mutate(scenario = "Until 7d no cases") }) %>%
    mutate(seed = "Healthcare personnel")


## Plot output ---------------

tmp <- rbind(ssr, lsr, hcp) %>% 
    group_by(sim, scenario, seed) %>% 
    filter(step == max(step)) %>%
    ungroup() %>%
    mutate(symp_R  = symptomatic_R,
           hosp_R  = H_SR + HP_SR + H_LR + HP_LR)


cols <- RColorBrewer::brewer.pal(3, "Spectral")


# Final figure S1
pS1 <- 
  tmp %>% filter(scenario == "Until 7d no cases") %>%
  group_by(sim, seed) %>% 
  filter(step == max(step)) %>%
  filter(symptomatic_R > 2) %>%
  ungroup() %>%
  mutate(hosp_R  = H_SR + HP_SR + H_LR + HP_LR) %>%
  ggplot(aes(x = seed, y = symptomatic_R, fill = seed)) +
  geom_violin(alpha = 0.5) + geom_boxplot(width = 0.15) + geom_jitter(width = 0.01, alpha = 0.25) + 
  scale_fill_manual(guide = "none", values = cols) +
  scale_y_continuous(limits = c(0, 80) ) +
  get_theme(txt = 9) + labs(x = "Initial infected individual", y = "Total symptomatic illnesses")

