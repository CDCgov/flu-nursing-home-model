## Set-up --------------------------------

# load packages
pacman::p_load(odin, dde, magrittr, tidyverse, here, patchwork)

# load eqns etc
source(here("rscripts", "helper_funs.R"))
source(here("rscripts", "eqns.R"))
source(here("rscripts", "parameters.R"))


## Baseline (i): Simulate model with and without prophylaxis ---------------------

# 1. current guidelines
sims_df <- simulate_model(sir, n_sims = n_sims)
  
# 2. no prophy
sir$set_user(uptake_R = 0, uptake_HC = 0, uptake_R_prop = 0, uptake_HC_prop = 0)
sims_df_noP <- simulate_model(sir, n_sims = n_sims)
  
combined <- sims_df_noP %>% mutate(scenario = "None") %>%
    rbind(., { sims_df %>% mutate(scenario = "Until 7d no cases") }) %>%
    mutate(delay = "No delay")


## Sensitivity (ii): Simulate model with and without prophylaxis ---------------------

source(here("rscripts", "eqns_HCP_isolation_delay.R"))
source(here("rscripts", "parameters.R"))

# 1. current guidelines
sims_df <- simulate_model(sir, n_sims = n_sims)

# 2. no prophy
sir$set_user(uptake_R = 0, uptake_HC = 0, uptake_R_prop = 0, uptake_HC_prop = 0)
sims_df_noP <- simulate_model(sir, n_sims = n_sims)

combined_delay <- sims_df_noP %>% mutate(scenario = "None") %>%
  rbind(., { sims_df %>% mutate(scenario = "Until 7d no cases") }) %>%
  mutate(delay = "One day delay")


## Plot output ---------------

tmp <- rbind(combined, combined_delay) %>% 
  group_by(sim, scenario, delay) %>% 
  filter(step == max(step)) %>%
  ungroup() %>%
  mutate(symp_R  = symptomatic_R,
         hosp_R  = H_SR + HP_SR + H_LR + HP_LR)

# keep colours consistent with scenarios in main text figures
mycols <- RColorBrewer::brewer.pal(8, "Spectral")[c(1, 7)]

wdth <- 0.1

p_average <- tmp %>% filter(scenario %in% c("None", "Until 7d no cases")) %>%
  filter(symp_R != 0) %>%
  group_by(scenario, delay) %>%
  summarize(med = median(symp_R),
            lo = quantile(symp_R, 0.025),
            hi = quantile(symp_R, 0.975)) %>%
  ggplot(aes(x = delay, y = med, colour = scenario)) +
  geom_point(size = 2.5, position = position_dodge(width = wdth)) + 
  geom_errorbar(aes(ymin = lo, ymax = hi), width = wdth * 0.5,
                position = position_dodge(width = wdth)) + 
  scale_color_manual("Prophylaxis", values = mycols) +
  get_theme(legend.position = "top") + 
  labs(y = "Symptomatic illnesses", x = NULL)


# plot symptomatic cases outbreak size for all scenarios
p_cases <- lapply(seq(10, 80, by = 5), function(x) {
  tmp %>% group_by(scenario, delay) %>%
    summarize(prop = round(length(symp_R[symp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>% 
  filter(scenario %in% c("None", "Until 7d no cases")) %>%
  ggplot(aes(x = size, y = prop, color = scenario)) + 
  geom_point(size = 1.5) + geom_line() + 
  get_theme(legend.position = "top") +
  facet_wrap(~ delay) +
  scale_color_manual("Prophylaxis", values = mycols) +
  labs(x = "Total illnesses > n", y = "% Simulations") 

# plot hospitalization outbreak size for all scenarios
p_hosps <- lapply(seq(0, 20, by = 1), function(x) {
  tmp %>% group_by(scenario, delay) %>%
    summarize(prop = round(length(hosp_R[hosp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>% 
  filter(scenario %in% c("None", "Until 7d no cases")) %>%
  ggplot(aes(x = size, y = prop, color = scenario)) + 
  geom_point(size = 1.5) + geom_line() + 
  get_theme() +
  facet_wrap(~ delay) +
  scale_color_manual(guide = "none", values = mycols) +
  labs(x = "Total hospitalizations > n", y = "% Simulations") 

# Final figure S6
pS6 <- p_average / p_cases / p_hosps + plot_annotation(tag_levels = "A")

