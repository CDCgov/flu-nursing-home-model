## Load packages / eqns / helper functions -------------------------------

pacman::p_load(odin, dde, magrittr, tidyverse, patchwork, here)

# load helper fns (pre / post-hoc wrangling and model simulation)
source(here("rscripts", "helper_funs.R"))
source(here("rscripts", "eqns.R"))
source(here("rscripts", "parameters.R"))


## Run simulations: different initiation scenarios -------------------------

# Prophy according to current guidelines
sims_df <- simulate_model(sir, n_sims = n_sims)

# Prophy until one day after last case
sir$set_user(stop_lag = 1)

sims_df1 <- simulate_model(sir, n_sims = n_sims)

# Prophy until 3 days after last case
sir$set_user(stop_lag = 3)

sims_df3 <- simulate_model(sir, n_sims = n_sims)


# Prophy until 5 days after last case
sir$set_user(stop_lag = 5)

sims_df5 <- simulate_model(sir, n_sims = n_sims)


# Prophy until 10 days after last case
sir$set_user(stop_lag = 10)

sims_df_extra <- simulate_model(sir, n_sims = n_sims)


# Prophy for max 10d
sir$set_user(stop_lag = 7, max_prophy_duration = 10, min_prophy_duration = 7)

sims_df10 <- simulate_model(sir, n_sims = n_sims)


# Prophy for max 14d
sir$set_user(stop_lag = 7, max_prophy_duration = 14, min_prophy_duration = 7)

sims_df14 <- simulate_model(sir, n_sims = n_sims)


# No prophy
sir$set_user(uptake_R = 0, uptake_HC = 0, uptake_R_prop = 0, uptake_HC_prop = 0)

sims_df_noP <- simulate_model(sir, n_sims = n_sims)


# Combine output for later analysis
combined_df <- sims_df_noP %>% mutate(scenario = "None") %>%
  rbind(., { sims_df   %>% mutate(scenario = "Until 7d no cases") }) %>% 
  rbind(., { sims_df1  %>% mutate(scenario = "Until 1d no cases") }) %>% 
  rbind(., { sims_df3  %>% mutate(scenario = "Until 3d no cases") }) %>% 
  rbind(., { sims_df5  %>% mutate(scenario = "Until 5d no cases") }) %>% 
  rbind(., { sims_df_extra  %>% mutate(scenario = "Until 10d no cases") }) %>%
  rbind(., { sims_df10 %>% mutate(scenario = "Fixed 10d") }) %>% 
  rbind(., { sims_df14 %>% mutate(scenario = "Fixed 14d") }) %>% 
  mutate(scenario = factor(scenario, levels = c("None", paste0("Fixed ", c(10, 14), "d"), 
                                                paste0("Until ", c(1, 3, 5, 7, 10), "d no cases"))))


## Example comparison of outputs using summary functions ----------------

# colour scale chosen to keep scenario colors consistent with later plots

# get total # scenarios (for total # levels in scale)
n_scenarios <- length(unique(combined_df$scenario))
# find which positions correspond to the no prophy / full prophy scenarios
whichlevels <- which(sort(unique(combined_df$scenario)) %in% c("None", 
                                                               "Until 7d no cases"))
# get the corresponding color scale
comp_cols <- RColorBrewer::brewer.pal(n_scenarios, "Spectral")[whichlevels]


# summary stats: no prophy vs current guidance prophy
comp <- compare_stats(sims_df_noP, sims_df, cols = comp_cols)

comp$table

comp$change

p_violin <- comp$plot 

p_violin2 <- p_violin + scale_x_discrete(labels = c("None", "Until 7d\nno cases\n(current)"))


# symptomatic infection time-series
set.seed(123)

p_symp_mean <- rbind(sims_df_noP %>% mutate(scenario = "None"), 
                     sims_df     %>% mutate(scenario = "Until 7d\nnocases")) %>%
  mutate(Residents = I_SR + IP_SR + Tr_SR +  I_LR + IP_LR + Tr_LR,
         HCP = I_HC + IP_HC + Tr_HC) %>%
  group_by(step, scenario) %>% summarize(Residents = mean(Residents), HCP = mean(HCP)) %>%
  ungroup() %>% gather(var, val, Residents, HCP) %>%
  mutate(sim = 1)

p_symp <- rbind(sims_df_noP %>% mutate(scenario = "None"), 
                sims_df     %>% mutate(scenario = "Until 7d\nnocases")) %>%
  filter(sim %in% sample(1:n_sims, replace = FALSE, size = n_sample_sims)) %>%
  mutate(Residents = I_SR + IP_SR + Tr_SR +  I_LR + IP_LR + Tr_LR,
         HCP = I_HC + IP_HC + Tr_HC) %>%
  gather(var, val, Residents, HCP) %>%
  ggplot(aes(x = step, y = val, color = scenario, group = interaction(sim, scenario))) + 
  geom_line(alpha = 0.15) + 
  geom_line(data = p_symp_mean, size = 1.1) +
  scale_color_manual(guide = "none", values = comp_cols) +
  facet_wrap(~ var, ncol = 1) + get_theme() +
  labs(x = "Days", y = "Number symptomatic")


# % simulations above a certain size
tmp <- combined_df %>%  
  mutate(symp_R  = symptomatic_R,
         hosp_R  = H_SR + HP_SR + H_LR + HP_LR) %>%
  group_by(sim, scenario) %>% 
  
  # filter to the day outbreak ends (so cases etc are only counted up until that point)
  filter(step == max(step)) %>%
  ungroup() 

p_cases <- lapply(seq(10, 80, by = 5), function(x) {
  tmp %>% group_by(scenario) %>%
    summarize(prop = round(length(symp_R[symp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = scenario)) + 
  geom_point(size = 1.5) + geom_line(size = 1.2) + get_theme() +
  scale_color_brewer(guide = "none", palette = "Spectral") +
  labs(x = "Total illnesses > n", y = "% Simulations")


p_hosps <- lapply(seq(0, 20, by = 1), function(x) {
  tmp %>% group_by(scenario) %>%
    summarize(prop = round(length(hosp_R[hosp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = scenario)) + 
  geom_point(size = 1.5) + geom_line(size = 1.2) + get_theme() +
  scale_color_brewer("Prophylaxis\nduration", palette = "Spectral") +
  labs(x = "Total hospitalizations > n", y = "% Simulations")


# Final figure 2
p2 <- (p_symp + p_violin2) / (p_cases + p_hosps) + 
  plot_annotation(tag_levels = "A") + plot_layout(guides = "collect")

