## Load packages / eqns / helper functions -------------------------------

pacman::p_load(odin, dde, magrittr, tidyverse, patchwork, here)

# load eqns, default parameters, and post-hoc analysis fns
source(here("rscripts", "helper_funs.R"))
source(here("rscripts", "eqns.R"))
source(here("rscripts", "parameters.R"))


## Run simulations: different duration scenarios -------------------------

# Prophy according to current guidelines (2 cases identified)
# relax 'within 72hrs' rule for fairer comparison at higher initiation thresholds
sir$set_user(start_lag = 8)

sims_df <- simulate_model(sir, n_sims = n_sims)

# Prophy when 1 case identified
sir$set_user(case_threshold = 1)

sims_df1 <- simulate_model(sir, n_sims = n_sims)

# Prophy when 3 cases identified
sir$set_user(case_threshold = 3)

sims_df3 <- simulate_model(sir, n_sims = n_sims)


# Prophy when 4 cases identified
sir$set_user(case_threshold = 4)

sims_df4 <- simulate_model(sir, n_sims = n_sims)


# Prophy when 5 cases identified
sir$set_user(case_threshold = 5)

sims_df5 <- simulate_model(sir, n_sims = n_sims)


# Prophy when 6 cases identified
sir$set_user(case_threshold = 6)

sims_df6 <- simulate_model(sir, n_sims = n_sims)


# Prophy when 7 cases identified
sir$set_user(case_threshold = 7)

sims_df7 <- simulate_model(sir, n_sims = n_sims)



# Combine output for later analysis
combined_df <- sims_df %>% mutate(scenario = "2 cases") %>% 
  rbind(., { sims_df1  %>% mutate(scenario = "1 case") }) %>% 
  rbind(., { sims_df3  %>% mutate(scenario = "3 cases") }) %>% 
  #rbind(., { sims_df4  %>% mutate(scenario = "4 cases") }) %>% 
  rbind(., { sims_df5  %>% mutate(scenario = "5 cases") }) %>% 
  #rbind(., { sims_df6  %>% mutate(scenario = "6 cases") }) %>% 
  rbind(., { sims_df7  %>% mutate(scenario = "7 cases") }) %>% 
  mutate(scenario = factor(scenario, levels = c("1 case", paste0(c(2:7), " cases"))))


## Example comparison of outputs using summary functions ----------------

# % simulations above a certain size
tmp <- combined_df %>%  
          mutate(symp_R  = symptomatic_R,
                 hosp_R  = H_SR + HP_SR + H_LR + HP_LR) %>%
          group_by(sim, scenario) %>% 
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
  geom_point(size = 1.5) + geom_line(size = 1.2) + 
  get_theme(legend.position = c(0.8, 0.58)) +
  scale_color_brewer("Initiation\nthreshold", palette = "Spectral") +
  labs(x = "Total hospitalizations > n", y = "% Simulations")


# peak cases
peak <- combined_df %>%  
  mutate(R = I_SR + IP_SR + Tr_SR + I_LR + IP_LR + Tr_LR,
         H = I_HC + IP_HC + Tr_HC) 
  
peakR <- peak %>% group_by(sim, scenario) %>%
  # consider only outbreaks that take off (~ e.g. with > 2 resident illnesses)
  filter(max(symptomatic_R) > 2) %>%
  filter(R == max(R)) %>%
  ungroup() %>% group_by(scenario) %>%
  summarise(md_day = median(step),
            lo_day = quantile(step, 0.025),
            hi_day = quantile(step, 0.975),
            md_sym = median(R),
            lo_sym = quantile(R, 0.025),
            hi_sym = quantile(R, 0.975)) %>%
  ungroup()

peak_day <- ggplot(peakR, aes(x = scenario, y = md_day, fill = scenario)) + 
  geom_bar(stat = "identity", width = 0.75, color = "black", size = 0.1) + 
  geom_errorbar(aes(ymin = lo_day, ymax = hi_day), width = 0.2) +
  scale_fill_brewer(guide = "none", palette = "Spectral") + get_theme() +
  labs(x = "Initiation threshold", y = "Day of peak symptomatic illness")

peak_sym <- ggplot(peakR, aes(x = scenario, y = md_sym, fill = scenario)) + 
  geom_bar(stat = "identity", width = 0.75, color = "black", size = 0.1) + 
  geom_errorbar(aes(ymin = lo_sym, ymax = hi_sym), width = 0.2) +
  scale_fill_brewer(guide = "none", palette = "Spectral") + get_theme() +
  labs(x = "Initiation threshold", y = "Peak symptomatic illness")


# Final figure 3
p3 <- (p_cases + p_hosps) / (peak_day + peak_sym) + plot_annotation(tag_levels = "A")


