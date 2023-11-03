## Set-up --------------------------------

# load packages
pacman::p_load(odin, dde, magrittr, tidyverse, patchwork, here)

# load ODES, post-simulation helper functions, and baseline parameter settings
source(here("rscripts", "helper_funs.R"))
source(here("rscripts", "eqns.R"))
source(here("rscripts", "parameters.R"))


## Parameter values to explore ---------------------------

# will express as a percentage of the within-NH transmission rate
v_import <- seq(0, 0.3, by = 0.05)


## Run sims across different parameter vectors --------------------------------

tmp <- list()


for (i in 1:length(v_import)) {
  sir$set_user(import_rate = v_import[i] * beta)
  
  tmp[[i]] <- simulate_model(sir) %>% mutate(import = v_import[i])
  
}

sims <- bind_rows(tmp)


## Plot results -------------------------------

cols <- RColorBrewer::brewer.pal(length(v_import), "Spectral")

# Final figure S2
pS2 <- sims %>% 
  group_by(sim, import) %>% 
  filter(step == max(step)) %>%
  filter(symptomatic_R > 2) %>%
  ungroup() %>%
  mutate(hosp_R  = H_SR + HP_SR + H_LR + HP_LR) %>%
  ggplot(aes(x = import, y = symptomatic_R, fill = as.factor(import))) +
  geom_violin(alpha = 0.5) + geom_boxplot(width = 0.02) + geom_jitter(width = 0.001, alpha = 0.25) + 
  scale_fill_manual(guide = "none", values = cols) +
  get_theme() + labs(x = "Fraction of within-facility transmission", y = "Total symptomatic illnesses")

