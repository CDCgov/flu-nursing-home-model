## Load packages / eqns / helper functions -------------------------------

pacman::p_load(odin, dde, magrittr, tidyverse, patchwork, here)

# load eqns, default parameters, and post-hoc analysis fns
source(here("rscripts", "helper_funs.R"))
source(here("rscripts", "eqns.R"))
source(here("rscripts", "parameters.R"))


## Run simulations: different uptake / adherence values ------------------------

av_grid <- expand.grid(uptake    = seq(0.05, 0.95, by = 0.1),
                       adherence = c(5, 1e10))

nvals  <- length(unique(av_grid$uptake))


## HCP ------

print("Varying HCP uptake/adherence...")

res <- list()

for (r in 1:nrow(av_grid)){
  
  sir$set_user(uptake_HC = prop_to_rate(av_grid$uptake[r]), 
               uptake_HC_prop = av_grid$uptake[r],
               c_stop_HC = 1/av_grid$adherence[r])
  
  res[[r]] <- simulate_model(sir, n_sims = n_sims) %>% 
    mutate(uptake = av_grid$uptake[r], adherence = av_grid$adherence[r])
}

res <- bind_rows(res)

tmp_HC <- res %>% 
  group_by(sim, uptake, adherence) %>% 
  filter(step == max(step)) %>%
  ungroup() %>%
  mutate(symp_R  = symptomatic_R,
         hosp_R  = H_SR + HP_SR + H_LR + HP_LR) 

rm(res)

## Residents ----

print("Varying resident uptake/adherence...")

# reset HC values
sir$set_user(uptake_HC = prop_to_rate(0.65), uptake_HC_prop = 0.65, c_stop_HC = 1/5)

res <- list()

for (r in 1:nrow(av_grid)){
  sir$set_user(uptake_R = prop_to_rate(av_grid$uptake[r]), 
               uptake_R_prop = av_grid$uptake[r],
               c_stop_R = 1/av_grid$adherence[r])
  
  res[[r]] <- simulate_model(sir, n_sims = n_sims) %>% 
    mutate(uptake = av_grid$uptake[r], adherence = av_grid$adherence[r])
}

res <- bind_rows(res)

tmp_R <- res %>% 
  group_by(sim, uptake, adherence) %>% 
  filter(step == max(step)) %>%
  ungroup() %>%
  mutate(symp_R  = symptomatic_R,
         hosp_R  = H_SR + HP_SR + H_LR + HP_LR) 



## Plot output ----------------

adh_label <- function(label) {
  paste0(label, "d adherence")
}

mycols <- RColorBrewer::brewer.pal(nvals, "Spectral")

linesize <- 1
ptsize   <- 2.2


## HCP (short adherence) ----

tmp_HC_short <- tmp_HC %>% filter(adherence == 5)

p_cases <- lapply(seq(10, 80, by = 5), function(x) {
  tmp_HC_short %>% group_by(uptake, adherence) %>%
    summarize(prop = round(length(symp_R[symp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = as.factor(uptake * 100))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + get_theme() +
  scale_color_manual("HCP\nuptake (%)", values = mycols) +
  #facet_wrap(~ adherence, labeller = labeller(adherence = adh_label)) +
  labs(x = "Total illnesses > n", y = "% Simulations")


p_hosps <- lapply(seq(0, 20, by = 1), function(x) {
  tmp_HC_short %>% group_by(uptake, adherence) %>%
    summarize(prop = round(length(hosp_R[hosp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = as.factor(uptake * 100))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + get_theme() +
  scale_color_manual("HCP\nuptake (%)", values = mycols) +
  #facet_wrap(~ adherence, labeller = labeller(adherence = adh_label)) +
  labs(x = "Total hospitalizations > n", y = "% Simulations") 


pH_short <- p_cases + p_hosps


## HCP (long adherence) ----------- 

tmp_HC_long <- tmp_HC %>% filter(adherence == 1e10)

p_cases <- lapply(seq(10, 80, by = 5), function(x) {
  tmp_HC_long %>% group_by(uptake, adherence) %>%
    summarize(prop = round(length(symp_R[symp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = as.factor(uptake * 100))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + get_theme() +
  scale_color_manual("HCP\nuptake (%)", values = mycols) +
  #facet_wrap(~ adherence, labeller = labeller(adherence = adh_label)) +
  labs(x = "Total illnesses > n", y = "% Simulations")


p_hosps <- lapply(seq(0, 20, by = 1), function(x) {
  tmp_HC_long %>% group_by(uptake, adherence) %>%
    summarize(prop = round(length(hosp_R[hosp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = as.factor(uptake * 100))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + get_theme() +
  scale_color_manual("HCP\nuptake (%)", values = mycols) +
  #facet_wrap(~ adherence, labeller = labeller(adherence = adh_label)) +
  labs(x = "Total hospitalizations > n", y = "% Simulations") 


pH_long <- p_cases + p_hosps


## Residents (long adherence) ------

tmp_R_long <- tmp_R %>% filter(adherence == 1e10)

p_cases <- lapply(seq(10, 80, by = 5), function(x) {
  tmp_R_long %>% group_by(uptake, adherence) %>%
    summarize(prop = round(length(symp_R[symp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = as.factor(uptake * 100))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + get_theme() +
  scale_color_manual("Resident\nuptake (%)", values = mycols) +
  #facet_wrap(~ adherence, labeller = labeller(adherence = adh_label)) +
  labs(x = "Total illnesses > n", y = "% Simulations")

p_hosps <- lapply(seq(0, 20, by = 1), function(x) {
  tmp_R_long %>% group_by(uptake, adherence) %>%
    summarize(prop = round(length(hosp_R[hosp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = as.factor(uptake * 100))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + get_theme() +
  scale_color_manual("Resident\nuptake (%)", values = mycols) +
  #facet_wrap(~ adherence, labeller = labeller(adherence = adh_label)) +
  labs(x = "Total hospitalizations > n", y = "% Simulations") 

pR_long <- p_cases + p_hosps + plot_layout(guides = "collect") + plot_annotation("Resident uptake, perfect compliance")


# Final figure 4
p4 <- pR_long / (pH_short / pH_long  + plot_layout(guides = "collect")) + 
  plot_annotation(tag_levels = 'A') + plot_layout(heights = c(1, 2))


