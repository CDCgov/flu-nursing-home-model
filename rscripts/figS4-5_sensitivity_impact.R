## Set-up --------------------------------

# load packages
pacman::p_load(odin, dde, magrittr, tidyverse, patchwork, here)

# load eqns
source(here("rscripts", "helper_funs.R"))
source(here("rscripts", "eqns.R"))


## Set variables / values to investigate ---------------------

# CHOOSE
# (1)
cname  <- "R0" # option 1: R0 ( = Figure S4)
cvals  <- 3:8  # possible R0 values

# OR
# (2)
cname <- "Pe" # option 2: antiviral effect on susceptibility (= Figure S5)
cvals <- seq(0.3, 0.8, by = 0.1)


combined_list <- list()

for(r in 1:length(cvals)) {
  
  print(paste(cname, "=", cvals[r]))
  
  # reset to default values
  source(here("rscripts", "parameters.R"))
  
  if (cname == "R0") {
    # reset beta based on new R0 val
    beta <- cvals[r] / max( eigen(contact_matrix * pop_matrix / gamma)$values )
    
    sir$set_user(beta = beta)
  }
  
  if (cname == "SSR") {
    sir$set_user(S_ini_SR = cvals[r] * N_factor - 1, 
                 S_ini_LR = (100 - cvals[r])  * N_factor)
  }

  if (cname == "Pe") {
    sir$set_user(Pe_on = cvals[r])
  }
  
  
  # simulate different prophylaxis scenarios
  # 1. current guidelines
  sims_df <- simulate_model(sir, n_sims = n_sims)
  
  # 2. no prophy
  sir$set_user(uptake_R = 0, uptake_HC = 0, uptake_R_prop = 0, uptake_HC_prop = 0)
  sims_df_noP <- simulate_model(sir, n_sims = n_sims)
  
  combined_list[[r]] <- sims_df_noP %>% mutate(scenario = "None") %>%
                        rbind(., { sims_df %>% mutate(scenario = "Until 7d no cases") }) %>% 
                        mutate(parval  = cvals[r],
                               parname = paste(cname, "=", cvals[r]))
}

combined <- bind_rows(combined_list)


## Plot output ---------------

# set color scale to match scenarios in main figure 1
mycols <- RColorBrewer::brewer.pal(8, "Spectral")[c(1, 7)]


tmp <- combined %>% 
  group_by(sim, scenario, parval, parname) %>% 
  filter(step == max(step))  %>%
  ungroup() %>%
  mutate(symp_R  = symptomatic_R,
         hosp_R  = H_SR + HP_SR + H_LR + HP_LR) 

# change Pe variables for nicer plotting
if (cname == "Pe") {
  tmp <- tmp %>% 
    mutate(parname = paste0("Reduction = ", parval * 100, "%"),
           parval  = parval * 100)
  xtitle = "Antiviral reduction in susceptibility (%)"
} else if (cname == "R0") xtitle = "Basic reproduction number, R0"

# customize position_dodge width depending on parameter magnitudes
wdth <- unique(tmp$parval)[1] * 0.1

# plot average impact of prophy vs no prophy
p_average <- tmp %>% filter(scenario %in% c("None", "Until 7d no cases")) %>%
  filter(symp_R != 0) %>%
  group_by(scenario, parval) %>%
  summarize(med = median(symp_R),
            lo = quantile(symp_R, 0.025),
            hi = quantile(symp_R, 0.975)) %>%
  ggplot(aes(x = parval, y = med, colour = scenario)) +
  geom_point(size = 2.5, position = position_dodge(width = wdth)) + 
  geom_errorbar(aes(ymin = lo, ymax = hi), width = wdth * 0.5,
                position = position_dodge(width = wdth)) + 
  scale_color_manual("Prophylaxis", values = mycols) +
  get_theme(legend.position = "top") + 
  labs(y = "Symptomatic illnesses", x = xtitle)


# plot symptomatic cases outbreak size
p_cases <- 
  lapply(seq(10, 80, by = 5), function(x) {
  tmp %>% group_by(scenario, parname) %>%
    summarize(prop = round(length(symp_R[symp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = scenario)) + 
  geom_point(size = 1.5) + geom_line() + 
  facet_wrap(~ parname) + get_theme(legend.position = "top") +
  scale_color_manual(guide = "none", values = mycols) +
  labs(x = "Outbreak size > n", y = "% Simulations")


# plot hospitalizations outbreak size
p_hosps <- 
  lapply(seq(0, 20, by = 1), function(x) {
  tmp %>% group_by(scenario, parname) %>%
    summarize(prop = round(length(hosp_R[hosp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() %>%
  ggplot(aes(x = size, y = prop, color = scenario)) + 
  geom_point(size = 1.5) + geom_line() + 
  facet_wrap(~ parname) + get_theme() +
  scale_color_manual(guide = "none", values = mycols) +
  labs(x = "Total hospitalizations > n", y = "% Simulations") 

# Final figure (S4 for R0, S5 for Pe)
pS4_5 <- p_average / p_cases / p_hosps + 
            plot_annotation(tag_levels = "A") + plot_layout(heights = c(0.7, 1, 1))

