## Set-up --------------------------------

# load packages
pacman::p_load(odin, dde, magrittr, tidyverse, patchwork, here)

# load ODES, post-simulation helper functions, and baseline parameter settings
source(here("rscripts", "helper_funs.R"))
source(here("rscripts", "eqns.R"))
source(here("rscripts", "parameters.R"))


## Parameter values to explore ---------------------------

nvalues <- 6

vR0  <- seq(3, 8,     length.out = nvalues)    # R0
vCSS <- seq(1, 6,     length.out = nvalues)    # HCW-HCW contact rates
vCSR <- seq(6, 11,    length.out = nvalues)    # HCW-Resident contact rates
vIPC <- seq(0, 0.5,   length.out = nvalues)    # effect of IPC
vPe  <- seq(0.3, 0.8, length.out = nvalues)    # antiviral effect on susceptibility
vPei <- seq(0, 0.5,   length.out = nvalues)    # antiviral effect on transmission
vSSR <- seq(30, 80,   length.out = nvalues)    # No. short stay residents


## Run sims across different parameter vectors --------------------------------

tmp <- list()


# R0 ---------------------

print("R0")

for (i in 1:nvalues) {
  beta <- vR0[i] / max( eigen(contact_matrix * pop_matrix / gamma)$values )
  
  sir$set_user(beta = beta)
  
  tmp[[i]] <- simulate_model(sir) %>% mutate(param = "R0", value = vR0[i])
  
}

R0sims <- bind_rows(tmp)

# reset original values
beta <- R0 / max( eigen(contact_matrix * pop_matrix / gamma)$values )
sir$set_user(beta = beta)


# C_SS ---------------------

print("C_SS")

for (i in 1:nvalues) {
  contact_matrix <- rbind(c(C_RR, C_RS), c(C_SR, vCSS[i]))
  
  beta <- R0 / max( eigen(contact_matrix * pop_matrix / gamma)$values )
  
  sir$set_user(beta = beta, C_SS = vCSS[i])
  
  tmp[[i]] <- simulate_model(sir) %>% mutate(param = "C_SS", value = vCSS[i])
  
}

CSSsims <- bind_rows(tmp)

# reset original values
contact_matrix <- rbind(c(C_RR, C_RS), c(C_SR, C_SS))
beta <- R0 / max( eigen(contact_matrix * pop_matrix / gamma)$values )
sir$set_user(beta = beta, C_SS = C_SS)


# C_SR ---------------------

print("C_SR")

for (i in 1:nvalues) {
  contact_matrix <- rbind(c(C_RR, vCSR[i]), c(vCSR[i], C_SS))
  
  beta <- R0 / max( eigen(contact_matrix * pop_matrix / gamma)$values )
  
  sir$set_user(beta = beta, C_RS = vCSR[i], C_SR = vCSR[i])
  
  tmp[[i]] <- simulate_model(sir) %>% mutate(param = "C_SR", value = vCSR[i])
  
}

CSRsims <- bind_rows(tmp)

# reset original values
contact_matrix <- rbind(c(C_RR, C_RS), c(C_SR, C_SS))
beta <- R0 / max( eigen(contact_matrix * pop_matrix / gamma)$values )
sir$set_user(beta = beta, C_RS = C_RS, C_SR = C_SR)


# IPC effect ---------------------

print("Case isolation")

for (i in 1:nvalues) {
  sir$set_user(ipc = vIPC[i])
  
  tmp[[i]] <- simulate_model(sir) %>% mutate(param = "IPC", value = vIPC[i])
  
}

IPCsims <- bind_rows(tmp)

# reset original values
sir$set_user(ipc = 0.35)


# Pe effect ---------------------

print("Antiviral susceptibility effect")

for (i in 1:nvalues) {
  sir$set_user(Pe_on = vPe[i])
  
  tmp[[i]] <- simulate_model(sir) %>% mutate(param = "Pe", value = vPe[i])
  
}

PEsims <- bind_rows(tmp)

# reset original values
sir$set_user(Pe_on = 0.5)


# Pei effect ---------------------

print("Antiviral transmission effect")

for (i in 1:nvalues) {
  sir$set_user(Pei_on = vPei[i], Pei_Tr = vPei[i])
  
  tmp[[i]] <- simulate_model(sir) %>% mutate(param = "Pei", value = vPei[i])
  
}

PEIsims <- bind_rows(tmp)

# reset original values
sir$set_user(Pei_on = 0.1, Pei_Tr = 0.1)


# Short-stay residents ---------------------

print("Short-stay residents")

for (i in 1:nvalues) {
  sir$set_user(S_ini_SR = vSSR[i] - 1, 
               S_ini_LR = 100 - vSSR[i])
  
  tmp[[i]] <- simulate_model(sir) %>% mutate(param = "SSR", value = vSSR[i])
  
}

SSRsims <- bind_rows(tmp)

# reset original values
sir$set_user(S_ini_SR = 78  * N_factor - 1, S_ini_LR = 22  * N_factor)



## Combine output for analysis ---------------------

print("Combining output...")

# combine output
combined <- rbind(R0sims, CSSsims, CSRsims, IPCsims, PEsims, PEIsims, SSRsims) %>%
  group_by(param, value, sim) %>%
  filter(step == max(step))  %>%
  ungroup() %>%
  
  # get burden measures
  mutate(symp_R  = symptomatic_R,
         hosp_R  = H_SR + HP_SR + H_LR + HP_LR)

# get % outbreaks above a certain # cases
pr_case <- lapply(seq(10, 80, by = 5), function(x) {
  combined %>% group_by(param, value) %>%
    summarize(prop = round(length(symp_R[symp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows() 

# get % outbreaks above a certain # hospitalizations
pr_hosp <- lapply(seq(0, 20, by = 1), function(x) {
  combined %>% group_by(param, value) %>%
    summarize(prop = round(length(hosp_R[hosp_R > x])/n() * 100, 1)) %>%
    mutate(size = x)
}) %>%
  bind_rows()


pr_all <- pr_case %>% mutate(metric = "Illnesses") %>%
  rbind(., pr_hosp %>% mutate(metric = "Hospitalizations"))


## Plot results -------------------------------

cols <- RColorBrewer::brewer.pal(nvalues, "Spectral")

ptsize <- 1.5
linesize <- 1.2

xtitle <- "Total number > n"
ytitle <- "% Simulations"

p_R0 <- pr_all %>% filter(param == "R0") %>%
  ggplot(aes(x = size, y = prop, color = as.factor(value))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + 
  facet_wrap(~ metric, scales = "free_x") +
  get_theme() +
  scale_color_manual("R0", values = cols) +
  labs(x = xtitle, y = ytitle)

p_CSS <- pr_all %>% filter(param == "C_SS") %>%
  ggplot(aes(x = size, y = prop, color = as.factor(value))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + 
  facet_wrap(~ metric, scales = "free_x") +
  get_theme() +
  scale_color_manual("HCP-HCP\ncontact", values = cols) +
  labs(x = xtitle, y = ytitle)

p_CSR <- pr_all %>% filter(param == "C_SR") %>%
  ggplot(aes(x = size, y = prop, color = as.factor(value))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + 
  facet_wrap(~ metric, scales = "free_x") +
  get_theme() +
  scale_color_manual("HCP-resident\ncontact", values = cols) +
  labs(x = xtitle, y = ytitle)

p_IPC <- pr_all %>% filter(param == "IPC") %>%
  ggplot(aes(x = size, y = prop, color = as.factor(value))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + 
  facet_wrap(~ metric, scales = "free_x") +
  get_theme() +
  scale_color_manual("Case\nisolation\neffect", values = cols) +
  labs(x = xtitle, y = ytitle)

p_Pe <- pr_all %>% filter(param == "Pe") %>%
  ggplot(aes(x = size, y = prop, color = as.factor(value))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + 
  facet_wrap(~ metric, scales = "free_x") +
  get_theme() +
  scale_color_manual("Antiviral\nsusceptibility", values = cols) +
  labs(x = xtitle, y = ytitle)

p_Pei <- pr_all %>% filter(param == "Pei") %>%
  ggplot(aes(x = size, y = prop, color = as.factor(value))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + 
  facet_wrap(~ metric, scales = "free_x") +
  get_theme() +
  scale_color_manual("Antiviral\ntransmission", values = cols) +
  labs(x = xtitle, y = ytitle)

p_SSR <- pr_all %>% filter(param == "SSR") %>%
  ggplot(aes(x = size, y = prop, color = as.factor(value))) + 
  geom_point(size = ptsize) + geom_line(size = linesize) + 
  facet_wrap(~ metric, scales = "free_x") +
  get_theme() +
  scale_color_manual("Short-stay\nresidents", values = cols) +
  labs(x = xtitle, y = ytitle)

# Final figure S3
pS3 <- p_R0 / p_CSS / p_CSR / p_IPC / p_Pe / p_Pei / p_SSR

