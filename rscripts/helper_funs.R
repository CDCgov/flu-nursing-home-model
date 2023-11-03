# Helper functions for running and analysing results

#' Convert proportions / probabilities to rate parameters
#'
#' @param prop proportion to be transformed between 0 and 1 (must not be equal to 1)
#' @param unit time interval; default = 1 time step
#'
#' @return
#' @export
#'
#' @examples
prop_to_rate <- function(prop, unit = 1) {
  -log(1 - prop)/unit
}



#' Derive the days an outbreak started / finished (used for calculating e.g. # cases up until that point)
#'
#' @param out
#' @param start_lag integer value for period of days to consider when defining outbreak start
#' @param start_val integer value for number of new cases above which an outbreak is defined to have started
#' @param stop_lag integer value for period of days to consider when defining outbreak end
#' @param stop_val integer value for number of new cases below which an outbreak is defined to have ended
#'
#' @return
#' @export
#'
#' @examples
epi_duration <- function(out, start_lag = 1, start_val = 1, stop_lag = 7, stop_val = 1, sim = TRUE) {

  if(sim) {
    tmp <- out %>%  mutate(new_cases_stop  = symptomatic_R - lag(symptomatic_R, stop_lag),
                           new_cases_start = symptomatic_R - lag(symptomatic_R, start_lag),
                           epi_start = step[which(new_cases_start > start_val)[1]],
                           total_cases = I_SR + I_LR + Tr_SR + Tr_LR + IP_SR + IP_LR) %>% 
      group_by(sim) %>%
      filter(step > step[which.max(total_cases)]) %>%
      mutate(epi_end = step[which(new_cases_stop < stop_val)[1]] - stop_lag) %>%
      ungroup() %>% distinct(sim, epi_start, epi_end)
    
    return(tmp)
    
  } else {
    tmp <- out %>%  mutate(new_cases_stop  = symptomatic_R - lag(symptomatic_R, stop_lag),
                           new_cases_start = symptomatic_R - lag(symptomatic_R, start_lag),
                           epi_start = step[which(new_cases_start > start_val)[1]],
                           total_cases = I_SR + I_LR + Tr_SR + Tr_LR + IP_SR + IP_LR) %>% 
      filter(step > step[which.max(total_cases)]) %>%
      mutate(epi_end = step[which(new_cases_stop < stop_val)[1]] - stop_lag) %>%
      ungroup() %>% distinct(epi_end)
    
    return(as.numeric(tmp$epi_end[1]))
  }
}


#' Simulate a model (already loaded into odin)
#'
#' @param sir Model to run. Should have already been loaded using sir <- sir_generator$new()
#' @param n_sims Number of simulations to perform. Default = 200.
#'
#' @return
#' @export
#'
#' @examples
simulate_model <- function(sir, n_sims = 200) {
  set.seed(1)
  
  sims_df = as.data.frame(NULL)
  
  for (each_sim in seq(1,n_sims)){
    each_res <- sir$run(0:ndays)
    each_res_df = as.data.frame(each_res)
    each_res_df['sim'] = each_sim
    sims_df = rbind(sims_df, each_res_df)
  }
  
  sims_df
}




#' Get summary stats / plots from simulation outputs
#'
#' @param df_base data.frame of simulation outputs from odin without prophylaxis
#' @param df_new data.frame of simulation outputs from odin with prophylaxis
#' @param only_outbreaks default = TRUE. Should simulations that never turned into outbreaks be omitted from summary analyses?
#' @param entire_sim default = TRUE. Calculate burden based on the entire simulation timeframe? If FALSe, calculates burden up to the end of the main outbreak.
#' @param cols vector of hexcodes for chosen color scale. Default is 'Spectral' palette from the Rcolorbrewer package. Length must be >= 2.
#'
#' @return
#' @export
#'
#' @examples
compare_stats <- function(df_base, df_new, only_outbreaks = TRUE, entire_sim = TRUE,
                          cols = RColorBrewer::brewer.pal(2, "Spectral")) {

  # calculate the day each simulated outbreak started / ended 
  # (will used later as a marker for when to stop counting new cases etc)
  epi_end_base <- epi_duration(df_base, start_lag = 3, start_val = 2, stop_lag = 7)
  epi_end_new  <- epi_duration(df_new , start_lag = 3, start_val = 2, stop_lag = 7)
  
  # join the prophy / no prophy outputs and the estimates of outbreak start / end
  tmp <- df_base %>% mutate(scenario = "None") %>% left_join(epi_end_base) %>%
    rbind(., { df_new %>% mutate(scenario = "Until 7d\nno cases") %>% left_join(epi_end_new) }) %>% 
    
    # get burden estimates (symptomatic cases, hospitalizations, case-hospitalization ratio)
    # & get # doses used for prophy vs treatment
    mutate(symp_R  = symptomatic_R,
           symp_HC = symptomatic_HC,
           hosp_R  = H_SR + HP_SR + H_LR + HP_LR,
           chr_R = (hosp_R / symp_R) * 100,
           tx_doses_R     = Av_R_T/10,
           tx_doses_H     = Av_H_T/10,
           prophy_doses_R = ifelse(scenario == "None", 0, ceiling(Av_R_P/10)),
           prophy_doses_H = ifelse(scenario == "None", 0, ceiling(Av_H_P/10)),
           all_doses = tx_doses_R + tx_doses_H + prophy_doses_R + prophy_doses_H
    ) %>%
    group_by(sim, scenario) 
  
  # ignore outbreaks that never took off (zeros could bias estimates downwards)
  if (only_outbreaks) {
    tmp2 <- tmp %>% filter(max(symp_R) > 2) 
  } else {
    tmp2 <- tmp
  }
  

  # calculate peak burden
  peak <- tmp2 %>% group_by(sim, scenario) %>%
    summarize(peak_R  = max(I_SR + IP_SR + Tr_SR + I_LR + IP_LR + Tr_LR),
              peak_HC = max(I_HC + IP_HC + Tr_HC)) %>%
    ungroup()
    
  # count # doses only up until prophylaxis stops (prophy_track_stop = 0)
  # or until the outbreak stops (epi_end = 0) in the case of the no prophy scenario
  doses <- tmp2 %>% 
    group_by(scenario, sim) %>%
    filter(prophy_track_stop == 0 | (scenario == "None" & step == epi_end)) %>%
    filter(step == min(step)) %>% ungroup() %>%
    select(sim, scenario, tx_doses_R, tx_doses_H, prophy_doses_R, prophy_doses_H, all_doses)
  
  # calculate duration of epidemic and then filter df to the last day of the outbreak
  # so that cases, hosps etc are only counted up until that day
  
  if (entire_sim) {
    tmp2 <- tmp2 %>% 
      mutate(prophy_duration = ifelse(scenario == "None", 0, epi_end - (epi_start + 1) )) %>%
      filter(step == max(step)) %>%
      ungroup() %>% select(sim, scenario, total_R, total_HC, symp_R, symp_HC, hosp_R, chr_R, prophy_duration)
  } else {
    tmp2 <- tmp2 %>% 
      mutate(prophy_duration = ifelse(scenario == "None", 0, epi_end - (epi_start + 1) )) %>%
      filter(step == epi_end) %>%
      ungroup()  %>% select(sim, scenario, total_R, total_HC,  symp_R, symp_HC, hosp_R, chr_R, prophy_duration)
  }

  
  # calculte mean, 95 CIs for each scenario (prophy / no prophy) and print as a table
  table <- tmp2 %>% left_join(doses) %>% left_join(peak) %>%
    group_by(scenario) %>%
    summarize(across(c(total_R, total_HC, symp_R:all_doses, peak_R, peak_HC), function(x) {
      tmp <- round(quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
      paste0(tmp[2], " (", tmp[1], " - ", tmp[3], ")")
    })) %>% ungroup()
  
  
  # !! doesn't really make sense to calculate % reduction on a per-sim basis since they're not paired
  change <- tmp2 %>% left_join(doses) %>% left_join(peak) %>%
    group_by(sim) %>%
    mutate(across(c(total_R, total_HC, symp_R:all_doses, peak_R, peak_HC), function(x) {
                    x <- (x[1] - x) / x[1] * 100
                    })) %>% 
    ungroup() %>% group_by(scenario) %>%
    summarize(across(c(total_R, total_HC, symp_R:all_doses, peak_R, peak_HC), function(x) {
                      tmp <- round(quantile(x, c(0.025, 0.5, 0.975), na.rm = T))
                      paste0(tmp[2], " (", tmp[1], " - ", tmp[3], ")")
                    })) %>% ungroup()
  
  
  # plot full distributions
  plot <- tmp2 %>% 
    rename(Illnesses = symp_R, Hospitalizations = hosp_R) %>%
    gather(var, val, Illnesses, Hospitalizations) %>%
    ggplot(aes(x = scenario,y = val, fill = scenario)) + 
    geom_violin(alpha = 0.3) + geom_boxplot(lwd = 0.7, width = 0.2) +
    geom_jitter(width = 0.1, alpha = 0.1) +
    scale_fill_manual(guide = "none", values = cols) +
    facet_wrap(~ var, scales = "free_y") +
    get_theme() + labs(x = "", y = "")
  
  # return raw results
  raw <- tmp2 %>% left_join(doses) %>% left_join(peak)
  
  return(list(table = table, plot = plot, raw = raw, change = change))
}



#' Get % outbreaks above a certain threshold
#'
#' @param df_base data.frame of simulation outputs from odin without prophylaxis
#' @param df_new data.frame of simulation outputs from odin with prophylaxis
#'
#' @return
#' @export
#'
#' @examples
compare_probs <- function(df_base, df_new) {
  
  # calculate outbreak start / end days (to know when to stop counting new cases)
  epi_end_base <- epi_duration(df_base, start_lag = 3, start_val = 2, stop_lag = 7)
  epi_end_new  <- epi_duration(df_new , start_lag = 3, start_val = 2, stop_lag = 7)
  
  # join prophy / no prophy scenarios and the estimates of outbreak start / end
  tmp <- df_base %>% mutate(scenario = "No prophy") %>% left_join(epi_end_base) %>%
    rbind(., { df_new %>% mutate(scenario = "Prophy") %>% left_join(epi_end_new) }) %>% 
    
    # calculate symptomatic cases / hosps
    mutate(symp_R  = symptomatic_R,
           hosp_R  = H_SR + HP_SR + H_LR + HP_LR) %>%
    group_by(sim, scenario) %>% 
    
    # filter to the day outbreak ends (so cases etc are only counted up until that point)
    filter(step == epi_end) %>%
    ungroup() 
  
  # loop over possible outbreak sizes (here 10-80) & calc % over that size
  cases <- lapply(seq(10, 80, by = 5), function(x) {
    tmp %>% group_by(scenario) %>%
      summarize(prop = round(length(symp_R[symp_R > x])/n() * 100, 1)) %>%
      mutate(size = x)
  }) %>%
    bind_rows()
  
  # loop over possible outbreak sizes (here 0-20 for hosps) & calc % over that size
  hosps <- lapply(seq(0, 20, by = 1), function(x) {
    tmp %>% group_by(scenario) %>%
      summarize(prop = round(length(hosp_R[hosp_R > x])/n() * 100, 1)) %>%
      mutate(size = x)
  }) %>%
    bind_rows()
  
  return(list(cases = cases, hosps = hosps))
}


#' Set default plotting theme
#'
#' @param txt integer value for text size. Default is 12.
#'
#' @return
#' @export
#'
#' @examples
get_theme <- function(txt = 12, ...){
  theme_light() + theme(axis.text    = element_text(size = txt),
                        axis.title   = element_text(size = txt),
                        strip.text   = element_text(size = txt),
                        legend.text  = element_text(size = txt - 1),
                        legend.title = element_text(size = txt - 1),
                        title        = element_text(size = txt), 
                        ...)
}
