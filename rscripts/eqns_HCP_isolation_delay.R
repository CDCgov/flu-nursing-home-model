## Stochastic compartmental model framework using the odin package

## Model code and equations -----

sir_generator <- odin::odin({
  ## Notes on equations for transitions between compartments -----
  # can only do one draw from a compartment at a timestep- otherwise will overdraw compartments and get negative numbers
  # SR = short stay resident, LR = long stay resident, HC = healthcare personel
  # no P = not on prophylaxis, P = prophylaxis
  # terms are labeled by were they are coming from and going to. i.e. SE = from susceptible to exposed
  # leaving facility is not explicitly written, but is included in the _tot draws 
  
  ## Equations short-stay -----
  update(S_SR) <- S_SR - n_SESP_tot_SR + n_SPS_SR + n_import_S_SR  # S short stay. Leaving: S-> E, S -> SP, exit facility. Entering: SP-> S, enter facility
  update(E_SR) <- E_SR + n_SE_SR - n_EIAEP_tot_SR + n_EPE_SR  # E short stay. Leaving: E-> I, E-> A, E-> EP, exit facility. Entering: S-> E, EP-> E
  update(A_SR) <- A_SR + n_EA_SR - n_ARAP_tot_SR + n_APA_SR+ n_import_A_SR # A short stay. Leaving: A-> R, A-> AP, exit facility. Entering: E-> A, AP-> A, enter facility
  update(I_SR) <- I_SR + n_EI_SR  - n_IHTrIP_tot_SR + n_import_I_SR# I short stay. Leaving: I-> H, I-> Tr, I-> IP, exit facility. Entering: E-> I, enter facility. All I are either treated or hospitalized
  update(Tr_SR) <- Tr_SR + n_ITr_SR - n_TrHR_tot_SR # Tr short stay. Leaving: Tr-> H, Tr-> R. Entering: I-> Tr
  update(H_SR) <- H_SR + n_IH_SR + n_TrH_SR # H short stay. Leaving: none. Entering: I-> H, Tr-> H. H accumulates. 
  update(R_SR) <- R_SR + n_TrR_SR + n_AR_SR + n_IR_SR - n_RRP_tot_SR + n_RPR_SR # R short stay. Leaving: R->RP, exit facility. Entering: Tr-> R, A-> R, RP-> R
  
  update(SP_SR) <- SP_SR - n_SPEPS_tot_SR + n_SSP_SR + n_import_SP_SR # SP short stay. Leaving: SP-> EP, SP -> S, exit facility. Entering: S-> SP, enter facility
  update(EP_SR) <- EP_SR + n_SPEP_SR - n_EPIPAPE_tot_SR + n_EEP_SR# EP short stay. Leaving: EP-> IP, EP-> AP, EP-> E, exit facility. Entering: SP-> EP, E-> EP
  update(AP_SR) <- AP_SR + n_EPAP_SR - n_APRPA_tot_SR + n_AAP_SR + n_import_AP_SR # AP short stay. Leaving: AP-> RP, AP-> A, exit facility. Entering: EP-> AP, A-> AP, enter facility
  update(IP_SR) <- IP_SR + n_EPIP_SR - n_IPHPRPI_tot_SR +n_IIP_SR + n_TrIP_SR + n_import_IP_SR# IP short stay. Leaving: IP-> HP, I-> RP, IP-> I, exit facility. Entering: EP-> IP, I-> IP, enter facility. IP can recover or hospitalize
  update(HP_SR) <- HP_SR + n_IPHP_SR # HP short stay. Leaving: none. Entering: IP-> HP. HP accumulates.
  update(RP_SR) <- RP_SR + n_APRP_SR + n_IPRP_SR + n_RRP_SR - n_RPR_tot_SR# RP short stay. Leaving: RP->R, exit facility. Entering: AP-> RP, IP-> RP, R-> RP
  
  ## Equations long-stay -----
  
  update(S_LR) <- S_LR - n_SESP_tot_LR  + n_SPS_LR+ n_import_S_LR # S long stay. Leaving: S-> E, S -> SP, exit facility. Entering: SP-> S, enter facility
  update(E_LR) <- E_LR + n_SE_LR - n_EIAEP_tot_LR + n_EPE_LR # E longstay. Leaving: E-> I, E-> A, E-> EP, exit facility. Entering: S-> E, EP-> E
  update(A_LR) <- A_LR + n_EA_LR - n_ARAP_tot_LR + n_APA_LR+ n_import_A_LR # A long stay. Leaving: A-> R, A-> AP, exit facility. Entering: E-> A, AP-> A, enter facility
  update(I_LR) <- I_LR + n_EI_LR  - n_IHTrIP_tot_LR + n_import_I_LR# I long stay. Leaving: I-> H, I-> Tr, I-> IP, exit facility. Entering: E-> I, enter facility. All I are either treated or hospitalized
  update(Tr_LR) <- Tr_LR + n_ITr_LR - n_TrHR_tot_LR # Tr long stay. Leaving: Tr-> H, Tr-> R. Entering: I-> Tr
  update(H_LR) <- H_LR + n_IH_LR + n_TrH_LR # H long stay. Leaving: none. Entering: I-> H, Tr-> H. H accumulates.
  update(R_LR) <- R_LR + n_TrR_LR + n_AR_LR + n_IR_LR - n_RRP_tot_LR + n_RPR_LR# R long stay. Leaving: R->RP, exit facility. Entering: Tr-> R, A-> R, RP-> R
  
  update(SP_LR) <- SP_LR - n_SPEPS_tot_LR + n_SSP_LR + n_import_SP_LR# SP long stay. Leaving: SP-> EP, SP -> S, exit facility. Entering: S-> SP, enter facility
  update(EP_LR) <- EP_LR + n_SPEP_LR - n_EPIPAPE_tot_LR + n_EEP_LR# EP long stay. Leaving: EP-> IP, EP-> AP, EP-> E, exit facility. Entering: SP-> EP, E-> EP
  update(AP_LR) <- AP_LR + n_EPAP_LR - n_APRPA_tot_LR + n_AAP_LR + n_import_AP_LR# AP long stay. Leaving: AP-> RP, AP-> A, exit facility. Entering: EP-> AP, A-> AP, enter facility
  update(IP_LR) <- IP_LR + n_EPIP_LR - n_IPHPRPI_tot_LR +n_IIP_LR + n_TrIP_LR + n_import_IP_LR# IP long stay. Leaving: IP-> HP, I-> RP, IP-> I, exit facility. Entering: EP-> IP, I-> IP, enter facility. IP can recover or hospitalize
  update(HP_LR) <- HP_LR + n_IPHP_LR # HP long stay. Leaving: none. Entering: IP-> HP. HP accumulates.
  update(RP_LR) <- RP_LR + n_APRP_LR + n_IPRP_LR + n_RRP_LR - n_RPR_tot_LR# RP long stay. Leaving: RP->R, exit facility. Entering: AP-> RP, IP-> RP, R-> RP
  
  ## Equations healthcare personnel -----
  
  update(S_HC) <- S_HC - n_SESP_tot_HC  + n_SPS_HC+ n_import_S_HC# S HCW. Leaving: S-> E, S -> SP, exit facility. Entering: SP-> S, enter facility
  update(E_HC) <- E_HC + n_SE_HC - n_EIAEP_tot_HC + n_EPE_HC# E HCW. Leaving: E-> I, E-> A, E-> EP, exit facility. Entering: S-> E, EP-> E
  update(A_HC) <- A_HC + n_EA_HC - n_ARAP_tot_HC + n_APA_HC+ n_import_A_HC# A HCW. Leaving: A-> R, A-> AP, exit facility. Entering: E-> A, AP-> A, enter facility
  update(I_HC) <- I_HC + n_EI_HC  - n_IHTrIP_tot_HC + n_import_I_HC# I HCW. Leaving: I-> H, I-> Tr, I-> IP, exit facility. Entering: E-> I, enter facility. All I are either treated or hospitalized
  update(Tr_HC) <- Tr_HC + n_ITr_HC - n_TrHR_tot_HC# Tr HCW. Leaving: Tr-> H, Tr-> R. Entering: I-> Tr
  update(H_HC) <- H_HC + n_IH_HC + n_TrH_HC # H HCW. Leaving: none. Entering: I-> H, Tr-> H. H accumulates.
  update(R_HC) <- R_HC + n_TrR_HC + n_AR_HC + n_IR_HC - n_RRP_tot_HC + n_RPR_HC# R HCW. Leaving: R->RP, exit facility. Entering: Tr-> R, A-> R, RP-> R
  
  update(SP_HC) <- SP_HC - n_SPEPS_tot_HC + n_SSP_HC+ n_import_SP_HC# SP HCW. Leaving: SP-> EP, SP -> S, exit facility. Entering: S-> SP, enter facility
  update(EP_HC) <- EP_HC + n_SPEP_HC - n_EPIPAPE_tot_HC + n_EEP_HC# EP HCW. Leaving: EP-> IP, EP-> AP, EP-> E, exit facility. Entering: SP-> EP, E-> EP
  update(AP_HC) <- AP_HC + n_EPAP_HC - n_APRPA_tot_HC + n_AAP_HC + n_import_AP_HC# AP HCW. Leaving: AP-> RP, AP-> A, exit facility. Entering: EP-> AP, A-> AP, enter facility
  update(IP_HC) <- IP_HC + n_EPIP_HC - n_IPHPRPI_tot_HC +n_IIP_HC + n_TrIP_HC + n_import_IP_HC# IP HCW. Leaving: IP-> HP, I-> RP, IP-> I, exit facility. Entering: EP-> IP, I-> IP, enter facility. IP can recover or hospitalize
  update(HP_HC) <- HP_HC + n_IPHP_HC  # HP HCW. Leaving: none. Entering: IP-> HP. HP accumulates.
  update(RP_HC) <- RP_HC + n_APRP_HC + n_IPRP_HC + n_RRP_HC - n_RPR_tot_HC# RP HCW. Leaving: RP->R, exit facility. Entering: AP-> RP, IP-> RP, R-> RP
  
  
  
  ## Entrance to the facility from the community -----
  
  # balance with total individuals leaving so that population size remains constant
  total_import_SR <- n_export_S_SR  + n_export_E_SR  + n_export_A_SR  + n_export_I_SR  + n_export_R_SR + n_export_Tr_SR +
    n_export_SP_SR + n_export_EP_SR + n_export_AP_SR + n_export_IP_SR + n_export_RP_SR
  
  total_import_LR <- n_export_S_LR  + n_export_E_LR  + n_export_A_LR  + n_export_I_LR  + n_export_R_LR + n_export_Tr_LR +
    n_export_SP_LR + n_export_EP_LR + n_export_AP_LR + n_export_IP_LR + n_export_RP_LR
  
  total_import_HC <- n_export_S_HC  + n_export_E_HC  + n_export_A_HC  + n_export_I_HC  + n_export_R_HC + n_export_Tr_HC +
    n_export_SP_HC + n_export_EP_HC + n_export_AP_HC + n_export_IP_HC + n_export_RP_HC
  
  # partition incomers among prophylaxis and non-prophylaxis
  # 1. find proportion who start on prophy vs those who don't
  update(prophy_incoming_R)  <- if ( new_cases_lag >= case_threshold && NP < 1) uptake_R_prop  else prophy_incoming_R
  update(prophy_incoming_HC) <- if ( new_cases_lag >= case_threshold && NP < 1) uptake_HC_prop else prophy_incoming_HC
  
  # 2. distribute among infection classes
  #    here we assume all incomers are susceptible, so distribute between S / SP
  n_import_S_SR <- round(total_import_SR * (1 - prophy_incoming_R))
  n_import_A_SR <- 0
  n_import_I_SR <- 0
  n_import_SP_SR <- total_import_SR - n_import_S_SR
  n_import_AP_SR <- 0
  n_import_IP_SR <- 0
  
  n_import_S_LR <- round(total_import_LR * (1 - prophy_incoming_R))
  n_import_A_LR <- 0
  n_import_I_LR <- 0
  n_import_SP_LR <- total_import_LR - n_import_S_LR
  n_import_AP_LR <- 0
  n_import_IP_LR <- 0
  
  n_import_S_HC <- round(total_import_HC * (1 - prophy_incoming_HC))
  n_import_A_HC <- 0
  n_import_I_HC <- 0
  n_import_SP_HC <- total_import_HC - n_import_S_HC
  n_import_AP_HC <- 0
  n_import_IP_HC <- 0
  
  
  ## Start of population draws ---------------------
  # all draws are based on: https://mrc-ide.github.io/odin/articles/discrete.html
  # When individuals from one population need to be distributed across multiple compartments,
  # need the multinomial distribution
  
  ## Calculate leaving probabilities for no prophylaxis -------
  ## Draw from S ----
  
  ipc0 <- 1 - ipc
  
  infectious_resid <- A_SR + A_LR + (1-Pei)*(AP_SR + AP_LR) + ipc0 * (I_SR + I_LR + (1-Pei_Tr)*(Tr_SR + Tr_LR + IP_SR + IP_LR) )
  
  # Note:'+ n_EI_HC' etc includes staff who have just become symptomatic (i.e. allows a day delay in staff isolation)
  infectious_staff <- A_HC + (1-Pei) * AP_HC + n_EI_HC + (1-Pei_Tr)* n_EPIP_HC
  
  # force of infection
  lambda_S_SR <- beta * (C_RR * infectious_resid + C_RS * infectious_staff) 
  
  calc_lambdaS_SR<-(lambda_S_SR + prophy_start_R + mu_SR) #those leaving S_SR due to transmission, start prophy, and leaving facility
  p_SESP_tot_SR <- 1 - exp(- calc_lambdaS_SR) # total probability of leaving S_SR
  n_SESP_tot_SR <- rbinom(S_SR, p_SESP_tot_SR)# draw total prob leaving S_SR from whole S_SR pop
  n_SESP_SR[] <- rmultinom(n_SESP_tot_SR, h_s_SR)# then, split the total leaving S_SR into the different comparements
  h_s_SR[1] <- lambda_S_SR / calc_lambdaS_SR # splits based on relative probabilities of going to the different compartments
  h_s_SR[2] <- prophy_start_R / calc_lambdaS_SR # go to a
  h_s_SR[3] <- mu_SR / calc_lambdaS_SR
  dim( h_s_SR) <- 3 # just a vector that hold the relative probs
  dim(n_SESP_SR) <- 3# a vector that holds the assigned number
  n_SE_SR <- n_SESP_SR[1] 
  n_SSP_SR <- n_SESP_SR[2]
  n_export_S_SR <- n_SESP_SR[3]
  
  lambda_S_LR <- lambda_S_SR
  calc_lambdaS_LR<-(lambda_S_LR + prophy_start_R +mu_LR)
  
  p_SESP_tot_LR <- 1 - exp(- calc_lambdaS_LR) # S to E
  n_SESP_tot_LR <- rbinom(S_LR, p_SESP_tot_LR)
  n_SESP_LR[] <- rmultinom(n_SESP_tot_LR, h_s_LR)
  h_s_LR[1] <- lambda_S_LR / calc_lambdaS_LR # go to rp
  h_s_LR[2] <-   prophy_start_R / calc_lambdaS_LR # go to a
  h_s_LR[3] <-   mu_LR / calc_lambdaS_LR # go to a
  dim( h_s_LR) <- 3
  dim(n_SESP_LR) <- 3
  n_SE_LR <- n_SESP_LR[1]
  n_SSP_LR <- n_SESP_LR[2]
  n_export_S_LR <- n_SESP_LR[3]
  
  lambda_S_HC <- beta * (C_SR * infectious_resid + C_SS * infectious_staff) 
  
  calc_lambdaS_HC<-(lambda_S_HC + prophy_start_HC + mu_HC)
  p_SESP_tot_HC <- 1 - exp(- calc_lambdaS_HC) # S to E
  n_SESP_tot_HC <- rbinom(S_HC, p_SESP_tot_HC)
  n_SESP_HC[] <- rmultinom(n_SESP_tot_HC, h_s_HC)
  h_s_HC[1] <- lambda_S_HC / calc_lambdaS_HC # go to rp
  h_s_HC[2] <-   prophy_start_HC / calc_lambdaS_HC # go to a
  h_s_HC[3] <-   mu_HC / calc_lambdaS_HC # go to a
  dim( h_s_HC) <- 3
  dim(n_SESP_HC) <- 3
  n_SE_HC <- n_SESP_HC[1]
  n_SSP_HC <- n_SESP_HC[2]
  n_export_S_HC <- n_SESP_HC[3]
  
  
  ## Draw from E ----
  
  p_EIAEP_tot_SR <- 1 - exp(-(sigma + prophy_start_R + mu_SR)) # I to R
  n_EIAEP_tot_SR <- rbinom(E_SR, p_EIAEP_tot_SR )
  n_EIAEP_SR[] <- rmultinom(n_EIAEP_tot_SR, h_e_SR)
  h_e_SR[1] <- (1-fv_res)*sigma  / (sigma + prophy_start_R + mu_SR)
  h_e_SR[2] <-  fv_res *sigma    / (sigma + prophy_start_R + mu_SR)
  h_e_SR[3] <-  (prophy_start_R) / (sigma + prophy_start_R + mu_SR) 
  h_e_SR[4] <-  (mu_SR)          / (sigma + prophy_start_R + mu_SR) 
  dim( h_e_SR) <- 4
  dim(n_EIAEP_SR) <- 4
  n_EA_SR <- n_EIAEP_SR[1]
  n_EI_SR <- n_EIAEP_SR[2]
  n_EEP_SR <- n_EIAEP_SR[3]
  n_export_E_SR <- n_EIAEP_SR[4]
  
  p_EIAEP_tot_LR <- 1 - exp(-(sigma + prophy_start_R + mu_LR)) # I to R
  n_EIAEP_tot_LR <- rbinom(E_LR, p_EIAEP_tot_LR )
  n_EIAEP_LR[] <- rmultinom(n_EIAEP_tot_LR, h_e_LR)
  h_e_LR[1] <- (1-fv_res)*sigma / (sigma + prophy_start_R + mu_LR)# go to rp
  h_e_LR[2] <-  fv_res*sigma / (sigma + prophy_start_R + mu_LR)# go to a
  h_e_LR[3] <-  (prophy_start_R) /(sigma + prophy_start_R + mu_LR) # go to a
  h_e_LR[4] <-  (mu_LR) /(sigma + prophy_start_R + mu_LR) # go to a
  dim( h_e_LR) <- 4
  dim(n_EIAEP_LR) <- 4
  n_EA_LR <- n_EIAEP_LR[1]
  n_EI_LR <- n_EIAEP_LR[2]
  n_EEP_LR <- n_EIAEP_LR[3]
  n_export_E_LR<- n_EIAEP_LR[4]
  
  p_EIAEP_tot_HC <- 1 - exp(-(sigma + prophy_start_HC + mu_HC)) # I to R
  n_EIAEP_tot_HC <- rbinom(E_HC, p_EIAEP_tot_HC )
  n_EIAEP_HC[] <- rmultinom(n_EIAEP_tot_HC, h_e_HC)
  h_e_HC[1] <- (1-fv_hc)*sigma / (sigma + prophy_start_HC+ mu_HC)# go to rp
  h_e_HC[2] <-  fv_hc*sigma / (sigma + prophy_start_HC+ mu_HC)# go to a
  h_e_HC[3] <-  (prophy_start_HC) /(sigma + prophy_start_HC+ mu_HC) # go to a
  h_e_HC[4] <-  (mu_HC) /(sigma + prophy_start_HC + mu_HC) # go to a
  dim( h_e_HC) <- 4
  dim(n_EIAEP_HC) <- 4
  n_EA_HC <- n_EIAEP_HC[1]
  n_EI_HC <- n_EIAEP_HC[2]
  n_EEP_HC <- n_EIAEP_HC[3]
  n_export_E_HC <- n_EIAEP_HC[4]
  
  
  ## Draw from A ---------
  
  p_ARAP_tot_SR <- 1 - exp(-(gamma + prophy_start_R + mu_SR)) # I to R
  n_ARAP_tot_SR <- rbinom(A_SR, p_ARAP_tot_SR)
  n_ARAP_SR[] <- rmultinom(n_ARAP_tot_SR, h_a_SR)
  h_a_SR[1] <- (gamma) / (gamma + prophy_start_R + mu_SR)# + go to rp
  h_a_SR[2] <-  (prophy_start_R) / (gamma + prophy_start_R+ mu_SR) # go to a
  h_a_SR[3] <-  (mu_SR) / (gamma + prophy_start_R + mu_SR) # go to a
  dim( h_a_SR) <- 3
  dim(n_ARAP_SR) <- 3
  n_AR_SR <- n_ARAP_SR[1]
  n_AAP_SR <- n_ARAP_SR[2]
  n_export_A_SR <- n_ARAP_SR[3]
  
  p_ARAP_tot_LR <- 1 - exp(-(gamma + prophy_start_R + mu_LR)) 
  n_ARAP_tot_LR <- rbinom(A_LR, p_ARAP_tot_LR)
  n_ARAP_LR[] <- rmultinom(n_ARAP_tot_LR, h_a_LR )
  h_a_LR[1] <- (gamma) / (gamma + prophy_start_R + mu_LR)
  h_a_LR[2] <-  (prophy_start_R) / (gamma + prophy_start_R + mu_LR) 
  h_a_LR[3] <-  (mu_LR) / (gamma + prophy_start_R + mu_LR) 
  dim( h_a_LR) <- 3
  dim(n_ARAP_LR) <- 3
  n_AR_LR <- n_ARAP_LR[1]
  n_AAP_LR <- n_ARAP_LR[2]
  n_export_A_LR<- n_ARAP_LR[3]
  
  p_ARAP_tot_HC <- 1 - exp(-(gamma + prophy_start_HC + mu_HC)) # I to R
  n_ARAP_tot_HC <- rbinom(A_HC, p_ARAP_tot_HC)
  n_ARAP_HC[] <- rmultinom(n_ARAP_tot_HC, h_a_HC)
  h_a_HC[1] <- (gamma) / (gamma + prophy_start_HC + mu_HC)# go to rp
  h_a_HC[2] <-  (prophy_start_HC) / (gamma + prophy_start_HC + mu_HC) # go to a
  h_a_HC[3] <-  (mu_HC) / (gamma + prophy_start_HC+ mu_HC) # go to a
  dim( h_a_HC) <- 3
  dim(n_ARAP_HC) <- 3
  n_AR_HC <- n_ARAP_HC[1]
  n_AAP_HC <- n_ARAP_HC[2]
  n_export_A_HC <- n_ARAP_HC[3]
  
  
  ## Draw from I ----
  
  p_IHTrIP_tot_SR <- 1 - exp(-(delta + gamma + prophy_start_R + mu_SR)) # I to R
  n_IHTrIP_tot_SR <- rbinom(I_SR, p_IHTrIP_tot_SR )
  n_IHTrIP_SR[] <- rmultinom(n_IHTrIP_tot_SR, h_i_SR)
  h_i_SR[1] <- chr_res * gamma / (delta + gamma + prophy_start_R+ mu_SR)# go to rp
  h_i_SR[2] <-  delta  / (delta + gamma + prophy_start_R+ mu_SR)# go to a
  h_i_SR[3] <-  (prophy_start_R) /(delta + gamma + prophy_start_R+ mu_SR) # go to a
  h_i_SR[4] <-  (mu_SR) /(delta + gamma + prophy_start_R+ mu_SR) # go to a
  h_i_SR[5] <-  (1 - chr_res) * gamma /(delta + gamma + prophy_start_R + mu_SR) # go to a
  dim( h_i_SR) <- 5
  dim(n_IHTrIP_SR) <- 5
  n_IH_SR <- n_IHTrIP_SR[1]
  n_ITr_SR <- n_IHTrIP_SR[2]
  n_IIP_SR <- n_IHTrIP_SR[3]
  n_export_I_SR <- n_IHTrIP_SR[4]
  n_IR_SR <- n_IHTrIP_SR[5]
  
  p_IHTrIP_tot_LR <- 1 - exp(-(delta + gamma + prophy_start_R + mu_LR)) # I to R
  n_IHTrIP_tot_LR <- rbinom(I_LR, p_IHTrIP_tot_LR )
  n_IHTrIP_LR[] <- rmultinom(n_IHTrIP_tot_LR, h_i_LR)
  h_i_LR[1] <- chr_res* gamma  / (delta + gamma + prophy_start_R + mu_LR)# go to rp
  h_i_LR[2] <-  delta  / (delta + gamma + prophy_start_R + mu_LR)# go to a
  h_i_LR[3] <-  (prophy_start_R) /(delta + gamma + prophy_start_R + mu_LR) # go to a
  h_i_LR[4] <-  (mu_LR) / (delta + gamma + prophy_start_R + mu_LR) # go to a
  h_i_LR[5] <-  (1 - chr_res) * gamma / (delta + gamma + prophy_start_R + mu_LR) # go to a
  dim( h_i_LR) <- 5
  dim(n_IHTrIP_LR) <- 5
  n_IH_LR <- n_IHTrIP_LR[1]
  n_ITr_LR <- n_IHTrIP_LR[2]
  n_IIP_LR <- n_IHTrIP_LR[3]
  n_export_I_LR <- n_IHTrIP_LR[4]
  n_IR_LR <- n_IHTrIP_LR[5]
  
  p_IHTrIP_tot_HC <- 1 - exp(-(delta + gamma + prophy_start_HC + mu_HC)) # I to R
  n_IHTrIP_tot_HC <- rbinom(I_HC, p_IHTrIP_tot_HC )
  n_IHTrIP_HC[] <- rmultinom(n_IHTrIP_tot_HC, h_i_HC)
  h_i_HC[1] <- chr_hc * gamma  / (delta + gamma + prophy_start_HC+ mu_HC)# go to rp
  h_i_HC[2] <-  delta  / (delta + gamma + prophy_start_HC+ mu_HC)# go to a
  h_i_HC[3] <-  (prophy_start_HC) / (delta + gamma + prophy_start_HC+ mu_HC) # go to a
  h_i_HC[4] <-  (mu_HC) / (delta + gamma + prophy_start_HC+ mu_HC) # go to a
  h_i_HC[5] <-  (1 - chr_hc) * gamma /(delta + gamma + prophy_start_HC + mu_HC) # go to a
  dim( h_i_HC) <- 5
  dim(n_IHTrIP_HC) <- 5
  n_IH_HC <- n_IHTrIP_HC[1]
  n_ITr_HC <- n_IHTrIP_HC[2]
  n_IIP_HC <- n_IHTrIP_HC[3]
  n_export_I_HC <- n_IHTrIP_HC[4]
  n_IR_HC <- n_IHTrIP_HC[5]
  
  
  ## Draw from Tr ----
  
  # gamma vs epsilon: 
  # assume faster recovery (epsilon) only applies to people who recover, 
  # not to ppl who will require hospitalization
  
  to_HR_res <- chr_res * (1 - Pes) * gamma + (1 - chr_res * (1 - Pes)) * epsilon
  to_HR_hc  <- chr_hc  * (1 - Pes) * gamma + (1 - chr_hc  * (1 - Pes)) * epsilon
  
  p_TrHR_tot_SR <- 1-exp(-(to_HR_res + mu_SR + prophy_start_R)) 
  n_TrHR_tot_SR <- rbinom(Tr_SR, p_TrHR_tot_SR)
  n_TrHR_SR[] <- rmultinom(n_TrHR_tot_SR, g_SR)
  g_SR[1] <-      chr_res * (1 - Pes) * gamma  / (to_HR_res + mu_SR + prophy_start_R)
  g_SR[2] <- (1 - chr_res * (1 - Pes))* epsilon/ (to_HR_res + mu_SR + prophy_start_R)
  g_SR[3] <-                           (mu_SR) / (to_HR_res + mu_SR + prophy_start_R)
  g_SR[4] <-                  (prophy_start_R) / (to_HR_res + mu_SR + prophy_start_R)
  dim(g_SR) <- 4
  dim( n_TrHR_SR) <- 4
  n_TrH_SR <- n_TrHR_SR[1]
  n_TrR_SR <- n_TrHR_SR[2]
  n_export_Tr_SR <- n_TrHR_SR[3]
  n_TrIP_SR <- n_TrHR_SR[4]
  
  p_TrHR_tot_LR <- 1-exp(-(to_HR_res + mu_LR + prophy_start_R))
  n_TrHR_tot_LR <- rbinom(Tr_LR, p_TrHR_tot_LR)
  n_TrHR_LR[] <- rmultinom(n_TrHR_tot_LR, g_LR)
  g_LR[1] <-      chr_res * (1 - Pes) * gamma  / (to_HR_res + mu_LR + prophy_start_R)
  g_LR[2] <- (1 - chr_res * (1 - Pes))* epsilon/ (to_HR_res + mu_LR + prophy_start_R)
  g_LR[3] <-                           (mu_LR) / (to_HR_res + mu_LR + prophy_start_R)
  g_LR[4] <-                  (prophy_start_R) / (to_HR_res + mu_LR + prophy_start_R)
  dim(g_LR) <- 4
  dim( n_TrHR_LR) <- 4
  n_TrH_LR <- n_TrHR_LR[1]
  n_TrR_LR <- n_TrHR_LR[2]
  n_export_Tr_LR <- n_TrHR_LR[3]
  n_TrIP_LR <- n_TrHR_LR[4]
  
  p_TrHR_tot_HC <- 1-exp(-(to_HR_hc + mu_HC + prophy_start_HC))
  n_TrHR_tot_HC <- rbinom(Tr_HC, p_TrHR_tot_HC)
  n_TrHR_HC[] <- rmultinom(n_TrHR_tot_HC, g_HC)
  g_HC[1] <-      chr_hc * (1 - Pes) * gamma  / (to_HR_hc + mu_HC + prophy_start_HC)
  g_HC[2] <- (1 - chr_hc * (1 - Pes))* epsilon/ (to_HR_hc + mu_HC + prophy_start_HC)
  g_HC[3] <-                          (mu_HC) / (to_HR_hc + mu_HC + prophy_start_HC)
  g_HC[4] <-                (prophy_start_HC) / (to_HR_hc + mu_HC + prophy_start_HC)
  dim(g_HC) <- 4
  dim( n_TrHR_HC) <- 4
  n_TrH_HC <- n_TrHR_HC[1]
  n_TrR_HC <- n_TrHR_HC[2]
  n_export_Tr_HC <- n_TrHR_HC[3]
  n_TrIP_HC <- n_TrHR_HC[4]
  
  ## Draw from R ----
  
  p_RRP_tot_SR <- 1-exp(-(prophy_start_R+ mu_SR))
  n_RRP_tot_SR<- rbinom(R_SR, p_RRP_tot_SR)
  n_RRP_SR_vec[] <- rmultinom(n_RRP_tot_SR, h_r_SR)
  h_r_SR[1] <- prophy_start_R/ (prophy_start_R+ mu_SR)# + go to rp
  h_r_SR[2] <-   mu_SR/ (prophy_start_R+ mu_SR) # go to a
  dim( h_r_SR) <- 2
  dim(n_RRP_SR_vec) <- 2
  n_RRP_SR <- n_RRP_SR_vec[1]
  n_export_R_SR <- n_RRP_SR_vec[2]
  
  
  p_RRP_tot_LR <- 1-exp(-(prophy_start_R+ mu_LR))
  n_RRP_tot_LR<- rbinom(R_LR, p_RRP_tot_LR)
  n_RRP_LR_vec[] <- rmultinom(n_RRP_tot_LR, h_r_LR)
  h_r_LR[1] <- prophy_start_R/ (prophy_start_R+ mu_LR)# + go to rp
  h_r_LR[2] <-   mu_LR/ (prophy_start_R+ mu_LR) # go to a
  dim( h_r_LR) <- 2
  dim(n_RRP_LR_vec) <- 2
  n_RRP_LR <- n_RRP_LR_vec[1]
  n_export_R_LR <- n_RRP_LR_vec[2]
  
  p_RRP_tot_HC <- 1-exp(-(prophy_start_HC+ mu_HC))
  n_RRP_tot_HC<- rbinom(R_HC, p_RRP_tot_HC)
  n_RRP_HC_vec[] <- rmultinom(n_RRP_tot_HC, h_r_HC)
  h_r_HC[1] <- prophy_start_HC/ (prophy_start_HC+ mu_HC)# + go to rp
  h_r_HC[2] <-   mu_HC/ (prophy_start_HC+ mu_HC) # go to a
  dim( h_r_HC) <- 2
  dim(n_RRP_HC_vec) <- 2
  n_RRP_HC <- n_RRP_HC_vec[1]
  n_export_R_HC <- n_RRP_HC_vec[2]
  
  
  ## Calculate leaving probabilities for prophylaxis -------
  ## Draw from SP -----
  
  lambda_SP_SR <- lambda_S_SR
  calc_lambdaP_SR <- ( (1 - Pe) * lambda_SP_SR + c_stop_R + mu_SR)
  
  p_SPEPS_tot_SR <- 1 - exp(-calc_lambdaP_SR) # S to E
  n_SPEPS_tot_SR <- rbinom(SP_SR, p_SPEPS_tot_SR)
  n_SPEPS_SR[] <- rmultinom(n_SPEPS_tot_SR, q_s_SR)
  q_s_SR[1] <- (1 - Pe) * lambda_SP_SR / calc_lambdaP_SR# go to rp
  q_s_SR[2] <-  c_stop_R / calc_lambdaP_SR# go to a
  q_s_SR[3] <-  mu_SR / calc_lambdaP_SR# go to a
  dim( q_s_SR) <- 3
  dim(n_SPEPS_SR) <- 3
  n_SPEP_SR <- n_SPEPS_SR[1]
  n_SPS_SR <- n_SPEPS_SR[2]
  n_export_SP_SR <- n_SPEPS_SR[3]              
  
  lambda_SP_LR <- lambda_S_LR
  calc_lambdaP_LR<-( (1 - Pe) * lambda_SP_LR + c_stop_R + mu_LR)
  p_SPEPS_tot_LR <- 1 - exp(-calc_lambdaP_LR) # S to E
  n_SPEPS_tot_LR <- rbinom(SP_LR, p_SPEPS_tot_LR)
  n_SPEPS_LR[] <- rmultinom(n_SPEPS_tot_LR, q_s_LR)
  q_s_LR[1] <- (1 - Pe) * lambda_SP_LR / calc_lambdaP_LR# go to rp
  q_s_LR[2] <-  c_stop_R / calc_lambdaP_LR# go to a
  q_s_LR[3] <-  mu_LR / calc_lambdaP_LR# go to a
  dim( q_s_LR) <- 3
  dim(n_SPEPS_LR) <- 3
  n_SPEP_LR <- n_SPEPS_LR[1]
  n_SPS_LR <- n_SPEPS_LR[2]
  n_export_SP_LR <- n_SPEPS_LR[3]
  
  lambda_SP_HC <- lambda_S_HC
  calc_lambdaP_HC<-( (1 - Pe) * lambda_SP_HC + c_stop_HC + mu_HC)
  p_SPEPS_tot_HC <- 1 - exp(-calc_lambdaP_HC) # S to E
  n_SPEPS_tot_HC <- rbinom(SP_HC, p_SPEPS_tot_HC)
  n_SPEPS_HC[] <- rmultinom(n_SPEPS_tot_HC, q_s_HC)
  q_s_HC[1] <- (1 - Pe) * lambda_SP_HC / calc_lambdaP_HC# go to rp
  q_s_HC[2] <-  c_stop_HC / calc_lambdaP_HC# go to a
  q_s_HC[3] <-  mu_HC / calc_lambdaP_HC# go to a
  dim( q_s_HC) <- 3
  dim(n_SPEPS_HC) <- 3
  n_SPEP_HC <- n_SPEPS_HC[1]
  n_SPS_HC <- n_SPEPS_HC[2]
  n_export_SP_HC <- n_SPEPS_HC[3]
  
  ## Draw from EP -----
  
  p_EPIPAPE_tot_SR <- 1 - exp(-(sigma + c_stop_R + mu_SR)) # I to R
  n_EPIPAPE_tot_SR <- rbinom(EP_SR, p_EPIPAPE_tot_SR )
  n_EPIPAPE_SR[] <- rmultinom(n_EPIPAPE_tot_SR, q_e_SR)
  q_e_SR[1] <- (1-fv_res)*sigma / (sigma + c_stop_R + mu_SR)# go to rp
  q_e_SR[2] <-     fv_res*sigma / (sigma + c_stop_R + mu_SR)# go to a
  q_e_SR[3] <-     (c_stop_R) /(sigma + c_stop_R + mu_SR) # go to a
  q_e_SR[4] <-      (mu_SR) /(sigma + c_stop_R + mu_SR) # go to a
  dim( q_e_SR) <- 4
  dim(n_EPIPAPE_SR) <- 4
  n_EPAP_SR <- n_EPIPAPE_SR[1]
  n_EPIP_SR <- n_EPIPAPE_SR[2]
  n_EPE_SR <- n_EPIPAPE_SR[3]
  n_export_EP_SR <- n_EPIPAPE_SR[4]
  
  p_EPIPAPE_tot_LR <- 1 - exp(-(sigma + c_stop_R + mu_LR)) # I to R
  n_EPIPAPE_tot_LR <- rbinom(EP_LR, p_EPIPAPE_tot_LR )
  n_EPIPAPE_LR[] <- rmultinom(n_EPIPAPE_tot_LR, q_e_LR)
  q_e_LR[1] <- (1-fv_res)*sigma / (sigma + c_stop_R + mu_LR)# go to rp
  q_e_LR[2] <-    fv_res*sigma / (sigma + c_stop_R + mu_LR)# go to a
  q_e_LR[3] <-     (c_stop_R) /(sigma + c_stop_R + mu_LR) # go to a
  q_e_LR[4] <-      (mu_LR) /(sigma + c_stop_R + mu_LR) # go to a
  dim( q_e_LR) <- 4
  dim(n_EPIPAPE_LR) <- 4
  n_EPAP_LR <- n_EPIPAPE_LR[1]
  n_EPIP_LR <- n_EPIPAPE_LR[2]
  n_EPE_LR <- n_EPIPAPE_LR[3]
  n_export_EP_LR <- n_EPIPAPE_LR[4]
  
  p_EPIPAPE_tot_HC <- 1 - exp(-(sigma + c_stop_HC + mu_HC)) # I to R
  n_EPIPAPE_tot_HC <- rbinom(EP_HC, p_EPIPAPE_tot_HC )
  n_EPIPAPE_HC[] <- rmultinom(n_EPIPAPE_tot_HC, q_e_HC)
  q_e_HC[1] <- (1-fv_hc)*sigma / (sigma + c_stop_HC + mu_HC)# go to rp
  q_e_HC[2] <-     fv_hc*sigma / (sigma + c_stop_HC + mu_HC)# go to a
  q_e_HC[3] <-     (c_stop_HC) /(sigma + c_stop_HC + mu_HC) # go to a
  q_e_HC[4] <-      (mu_HC) /(sigma + c_stop_HC + mu_HC) # go to a
  dim( q_e_HC) <- 4
  dim(n_EPIPAPE_HC) <- 4
  n_EPAP_HC <- n_EPIPAPE_HC[1]
  n_EPIP_HC <- n_EPIPAPE_HC[2]
  n_EPE_HC <- n_EPIPAPE_HC[3]
  n_export_EP_HC <- n_EPIPAPE_HC[4]
  
  
  ## Draw from AP -----
  
  p_APRPA_tot_SR <- 1 - exp(-(epsilon0 + c_stop_R + mu_SR)) # I to R
  n_APRPA_tot_SR <- rbinom(AP_SR, p_APRPA_tot_SR)
  n_APRPA_SR[] <- rmultinom(n_APRPA_tot_SR, q_a_SR)
  q_a_SR[1] <-    epsilon0 / (epsilon0 + c_stop_R + mu_SR)# go to rp
  q_a_SR[2] <-  (c_stop_R) / (epsilon0 + c_stop_R + mu_SR) # go to a
  q_a_SR[3] <-     (mu_SR) / (epsilon0 + c_stop_R + mu_SR) # go to a
  dim( q_a_SR) <- 3
  dim(n_APRPA_SR) <- 3
  n_APRP_SR <- n_APRPA_SR[1]
  n_APA_SR <- n_APRPA_SR[2]
  n_export_AP_SR <- n_APRPA_SR[3]
  
  p_APRPA_tot_LR <- 1 - exp(-(epsilon0 + c_stop_R + mu_LR)) # I to R
  n_APRPA_tot_LR <- rbinom(AP_LR, p_APRPA_tot_LR)
  n_APRPA_LR[] <- rmultinom(n_APRPA_tot_LR, q_a_LR)
  q_a_LR[1] <-    epsilon0 / (epsilon0 + c_stop_R + mu_LR)# go to rp
  q_a_LR[2] <-  (c_stop_R) / (epsilon0 + c_stop_R + mu_LR) # go to a
  q_a_LR[3] <-     (mu_LR) / (epsilon0 + c_stop_R + mu_LR) # go to a
  dim( q_a_LR) <- 3
  dim(n_APRPA_LR) <- 3
  n_APRP_LR <- n_APRPA_LR[1]
  n_APA_LR <- n_APRPA_LR[2]
  n_export_AP_LR <- n_APRPA_LR[3]
  
  p_APRPA_tot_HC <- 1 - exp(-(epsilon0 + c_stop_HC + mu_HC)) # I to R
  n_APRPA_tot_HC <- rbinom(AP_HC, p_APRPA_tot_HC)
  n_APRPA_HC[] <- rmultinom(n_APRPA_tot_HC, q_a_HC)
  q_a_HC[1] <-     epsilon0 / (epsilon0 + c_stop_HC + mu_HC)# go to rp
  q_a_HC[2] <-  (c_stop_HC) / (epsilon0 + c_stop_HC + mu_HC) # go to a
  q_a_HC[3] <-      (mu_HC) / (epsilon0 + c_stop_HC + mu_HC) # go to a
  dim( q_a_HC) <- 3
  dim(n_APRPA_HC) <- 3
  n_APRP_HC <- n_APRPA_HC[1]
  n_APA_HC <- n_APRPA_HC[2]
  n_export_AP_HC <- n_APRPA_HC[3]
  
  ## Draw from IP -----
  
  # gamma vs epsilon: 
  # assume faster recovery (epsilon) only applies to people who recover, 
  # not to ppl who will require hospitalization
  # to_HR_res and to_HR_hc defined above (~ line 300)
  
  # assume no loss to compliance while people are symptomatic (i.e. c_stop = 0)
  # (sending back to I class would have prob of starting Av again by --> Tr)
  
  p_IPHPRPI_tot_SR <- 1 - exp(-(to_HR_res + mu_SR)) # I to R
  n_IPHPRPI_tot_SR <- rbinom(IP_SR, p_IPHPRPI_tot_SR )
  n_IPHPRPI_SR[] <- rmultinom(n_IPHPRPI_tot_SR, q_i_SR)
  q_i_SR[1] <-       chr_res * (1 - Pes)  * gamma  / (to_HR_res + mu_SR)# go to rp
  q_i_SR[2] <-  (1 - chr_res * (1 - Pes)) * epsilon/ (to_HR_res + mu_SR)# go to a
  q_i_SR[3] <-         (mu_SR)   /(to_HR_res + mu_SR) # go to a
  dim( q_i_SR) <- 3
  dim(n_IPHPRPI_SR) <- 3
  n_IPHP_SR <- n_IPHPRPI_SR[1]
  n_IPRP_SR <- n_IPHPRPI_SR[2]
  n_export_IP_SR <- n_IPHPRPI_SR[3]
  
  p_IPHPRPI_tot_LR <- 1 - exp(-(to_HR_res + mu_LR)) # I to R
  n_IPHPRPI_tot_LR <- rbinom(IP_LR, p_IPHPRPI_tot_LR )
  n_IPHPRPI_LR[] <- rmultinom(n_IPHPRPI_tot_LR, q_i_LR)
  q_i_LR[1] <-       chr_res * (1 - Pes) * gamma   / (to_HR_res + mu_LR)# go to rp
  q_i_LR[2] <-  (1 - chr_res * (1 - Pes)) * epsilon/ (to_HR_res + mu_LR)# go to a
  q_i_LR[3] <-  (mu_LR) /(to_HR_res + mu_LR) # go to a
  dim( q_i_LR) <- 3
  dim(n_IPHPRPI_LR) <- 3
  n_IPHP_LR <- n_IPHPRPI_LR[1]
  n_IPRP_LR <- n_IPHPRPI_LR[2]
  n_export_IP_LR <- n_IPHPRPI_LR[3]
  
  p_IPHPRPI_tot_HC <- 1 - exp(-(to_HR_hc + mu_HC)) # I to R
  n_IPHPRPI_tot_HC <- rbinom(IP_HC, p_IPHPRPI_tot_HC )
  n_IPHPRPI_HC[] <- rmultinom(n_IPHPRPI_tot_HC, q_i_HC)
  q_i_HC[1] <- chr_hc * (1 - Pes) * gamma / (to_HR_hc + mu_HC)# go to rp
  q_i_HC[2] <-  (1 - chr_hc * (1 - Pes)) * epsilon / (to_HR_hc + mu_HC)# go to a
  q_i_HC[3] <-  (mu_HC) /(to_HR_hc + mu_HC) # go to a
  dim( q_i_HC) <- 3
  dim(n_IPHPRPI_HC) <- 3
  n_IPHP_HC <- n_IPHPRPI_HC[1]
  n_IPRP_HC <- n_IPHPRPI_HC[2]
  n_export_IP_HC <- n_IPHPRPI_HC[3]
  
  ## Draw from RP -----
  
  p_RPR_tot_SR <- 1-exp(-(c_stop_R + mu_SR))
  n_RPR_tot_SR <- rbinom(RP_SR, p_RPR_tot_SR)
  n_RPR_SR_vec[] <- rmultinom(n_RPR_tot_SR, q_r_SR)
  q_r_SR[1] <- c_stop_R/ (c_stop_R+ mu_SR)# + go to rp
  q_r_SR[2] <-   mu_SR/ (c_stop_R+ mu_SR) # go to a
  dim( q_r_SR) <- 2
  dim(n_RPR_SR_vec) <- 2
  n_RPR_SR <- n_RPR_SR_vec[1]
  n_export_RP_SR <- n_RPR_SR_vec[2]
  
  p_RPR_tot_LR <- 1-exp(-(c_stop_R + mu_LR))
  n_RPR_tot_LR <- rbinom(RP_LR, p_RPR_tot_LR)
  n_RPR_LR_vec[] <- rmultinom(n_RPR_tot_LR, q_r_LR)
  q_r_LR[1] <- c_stop_R/ (c_stop_R+ mu_LR)# + go to rp
  q_r_LR[2] <-   mu_LR/ (c_stop_R+ mu_LR) # go to a
  dim( q_r_LR) <- 2
  dim(n_RPR_LR_vec) <- 2
  n_RPR_LR <- n_RPR_LR_vec[1]
  n_export_RP_LR <- n_RPR_LR_vec[2]
  
  p_RPR_tot_HC <- 1-exp(-(c_stop_HC + mu_HC))
  n_RPR_tot_HC <- rbinom(RP_HC, p_RPR_tot_HC)
  n_RPR_HC_vec[] <- rmultinom(n_RPR_tot_HC, q_r_HC)
  q_r_HC[1] <- c_stop_HC/ (c_stop_HC+ mu_HC)# + go to rp
  q_r_HC[2] <-   mu_HC/ (c_stop_HC+ mu_HC) # go to a
  dim( q_r_HC) <- 2
  dim(n_RPR_HC_vec) <- 2
  n_RPR_HC <- n_RPR_HC_vec[1]
  n_export_RP_HC <- n_RPR_HC_vec[2]
  
  
  ## Total population size ------
  
  N <- S_SR  + E_SR  + A_SR  + I_SR  + Tr_SR + H_SR  + R_SR + 
       SP_SR + EP_SR + AP_SR + IP_SR + HP_SR + RP_SR + 
       S_LR  + E_LR  + A_LR  + I_LR  + Tr_LR + H_LR + R_LR + 
       SP_LR + EP_LR + AP_LR + IP_LR + HP_LR + RP_LR +
       S_HC  + E_HC  + A_HC  + I_HC  + Tr_HC + H_HC + R_HC + 
       SP_HC + EP_HC + AP_HC + IP_HC + HP_HC + RP_HC 
  
  # Total on prophylaxis
  NP <- SP_SR + EP_SR + AP_SR + IP_SR + HP_SR + RP_SR + 
        SP_LR + EP_LR + AP_LR + IP_LR + HP_LR + RP_LR +
        SP_HC + EP_HC + AP_HC + IP_HC + HP_HC + RP_HC 
  
  ## Prophylaxis initiation (varies over time) ----- 
  # When # cases passes a certain threshold, the prophy param switches from 0 (none go to prophy before this) to the % uptaking prophy. 
  # Then one timestep later, switches back to 0. 
  
  # New identified (symptomatic) cases are those going from E -> I (e.g. n_EI_SR, n_EPIP_SR etc)
  # Calculate new cases at all possible lags
  new_cases_7d <- delay(n_EI_SR + n_EI_LR + n_EPIP_SR + n_EPIP_LR, 7)
  new_cases_6d <- delay(n_EI_SR + n_EI_LR + n_EPIP_SR + n_EPIP_LR, 6)
  new_cases_5d <- delay(n_EI_SR + n_EI_LR + n_EPIP_SR + n_EPIP_LR, 5)
  new_cases_4d <- delay(n_EI_SR + n_EI_LR + n_EPIP_SR + n_EPIP_LR, 4)
  new_cases_3d <- delay(n_EI_SR + n_EI_LR + n_EPIP_SR + n_EPIP_LR, 3)
  new_cases_2d <- delay(n_EI_SR + n_EI_LR + n_EPIP_SR + n_EPIP_LR, 2)
  new_cases_1d <- delay(n_EI_SR + n_EI_LR + n_EPIP_SR + n_EPIP_LR, 1)
  new_cases_0d <-       n_EI_SR + n_EI_LR + n_EPIP_SR + n_EPIP_LR
  
  # Sum only lags relevant for the current scenario 
  update(new_cases_lag) <- (start_lag >= 8 && step >= 8) * new_cases_7d + (start_lag >= 7 && step >= 7) * new_cases_6d + 
                           (start_lag >= 6 && step >= 6) * new_cases_5d + (start_lag >= 5 && step >= 5) * new_cases_4d + 
                           (start_lag >= 4 && step >= 4) * new_cases_3d + (start_lag >= 3 && step >= 3) * new_cases_2d + 
                           (start_lag >= 2 && step >= 2) * new_cases_1d + (start_lag >= 1 && step >= 1) * new_cases_0d
  
  prophy_start_R  <- if ( new_cases_lag >= case_threshold && NP < 1) uptake_R  else 0
  prophy_start_HC <- if ( new_cases_lag >= case_threshold && NP < 1) uptake_HC else 0
  
  # Track when prophy is turned on
  update(prophy_track_start) <- if (new_cases_lag >= case_threshold && NP < 1) 1 else 0
  update(prophy_start_time)  <- if (new_cases_lag >= case_threshold && NP < 1) step else prophy_start_time
  
  
  ## Prophylaxis stopping ----- 
  
  # Calculate new cases in last Xd ago for prophy stopping condition (include only one timepts defined by stop_lag)
  new_cases_lag_stop <- (stop_lag >= 8 && step >= 8) * new_cases_7d + (stop_lag >= 7 && step >= 7) * new_cases_6d + 
                        (stop_lag >= 6 && step >= 6) * new_cases_5d + (stop_lag >= 5 && step >= 5) * new_cases_4d + 
                        (stop_lag >= 4 && step >= 4) * new_cases_3d + (stop_lag >= 3 && step >= 3) * new_cases_2d + 
                        (stop_lag >= 2 && step >= 2) * new_cases_1d + (stop_lag >= 1 && step >= 1) * new_cases_0d
  
  
  ## Define case- and time-based definitions for when prophy should be turned off:
  # First check prophy is being administered
  update(prophy_track_stop) <- if (( NP > 0) &&
                                   # Case-based definition for stopping
                                   (step >= (prophy_start_time + min_prophy_duration) && new_cases_lag_stop < 1) || 
                                   # Time-based definition for stopping
                                   (step >= prophy_start_time + max_prophy_duration)) 
    0 else prophy_track_stop
  
  # Use the tracker to turn off the prophylaxis efficacy parameters
  Pei <- Pei_on * prophy_track_stop
  Pe  <- Pe_on  * prophy_track_stop
  
  # if prophy has stopped, there is no faster recovery for AP class                       
  epsilon0 <- if (prophy_track_stop == 1) epsilon else gamma
  
  
  ## Tracking prophylaxis use -----
  
  # Doses for treatment = ppl entering Tr & IP classes x 10d course
  update(Av_R_T) <- Av_R_T + (n_ITr_SR + n_EPIP_SR + n_IIP_SR + n_ITr_LR + n_EPIP_LR + n_IIP_LR) * 10
  update(Av_H_T) <- Av_H_T + (n_ITr_HC + n_EPIP_HC + n_IIP_HC) * 10
  
  # Doses for prophy = 1x daily for everyone currently on PEP except IP (- get Tx dose) 
  # Ignore hospitalized cases for the time being
  update(Av_R_P) <- Av_R_P + SP_SR + EP_SR + AP_SR + RP_SR + SP_LR + EP_LR + AP_LR + RP_LR
  update(Av_H_P) <- Av_H_P + SP_HC + EP_HC + AP_HC + RP_HC

  initial(Av_R_T) <- 0
  initial(Av_H_T) <- 0
  initial(Av_R_P) <- 0
  initial(Av_H_P) <- 0
  
  ## Tracking total burden -------------------------------
  
  initial(symptomatic_R)  <- 0
  initial(symptomatic_HC) <- 0
  
  update(symptomatic_R)  <- symptomatic_R  + new_cases_0d
  update(symptomatic_HC) <- symptomatic_HC + n_EI_HC + n_EPIP_HC + n_import_I_HC + n_import_IP_HC
  
  initial(total_R)  <- 0 
  initial(total_HC) <- 0
  
  # Note: could include imported cases (n_import_I / IP / A / AP) in below eqns 
  # but since these are = 0, they are omitted for simplicity
  update(total_R)  <- total_R  + n_SE_SR + n_SPEP_SR + n_SE_LR + n_SPEP_LR
  update(total_HC) <- total_HC + n_SE_HC + n_SPEP_HC
  
  initial(foi_SR) <- 0
  initial(foi_LR) <- 0
  initial(foi_HC) <- 0
  
  update(foi_SR) <- lambda_S_SR
  update(foi_LR) <- lambda_S_LR
  update(foi_HC) <- lambda_S_HC
  
  
  ## Initial states ----
  
  initial(S_SR) <- S_ini_SR
  initial(E_SR) <- E_ini_SR
  initial(A_SR) <- A_ini_SR
  initial(I_SR) <- I_ini_SR
  initial(Tr_SR) <- 0
  initial(R_SR) <- 0
  initial(H_SR) <- 0
  initial(SP_SR) <- SP_ini_SR
  initial(EP_SR) <- EP_ini_SR
  initial(AP_SR) <- AP_ini_SR
  initial(IP_SR) <- IP_ini_SR
  initial(RP_SR) <- 0
  initial(HP_SR) <- 0
  
  initial(S_LR) <- S_ini_LR
  initial(E_LR) <- E_ini_LR
  initial(A_LR) <- A_ini_LR
  initial(I_LR) <- I_ini_LR
  initial(Tr_LR) <- 0
  initial(R_LR) <- 0
  initial(H_LR) <- 0
  initial(SP_LR) <- SP_ini_LR
  initial(EP_LR) <- EP_ini_LR
  initial(AP_LR) <- AP_ini_LR
  initial(IP_LR) <- IP_ini_LR
  initial(RP_LR) <- 0
  initial(HP_LR) <- 0
  
  initial(S_HC) <- S_ini_HC
  initial(E_HC) <- E_ini_HC
  initial(A_HC) <- A_ini_HC
  initial(I_HC) <- I_ini_HC
  initial(Tr_HC) <- 0
  initial(R_HC) <- 0
  initial(H_HC) <- 0
  initial(SP_HC) <- SP_ini_HC
  initial(EP_HC) <- EP_ini_HC
  initial(AP_HC) <- AP_ini_HC
  initial(IP_HC) <- IP_ini_HC
  initial(RP_HC) <- 0
  initial(HP_HC) <- 0

  
  ## User defined Initial conditions ------------------------
  
  S_ini_SR <- user(78)
  E_ini_SR <- user(0)
  A_ini_SR <- user(1)
  I_ini_SR <- user(0)
  
  SP_ini_SR <- user(0)
  EP_ini_SR <- user(0)
  AP_ini_SR<- user(0)
  IP_ini_SR <- user(0)
  
  S_ini_LR <- user(22)
  E_ini_LR <- user(0)
  A_ini_LR <- user(0)
  I_ini_LR <- user(0)
  
  SP_ini_LR <- user(0)
  EP_ini_LR <- user(0)
  AP_ini_LR<- user(0)
  IP_ini_LR <- user(0)
  
  S_ini_HC <- user(100)
  E_ini_HC <- user(0)
  A_ini_HC <- user(0)
  I_ini_HC <- user(0)
  
  SP_ini_HC <- user(0)
  EP_ini_HC <- user(0)
  AP_ini_HC<- user(0)
  IP_ini_HC <- user(0)
  
  
  ## User defined parameters (not related to prophylaxis)------------------------
  
  beta  <- user(0.05)
  gamma <- user(1/4)
  sigma <- user(1/1.5)
  delta <- user(1/1)
  
  ipc <- user(0.3)
  
  fv_res  <- user()
  fv_hc   <- user()
  chr_res <- user(0.45)
  chr_hc  <- user(0.006)
  
  theta   <- user(3/2)
  epsilon <- gamma * theta
  
  C_RR <- user(2)
  C_RS <- user(6)
  C_SR <- user(6)
  C_SS <- user(2)
  
  mu_SR <- user(1/25)
  mu_LR <- user(1/406)
  mu_HC <- user(1/1000000000)
  
  
  ## Prophylaxis initiation/cessation -------------------
  
  # variables to track when prophylaxis starts / stops
  initial(prophy_track_start) <- 0
  initial(prophy_start_time)  <- 0
  initial(prophy_track_stop)  <- 1
  
  # variables to partition incomers to prophylaxis or not
  initial(prophy_incoming_R)  <- 0
  initial(prophy_incoming_HC) <- 0
  
  # tracking # cases for when to initiate prophylaxis
  initial(new_cases_lag) <- 0
  
  case_threshold <- user(2)
  
  # max and min duration of prophy
  max_prophy_duration <- user(365)
  min_prophy_duration <- user(14)
  
  # days of over which cases are summed for prophy initiation (e.g. 2 cases in 3 days)
  start_lag <- user(3)
  
  # days that prophy continues after outbreak over (for conditions without max duration)
  stop_lag <- user(7)
  
  # prophylaxis uptake parameters (proportions)
  uptake_R_prop  <- user(0.9)
  uptake_HC_prop <- user(0.65)
  
  # prophylaxis uptake parameters (rates)
  uptake_R  <- user()
  uptake_HC <- user()
  
  # prophylaxis loss to compliance
  c_stop_R  <- user(1e-10)
  c_stop_HC <- user(1/5)
  
  # antiviral efficacy parameters
  # against onward transmission
  Pei_on <- user(0.1)
  Pei_Tr <- user(0.1)
  
  # against hospitalization (assume everyone w/ symptoms is on Tx dose, so continues no matter what)
  Pes <- user(0.6)
  
  # against infection (only relevant for prophylaxis)
  Pe_on <- user(0.5)
  
  
  
  ### End -----
}, verbose = FALSE)
