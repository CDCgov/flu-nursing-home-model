# Baseline parameter setting for stochastic model

f_res  <- 0.55
f_hc   <- 0.4
VC_res <- 0.75
VC_hc  <- 0.65
VE_res <- 0.25
VE_hc  <- 0.4

C_RR <- 2  
C_RS <- C_SR <- 6
C_SS <- 2

# if investigating larger population sizes, set N_factor > 1
N_factor <- 1

contact_matrix <- rbind(c(C_RR, C_RS), c(C_SR, C_SS))
pop_matrix     <- rbind(c(100,  100),  c(100,  100)) * N_factor

R0    <- 7
gamma <- 1/4


# R0 = beta * (max eigenvalue of NGM) <=> beta = R0 / (max eigenvalue)
beta <- R0 / max( eigen(contact_matrix * pop_matrix / gamma)$values )

sir <- sir_generator$new( S_ini_SR = 78  * N_factor - 1, # one SSR is initially infected
                          S_ini_LR = 22  * N_factor, 
                          S_ini_HC = 100 * N_factor,
                          beta = beta, gamma = gamma, ipc = 0.35, chr_res = 0.5,
                          C_RR = C_RR, C_RS = C_RS, C_SR = C_SR, C_SS = C_SS, 
                          fv_res = f_res * (1 - VC_res * VE_res), 
                          fv_hc  = f_hc  * (1 - VC_hc * VE_hc),
                          uptake_R = prop_to_rate(0.9), uptake_HC = prop_to_rate(0.65),
                          uptake_R_prop = 0.9, uptake_HC_prop = 0.65
)

ndays  <- 90
n_sims <- 500
n_sample_sims <- 40