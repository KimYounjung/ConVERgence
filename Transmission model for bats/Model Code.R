###############################################################################

# SEIR_generator_1_C1 is for the first year for Colony 1
SEIR_generator_1_C1 = odin::odin( {
  
  beta = beta0 + beta1*((sin(2*3.141593*(tunit-beta2))+1)/2)

  # juvenile 
  deriv(S_j) = - death_rate*S_j - beta*(I_j+I_a)*S_j + sigma_j*R_j + sigma_m*M_j + birth_rate_S
  deriv(E_j) = - death_rate*E_j + beta*(I_j+I_a)*S_j - tau_j*E_j     
  deriv(I_j) = - death_rate*I_j + tau_j*E_j - gamma_j*I_j                     
  deriv(R_j) = - death_rate*R_j + gamma_j*I_j - sigma_j*R_j  
  deriv(M_j) = - death_rate*M_j - sigma_m*M_j + birth_rate_M
  
  # adults 
  deriv(S_a) = - death_rate*S_a - beta*(I_j+I_a)*S_a + sigma_a*R_a    
  deriv(E_a) = - death_rate*E_a + beta*(I_j+I_a)*S_a - tau_a*E_a                 
  deriv(I_a) = - death_rate*I_a + tau_a*E_a - gamma_a*I_a                  
  deriv(R_a) = - death_rate*R_a + gamma_a*I_a - sigma_a*R_a    
  
  # cumulative number of infections
  deriv(cum_inf_j) = beta*(I_j+I_a)*S_j
  deriv(cum_inf_a) = beta*(I_j+I_a)*S_a
  deriv(tunit) = 1/365  
  
  ## Initial states:
  initial(S_j) = 0 
  initial(E_j) = 0 
  initial(I_j) = 0 
  initial(R_j) = 0 
  initial(M_j) = 0 
  
  # 659.0835 is the mean population size of Colony 1
  
  initial(S_a) = 659.0835 - 1 -1 
  #initial(S_a) = 1042.75 - 1 -1 - 104.275 
  initial(E_a) = 1 
  initial(I_a) = 1 
  initial(R_a) = 0
  #initial(R_a) = 104.275
  
  initial(cum_inf_j) = 0
  initial(cum_inf_a) = 0
  initial(tunit) = 1/365  
  
  times[] = user()
  dim(times) = user()
  
  # Parameters to be varied over time   
  beta0 = interpolate(times, beta0_0, "constant")
  beta0_0[] = user()
  dim(beta0_0) = length(times)
  
  beta1 = interpolate(times, beta1_0, "constant")
  beta1_0[] = user()
  dim(beta1_0) = length(times)
  
  beta2 = interpolate(times, beta2_0, "constant")
  beta2_0[] = user()
  dim(beta2_0) = length(times)
  
  tau_j = interpolate(times, tau_j_0, "constant")
  tau_j_0[] = user()
  dim(tau_j_0) = length(times)
  
  tau_a = interpolate(times, tau_a_0, "constant")
  tau_a_0[] = user()
  dim(tau_a_0) = length(times)
  
  gamma_j = interpolate(times, gamma_j_0, "constant")
  gamma_j_0[] = user()
  dim(gamma_j_0) = length(times)
  
  gamma_a = interpolate(times, gamma_a_0, "constant")
  gamma_a_0[] = user()
  dim(gamma_a_0) = length(times)
  
  sigma_j = interpolate(times, sigma_j_0, "constant")
  sigma_j_0[] = user()
  dim(sigma_j_0) = length(times)

  sigma_a = interpolate(times, sigma_a_0, "constant")
  sigma_a_0[] = user()
  dim(sigma_a_0) = length(times)
  
  sigma_m = interpolate(times, sigma_m_0, "constant")
  sigma_m_0[] = user()
  dim(sigma_m_0) = length(times)
  
  death_rate = interpolate(times, death_rate_0, "constant")
  death_rate_0[] = user()
  dim(death_rate_0) = length(times)
  
  birth_rate_S = interpolate(times, birth_rate_S_0, "constant")
  birth_rate_S_0[] = user()
  dim(birth_rate_S_0) = length(times)
  
  birth_rate_M = interpolate(times, birth_rate_M_0, "constant")
  birth_rate_M_0[] = user()
  dim(birth_rate_M_0) = length(times)
  
})
# SEIR_generator_1_C2 is for the first year for Colony 2
SEIR_generator_1_C2 = odin::odin( {
  
  beta = beta0 + beta1*((sin(2*3.141593*(tunit-beta2))+1)/2)

  deriv(S_j) = - death_rate*S_j - beta*(I_j+I_a)*S_j + sigma_j*R_j + sigma_m*M_j + birth_rate_S
  deriv(E_j) = - death_rate*E_j + beta*(I_j+I_a)*S_j - tau_j*E_j     
  deriv(I_j) = - death_rate*I_j + tau_j*E_j - gamma_j*I_j                     
  deriv(R_j) = - death_rate*R_j + gamma_j*I_j - sigma_j*R_j  
  deriv(M_j) = - death_rate*M_j - sigma_m*M_j + birth_rate_M
  
  # adults 
  deriv(S_a) = - death_rate*S_a - beta*(I_j+I_a)*S_a + sigma_a*R_a    
  deriv(E_a) = - death_rate*E_a + beta*(I_j+I_a)*S_a - tau_a*E_a                 
  deriv(I_a) = - death_rate*I_a + tau_a*E_a - gamma_a*I_a                  
  deriv(R_a) = - death_rate*R_a + gamma_a*I_a - sigma_a*R_a    
  
  # cumulative number of infections
  deriv(cum_inf_j) = beta*(I_j+I_a)*S_j
  deriv(cum_inf_a) = beta*(I_j+I_a)*S_a
  deriv(tunit) = 1/365  
  
  ## Initial states:
  initial(S_j) = 0 
  initial(E_j) = 0 
  initial(I_j) = 0 
  initial(R_j) = 0 
  initial(M_j) = 0 
  
  # 668.75 is the mean population size of Colony 2
  initial(S_a) = 668.75 - 1 -1
  #initial(S_a) = 957.5 - 1 -1 - 95.75
  initial(E_a) = 1 
  initial(I_a) = 1 
  initial(R_a) = 0
  #initial(R_a) = 95.75
  
  initial(cum_inf_j) = 0
  initial(cum_inf_a) = 0
  initial(tunit) = 1/365
  
  times[] = user()
  dim(times) = user()
  
  # Parameters to be varied over time   
  beta0 = interpolate(times, beta0_0, "constant")
  beta0_0[] = user()
  dim(beta0_0) = length(times)
  
  beta1 = interpolate(times, beta1_0, "constant")
  beta1_0[] = user()
  dim(beta1_0) = length(times)
  
  beta2 = interpolate(times, beta2_0, "constant")
  beta2_0[] = user()
  dim(beta2_0) = length(times)
  
  tau_j = interpolate(times, tau_j_0, "constant")
  tau_j_0[] = user()
  dim(tau_j_0) = length(times)
  
  tau_a = interpolate(times, tau_a_0, "constant")
  tau_a_0[] = user()
  dim(tau_a_0) = length(times)
  
  gamma_j = interpolate(times, gamma_j_0, "constant")
  gamma_j_0[] = user()
  dim(gamma_j_0) = length(times)
  
  gamma_a = interpolate(times, gamma_a_0, "constant")
  gamma_a_0[] = user()
  dim(gamma_a_0) = length(times)
  
  sigma_j = interpolate(times, sigma_j_0, "constant")
  sigma_j_0[] = user()
  dim(sigma_j_0) = length(times)
  
  sigma_a = interpolate(times, sigma_a_0, "constant")
  sigma_a_0[] = user()
  dim(sigma_a_0) = length(times)
  
  sigma_m = interpolate(times, sigma_m_0, "constant")
  sigma_m_0[] = user()
  dim(sigma_m_0) = length(times)
  
  death_rate = interpolate(times, death_rate_0, "constant")
  death_rate_0[] = user()
  dim(death_rate_0) = length(times)
  
  birth_rate_S = interpolate(times, birth_rate_S_0, "constant")
  birth_rate_S_0[] = user()
  dim(birth_rate_S_0) = length(times)
  
  birth_rate_M = interpolate(times, birth_rate_M_0, "constant")
  birth_rate_M_0[] = user()
  dim(birth_rate_M_0) = length(times)
  
})
# SEIR_generator_2 is for the remaining years 
SEIR_generator_2_1 = odin::odin( {
  
  beta = beta0 + beta1*((sin(2*3.141593*(tunit-beta2))+1)/2)

  deriv(S_j) = - death_rate*S_j - beta*(I_j+I_a)*S_j + sigma_j*R_j + sigma_m*M_j + birth_rate_S
  deriv(E_j) = - death_rate*E_j + beta*(I_j+I_a)*S_j - tau_j*E_j     
  deriv(I_j) = - death_rate*I_j + tau_j*E_j - gamma_j*I_j                     
  deriv(R_j) = - death_rate*R_j + gamma_j*I_j - sigma_j*R_j  
  deriv(M_j) = - death_rate*M_j - sigma_m*M_j + birth_rate_M
  
  # adults during non-hibernation (note: migration is included here)
  deriv(S_a) = migration_rate_S - death_rate*S_a - beta*(I_j+I_a)*S_a + sigma_a*R_a 
  deriv(E_a) = migration_rate_E - death_rate*E_a + beta*(I_j+I_a)*S_a - tau_a*E_a                 
  deriv(I_a) = migration_rate_I - death_rate*I_a + tau_a*E_a - gamma_a*I_a                  
  deriv(R_a) = migration_rate_R - death_rate*R_a + gamma_a*I_a - sigma_a*R_a 
  
  # cumulative number of infections
  deriv(cum_inf_j) = beta*(I_j+I_a)*S_j
  deriv(cum_inf_a) = beta*(I_j+I_a)*S_a
  deriv(tunit) = 1/365  
  
  ## Initial states:
  initial(S_j) = 0 
  initial(E_j) = 0 
  initial(I_j) = 0 
  initial(R_j) = 0 
  initial(M_j) = 0 
  
  initial(S_a) = 0 
  initial(E_a) = 0  
  initial(I_a) = 0 
  initial(R_a) = 0
  
  initial(cum_inf_j) = 0
  initial(cum_inf_a) = 0
  initial(tunit) = 1/365
  
  times[] = user()
  dim(times) = user()
  
  # Parameters to be varied over time 
  beta0 = interpolate(times, beta0_0, "constant")
  beta0_0[] = user()
  dim(beta0_0) = length(times)
  
  beta1 = interpolate(times, beta1_0, "constant")
  beta1_0[] = user()
  dim(beta1_0) = length(times)
  
  beta2 = interpolate(times, beta2_0, "constant")
  beta2_0[] = user()
  dim(beta2_0) = length(times)
  
  tau_j = interpolate(times, tau_j_0, "constant")
  tau_j_0[] = user()
  dim(tau_j_0) = length(times)
  
  tau_a = interpolate(times, tau_a_0, "constant")
  tau_a_0[] = user()
  dim(tau_a_0) = length(times)
  
  gamma_j = interpolate(times, gamma_j_0, "constant")
  gamma_j_0[] = user()
  dim(gamma_j_0) = length(times)
  
  gamma_a = interpolate(times, gamma_a_0, "constant")
  gamma_a_0[] = user()
  dim(gamma_a_0) = length(times)
  
  sigma_j = interpolate(times, sigma_j_0, "constant")
  sigma_j_0[] = user()
  dim(sigma_j_0) = length(times)
  
  sigma_a = interpolate(times, sigma_a_0, "constant")
  sigma_a_0[] = user()
  dim(sigma_a_0) = length(times)
  
  sigma_m = interpolate(times, sigma_m_0, "constant")
  sigma_m_0[] = user()
  dim(sigma_m_0) = length(times)
  
  death_rate = interpolate(times, death_rate_0, "constant")
  death_rate_0[] = user()
  dim(death_rate_0) = length(times)
  
  birth_rate_S = interpolate(times, birth_rate_S_0, "constant")
  birth_rate_S_0[] = user()
  dim(birth_rate_S_0) = length(times)
  
  birth_rate_M = interpolate(times, birth_rate_M_0, "constant")
  birth_rate_M_0[] = user()
  dim(birth_rate_M_0) = length(times)
  
  migration_rate_S = interpolate(times, migration_rate_S_0, "constant")
  migration_rate_E = interpolate(times, migration_rate_E_0, "constant")
  migration_rate_I = interpolate(times, migration_rate_I_0, "constant")
  migration_rate_R = interpolate(times, migration_rate_R_0, "constant")
  migration_rate_S_0[] = user()
  migration_rate_E_0[] = user()
  migration_rate_I_0[] = user()
  migration_rate_R_0[] = user()
  dim(migration_rate_S_0) = length(times)
  dim(migration_rate_E_0) = length(times)
  dim(migration_rate_I_0) = length(times)
  dim(migration_rate_R_0) = length(times)
  
})
SEIR_generator_2_2 = odin::odin( {

  beta = beta0 + beta1*((sin(2*3.141593*(tunit-beta2))+1)/2)

  deriv(S_j) = - death_rate*S_j - beta*(I_j+I_a)*S_j + sigma_j*R_j + sigma_m*M_j + birth_rate_S
  deriv(E_j) = - death_rate*E_j + beta*(I_j+I_a)*S_j - tau_j*E_j     
  deriv(I_j) = - death_rate*I_j + tau_j*E_j - gamma_j*I_j                     
  deriv(R_j) = - death_rate*R_j + gamma_j*I_j - sigma_j*R_j  
  deriv(M_j) = - death_rate*M_j - sigma_m*M_j + birth_rate_M
  
  # adults during non-hibernation (note: migration is included here)
  deriv(S_a) = migration_rate_S - death_rate*S_a - beta*(I_j+I_a)*S_a + sigma_a*R_a 
  deriv(E_a) = migration_rate_E - death_rate*E_a + beta*(I_j+I_a)*S_a - tau_a*E_a                 
  deriv(I_a) = migration_rate_I - death_rate*I_a + tau_a*E_a - gamma_a*I_a                  
  deriv(R_a) = migration_rate_R - death_rate*R_a + gamma_a*I_a - sigma_a*R_a 
  
  # cumulative number of infections
  deriv(cum_inf_j) = beta*(I_j+I_a)*S_j
  deriv(cum_inf_a) = beta*(I_j+I_a)*S_a
  deriv(tunit) = 1/366  
  
  ## Initial states:
  initial(S_j) = 0 
  initial(E_j) = 0 
  initial(I_j) = 0 
  initial(R_j) = 0 
  initial(M_j) = 0 
  
  initial(S_a) = 1 
  initial(E_a) = 0  
  initial(I_a) = 0 
  initial(R_a) = 0
  
  initial(cum_inf_j) = 0
  initial(cum_inf_a) = 0
  initial(tunit) = 1/366
  
  times[] = user()
  dim(times) = user()
  
  # Parameters to be varied over time 
  beta0 = interpolate(times, beta0_0, "constant")
  beta0_0[] = user()
  dim(beta0_0) = length(times)
  
  beta1 = interpolate(times, beta1_0, "constant")
  beta1_0[] = user()
  dim(beta1_0) = length(times)
  
  beta2 = interpolate(times, beta2_0, "constant")
  beta2_0[] = user()
  dim(beta2_0) = length(times)
  
  tau_j = interpolate(times, tau_j_0, "constant")
  tau_j_0[] = user()
  dim(tau_j_0) = length(times)
  
  tau_a = interpolate(times, tau_a_0, "constant")
  tau_a_0[] = user()
  dim(tau_a_0) = length(times)
  
  gamma_j = interpolate(times, gamma_j_0, "constant")
  gamma_j_0[] = user()
  dim(gamma_j_0) = length(times)
  
  gamma_a = interpolate(times, gamma_a_0, "constant")
  gamma_a_0[] = user()
  dim(gamma_a_0) = length(times)
  
  sigma_j = interpolate(times, sigma_j_0, "constant")
  sigma_j_0[] = user()
  dim(sigma_j_0) = length(times)
  
  sigma_a = interpolate(times, sigma_a_0, "constant")
  sigma_a_0[] = user()
  dim(sigma_a_0) = length(times)
  
  sigma_m = interpolate(times, sigma_m_0, "constant")
  sigma_m_0[] = user()
  dim(sigma_m_0) = length(times)
  
  death_rate = interpolate(times, death_rate_0, "constant")
  death_rate_0[] = user()
  dim(death_rate_0) = length(times)
  
  birth_rate_S = interpolate(times, birth_rate_S_0, "constant")
  birth_rate_S_0[] = user()
  dim(birth_rate_S_0) = length(times)
  
  birth_rate_M = interpolate(times, birth_rate_M_0, "constant")
  birth_rate_M_0[] = user()
  dim(birth_rate_M_0) = length(times)
  
  migration_rate_S = interpolate(times, migration_rate_S_0, "constant")
  migration_rate_E = interpolate(times, migration_rate_E_0, "constant")
  migration_rate_I = interpolate(times, migration_rate_I_0, "constant")
  migration_rate_R = interpolate(times, migration_rate_R_0, "constant")
  migration_rate_S_0[] = user()
  migration_rate_E_0[] = user()
  migration_rate_I_0[] = user()
  migration_rate_R_0[] = user()
  dim(migration_rate_S_0) = length(times)
  dim(migration_rate_E_0) = length(times)
  dim(migration_rate_I_0) = length(times)
  dim(migration_rate_R_0) = length(times)
  
})

# ODE for Colony 1 
runSim_C1 = function(theta){
  
  out = vector(mode = "list", length = 7)
  SEIR_compiled_1 = SEIR_generator_1_C1$new(times = seq(1,365), 
                                            tau_j_0 = c(rep(theta[1],197), rep(1e-10, 168)),  
                                            tau_a_0 = c(rep(theta[1],197), rep(1e-10, 168)),  
                                            gamma_j_0 = c(rep(theta[2],197), rep(1e-10, 168)),    
                                            gamma_a_0 = c(rep(theta[2],197), rep(1e-10, 168)),    
                                            sigma_j_0 = c(rep(theta[3],197), rep(theta[4], 168)),
                                            sigma_a_0 = c(rep(theta[3],197), rep(theta[4], 168)),  
                                            sigma_m_0 = c(rep(theta[5],197), rep(theta[6], 168)),    
                                            beta0_0 = c(rep(theta[7],197), rep(0, 168)),  
                                            beta1_0 = c(rep(theta[8],197), rep(0, 168)), 
                                            beta2_0 = c(rep(theta[9],365)),  
                                            death_rate_0 = c(rep(death_cons, 365)),                
                                            birth_rate_S_0 = c(rep(0,55), pop_C1[1], rep(0,309)),
                                            birth_rate_M_0 = c(rep(0,55), 0, rep(0,309))) 
  out[[1]] = as.data.frame(SEIR_compiled_1$run(seq(1, 365), maxsteps = 1e7))
  
  # compute the proportion of bats in each compartment. It will be used to compute the number of bats to migrate in the next year 
  S_prop = sum(out[[1]][nrow(out[[1]]), c("S_j", "S_a")])/sum(out[[1]][nrow(out[[1]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
  E_prop = sum(out[[1]][nrow(out[[1]]), c("E_j", "E_a")])/sum(out[[1]][nrow(out[[1]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
  I_prop = sum(out[[1]][nrow(out[[1]]), c("I_j", "I_a")])/sum(out[[1]][nrow(out[[1]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
  R_prop = sum(out[[1]][nrow(out[[1]]), c("R_j", "M_j", "R_a")])/sum(out[[1]][nrow(out[[1]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
  
  for(i in 2:nyear) {
    
    # SEIR_pop: a vector representing the target number of bats to migrate for each compartment 
    SEIR_pop = vector(length = 4)
    SEIR_pop[1] = pop_C1[i]*S_prop
    SEIR_pop[2] = pop_C1[i]*E_prop
    SEIR_pop[3] = pop_C1[i]*I_prop
    SEIR_pop[4] = pop_C1[i]*R_prop
    
    # note that leap years (2012, 2016, 2020) have one extra date (Feb 29) 
    # for leap years (2012 a leap year, but the study starts from April)
    if(study_year[i] %in% c(2016, 2020)) { 
      SEIR_compiled_2 = SEIR_generator_2_2$new(times = seq(1,366), 
                                               tau_j_0 = c(rep(theta[1],197), rep(1e-10, 169)),  
                                               tau_a_0 = c(rep(theta[1],197), rep(1e-10, 169)),  
                                               gamma_j_0 = c(rep(theta[2],197), rep(1e-10, 169)),    
                                               gamma_a_0 = c(rep(theta[2],197), rep(1e-10, 169)),    
                                               sigma_j_0 = c(rep(theta[3],197), rep(theta[4], 169)),
                                               sigma_a_0 = c(rep(theta[3],197), rep(theta[4], 169)),  
                                               sigma_m_0 = c(rep(theta[5],197), rep(theta[6], 169)),    
                                               beta0_0 = c(rep(theta[7],197), rep(0, 169)),  
                                               beta1_0 = c(rep(theta[8],197), rep(0, 169)), 
                                               beta2_0 = c(rep(theta[9],366)),    
                                               death_rate_0 = c(rep(death_cons, 366)),                
                                               birth_rate_S_0 = c(rep(0,55), sum(SEIR_pop[1:3]), rep(0,310)),
                                               birth_rate_M_0 = c(rep(0,55), SEIR_pop[4], rep(0,310)), 
                                               migration_rate_S_0 = c(SEIR_pop[1], rep(0,365)),
                                               migration_rate_E_0 = c(SEIR_pop[2], rep(0,365)),
                                               migration_rate_I_0 = c(SEIR_pop[3], rep(0,365)),
                                               migration_rate_R_0 = c(SEIR_pop[4], rep(0,365)))
      out[[i]] = as.data.frame(SEIR_compiled_2$run(seq(1,366,1), maxsteps = 1e7))}
    
    # for none-leap years
    if(!study_year[i] %in% c(2016, 2020)) { 
      SEIR_compiled_2 = SEIR_generator_2_1$new(times = seq(1,365), 
                                               tau_j_0 = c(rep(theta[1],197), rep(1e-10, 168)),  
                                               tau_a_0 = c(rep(theta[1],197), rep(1e-10, 168)),  
                                               gamma_j_0 = c(rep(theta[2],197), rep(1e-10, 168)),    
                                               gamma_a_0 = c(rep(theta[2],197), rep(1e-10, 168)),    
                                               sigma_j_0 = c(rep(theta[3],197), rep(theta[4], 168)),
                                               sigma_a_0 = c(rep(theta[3],197), rep(theta[4], 168)),  
                                               sigma_m_0 = c(rep(theta[5],197), rep(theta[6], 168)),    
                                               beta0_0 = c(rep(theta[7],197), rep(0, 168)),  
                                               beta1_0 = c(rep(theta[8],197), rep(0, 168)), 
                                               beta2_0 = c(rep(theta[9],365)),  
                                               death_rate_0 = c(rep(death_cons, 365)),                
                                               birth_rate_S_0 = c(rep(0,55), sum(SEIR_pop[1:3]), rep(0,309)),
                                               birth_rate_M_0 = c(rep(0,55), SEIR_pop[4], rep(0,309)), 
                                               migration_rate_S_0 = c(SEIR_pop[1], rep(0,364)),
                                               migration_rate_E_0 = c(SEIR_pop[2], rep(0,364)),
                                               migration_rate_I_0 = c(SEIR_pop[3], rep(0,364)),
                                               migration_rate_R_0 = c(SEIR_pop[4], rep(0,364)))
      out[[i]] = as.data.frame(SEIR_compiled_2$run(seq(1,365,1), maxsteps = 1e7))}
    
    S_prop = sum(out[[i]][nrow(out[[i]]), c("S_j", "S_a")])/sum(out[[i]][nrow(out[[i]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
    E_prop = sum(out[[i]][nrow(out[[i]]), c("E_j", "E_a")])/sum(out[[i]][nrow(out[[i]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
    I_prop = sum(out[[i]][nrow(out[[i]]), c("I_j", "I_a")])/sum(out[[i]][nrow(out[[i]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
    R_prop = sum(out[[i]][nrow(out[[i]]), c("R_j", "M_j", "R_a")])/sum(out[[i]][nrow(out[[i]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
    
  }
  
  out = do.call("rbind", out)
  colnames(out)[1] = "time" 
  out$time = seq(1, (365*5+366*2), 1)
  
  return(out)
}


# ODE for Colony 2 
runSim_C2 = function(theta){
  
  out = vector(mode = "list", length = 7)
  SEIR_compiled_1 = SEIR_generator_1_C2$new(times = seq(1,365), 
                                            tau_j_0 = c(rep(theta[1],197), rep(1e-10, 168)),  
                                            tau_a_0 = c(rep(theta[1],197), rep(1e-10, 168)),  
                                            gamma_j_0 = c(rep(theta[2],197), rep(1e-10, 168)),    
                                            gamma_a_0 = c(rep(theta[2],197), rep(1e-10, 168)),    
                                            sigma_j_0 = c(rep(theta[3],197), rep(theta[4], 168)),
                                            sigma_a_0 = c(rep(theta[3],197), rep(theta[4], 168)),  
                                            sigma_m_0 = c(rep(theta[5],197), rep(theta[6], 168)),    
                                            beta0_0 = c(rep(theta[7],197), rep(0, 168)),  
                                            beta1_0 = c(rep(theta[8],197), rep(0, 168)), 
                                            beta2_0 = c(rep(theta[9],365)),   
                                            death_rate_0 = c(rep(death_cons, 365)),                
                                            birth_rate_S_0 = c(rep(0,55), pop_C2[1], rep(0,309)),
                                            birth_rate_M_0 = c(rep(0,55), 0, rep(0,309))) 
  out[[1]] = as.data.frame(SEIR_compiled_1$run(seq(1, 365), maxsteps = 1e7))
  
  # compute the proportion of bats in each compartment. It will be used to compute the number of bats to migrate in the next year 
  S_prop = sum(out[[1]][nrow(out[[1]]), c("S_j", "S_a")])/sum(out[[1]][nrow(out[[1]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
  E_prop = sum(out[[1]][nrow(out[[1]]), c("E_j", "E_a")])/sum(out[[1]][nrow(out[[1]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
  I_prop = sum(out[[1]][nrow(out[[1]]), c("I_j", "I_a")])/sum(out[[1]][nrow(out[[1]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
  R_prop = sum(out[[1]][nrow(out[[1]]), c("R_j", "M_j", "R_a")])/sum(out[[1]][nrow(out[[1]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
  
  for(i in 2:nyear) {
    
    # SEIR_pop: a vector representing the target number of bats to migrate for each compartment 
    SEIR_pop = vector(length = 4)
    SEIR_pop[1] = pop_C2[i]*S_prop
    SEIR_pop[2] = pop_C2[i]*E_prop
    SEIR_pop[3] = pop_C2[i]*I_prop
    SEIR_pop[4] = pop_C2[i]*R_prop
    
    # note that leap years (2012, 2016, 2020) have one extra date (Feb 29) 
    # for leap years (2012 a leap year, but the study starts from April)
    if(study_year[i] %in% c(2016, 2020)) { 
      SEIR_compiled_2 = SEIR_generator_2_2$new(times = seq(1,366), 
                                               tau_j_0 = c(rep(theta[1],197), rep(1e-10, 169)),  
                                               tau_a_0 = c(rep(theta[1],197), rep(1e-10, 169)),  
                                               gamma_j_0 = c(rep(theta[2],197), rep(1e-10, 169)),    
                                               gamma_a_0 = c(rep(theta[2],197), rep(1e-10, 169)),    
                                               sigma_j_0 = c(rep(theta[3],197), rep(theta[4], 169)),
                                               sigma_a_0 = c(rep(theta[3],197), rep(theta[4], 169)),  
                                               sigma_m_0 = c(rep(theta[5],197), rep(theta[6], 169)),    
                                               beta0_0 = c(rep(theta[7],197), rep(0, 169)),  
                                               beta1_0 = c(rep(theta[8],197), rep(0, 169)), 
                                               beta2_0 = c(rep(theta[9],366)),     
                                               death_rate_0 = c(rep(death_cons, 366)),                
                                               birth_rate_S_0 = c(rep(0,55), sum(SEIR_pop[1:3]), rep(0,310)),
                                               birth_rate_M_0 = c(rep(0,55), SEIR_pop[4], rep(0,310)), 
                                               migration_rate_S_0 = c(SEIR_pop[1], rep(0,365)),
                                               migration_rate_E_0 = c(SEIR_pop[2], rep(0,365)),
                                               migration_rate_I_0 = c(SEIR_pop[3], rep(0,365)),
                                               migration_rate_R_0 = c(SEIR_pop[4], rep(0,365)))
      out[[i]] = as.data.frame(SEIR_compiled_2$run(seq(1,366,1), maxsteps = 1e7))}
    
    # for none-leap years
    if(!study_year[i] %in% c(2016, 2020)) { 
      SEIR_compiled_2 = SEIR_generator_2_1$new(times = seq(1,365), 
                                               tau_j_0 = c(rep(theta[1],197), rep(1e-10, 168)),  
                                               tau_a_0 = c(rep(theta[1],197), rep(1e-10, 168)),  
                                               gamma_j_0 = c(rep(theta[2],197), rep(1e-10, 168)),    
                                               gamma_a_0 = c(rep(theta[2],197), rep(1e-10, 168)),    
                                               sigma_j_0 = c(rep(theta[3],197), rep(theta[4], 168)),
                                               sigma_a_0 = c(rep(theta[3],197), rep(theta[4], 168)),  
                                               sigma_m_0 = c(rep(theta[5],197), rep(theta[6], 168)),    
                                               beta0_0 = c(rep(theta[7],197), rep(0, 168)),  
                                               beta1_0 = c(rep(theta[8],197), rep(0, 168)), 
                                               beta2_0 = c(rep(theta[9],365)),     
                                               death_rate_0 = c(rep(death_cons, 365)),                
                                               birth_rate_S_0 = c(rep(0,55), sum(SEIR_pop[1:3]), rep(0,309)),
                                               birth_rate_M_0 = c(rep(0,55), SEIR_pop[4], rep(0,309)), 
                                               migration_rate_S_0 = c(SEIR_pop[1], rep(0,364)),
                                               migration_rate_E_0 = c(SEIR_pop[2], rep(0,364)),
                                               migration_rate_I_0 = c(SEIR_pop[3], rep(0,364)),
                                               migration_rate_R_0 = c(SEIR_pop[4], rep(0,364)))
      out[[i]] = as.data.frame(SEIR_compiled_2$run(seq(1,365,1), maxsteps = 1e7))}
    
    S_prop = sum(out[[i]][nrow(out[[i]]), c("S_j", "S_a")])/sum(out[[i]][nrow(out[[i]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
    E_prop = sum(out[[i]][nrow(out[[i]]), c("E_j", "E_a")])/sum(out[[i]][nrow(out[[i]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
    I_prop = sum(out[[i]][nrow(out[[i]]), c("I_j", "I_a")])/sum(out[[i]][nrow(out[[i]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
    R_prop = sum(out[[i]][nrow(out[[i]]), c("R_j", "M_j", "R_a")])/sum(out[[i]][nrow(out[[i]]), c("S_j","E_j","I_j","R_j","M_j","S_a","E_a","I_a","R_a")])
    
  }
  
  out = do.call("rbind", out)
  colnames(out)[1] = "time" 
  out$time = seq(1, (365*5+366*2), 1)
  
  return(out)
}
Like1 = function(C1A1_sero, C1A2_sero, C2A1_sero, C2A2_sero, 
                 out_C1, out_C2, 
                 sample_date_C1A1_sero, sample_date_C1A2_sero, sample_date_C2A1_sero, sample_date_C2A2_sero, 
                 Ts){
  
  p_C1A1_sero = rowSums(out_C1[sample_date_C1A1_sero, c("R_j","M_j")])/rowSums(out_C1[sample_date_C1A1_sero,c("S_j", "E_j", "I_j", "R_j", "M_j")])  
  p_C1A2_sero = out_C1[sample_date_C1A2_sero, c("R_a")]/rowSums(out_C1[sample_date_C1A2_sero,c("S_a", "E_a", "I_a", "R_a")]) 
  p_C2A1_sero = rowSums(out_C2[sample_date_C2A1_sero, c("R_j","M_j")])/rowSums(out_C2[sample_date_C2A1_sero,c("S_j", "E_j", "I_j", "R_j", "M_j")])  
  p_C2A2_sero = out_C2[sample_date_C2A2_sero, c("R_a")]/rowSums(out_C2[sample_date_C2A2_sero,c("S_a", "E_a", "I_a", "R_a")]) 
  
  #D = 5.077, 95% CI, 4.033 - 8.968.
  prior_gamma = dgamma(Ts[2], 15, 15*5.077, log = TRUE)
  
  L_sero = sum(dbinom(x = c(C1A1_sero$pos, C1A2_sero$pos, C2A1_sero$pos, C2A2_sero$pos),
                      size = c(C1A1_sero$total, C1A2_sero$total, C2A1_sero$total, C2A2_sero$total),
                      prob = c(p_C1A1_sero, p_C1A2_sero, p_C2A1_sero, p_C2A2_sero), log = TRUE), na.rm=TRUE)
  
  return(L_sero + prior_gamma)
} 

Like2 = function(C1A1_sero, C1A2_sero, C2A1_sero, C2A2_sero, 
                 out_C1, out_C2, 
                 sample_date_C1A1_sero, sample_date_C1A2_sero, sample_date_C2A1_sero, sample_date_C2A2_sero){
  
  p_C1A1_sero = rowSums(out_C1[sample_date_C1A1_sero, c("R_j","M_j")])/rowSums(out_C1[sample_date_C1A1_sero,c("S_j", "E_j", "I_j", "R_j", "M_j")])  
  p_C1A2_sero = out_C1[sample_date_C1A2_sero, c("R_a")]/rowSums(out_C1[sample_date_C1A2_sero,c("S_a", "E_a", "I_a", "R_a")]) 
  p_C2A1_sero = rowSums(out_C2[sample_date_C2A1_sero, c("R_j","M_j")])/rowSums(out_C2[sample_date_C2A1_sero,c("S_j", "E_j", "I_j", "R_j", "M_j")])  
  p_C2A2_sero = out_C2[sample_date_C2A2_sero, c("R_a")]/rowSums(out_C2[sample_date_C2A2_sero,c("S_a", "E_a", "I_a", "R_a")]) 
  
  L_sero = sum(dbinom(x = c(C1A1_sero$pos, C1A2_sero$pos, C2A1_sero$pos, C2A2_sero$pos),
                      size = c(C1A1_sero$total, C1A2_sero$total, C2A1_sero$total, C2A2_sero$total),
                      prob = c(p_C1A1_sero, p_C1A2_sero, p_C2A1_sero, p_C2A2_sero), log = TRUE), na.rm=TRUE)
  
  return(L_sero)
} 

################################################################################

data = read.csv(file = "Data.csv", header = TRUE) # import data
data$date = as.Date(data$date, "%Y-%m-%d")                        # make sure this is Date object

### source: Lyssavirus_summary.csv

# max number 
sum_max_C1 = data.frame(year = seq(2015, 2021, 1), number = c(400, 1720, 2700, 1153, 1100, NA, 836))
sum_max_C2 = data.frame(year = seq(2015, 2021, 1), number = c(NA, 1400, 1850, NA, 1850, NA, 250))

# inference number 
#sum_max_C1 = data.frame(year = seq(2012, 2021, 1), number = c(NA, NA, NA, 7.2*400, 7.2*350, 2700, 1153, 1100, NA, 7.2*300))
#sum_max_C2 = data.frame(year = seq(2012, 2021, 1), number = c(NA, NA, NA, NA, 7.2*300, 1850, NA, 1850, NA, 7.2*250))

sum_max_C1$mean = mean(sum_max_C1$number, na.rm = TRUE)
sum_max_C1$number[is.na(sum_max_C1$number)] = unique(sum_max_C1$mean) # assign the mean population size to the year without known population size
sum_max_C2$mean = mean(sum_max_C2$number, na.rm = TRUE)
sum_max_C2$number[is.na(sum_max_C2$number)] = unique(sum_max_C2$mean) # assign the mean population size to the year without known population size

nyear = 7
death_cons = 1/(15*365)
study_year = seq(2015,2021,1)

pop_C1 = sum_max_C1$number/2 
pop_C2 = sum_max_C2$number/2 

C1A1_sero = data[data$colony == "Colony 1" & data$age == "Young",]
C1A2_sero = data[data$colony == "Colony 1" & data$age == "Adult",]
C2A1_sero = data[data$colony == "Colony 2" & data$age == "Young",]
C2A2_sero = data[data$colony == "Colony 2" & data$age == "Adult",]

sum(C1A1_sero$total, C1A2_sero$total, C2A1_sero$total, C2A2_sero$total)

out = data.frame(time = seq(1, (365*5+366*2), 1)) # two study years had 366 days 
out$time = (1:nrow(out)) + as.Date("2015-03-31",'%Y-%m-%d')

sample_date_C1A1_sero = which(out$time %in% C1A1_sero$date) 
sample_date_C1A2_sero = which(out$time %in% C1A2_sero$date) 
sample_date_C2A1_sero = which(out$time %in% C2A1_sero$date) 
sample_date_C2A2_sero = which(out$time %in% C2A2_sero$date) 

################################################################################
### MCMC run ### 

for(c in 1:4) {
  
  prior_tau = c(1/58,1/7)            
  prior_gamma = c(0, 1)                 
  prior_sigma = c(1/(60*30),1/(2*30))  
  prior_beta0 = c(0, 10)   
  prior_beta1 = c(0, 10)   
  prior_beta2 = c(0, 1)   
  theta0 = c(runif(1, min=prior_tau[1],   max=prior_tau[2]), 
             runif(1, min=prior_gamma[1], max=prior_gamma[2]), 
             runif(4, min=prior_sigma[1], max=prior_sigma[2]), 
             runif(1, min=0,  max=0.001),
             runif(1, min=0,  max=0.001),
             runif(1, min=prior_beta2[1],  max=prior_beta2[2]))
  
  check = 100                         
  att = rep(0,length(theta0))
  acc = rep(0,length(theta0))
  
  prior = matrix(c(prior_tau, 
                   prior_gamma, 
                   rep(prior_sigma,4), 
                   prior_beta0, 
                   prior_beta1,
                   prior_beta2), 
                 nrow = 2, ncol = length(theta0))
  s = rep(1,length(theta0))
  
  iter = 1e5
  
  out_C1 = runSim_C1(theta = theta0)
  out_C2 = runSim_C2(theta = theta0)
  
  L0 = Like1(C1A1_sero = C1A1_sero, C1A2_sero = C1A2_sero, C2A1_sero = C2A1_sero, C2A2_sero = C2A2_sero, 
             out_C1 = out_C1, 
             out_C2 = out_C2, 
             sample_date_C1A1_sero = sample_date_C1A1_sero, sample_date_C1A2_sero = sample_date_C1A2_sero,
             sample_date_C2A1_sero = sample_date_C2A1_sero, sample_date_C2A2_sero = sample_date_C2A2_sero,
             Ts = theta0)
  
  theta = matrix(NA, nrow = iter, ncol = length(theta0)) 
  theta[1,] = theta0                                    
  logL = matrix(NA, nrow = iter, ncol = 1)             
  logL[1] = L0      
  
  for (i in 2:iter){
    
    for (j in 1:length(theta0)){
      
      att[j] = att[j] + 1  
      
      Ts = theta0                            
      Ts[j] = Ts[j]*exp(s[j]*rnorm(1,0,1))
      
      if (Ts[j] > prior[1,j] && Ts[j] < prior[2,j]){ 
        
        out_C1 = runSim_C1(theta = Ts)               
        out_C2 = runSim_C2(theta = Ts)               
        
        Ls1 = Like1(C1A1_sero = C1A1_sero, C1A2_sero = C1A2_sero, C2A1_sero = C2A1_sero, C2A2_sero = C2A2_sero, 
                    out_C1 = out_C1, 
                    out_C2 = out_C2, 
                    sample_date_C1A1_sero = sample_date_C1A1_sero, sample_date_C1A2_sero = sample_date_C1A2_sero,
                    sample_date_C2A1_sero = sample_date_C2A1_sero, sample_date_C2A2_sero = sample_date_C2A2_sero,
                    Ts = Ts)
        
        Ls2 = Like2(C1A1_sero = C1A1_sero, C1A2_sero = C1A2_sero, C2A1_sero = C2A1_sero, C2A2_sero = C2A2_sero, 
                    out_C1 = out_C1, 
                    out_C2 = out_C2, 
                    sample_date_C1A1_sero = sample_date_C1A1_sero, sample_date_C1A2_sero = sample_date_C1A2_sero,
                    sample_date_C2A1_sero = sample_date_C2A1_sero, sample_date_C2A2_sero = sample_date_C2A2_sero)
        
        r = exp(Ls1-L0)*Ts[j]/theta0[j] 
        
        if (runif(1,0,1) <= r) { 
          theta0 = Ts 
          L0 = Ls1
          LL = Ls2
          acc[j] = acc[j] + 1
        }}
      
      if(att[j] == check){
        print(paste0("Can sd of ", round(s[j], 3),
                     " for par[",j,"] gave acc rate ", acc[j]/att[j], " at ",i," iterations")) 
        if(acc[j]/att[j] < 0.2){ s[j] = s[j]*0.8 }
        if(acc[j]/att[j] > 0.4){ s[j] = s[j]*1.2 }
        acc[j] = att[j] = 0 
      }
      
    }
    
    logL[i] = LL
    theta[i,] = theta0
    
    t = i %% 301 
    if(t == 0) {
      graphics.off()
      par(mfrow=c(3,3))
      for(k in 1:length(theta0)) {
            plot(theta[,k][c(1:which(is.na(theta[,k]))[1]-1)], type = "l", ylab = "", xlab = theta[,k][which(is.na(theta[,k]))[1]-1])
            title(main= paste("par", k, sep = ":")) }
      plot(logL[c(1:which(is.na(logL))[1]-1)], type = "l", ylab = "LLike", xlab = logL[which(is.na(logL))[1]-1])
    }}
  write.csv(theta, file = paste("Model H2_", c, ".csv" , sep=""), row.names = FALSE)
}

################################################################################

post = vector(mode = "list", length = 4)
post[[1]] = read.csv(file = "Model H2_1.csv", header = TRUE)
post[[2]] = read.csv(file = "Model H2_2.csv", header = TRUE)
post[[3]] = read.csv(file = "Model H2_3.csv", header = TRUE)
post[[4]] = read.csv(file = "Model H2_4.csv", header = TRUE)

library(coda)
source("DBDA2E-utilities.R")

burn = 5e3
iter = 5e4
post_mcmc = vector(mode = "list", length = 4)
for(i in 1:4) {
  post_mcmc[[i]] = coda::mcmc(post[[i]][burn:iter,])}
post_mcmc_list = coda::mcmc.list(post_mcmc)

coda::effectiveSize(post_mcmc_list)
coda::autocorr.plot(post_mcmc_list, lag.max = 150)
gelman = coda::gelman.diag(post_mcmc_list)
gelman$psrf

