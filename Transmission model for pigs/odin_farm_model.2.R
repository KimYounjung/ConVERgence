################################################################################
# base model 
SEIR_cont <- odin::odin({
  
  n_pens <- user()
  n <- n_pens

  update(S[]) <- S[i] + n_RS[i] - n_SE[i] 
  update(E[]) <- E[i] + n_SE[i] - n_EI[i]
  update(I[]) <- I[i] + n_EI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i] - n_RS[i]
  
  lambda_prod[ , ] <- (R0_pen[i, j]*gamma/n_animals[i]) * I[j]
  lambda[] <- sum(lambda_prod[i, ]) 
  
  p_SE[] <- 1 - exp(-lambda[i]) # S to E
  p_EI   <- 1 - exp(-sigma)     # E to I
  p_IR   <- 1 - exp(-gamma)     # I to R
  p_RS   <- 1 - exp(-tau)       # R to S
  
  n_SE[] <- rbinom(S[i], p_SE[i])
  n_EI[] <- rbinom(E[i], p_EI)
  n_IR[] <- rbinom(I[i], p_IR)
  n_RS[] <- rbinom(R[i], p_RS)
  
  initial(S[]) <- S_ini[i]
  initial(E[]) <- E_ini[i]
  initial(I[]) <- I_ini[i]
  initial(R[]) <- R_ini[i]
  
  n_animals[] <- user()
  R0_pen[,]   <- user()   # pen-to-pen basic reproduction number, i.e. next generation matrix 
  sigma       <- user()   # rate of breakdown to active disease
  gamma       <- user()   # rate of recovery from active disease
  tau         <- user()   # background mortality
  S_ini[]     <- user()
  E_ini[]     <- user()
  I_ini[]     <- user()
  R_ini[]     <- user()
  
  ## dimensions
  dim(n_animals)   <- n
  dim(R0_pen)      <- c(n, n)
  dim(lambda_prod) <- c(n, n)
  dim(lambda)      <- n
  dim(S)           <- n
  dim(E)           <- n
  dim(I)           <- n
  dim(R)           <- n
  dim(n_SE)        <- n
  dim(n_EI)        <- n
  dim(n_IR)        <- n
  dim(n_RS)        <- n
  dim(p_SE)        <- n
  dim(S_ini)       <- n
  dim(E_ini)       <- n
  dim(I_ini)       <- n
  dim(R_ini)       <- n
  
})

################################################################################
### transmission paramters ### 
R0_int = seq(0.1, 0.6, 0.1)                  # R0 within pens to be explored
R0_btw = R0_int/10                           # R0 between neighbouring pens 
R0_house = R0_int/100                        # R0 between houses to be explored 
tau_range = rev(1/seq(30, 360, 30))          # immune duration to be explored
exponent = 0.2                               # exponent to model the extent of decrease of R0 by diatance between pens 
iter = 500                                   # no. of iteration for given parameter values
times = 36500                                # no. of days to be simulated in each iteration 
sigma = 1/2                                  # rate of becoming infectious after infection 
gamma = 1/5                                  # rate of recovery and becoming immune 
overall_R0 = vector(length = length(R0_int))

################################################################################
### farm structure ### 
# no. of houses for a given production stage 
S_nP_n_house = 10 # Sows in S_nP houses 1-9 are dry sows, and sows in S_nP house 10 are lactating sows mixed with their piglets 
S_yP_n_house = 1  # Those piglets are actually housed with sows in S_nP house 10 (labelled S_yP here to distinguish). 
W_n_house = 4     # Piglets move to W house 1 as weaners. Every 4 week, they move to the next W house (e.g. W house 2, 3, and 4)
F_n_house = 4     # Weaners in W house 4 move to F house 1 as fatteners. Every 4 week, they move to the next F house (e.g. F house 2, 3, and 4). They are slaughtered after F house 4. 

# no. of pens per house for a given production stage
S_nP_n_pen = 40   
S_yP_n_pen = 40
F_n_pen = 20
W_n_pen = 20

# no. of animals per pen for a given production stage
S_nP_n_animal =  1 # One sow per pen, but note that sows in S_nP house 10 are housed with their piglets in each pen. 
S_yP_n_animal = 10 # Note that piglets are housed with their sows in S_nP house 10 
W_n_animal = 20 
F_n_animal = 20 

# total no. of pens 
n_animals = c(rep(S_nP_n_animal, S_nP_n_pen*S_nP_n_house),
              rep(S_yP_n_animal, S_yP_n_pen*S_yP_n_house),
              rep(W_n_animal, W_n_pen*W_n_house),
              rep(F_n_animal, F_n_pen*F_n_house))
n_pens = length(n_animals)
n_newborn = 10

# run for the 1st production cycle (28 days)
d_cycle = 28
n_cycle = round(times/d_cycle)

pen_index = c(rep("S1_S", 40),rep("S2_S", 40),rep("S3_S", 40),rep("S4_S", 40),rep("S5_S", 40),
              rep("S6_S", 40),rep("S7_S", 40),rep("S8_S", 40),rep("S9_S", 40),rep("S10_S", 40), 
              rep("P_S", 40), 
              rep("W1_S", 20),rep("W2_S", 20),rep("W3_S", 20),rep("W4_S", 20), 
              rep("F1_S", 20),rep("F2_S", 20),rep("F3_S", 20),rep("F4_S", 20), 
              rep("S1_E", 40),rep("S2_E", 40),rep("S3_E", 40),rep("S4_E", 40),rep("S5_E", 40),
              rep("S6_E", 40),rep("S7_E", 40),rep("S8_E", 40),rep("S9_E", 40),rep("S10_E", 40), 
              rep("P_E", 40), 
              rep("W1_E", 20),rep("W2_E", 20),rep("W3_E", 20),rep("W4_E", 20), 
              rep("F1_E", 20),rep("F2_E", 20),rep("F3_E", 20),rep("F4_E", 20), 
              rep("S1_I", 40),rep("S2_I", 40),rep("S3_I", 40),rep("S4_I", 40),rep("S5_I", 40),
              rep("S6_I", 40),rep("S7_I", 40),rep("S8_I", 40),rep("S9_I", 40),rep("S10_I", 40), 
              rep("P_I", 40), 
              rep("W1_I", 20),rep("W2_I", 20),rep("W3_I", 20),rep("W4_I", 20), 
              rep("F1_I", 20),rep("F2_I", 20),rep("F3_I", 20),rep("F4_I", 20), 
              rep("S1_R", 40),rep("S2_R", 40),rep("S3_R", 40),rep("S4_R", 40),rep("S5_R", 40),
              rep("S6_R", 40),rep("S7_R", 40),rep("S8_R", 40),rep("S9_R", 40),rep("S10_R", 40), 
              rep("P_R", 40), 
              rep("W1_R", 20),rep("W2_R", 20),rep("W3_R", 20),rep("W4_R", 20), 
              rep("F1_R", 20),rep("F2_R", 20),rep("F3_R", 20),rep("F4_R", 20))

extinct = vector(mode = "list", length = length(R0_int))
for(i in 1:length(R0_int)) {
  extinct[[i]] = vector(mode = "list", length = length(tau_range))
  for(j in 1:length(tau_range)) {
    extinct[[i]][[j]] = vector(length = iter) }}

extinct.median = matrix(nrow = length(R0_int), ncol = length(tau_range))
extinct.0.025 = matrix(nrow = length(R0_int), ncol = length(tau_range))
extinct.0.975 = matrix(nrow = length(R0_int), ncol = length(tau_range))

seed = vector(mode = "list", length = length(R0_int))
for(i in 1:length(R0_int)) {
  seed[[i]] = vector(mode = "list", length = length(tau_range))
  for(j in 1:length(tau_range)) {
    seed[[i]][[j]] = vector(length = iter)}}

################################################################################
### Simulation run ###  
# For a given R0, 

for(R0_i in 1:length(R0_int)) {     
  
  # First, construct R0 arising from pens in the same house. 
  # For a given pair of pens, R0 is assumed to decrease exponentially by their distance.
  S_nP_R0 = matrix(R0_int[R0_i], nrow=S_nP_n_pen, ncol=S_nP_n_pen) 
  S_yP_R0 = matrix(R0_int[R0_i], nrow=S_yP_n_pen, ncol=S_yP_n_pen)
  W_R0 = matrix(R0_int[R0_i], nrow=W_n_pen, ncol=W_n_pen)
  F_R0 = matrix(R0_int[R0_i], nrow=F_n_pen, ncol=F_n_pen)
  
  for(pen_i in 1:S_nP_n_pen) {
    for(pen_j in 1:S_nP_n_pen) {
      if(pen_i > pen_j) {
        t = pen_i - pen_j 
        if(t%%2 == 0) S_nP_R0[pen_i,pen_j] = S_nP_R0[pen_j,pen_i] = R0_btw[R0_i]*exp(-exponent*abs(t - 1))       # R0 for pigs in different pens in the same line (R0 decreases exponentially by distance)
        if(t%%2 == 1) S_nP_R0[pen_i,pen_j] = S_nP_R0[pen_j,pen_i] = R0_btw[R0_i]*exp(-exponent*abs(t + 1))  }}}  # R0 for pigs in different pens in different lines (R0 decreases exponentially by distance)
  
  for(pen_i in 1:S_yP_n_pen) {
    for(pen_j in 1:S_yP_n_pen) {
      if(pen_i > pen_j) {
        t = pen_i - pen_j 
        if(t%%2 == 0) S_yP_R0[pen_i,pen_j] = S_yP_R0[pen_j,pen_i] = R0_btw[R0_i]*exp(-exponent*abs(t - 1))      # R0 for pigs in different pens in the same line (R0 decreases exponentially by distance)
        if(t%%2 == 1) S_yP_R0[pen_i,pen_j] = S_yP_R0[pen_j,pen_i] = R0_btw[R0_i]*exp(-exponent*abs(t + 1))  }}} # R0 for pigs in different pens in different lines (R0 decreases exponentially by distance)
  
  for(pen_i in 1:W_n_pen) {
    for(pen_j in 1:W_n_pen) {
      if(pen_i > pen_j) {
        t = pen_i - pen_j 
        if(t%%2 == 0) W_R0[pen_i,pen_j] = W_R0[pen_j,pen_i] = R0_btw[R0_i]*exp(-exponent*abs(t - 1))       # R0 for pigs in different pens in the same line (R0 decreases exponentially by distance)
        if(t%%2 == 1) W_R0[pen_i,pen_j] = W_R0[pen_j,pen_i] = R0_btw[R0_i]*exp(-exponent*abs(t + 1))   }}} # R0 for pigs in different pens in different lines (R0 decreases exponentially by distance)
  
  for(pen_i in 1:F_n_pen) {
    for(pen_j in 1:F_n_pen) {
      if(pen_i > pen_j) {
        t = pen_i - pen_j  
        if(t%%2 == 0) F_R0[pen_i,pen_j] = F_R0[pen_j,pen_i] = R0_btw[R0_i]*exp(-exponent*abs(t - 1))       # R0 for pigs in different pens in the same line (R0 decreases exponentially by distance)
        if(t%%2 == 1) F_R0[pen_i,pen_j] = F_R0[pen_j,pen_i] = R0_btw[R0_i]*exp(-exponent*abs(t + 1))   }}} # R0 for pigs in different pens in different lines (R0 decreases exponentially by distance)
  
  # Create a NGM template including all pens on the farm (n x n matrix, n=total number of pens on the farm).
  # First, fill the NGM with between-pen R0 values. 
  # Note that here we create a NGM for each between-R0 house value.
  R0_pen = matrix(R0_house[R0_i], nrow = n_pens, ncol = n_pens)
  
  # Second, fill NGM entries for pens in the same house with R0 arising from pens in the same house  
  for(house in 1:S_nP_n_house) {
    a = S_nP_n_pen*(house - 1) + 1
    b = S_nP_n_pen*house
    R0_pen[a:b, a:b] = S_nP_R0  }
  
  for(house in 1:S_yP_n_house) {
    c = S_yP_n_pen*(house - 1) + 1 + b
    d = S_yP_n_pen*house + b
    R0_pen[c:d, c:d] = S_yP_R0  }
  
  # Note that we assume lactating sows are housed with piglets. 
  # Therefore, lactating sows and piglets have R0 arising from pens in the same house.   
  R0_pen[a:b, c:d] = S_yP_R0
  R0_pen[c:d, a:b] = S_yP_R0
  
  for(house in 1:W_n_house) {
    e = W_n_pen*(house - 1) + 1 + d
    f = W_n_pen*house + d
    R0_pen[e:f, e:f] = W_R0 }
  
  for(house in 1:F_n_house) {
    g = F_n_pen*(house - 1) + 1 + f
    h = F_n_pen*house + f
    R0_pen[g:h, g:h] = F_R0 }
  
  # Now compute the farm-level R0, i.e. the dominant eigenvalue of the NGM 
  R0_eigen = eigen(R0_pen)
  overall_R0[R0_i] = as.numeric(R0_eigen$values[1]) 
  
  # For a given tau (i.e. rate of becoming susceptible, R to S)
  for(tau_i in 1:length(tau_range)) {     
    
    for(iteration in 1:iter) {
      
      S_ini = n_animals
      E_ini = rep(0, n_pens)
      I_ini = rep(0, n_pens)
      R_ini = rep(0, n_pens)
      
      # Randomly seed an infection
      
      r = sample(1:n_pens, 1)
      seed[[R0_i]][[tau_i]][iteration] = r
      S_ini[r] = S_ini[r] - 1
      I_ini[r] = 1
      
      model = SEIR_cont$new(n_pens = n_pens, 
                            n_animals = n_animals,
                            R0_pen = R0_pen, 
                            sigma = 1/2, 
                            gamma = 1/5,
                            tau = tau_range[tau_i],
                            S_ini = S_ini,
                            E_ini = E_ini,
                            I_ini = I_ini,
                            R_ini = R_ini)
      out = model$run(0:28)
      colnames(out) = c("time", pen_index)
      
      inspect.1 = out[,c(which(colnames(out) == "S1_E"), which(colnames(out) == "S2_E"), which(colnames(out) == "S3_E"), which(colnames(out) == "S4_E"),
                         which(colnames(out) == "S5_E"), which(colnames(out) == "S6_E"), which(colnames(out) == "S7_E"), which(colnames(out) == "S8_E"),
                         which(colnames(out) == "S9_E"), which(colnames(out) == "S10_E"),
                         which(colnames(out) == "P_E"),
                         which(colnames(out) == "W1_E"), which(colnames(out) == "W2_E"), which(colnames(out) == "W3_E"), which(colnames(out) == "W4_E"),     
                         which(colnames(out) == "F1_E"), which(colnames(out) == "F2_E"), which(colnames(out) == "F3_E"), which(colnames(out) == "F4_E"), 
                         which(colnames(out) == "S1_I"), which(colnames(out) == "S2_I"), which(colnames(out) == "S3_I"), which(colnames(out) == "S4_I"),
                         which(colnames(out) == "S5_I"), which(colnames(out) == "S6_I"), which(colnames(out) == "S7_I"), which(colnames(out) == "S8_I"),
                         which(colnames(out) == "S9_I"), which(colnames(out) == "S10_I"),
                         which(colnames(out) == "P_I"),
                         which(colnames(out) == "W1_I"), which(colnames(out) == "W2_I"), which(colnames(out) == "W3_I"), which(colnames(out) == "W4_I"),     
                         which(colnames(out) == "F1_I"), which(colnames(out) == "F2_I"), which(colnames(out) == "F3_I"), which(colnames(out) == "F4_I"))]
      
      inspect.1 = rowSums(inspect.1) 
      inspect.2 = which(inspect.1 == 0) 
      
      if(length(inspect.2) > 0) {
        extinct[[R0_i]][[tau_i]][iteration] = inspect.2[1] 
        print(inspect.2[1]) }
      
      # if not extinct, proceed to the next production cycle
      if(length(inspect.2) == 0) {
        
        for(cycle_i in 2:n_cycle) {
          
          # determine the proportion of 'exposed' among piglets, which depends on the proportion of 'exposed'/'infected' 
          # piglets move to weaner group 1
          new_W1_S = out[nrow(out), which(colnames(out) == "P_S")] 
          new_W1_E = out[nrow(out), which(colnames(out) == "P_E")] 
          new_W1_I = out[nrow(out), which(colnames(out) == "P_I")] 
          new_W1_R = out[nrow(out), which(colnames(out) == "P_R")] 
          
          new_W1_S = c(sum(new_W1_S[1:2]),  sum(new_W1_S[3:4]),  sum(new_W1_S[5:6]),  sum(new_W1_S[7:8]),  sum(new_W1_S[9:10]),
                       sum(new_W1_S[11:12]),sum(new_W1_S[13:14]),sum(new_W1_S[15:16]),sum(new_W1_S[17:18]),sum(new_W1_S[19:20]),
                       sum(new_W1_S[21:22]),sum(new_W1_S[23:24]),sum(new_W1_S[25:26]),sum(new_W1_S[27:28]),sum(new_W1_S[29:30]),
                       sum(new_W1_S[31:32]),sum(new_W1_S[33:34]),sum(new_W1_S[35:36]),sum(new_W1_S[37:38]),sum(new_W1_S[39:40]))
          new_W1_E = c(sum(new_W1_E[1:2]),  sum(new_W1_E[3:4]),  sum(new_W1_E[5:6]),  sum(new_W1_E[7:8]),  sum(new_W1_E[9:10]),
                       sum(new_W1_E[11:12]),sum(new_W1_E[13:14]),sum(new_W1_E[15:16]),sum(new_W1_E[17:18]),sum(new_W1_E[19:20]),
                       sum(new_W1_E[21:22]),sum(new_W1_E[23:24]),sum(new_W1_E[25:26]),sum(new_W1_E[27:28]),sum(new_W1_E[29:30]),
                       sum(new_W1_E[31:32]),sum(new_W1_E[33:34]),sum(new_W1_E[35:36]),sum(new_W1_E[37:38]),sum(new_W1_E[39:40]))
          new_W1_I = c(sum(new_W1_I[1:2]),  sum(new_W1_I[3:4]),  sum(new_W1_I[5:6]),  sum(new_W1_I[7:8]),  sum(new_W1_I[9:10]),
                       sum(new_W1_I[11:12]),sum(new_W1_I[13:14]),sum(new_W1_I[15:16]),sum(new_W1_I[17:18]),sum(new_W1_I[19:20]),
                       sum(new_W1_I[21:22]),sum(new_W1_I[23:24]),sum(new_W1_I[25:26]),sum(new_W1_I[27:28]),sum(new_W1_I[29:30]),
                       sum(new_W1_I[31:32]),sum(new_W1_I[33:34]),sum(new_W1_I[35:36]),sum(new_W1_I[37:38]),sum(new_W1_I[39:40]))
          new_W1_R = c(sum(new_W1_R[1:2]),  sum(new_W1_R[3:4]),  sum(new_W1_R[5:6]),  sum(new_W1_R[7:8]),  sum(new_W1_R[9:10]),
                       sum(new_W1_R[11:12]),sum(new_W1_R[13:14]),sum(new_W1_R[15:16]),sum(new_W1_R[17:18]),sum(new_W1_R[19:20]),
                       sum(new_W1_R[21:22]),sum(new_W1_R[23:24]),sum(new_W1_R[25:26]),sum(new_W1_R[27:28]),sum(new_W1_R[29:30]),
                       sum(new_W1_R[31:32]),sum(new_W1_R[33:34]),sum(new_W1_R[35:36]),sum(new_W1_R[37:38]),sum(new_W1_R[39:40]))
          
          # prepare new initial groups for the next run by shifting pig groups to the next production cycle. Note that new piglets' infection status depends on their mothers' infection status. 
          new_P_S = out[nrow(out), c(which(colnames(out) == "S9_S"))]*n_newborn
          new_P_E = (out[nrow(out), c(which(colnames(out) == "S9_E"))] + out[nrow(out), c(which(colnames(out) == "S9_I"))])*n_newborn
          new_P_I = rep(0, S_yP_n_pen*S_yP_n_house)
          new_P_R = out[nrow(out), c(which(colnames(out) == "S9_R"))]*n_newborn
          
          new_S_S = out[nrow(out), c(which(colnames(out) == "S10_S"), which(colnames(out) == "S1_S"),
                                     which(colnames(out) == "S2_S"),  which(colnames(out) == "S3_S"),
                                     which(colnames(out) == "S4_S"),  which(colnames(out) == "S5_S"),
                                     which(colnames(out) == "S6_S"),  which(colnames(out) == "S7_S"),
                                     which(colnames(out) == "S8_S"),  which(colnames(out) == "S9_S"))]
          new_W_S = out[nrow(out), c(which(colnames(out) == "W1_S"),  which(colnames(out) == "W2_S"), which(colnames(out) == "W3_S"))]
          new_F_S = out[nrow(out), c(which(colnames(out) == "W4_S"),  which(colnames(out) == "F1_S"),
                                     which(colnames(out) == "F2_S"),  which(colnames(out) == "F3_S"))]
          
          new_S_E = out[nrow(out), c(which(colnames(out) == "S10_E"), which(colnames(out) == "S1_E"),
                                     which(colnames(out) == "S2_E"),  which(colnames(out) == "S3_E"),
                                     which(colnames(out) == "S4_E"),  which(colnames(out) == "S5_E"),
                                     which(colnames(out) == "S6_E"),  which(colnames(out) == "S7_E"),
                                     which(colnames(out) == "S8_E"),  which(colnames(out) == "S9_E"))]
          new_W_E = out[nrow(out), c(which(colnames(out) == "W1_E"),  which(colnames(out) == "W2_E"), which(colnames(out) == "W3_E"))]
          new_F_E = out[nrow(out), c(which(colnames(out) == "W4_E"),  which(colnames(out) == "F1_E"),
                                     which(colnames(out) == "F2_E"),  which(colnames(out) == "F3_E"))]
          
          new_S_I = out[nrow(out), c(which(colnames(out) == "S10_I"), which(colnames(out) == "S1_I"),
                                     which(colnames(out) == "S2_I"),  which(colnames(out) == "S3_I"),
                                     which(colnames(out) == "S4_I"),  which(colnames(out) == "S5_I"),
                                     which(colnames(out) == "S6_I"),  which(colnames(out) == "S7_I"),
                                     which(colnames(out) == "S8_I"),  which(colnames(out) == "S9_I"))]
          new_W_I = out[nrow(out), c(which(colnames(out) == "W1_I"),  which(colnames(out) == "W2_I"), which(colnames(out) == "W3_I"))]
          new_F_I = out[nrow(out), c(which(colnames(out) == "W4_I"),  which(colnames(out) == "F1_I"),
                                     which(colnames(out) == "F2_I"),  which(colnames(out) == "F3_I"))]
          
          new_S_R = out[nrow(out), c(which(colnames(out) == "S10_R"), which(colnames(out) == "S1_R"),
                                     which(colnames(out) == "S2_R"),  which(colnames(out) == "S3_R"),
                                     which(colnames(out) == "S4_R"),  which(colnames(out) == "S5_R"),
                                     which(colnames(out) == "S6_R"),  which(colnames(out) == "S7_R"),
                                     which(colnames(out) == "S8_R"),  which(colnames(out) == "S9_R"))]
          new_W_R = out[nrow(out), c(which(colnames(out) == "W1_R"),  which(colnames(out) == "W2_R"), which(colnames(out) == "W3_R"))]
          new_F_R = out[nrow(out), c(which(colnames(out) == "W4_R"),  which(colnames(out) == "F1_R"),
                                     which(colnames(out) == "F2_R"),  which(colnames(out) == "F3_R"))]
          
          # run the current production cycle
          model = SEIR_cont$new(n_pens = n_pens, 
                                n_animals = n_animals,
                                R0_pen = R0_pen, 
                                sigma = 1/2, 
                                gamma = 1/5,
                                tau = tau_range[tau_i],
                                S_ini = as.vector(c(new_S_S, new_P_S, new_W1_S, new_W_S, new_F_S)),
                                E_ini = as.vector(c(new_S_E, new_P_E, new_W1_E, new_W_E, new_F_E)),
                                I_ini = as.vector(c(new_S_I, new_P_I, new_W1_I, new_W_I, new_F_I)),
                                R_ini = as.vector(c(new_S_R, new_P_R, new_W1_R, new_W_R, new_F_R)))
          
          out = model$run(0:28)
          colnames(out) = c("time", pen_index)
          
          # check if the transmission has become extinct 
          inspect.3 = out[,c(which(colnames(out) == "S1_E"), which(colnames(out) == "S2_E"),
                             which(colnames(out) == "S3_E"), which(colnames(out) == "S4_E"),
                             which(colnames(out) == "S5_E"), which(colnames(out) == "S6_E"),
                             which(colnames(out) == "S7_E"), which(colnames(out) == "S8_E"),
                             which(colnames(out) == "S9_E"), which(colnames(out) == "S10_E"),
                             which(colnames(out) == "P_E"),
                             which(colnames(out) == "W1_E"), which(colnames(out) == "W2_E"),
                             which(colnames(out) == "W3_E"), which(colnames(out) == "W4_E"),     
                             which(colnames(out) == "F1_E"), which(colnames(out) == "F2_E"),
                             which(colnames(out) == "F3_E"), which(colnames(out) == "F4_E"), 
                             which(colnames(out) == "S1_I"), which(colnames(out) == "S2_I"),
                             which(colnames(out) == "S3_I"), which(colnames(out) == "S4_I"),
                             which(colnames(out) == "S5_I"), which(colnames(out) == "S6_I"),
                             which(colnames(out) == "S7_I"), which(colnames(out) == "S8_I"),
                             which(colnames(out) == "S9_I"), which(colnames(out) == "S10_I"),
                             which(colnames(out) == "P_I"),
                             which(colnames(out) == "W1_I"), which(colnames(out) == "W2_I"),
                             which(colnames(out) == "W3_I"), which(colnames(out) == "W4_I"),     
                             which(colnames(out) == "F1_I"), which(colnames(out) == "F2_I"),
                             which(colnames(out) == "F3_I"), which(colnames(out) == "F4_I"))]

          inspect.3 = rowSums(inspect.3) 
          inspect.4 = which(inspect.3 == 0) 
          
          # if extinct, stop the production cycle run, and save the day of extinction 
          if(length(inspect.4) > 0) {
            extinct[[R0_i]][[tau_i]][iteration] = inspect.4[1] + (cycle_i-1)*d_cycle
            print(inspect.4[1] + (cycle_i-1)*d_cycle)
            break    }
          
          # if not extinct until the end, save the last day of iteration 
          if(length(inspect.4) == 0 & cycle_i == n_cycle) {
            extinct[[R0_i]][[tau_i]][iteration] = n_cycle*d_cycle
            print(n_cycle*d_cycle)
            break    } 
          
        }}}
    
    extinct.summary = quantile(extinct[[R0_i]][[tau_i]], c(0.5, 0.025, 0.975))
    extinct.median[R0_i,tau_i] = extinct.summary[1]
    extinct.0.025[R0_i,tau_i] = extinct.summary[2]
    extinct.0.975[R0_i,tau_i] = extinct.summary[3]
    hist(extinct[[R0_i]][[tau_i]], 
         main = paste0("median extinction day of ", extinct.median[R0_i,tau_i],
                       " with pen_R0 of ", R0_int[R0_i],
                       " with overall_R0 of ", format(round(overall_R0[R0_i],2),2), 
                       " and Tau of ", 1/tau_range[tau_i], " days"),
         xlim = c(0, max(extinct[[R0_i]][[tau_i]])), breaks = 50)    
  }}
