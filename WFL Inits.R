# Define parameter initial values:
inits <- list(log_AvgR = 15,
              init_recDevs = rep(0, length(data$ages)),  # Domain = [-5,5].
              recDevs = rep(0, length(data$years)-1))    # Domain = [-5,5].
              
inits$log_M_juvenile <- rep(-0.7, data$n_M)  # Domain = [-3,0.69]. 
inits$log_M_adult <- rep(-1.2, data$n_M)     # Domain = [-3,0.69].

inits$log_S50 <- c(1.128240868, 1.045758089) # Domain = [0.6,2].
inits$log_S95_delta <- c(1.67E-09, 1.67E-09) # Domain = [0,2].
inits$log_q_rv <- -1.0              # Domain = [-5,1.5].
inits$log_sigma_rv <- c(-0.3, -0.8) # Domain = [-3.0,1.5].
inits$log_Finit <- -2.5             # Domain = [-3.0,0.69]
inits$logit_gamma_R <- -4.000

# Fixed values:
inits$sigma_R <- 0.5
inits$sigma_M <- c(0.05, 0.025)
inits$Minit_prior <- c(0.8, 0.3)
inits$log_q_prior_rv <- -1.14
inits$sigma_q_rv <- 0.2

inits$log_q_sen <- -1.0              # Domain = [-5,1.5].
inits$sigma_q_sen <- 0.2
inits$log_q_prior_sen <- -1.14
inits$log_sigma_sen <- c(-0.3) # Domain = [-3.0,1.5].
