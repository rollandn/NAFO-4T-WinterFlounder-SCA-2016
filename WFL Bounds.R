# Define parameter bounds:
bounds = data.frame(log_AvgR = c(-Inf, Inf),
                    sigma_R = c(-Inf, Inf),        # Fixed.
                    init_recDevs = c(-5, 5),
                    recDevs = c(-5, 5),
                    logit_gamma_R = c(-Inf, Inf),  # Fixed.
                    log_Finit = c(-3.0,0.69),
                    log_q_rv = c(-5, 1.5),
                    log_q_prior_rv = c(-Inf, Inf), # Fixed.
                    sigma_q_rv = c(-Inf, Inf),     # Fixed.
                    log_M_juvenile = c(-3, 0.69),
                    log_M_adult = c(-3, 0.69),
                    Minit_prior = c(-Inf, Inf),    # Fixed.  
                    sigma_M = c(-Inf, Inf),        # Fixed.
                    log_S50 = c(0.6, 2.5), 
                    log_S95_delta = c(0, 2.0),
                    log_sigma_rv = c(-3.0,1.5),
                    log_q_sen = c(-5, 1.5),
                    log_q_prior_sen = c(-Inf, Inf), # Fixed.
                    sigma_q_sen = c(-Inf, Inf),     # Fixed.
                    log_sigma_sen = c(-3.0,1.5))
   
bounds <- as.data.frame(t(bounds))
names(bounds) <- c("lower", "upper")

