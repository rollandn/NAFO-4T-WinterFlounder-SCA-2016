fit.model.sca <- function(data, inits, bounds, model, random = NULL, additional, optim = FALSE){
  
   # Phase 1 fitting:
   active.pars <- c("log_AvgR", "log_q_rv", "log_sigma_rv")
   obj <- MakeADFun(data = data, parameters = inits, map = control.parameters(inits, active = active.pars), random = random, DLL = model)
   obj$fn(obj$par)
   theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"])
   obj$par <- theta$par
   theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"])
   obj$par <- theta$par
   rep  <- sdreport(obj)
   fixed <- summary(rep, "fixed")
   parameters <- update.parameters(inits, fixed)

   #if ("log_q_correction" %in% names(inits)){
   #   active.pars <- c(active.pars, "log_q_correction")
   #   obj <- MakeADFun(data = data, parameters = inits, map = control.parameters(inits, active = active.pars), random = random, DLL = model)
   #   obj$fn(obj$par)
   #   theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"])
   #   obj$par <- theta$par
   #   theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"])
   #   obj$par <- theta$par
   #   rep  <- sdreport(obj)
   #   fixed <- summary(rep, "fixed")
   #   parameters <- update.parameters(inits, fixed)
   #}
   
   # Phase 1 fitting:
   if (all(c("log_q_sen", "log_sigma_sen") %in% names(inits))){
      active.pars <- c(active.pars, "log_q_sen", "log_sigma_sen")  
      obj <- MakeADFun(data = data, parameters = inits, map = control.parameters(inits, active = active.pars), random = random, DLL = model)
      obj$fn(obj$par)
      theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"])
      obj$par <- theta$par
      rep  <- sdreport(obj)
      fixed <- summary(rep, "fixed")
      parameters <- update.parameters(inits, fixed)
   }
   
   # Fit recruitment deviates:
   active.pars <- c(active.pars, "init_recDevs", "recDevs")
   obj <- MakeADFun(data = data, parameters = parameters, map = control.parameters(parameters, active = active.pars), random = random, DLL = model)
   obj$fn(obj$par)
   theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"], eval.max = 1000, iter.max = 1000)
   obj$par <- theta$par
   theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"], eval.max = 1000, iter.max = 1000)
   obj$par <- theta$par
   rep  <- sdreport(obj)
   fixed <- summary(rep, "fixed")
   parameters <- update.parameters(inits, fixed)
   
   # Phase 2:
   active.pars <- c(active.pars, "log_S50", "log_S95_delta")
   obj <- MakeADFun(data = data, parameters = parameters, map = control.parameters(parameters, active = active.pars), random = random, DLL = model)
   obj$fn(obj$par)
   theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"], eval.max = 1000, iter.max = 1000)
   obj$par <- theta$par
   rep  <- sdreport(obj)
   fixed <- summary(rep, "fixed")
   parameters <- update.parameters(inits, fixed)
   
   # Phase 4:
   save.pars <- active.pars 
   if ("log_M" %in% names(inits)) active.pars <- "log_M" else active.pars <- c("log_M_juvenile", "log_M_adult") 
   obj <- MakeADFun(data = data, parameters = parameters, map = control.parameters(parameters, active = active.pars), random = random, DLL = model)
   obj$fn(obj$par)
   theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"])
   obj$par <- theta$par
   rep  <- sdreport(obj)
   fixed <- summary(rep, "fixed")
   parameters <- update.parameters(inits, fixed)
   active.pars <- unique(c(save.pars, active.pars))
   
   # Additional parameters to be fit:
   if (!missing(additional)){
      additional <- additional[additional %in% names(inits)]
      if (length(additional) > 0){
         save.pars <- active.pars 
         active.pars <- additional
         print(additional)
         obj <- MakeADFun(data = data[dvars], parameters = parameters, map = control.parameters(parameters, active = active.pars), random = random, DLL = model)
         obj$fn(obj$par)
         theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"])
         obj$par <- theta$par
         rep  <- sdreport(obj)
         fixed <- summary(rep, "fixed")
         parameters <- update.parameters(inits, fixed)   
         active.pars <- unique(c(save.pars, active.pars))
      }
   }
   
   # Phase 5 - final fit:
   if ("log_M" %in% names(inits)) active.pars <- active.pars <- c(active.pars, "log_M") else active.pars <- c(active.pars, "log_M_juvenile", "log_M_adult") 
   if ("log_C" %in% names(inits)) active.pars <- active.pars <- c(active.pars, "log_C") 
   obj <- MakeADFun(data = data[dvars], parameters = parameters, map = control.parameters(parameters, active = active.pars), random = random, DLL = model)
   obj$par <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"], eval.max = 1000, iter.max = 1000)$par
   for (i in 1:20){
      theta <- nlminb(obj$par, obj$fn, obj$gr, lower = bounds[names(obj$par),"lower"], upper = bounds[names(obj$par),"upper"], eval.max = 1000, iter.max = 1000)
      obj$par <- theta$par
   }
   if (optim){
       for (i in 1:4){
          theta <- optim(obj$par, obj$fn, obj$gr, control = list(trace = 3, maxit = 5000))
          obj$par <- theta$par
       }
   }
   rep  <- sdreport(obj)
   fixed <- summary(rep, "fixed")
   parameters <- update.parameters(inits, fixed)
   
   res <- list()
   res$obj <- obj
   res$fixed <- fixed
   res$parameters <- parameters
   res$theta <- theta
   res$report <- obj$report(obj$par)
   res$inits <- inits
     
   return(res)
}   