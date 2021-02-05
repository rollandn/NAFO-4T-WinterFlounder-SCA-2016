library(gulf)
library(TMB)
rm(list = ls())
clg()
fp <- "U:/Projects/Winter Flounder/"

source("U:/TMB/TMB Utilities.R")
setwd("U:/Winter Flounder/Population Model")
source("WFL Data Inputs Revised.R")
source("WFL Inits.R")
source("WFL Bounds.R")
source("output.sca.R")
source("fit.model.sca.R")
source("plot.sca.R")
source("C:/R packages/gulf/R/clg.R")
source("solve_Baranov.R")
     
model <- "SCA_Model_3M_C_pred"
dvars <- data.cpp(paste0(model, ".cpp"))
pvars <- parameters.cpp(paste0(model, ".cpp"))

#compile(paste0(model, ".cpp"))
dyn.load(dynlib(model))

# Data plots: 
plot.sca(data, language = "french")
 
#===============================================================================
# Model:
#   Mortality - 2 age blocks (2-4) and (5+)
#             - 9 time blocks of 5 years (last block is four years).
#   Catches - uncertainty placed on 1984.

# Catch uncertainty:
data$C_obs <- data$C
data$sigma_C <- rep(0.1, length(data$years))
data$sigma_C[data$years == 1984] <- 2   # Large error for 1984.
inits$log_C <- log(data$C_obs)  # Initial estimates set to observed catches.
bounds["log_C", ] <- c(1, 10)   # Loose bounds.

# Mortality parameters:
data$M_block_year <- data$M_block
data$M_block_age <- 1*(data$ages %in% 1:4) + 2*(data$ages %in% 5:15)
inits$Minit_prior <- c(0.6, 0.44) # Priors for initial M.
inits$sigma_M <- c(0.1, 0.1)      # Log-scale error for initial M prior.
inits$log_M <- cbind(rep(-0.7, max(data$M_block_year)), rep(-0.5, max(data$M_block_year)))
bounds["log_M", ] <- c(-4.0, log(3))

# Survey index errors:
inits$log_sigma_rv <- -0.5  
inits$log_sigma_sen <- -0.5 

# Survey catchabilities:
inits$log_q_prior_rv <- -1.7  # Log-scale prior on mean RV catchability.
inits$log_q_prior_sen <- -2   # Log-scale prior on mean Sentinel catchability.
inits$sigma_q_rv <- 0.7       # Log-scale prior on RV catchability.
inits$sigma_q_sen <- 5        # Log-scale prior on Sentinel catchability.
inits$log_q_rv <- -1          # Initial value on log-scale RV catchability.
inits$log_q_sen <- -3         # Initial value on log-scale Sentinel catchability.
inits$log_q_correction <- 0

# Plot of q_rv prior:
clg(); windows(width = 6, height = 4); 
t <- seq(0, 1, len = 500); 
ty <- dlnorm(t, inits$log_q_prior_rv, inits$sigma_q_rv)
plot(t, ty, type = "l", lwd = 2, yaxs = "i", xaxs = "i", ylim = c(0, 1.05* max(ty)), 
     xlab = "RV survey catcability", cex.lab = 1.2, ylab = "Density", cex.axis = 1.0, col = "blue") 
qlnorm(c(0.025, 0.975), inits$log_q_prior_rv, inits$sigma_q_rv)

# Selectivity parameters:
inits$log_S50 <- c(1.128241, 1.045758)
inits$log_S95_delta <- c(0.5, 0.5)
bounds["log_S95_delta", ] <- c(-1.0, 2)
bounds["log_S50", ] <- c( 0.6, 2.5)  
bounds["log_q_correction", ] <- c(-Inf, Inf)  

#bounds["log_q_rv", ] <- c(log(0.04), log(0.7))


# Extend data into future:
data$C_obs   <- c(data$C_obs, 300, 300, 300, 300, 300)
#data$sigma_C <- c(data$sigma_C, 0.1, 0.1, 0.1, 0.1, 0.1) 
data$weight_at_age_rv <- rbind(data$weight_at_age_rv,
                               data$weight_at_age_rv["2016", ],
                               data$weight_at_age_rv["2016", ],
                               data$weight_at_age_rv["2016", ],
                               data$weight_at_age_rv["2016", ],
                               data$weight_at_age_rv["2016", ])
data$weight_at_age_fishery <- rbind(data$weight_at_age_fishery,
                               data$weight_at_age_fishery["2016", ],
                               data$weight_at_age_fishery["2016", ],
                               data$weight_at_age_fishery["2016", ],
                               data$weight_at_age_fishery["2016", ],
                               data$weight_at_age_fishery["2016", ])
data$maturity_at_age <- rbind(data$maturity_at_age,
                              data$maturity_at_age["2016", ],
                              data$maturity_at_age["2016", ],
                              data$maturity_at_age["2016", ],
                              data$maturity_at_age["2016", ],
                              data$maturity_at_age["2016", ])
data$S_rv_block <- c(data$S_rv_block, 1, 1, 1, 1, 1)
data$S_fishery_block <- c(data$S_fishery_block, 2, 2, 2, 2, 2)
data$M_block_year <- c(data$M_block_year, 9, 9, 9, 9, 9)

# Inits:
#inits$log_C <- c(inits$log_C)
#inits$recDevs <- c(inits$recDevs)


# Check that data and parameters are all available:
pvars %in% names(inits)
dvars %in% names(data)


m <- fit.model.sca(data = data[dvars], inits = inits[pvars], bounds = bounds, model = model, optim = FALSE)#, additional = "log_q_correction")
rep <- sdreport(m$obj)


# Fit with RV-2S model:
initnew <- update.parameters(inits[pvars], m$fixed)
initnew$log_q_rv <- -1.6
data$S_rv_block <- 1*(1973:2021 <= 1985) + 3*(1973:2021 > 1985)
initnew$log_S50 <- c(1, 1, 1)
initnew$log_S95_delta <- c(0.5, 0.5, 0.5)


names(data$C_obs) <- 1973:2021
C <- c(1, 100, 200, 300)

# Initialize result variable:
res <- list()
for (i in 1:length(C)){
   data$C_obs[as.character(2017:2021)] <- C[i]
   res[[i]] <- fit.model.sca(data = data[dvars], inits = initnew[pvars], bounds = bounds, model = model, optim = TRUE)#, additional = "log_q_correction")
   res[[i]]$data <- data  
}
for (i in 1:length(C)){
   res[[i]]$rep <- sdreport(res[[i]]$obj)
}
vars <- c("SSB_year", "N_young", "N_old", "recruitment_rate")
m <- list(NULL, NULL, NULL, NULL)
s <- list(NULL, NULL, NULL, NULL)
for (i in 1:length(vars)){
   for (j in 1:length(res)){
      tmp <- res[[j]]$rep$value
      index <- grep(vars[i], names(res[[j]]$rep$value))
      m[[i]] <- cbind(m[[i]], tmp[index])
      tmp <- res[[j]]$rep$sd
      index <- grep(vars[i], names(res[[j]]$rep$value))
      s[[i]] <- cbind(s[[i]], tmp[index])
   }
   colnames(m[[i]]) <- c(0, 100, 200, 300)
   colnames(s[[i]]) <- c(0, 100, 200, 300)
   if (i < 4){
      rownames(m[[i]]) <- 1973:2021
      rownames(s[[i]]) <- 1973:2021
   }else{
      rownames(m[[i]]) <- 1973:2014
      rownames(s[[i]]) <- 1973:2014
   }
}

# Recruitment rate:
windows(width = 10, height = 7); 
dbarplot(m[[4]][,1], xlab = "Year", ylab = "Recruitment rate (#/kg)", cex.lab = 1.4, cex.axis = 1.25, xaxt = "n")
axis(1, at = seq(1970, 2020, by = 5) )


names(m) <- vars
names(s) <- vars
years <- 2011:2021 
cols <- c("black", "blue", "green", "red")
lty = c("solid", "dashed", "dotted", "dashed")
scale <- c(1000, 1000000, 1000000)
labels <- c("SSB (x 1000 t)", "Age 2-4 abundance (millions)", "Age 5+ abundance (millions)")
labels <- c("BSR (x 1000 t)", "Age 2-4 abundance (millions)", "Age 5+ abundance (millions)")

clg()

#*************************************************************************************************************

output <- file.path(fp, "Projection.tiff")

tiff(output, width=3500, height=1500, compression="lzw", res=300)

par(mar = c(5, 6.5, 1.5, 1.5))    # c(bottom, left, top, right)

for (i in 1:length(vars)){

   sigma <- sqrt(log(1+(s[[i]]^2/m[[i]]^2)))
   mu <- log(m[[i]] / sigma)
   q2.5 <- qlnorm(0.10, mu, sigma)
   q97.5 <- qlnorm(0.90, mu, sigma)
   #plot(range(years), c(0, 1.25*max(m[[i]][as.character(years),])) / scale[i], 
   #     type = "n", xlab = "Year", ylab = labels[i], yaxs = "i", cex.lab = 1.6, cex.axis = 1.4)
   plot(range(years), c(0, 1.25*max(m[[i]][as.character(years),])) / scale[i], 
        type = "n", xlab = "Annee", ylab = labels[i], yaxs = "i", cex.lab = 1.6, cex.axis = 1.4)
   polygon(c(years[1:6], rev(years[1:6])), c(q2.5[as.character(years[1:6]), 1], rev(q97.5[as.character(years[1:6]), 1])) / scale[i], col = "grey85", border = NA)
   polygon(c(years[6:11], rev(years[6:11])), c(q2.5[as.character(years[6:11]), 1], rev(q97.5[as.character(years[6:11]), 1])) / scale[i], col = "lightblue", border = NA)
   for (j in 1:ncol(mu)) lines(years, m[[i]][as.character(years), j] / scale[i], col = cols[j], lwd = 2, lty = lty[j])
   box()
}
legend("topright", legend = paste((0:3)*100, "t"), lwd = 2, lty = lty, col = cols, bg = "white", cex = 1.6) 

dev.off()

#*************************************************************************************************************

   windows()
   sigma <- sqrt(log(1+(s[[1]]^2/m[[1]]^2)))
   mu <- log(m[[1]] / sigma)
   q2.5 <- qlnorm(0.025, mu, sigma)
   q97.5 <- qlnorm(0.975, mu, sigma)
   plot(range(years), c(0, 1.25*max(m[[1]][as.character(years),])) / scale[1], type = "n", xlab = "Year", 
        ylab = labels[1], yaxs = "i", ylim = c(70, 80), cex.lab = 1.6, cex.axis = 1.4)
   #polygon(c(years, rev(years)), c(q2.5[as.character(years), 1], rev(q97.5[as.character(years), 1])) / scale[1], col = "grey85", border = "grey60")
   for (j in 1:ncol(mu)) lines(years, m[[1]][as.character(years), j] / scale[1], col = cols[j], lwd = 2, lty = lty[j])
   legend("topleft", legend = paste((0:3)*100, "t"), lwd = 2, lty = lty, col = cols, bg = "white", cex = 1.6) 
   box()
   
rep  <- sdreport(res[[2]]$obj)

winbugs(rep$value[grep("SSB_year", names(rep$value))])
winbugs(apply(res[[2]]$report$N[,4:11], 1, sum))
winbugs(apply(res[[2]]$report$N[,6:11], 1, sum))

   
plot(2011:2021, v[as.character(2011:2021)] / 1000, xlab = "Year", ylab = "SSB(x 1000 tonnes)", type = "l", lwd = 2, col = "black")

dbarplot(apply(res[[2]]$report$N[,4:11] / 1000000, 1, sum))

# Fit with RV-1S and Fishery-3S model:
initnew <- update.parameters(inits, res[[1]]$fixed)
initnew$log_q_rv <- -1.6
data$S_rv_block <- rep(1, data$n_year)
data$S_fishery_block <- 2*(data$years %in% 1960:1985) + 3*(data$years %in% 1986:2005) + 4*(data$years %in% 2006:2016)
initnew$log_S50 <- c(1, 1, 1, 1)
initnew$log_S95_delta <- c(0.5, 0.5, 0.5, 0.5)
res[[3]] <- fit.model.sca(data = data[dvars], inits = initnew[pvars], bounds = bounds, model = model, optim = optim)#, additional = "log_q_correction")
res[[3]]$data <- data

# Fit with RV-1S and Fishery-3S model and 3M:
initnew <- update.parameters(initnew, res[[3]]$fixed)
initnew$log_q_rv <- -1.6
data$S_rv_block <- rep(1, data$n_year)
data$S_fishery_block <- 2*(data$years %in% 1960:1985) + 3*(data$years %in% 1986:2005) + 4*(data$years %in% 2006:2016)
initnew$log_S50 <- c(1, 1, 1, 1)
initnew$log_S95_delta <- c(0.5, 0.5, 0.5, 0.5)
data$M_block_age <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3)
initnew$Minit_prior <- c(0.6, 0.52, 0.44)
initnew$sigma_M <- c(0.1, 0.1, 0.1)
initnew$log_M <- cbind(rep(-0.7, data$n_M[1]), rep(-0.6, data$n_M[1]), rep(-0.5, data$n_M[1]))
res[[4]] <- fit.model.sca(data = data[dvars], inits = initnew[pvars], bounds = bounds, model = model, optim = optim)#, additional = "log_q_correction")
res[[4]]$data <- data

# Fit with RV-2S and Fishery-3S model:
initnew <- update.parameters(inits, res[[3]]$fixed)
initnew$log_q_rv <- -1.6
data <- res[[3]]$data
data$S_rv_block <- 1*(data$years <= 1985) + 2*(data$years > 1985)
data$S_fishery_block <- 3*(data$years %in% 1960:1985) + 4*(data$years %in% 1986:2005) + 5*(data$years %in% 2006:2016)
initnew$log_S50 <- c(1, 1, 1, 1, 1)
initnew$log_S95_delta <- c(0.5, 0.5, 0.5, 0.5, 0.5)
res[[5]] <- fit.model.sca(data = data[dvars], inits = initnew[pvars], bounds = bounds, model = model, optim = optim)#, additional = "log_q_correction")
res[[5]]$data <- data

# Fit with RV-2S and Fishery-3S model with 3M:
initnew <- update.parameters(inits, res[[5]]$fixed)
initnew$log_q_rv <- -1.6
data <- res[[5]]$data
data$M_block_age <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3)
initnew$Minit_prior <- c(0.6, 0.52, 0.44)
initnew$sigma_M <- c(0.1, 0.1, 0.1)
initnew$log_M <- cbind(rep(-0.7, data$n_M[1]), rep(-0.6, data$n_M[1]), rep(-0.5, data$n_M[1]))
res[[6]] <- fit.model.sca(data = data[dvars], inits = initnew[pvars], bounds = bounds, model = model, optim = optim)#, additional = "log_q_correction")
res[[6]]$data <- data

# Fit with RV-2S and Fishery-3S model with 3M:
initnew <- update.parameters(inits, res[[2]]$fixed)
initnew$log_q_rv <- -1.6
data <- res[[2]]$data
data$M_block_age <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3)
initnew$Minit_prior <- c(0.6, 0.52, 0.44)
initnew$sigma_M <- c(0.1, 0.1, 0.1)
initnew$log_M <- cbind(rep(-0.7, data$n_M[1]), rep(-0.6, data$n_M[1]), rep(-0.5, data$n_M[1]))
res[[7]] <- fit.model.sca(data = data[dvars], inits = initnew[pvars], bounds = bounds, model = model, optim = optim)#, additional = "log_q_correction")
res[[7]]$data <- data

res <- res[c(1,2,3,5,7,4,6)]

summary <- cbind(unlist(lapply(res, function(x) length(unique(x$data$S_rv_block)))),      # RV S
                 unlist(lapply(res, function(x) length(unique(x$data$S_fishery_block)))), # Fishery S
                 unlist(lapply(res, function(x) ncol(x$inits$log_M))))                    # Mortality
colnames(summary) <- c("RV_S", "Fishery_S", "Mortality")


# Constant selectivities, q-correction:
#for (i in 8:14){
#   res[[i]] <- fit.model.sca(data = res[[i-7]]$data[dvars], inits = res[[i-7]]$inits[pvars], bounds = bounds, model = model, optim = optim, additional = "log_q_correction")
#   res[[i]]$data <- res[[i-7]]$data
#}

# Model summary results:
SS <- unlist(lapply(res, function(x) x$theta$value))   # Sum-of-squares.
cbind(SS[1:5], SS[6:10])
p  <- unlist(lapply(res, function(x) length(x$theta$par)))  # Number of parameters.
cbind(p[1:5], p[6:10])
qrv <- unlist(lapply(res, function(x) exp(x$theta$par["log_q_rv"])))  
log_qrv <- unlist(lapply(res, function(x) x$theta$par["log_q_rv"]))  
cbind(qrv[1:5], qrv[6:10])
AIC <- 2 * p + 2 * SS
cbind(AIC[1:5], AIC[6:10])
qcor <- unlist(lapply(res, function(x) exp(x$theta$par["log_q_correction"])))[6:10]  

clg(); output.sca(res[[1]]$obj, res[[1]]$data) # RV 1S, Fishery 1S, 2M   No
clg(); output.sca(res[[2]]$obj, res[[2]]$data) # RV 2S, Fishery 1S, 2M   Doug Yes
clg(); output.sca(res[[3]]$obj, res[[3]]$data) # RV 1S, Fishery 3S, 2M   No
clg(); output.sca(res[[4]]$obj, res[[4]]$data) # RV 2S, Fishery 3S, 2M   Maybe, RV selectivity  is a bit weird...
clg(); output.sca(res[[5]]$obj, res[[5]]$data) # RV 2S, Fishery 1S, 3M   Maybe, RV selecetivity inverted
clg(); output.sca(res[[6]]$obj, res[[6]]$data) # RV 1S, Fishery 3S, 3M   No
clg(); output.sca(res[[7]]$obj, res[[7]]$data) # RV 2S, Fishery 3S, 3M   Doug Residuals of, RV S different.


S <- rbind(unique(res[[2]]$report$S_rv), unique(res[[2]]$report$S_fishery))
cols <- c("blue", "green", "red")
windows(width = 7, height = 4)
plot(range(data$ages), c(0, 1), yaxs = "i", xlab = "Age", cex.lab = 1.3, ylab = "Proportion", type = "n")
for (i in 1:nrow(S)){
   lines(data$ages, S[i,], col = cols[i], lwd = 2)
}
legend("topleft", 
       c("Survey 1973-1984", "Survey 1985-2016", "Fishery"),
       lwd = 2, col = cols)
       
2, 5, 7
# Model comparison::
report <- res[[2]]$report
   # Plot of age-grouped RV index:
   dimnames(report$N_rv) <- list(year = data$years, age = data$ages)
   age.group <- list(2:4, 5:7, 8:10, 11:12)
   m <- repvec(rep(c(1, 3), each = 4), nrow = 3)
   m <- rbind(m, m + 1)
   m <- cbind(0, m, 0)
   m <- rbind(0, m, 0)
   windows(width = 11, height = 7)
   par(mar = c(0, 3, 0, 0))    # c(bottom, left, top, right)
   layout(m)   
   for (i in 1:4){
      yy <- apply(report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum)
      zz <- apply(data$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum)
      plot(zz, data$years, type = "n", ylim = 1.05 * c(0, max(c(yy, zz))),
                      xlim = c(1972.5, 2018), xaxs = "i",  
                      xaxt = "n", cex.axis = 1.3)
      points(data$years, zz, pch = 21, bg = "grey", cex = 1.5)  
      
      yy <- apply(report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum)      
      for (j in 1:7) colnames(res[[j]]$report$N_rv) <- data$ages
      lines(data$years, apply(res[[1]]$report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum), col = "purple", lwd = 3)
      lines(data$years, apply(res[[2]]$report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum), col = "blue", lwd = 3)
      lines(data$years, apply(res[[5]]$report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum), col = "green", lwd = 3)
      lines(data$years, apply(res[[7]]$report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum), col = "red", lwd = 3)
                        
      col <- c("grey", "red")
      if (i == 3){ 
         legend("topright", c("Observed", "Model 1", "Model 2", "Model 5", "Model 7"), pch = c(21, NA, NA, NA, NA), pt.bg = col, 
                col = c("black", "purple", "blue", "green", "red"), pt.cex = c(2, 0, 0, 0, 0), cex = 1.5, bg = "white", lwd = c(0, 3, 3, 3, 3))
      }
      if (i %in% c(2, 4)) axis(1, at = seq(1970, 2020, by = 5),  cex.axis = 1.3)
      if (i %in% c(1)) mtext("Abundance (millions)", line = 3, side = 2, cex = 1.3, adj = -1.8)
      if (i %in% c(2,4)) mtext("Year", line = 3, side = 1, cex = 1.3)
      text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]),
           par("usr")[3] + 0.90 * diff(par("usr")[3:4]),
           paste("Age", paste(range(age.group[[i]]), collapse = "-")), cex = 1.75)
           
   }



clg(-23)
output.sca(res[[6]]$obj, res[[6]]$data)
clg(-c(23:24))
# 5-year projections:
k <- 3

# 
theta <- res[[k]]$theta$par
# Loop over remaining years:
n_age <- data$n_age
n_year <- data$n_year
C <- rep(0, 6) # Catch for the next five years.
C[1] <- data$C_obs[n_year]
M <- res[[k]]$report$M[n_year,]
N <- matrix(NA, 6, n_age)
N[1, ] <- res[[k]]$report$N[n_year, ]
Z <- matrix(NA, 6, n_age)
Z[1, ] <- res[[k]]$report$M[n_year, ] + res[[k]]$report$S_fishery[n_year, ] * res[[k]]$report$F[n_year]
F <- rep(0, 6)
F[1] <- res[[k]]$report$F[n_year]
weight_at_age_fishery <- data$weight_at_age_fishery[n_year, ]
S_fishery <- res[[k]]$report$S_fishery[n_year,]
SSB <- rep(0, 6)
SSB[1] <- res[[k]]$report$SSB_year[n_year]
for (t in 1:5){
   # Calculate population numbers for subsequent years:
   N[t+1,1] = mean(tail(res[[k]]$report$N[,1], 5),1)              # Recruitment for first age.
   for (a in 2:(n_age-1)) N[t+1,a] = N[t,a-1] * exp(-Z[t,a-1]);               # Abundance intermediate ages.
   N[t+1,n_age] = N[t,n_age-1] * exp(-Z[t,n_age-1]) + N[t,n_age] * exp(-Z[t,n_age]); # Abundance last age. 
   Bprime = S_fishery * N[t+1,] * weight_at_age_fishery;  
   F[t+1] <- solve_Baranov(Bprime, S_fishery, C[t+1], M, 10)
   Z[t+1,] =  M + S_fishery * F[t+1]; # Update total mortality.
   SSB[t+1] = sum(data$maturity_at_age[n_year, ] * N[t+1,] * data$weight_at_age_rv[n_year,])
}
   
k <- 2

# Perform retrospective analysis:
resp <- list()
rep  <- sdreport(res[[k]]$obj)
fixed <- summary(rep, "fixed")
n_year <- data$n_year
for (i in 1:5){
   # Create subset of data:
   datasub <- res[[k]]$data
   datasub$C_obs                 <- data$C_obs[1:(n_year-i)]               
   datasub$sigma_C               <- data$sigma_C[1:(n_year-i)]
   datasub$B_rv                  <- data$B_rv[1:(n_year-i)]
   datasub$P_rv                  <- data$P_rv[1:(n_year-i), ]
   datasub$P_fishery             <- data$P_fishery[1:(n_year-i), ]
   datasub$weight_at_age_rv      <- data$weight_at_age_rv[1:(n_year-i), ]
   datasub$weight_at_age_fishery <- data$weight_at_age_fishery[1:(n_year-i), ]  
   datasub$maturity_at_age       <- data$maturity_at_age[1:(n_year-i), ]    
   datasub$S_rv_block            <- data$S_rv_block[1:(n_year-i)]
   datasub$S_fishery_block       <- data$S_fishery_block[1:(n_year-i)] 
   datasub$M_block_year          <- data$M_block_year[1:(n_year-i)]        
   datasub$n_year                <- length(datasub$C_obs)
   datasub$years                 <- data$years[1:(n_year-i)]
   
   # Modify inits:
   initsnew <- update.parameters(inits, fixed)
   initsnew$log_C   <- initsnew$log_C[1:(length(inits$log_C)-i)]
   initsnew$recDevs <- inits$recDevs[1:(length(inits$recDevs)-i)]
   initsnew$log_M   <- inits$log_M[1:max(datasub$M_block_year), ]
   
   # Fit with RV-2S and Fishery-3S model:
 #               fit.model.sca(data = data[dvars], inits = initnew[pvars], bounds = bounds, model = model, optim = optim)
   resp[[i]] <- fit.model.sca(data = datasub[dvars], inits = initsnew[pvars], bounds = bounds, model = model)
   resp[[i]]$data <- data
   #resp[[i]]$report$N[,1]
   #clg(); output.sca(res[[3]]$obj, res[[3]]$data)
   #rep  <- sdreport(res[[3]]$obj)
   #fixed <- summary(rep, "fixed")
   #exp(res[[3]]$theta$par["log_q_rv"])

   # Loop over remaining years:
   #n_age <- datasub$n_age
   #n_year <- datasub$n_year
   #C <- rep(0, 6) # Catch for the next five years.
   #C[1] <- data$C_obs[n_year]
   #M <- res[[k]]$report$M[n_year,]
   #N <- matrix(NA, 6, n_age)
   #N[1, ] <- res[[k]]$report$N[n_year, ]
   #Z <- matrix(NA, 6, n_age)
   #Z[1, ] <- res[[k]]$report$M[n_year, ] + res[[k]]$report$S_fishery[n_year, ] * res[[k]]$report$F[n_year]
   #F <- rep(0, 6)
   #F[1] <- res[[k]]$report$F[n_year]
   #weight_at_age_fishery <- data$weight_at_age_fishery[n_year, ]
   #S_fishery <- res[[k]]$report$S_fishery[n_year,]
   #SSB <- rep(0, 6)
   #SSB[1] <- res[[k]]$report$SSB_year[n_year]
   #for (t in 1:5){
   #   # Calculate population numbers for subsequent years:
   #   N[t+1,1] = sample(tail(resp[[i]]$report$N[,1], 10),1)                                                # Recruitment for first age.
   #   for (a in 2:(n_age-1)) N[t+1,a] = N[t,a-1] * exp(-Z[t,a-1]);               # Abundance intermediate ages.
   #   N[t+1,n_age] = N[t,n_age-1] * exp(-Z[t,n_age-1]) + N[t,n_age] * exp(-Z[t,n_age]); # Abundance last age. 
   #   Bprime = S_fishery * N[t+1,] * weight_at_age_fishery;  
   #   F[t+1] <- solve_Baranov(Bprime, S_fishery, C[t+1], M, 10)
   #   Z[t+1,] =  M + S_fishery * F[t+1]; # Update total mortality.
   #   SSB[t+1] = sum(data$maturity_at_age[n_year, ] * N[t+1,] * data$weight_at_age_rv[n_year,])
   ##}
}

# Retrospective plots:
# F
m <- repvec(rep(c(1, 2), each = 4), nrow = 3)
m <- rbind(m, m + 2, m+4)
m <- cbind(0, m, 0)
m <- rbind(0, m, 0)

windows(width = 9, height = 10)
par(mar = c(0, 5, 0, 0))    # c(bottom, left, top, right)
layout(m) 
   
# Abundance 2-4:
plot(res[[k]]$data$year, apply(res[[k]]$report$N[, 1:3], 1, sum) / 1000000, lwd = 2, col = "black", type = "l", 
      xlab = "Year", cex.lab = 1.5, ylab = "Abundance age 2-4 (billions)", xaxt = "n", cex.axis = 1.2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], apply(resp[[i]]$report$N[, 1:3], 1, sum) / 1000000, lwd = 2, col = cols[i], lty = lty[i])
}
legend("topright", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# Abundance 5-12:
plot(res[[k]]$data$year, apply(res[[k]]$report$N[, 4:11], 1, sum) / 1000, lwd = 2, col = "black", type = "l", 
      xlab = "Year", cex.lab = 1.5, ylab = "Abundance age 5+ (millions)", cex.axis = 1.2, xaxt = "n",)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], apply(resp[[i]]$report$N[, 4:11], 1, sum) / 1000, lwd = 2, col = cols[i], lty = lty[i])
}
legend("topright", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# SSB
plot(res[[k]]$data$year, res[[k]]$report$SSB_year / 1000, lwd = 2, col = "black", type = "l", xlab = "Year", 
     cex.lab = 1.5, ylab = "SSB(x1000 tonnes)", cex.axis = 1.2, xaxt = "n")
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$SSB_year / 1000, lwd = 2, col = cols[i], lty = lty[i])
}
legend("topright", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# F
plot(res[[k]]$data$year, res[[k]]$report$F, lwd = 2, col = "black", type = "l", xlab = "Year", cex.lab = 1.5, ylab = "F", cex.axis = 1.2, xaxt = "n",)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$F, lwd = 2, col = cols[i], lty = lty[i])
}
legend("topright", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# M young
plot(res[[k]]$data$year, res[[k]]$report$M[,1], ylim = c(0, 1), lwd = 2, col = "black", type = "l", xlab = "Year", cex.lab = 1.5, ylab = "M ages 2-4", cex.axis = 1.2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$M[,1], lwd = 2, col = cols[i], lty = lty[i])
}
legend("topleft", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# M old
plot(res[[k]]$data$year, res[[k]]$report$M[,4], ylim = c(0.5, 1.5), lwd = 2, col = "black", type = "l", xlab = "Year", cex.lab = 1.5, ylab = "M ages 5+", cex.axis = 1.2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$M[,4], lwd = 2, col = cols[i], lty = lty[i])
}
legend("topleft", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

a <- run_mcmc(res[[k]]$obj, nsim = 1000, algorithm = "RWM", alpha = 0.5)

a <- run_mcmc(res[[k]]$obj, nsim = 1000, algorithm = "HMC", L = 8, diagnostic = TRUE)

run_mcmc.nuts(100, res[[1]]$obj$fn, res[[1]]$obj$gr, max_doublings = 4, Madapt = 20, delta = 0.5,  diagnostic = TRUE)
    
# Use this one!
sim <- run_mcmc(res[[k]]$obj, nsim = 1000, algorithm = "NUTS", max_doublings = 5, Madapt = 500, delta = 0.5,  diagnostic = TRUE)
sim.hmc <- run_mcmc(res[[k]]$obj, nsim = 1000, algorithm = "HMC", L = 8, diagnostic = TRUE)
sim.rwm <- run_mcmc(res[[k]]$obj, nsim = 1000, algorithm = "RWM", alpha = 0.5)

active.pars <- c("log_AvgR", "log_q_rv", "log_sigma_rv")
active.pars <- c("log_AvgR", "log_q_rv", "log_sigma_rv", "init_recDevs", "recDevs")
 active.pars <- c("log_AvgR", "log_q_rv", "log_sigma_rv", "log_M")
 
 "log_Finit", , "log_S50"          "log_S95_delta"    "log_sigma_rv"    


rep  <- sdreport(res[[k]]$obj)
fixed <- summary(rep, "fixed")
 
nsim <- 1000
obj <- res[[k]]$obj
# Get projections:
n_age <- data$n_age
n_year <- data$n_year
weight_at_age_fishery <- data$weight_at_age_fishery[n_year, ]
SSB <- matrix(NA, nrow = nsim, ncol = 6)
colnames(SSB) <- 2016:2021
N <- array(NA, dim = c(6, n_age, nsim))
dimnames(N) <- list(year = 2016:2021, age = data$ages, sim = 1:nsim)
M <- matrix(NA, nrow = nsim, ncol = n_age)
colnames(M) <- data$ages
F <- matrix(NA, nrow = nsim, ncol = 6)
colnames(F) <- 2016:2021
model <- obj # res[[k]]$obj
p <- update.parameters(res[[2]]$inits, res[[2]]$fixed) 
index <- c(109,118,119:124)
Sigma <- solve(obj$he(res[[2]]$obj$par))
Sigma <- Sigma[index,index]
Mu <- res[[2]]$fixed[index,1]
library(MASS)
report <- obj$report(obj$par)
for (i in 1:nsim){
   if ((i %% 100)==0) print(i)
   # Define catches:
   C <- rep(250, 6) # Catch for the next five years.
   C[1] <- data$C_obs[n_year]

   # Simulate parameters:
   sim <- mvrnorm(n = 1, Mu, Sigma)
   
   # Mortality:
   M[i,] <- exp(c(rep(sim[1], 3), rep(sim[2], 8)))
   
   # Selectivity:
   S50 <- exp(sim[4]) 
   S95 <- S50 + exp(sim[7])
   tmp <- log(19) / (S95 - S50)
   S_fishery <- 1 / (1 + exp(-tmp * (data$ages - S50))) 
   
   index <- grep("N.row(43)", names(rep$value), fixed = TRUE)
   m <- rep$value[index]
   s <- rep$cov[index,index]
   N[1, ,i] <- mvrnorm(n = 1, m, s)
   
   Z <- matrix(NA, 6, n_age)
   Z[1, ] <- M[i,] + S_fishery * report$F[n_year]
   F[i,1] <- report$F[n_year]
   
   index <- grep("SSB_year", names(rep$value), fixed = TRUE)
   m <- rep$value[index]
   s <- rep$cov[index,index]
   
   #SSB[i,1] <- report$SSB_year[n_year]
   SSB[i,1] <- sum(mvrnorm(n = 1, m, s))
   
   for (t in 1:5){
      # Calculate population numbers for subsequent years:
      N[t+1,1,i] = sample(tail(report$N[,1], 10),1)                                                # Recruitment for first age.
      for (a in 2:(n_age-1)) N[t+1,a,i] = N[t,a-1,i] * exp(-Z[t,a-1]);               # Abundance intermediate ages.
      N[t+1,n_age,i] = N[t,n_age-1,i] * exp(-Z[t,n_age-1]) + N[t,n_age,i] * exp(-Z[t,n_age]); # Abundance last age. 
      Bprime = S_fishery * N[t+1,,i] * weight_at_age_fishery;  
      F[i,t+1] <- solve_Baranov(Bprime, S_fishery, C[t+1], M[i,], 10)
      Z[t+1,] =  M[i,] + S_fishery * F[i,t+1]; # Update total mortality.
      SSB[i,t+1] = sum(data$maturity_at_age[n_year, ] * N[t+1,,i] * data$weight_at_age_rv[n_year,])
   }
}









