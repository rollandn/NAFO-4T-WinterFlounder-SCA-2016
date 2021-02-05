library(gulf)
library(TMB)
rm(list = ls())
clg()

source("U:/Projects/R/TMB/TMB Utilities.R")
setwd("U:/Projects/Winter Flounder/Tobie/Population Model")
source("WFL Data Inputs Revised.R")
source("WFL Inits.R")
source("WFL Bounds.R")
source("output.sca.R")
source("fit.model.sca.R")
source("plot.sca.R")
source("C:/R packages/gulf/R/clg.R")
source("solve_Baranov.R")

fp <- "U:/Projects/Winter Flounder/CSAS 2017/Figures"

model <- "SCA_Model_3M_C"
dvars <- data.cpp(paste0(model, ".cpp"))
pvars <- parameters.cpp(paste0(model, ".cpp"))

#compile(paste0(model, ".cpp"))
dyn.load(dynlib(model))

# Data plots: 
plot.sca(data)
 
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

# Check that data and parameters are all available:
pvars %in% names(inits)
dvars %in% names(data)

# Initialize result variable:
res <- list()
optim <- TRUE

# Initial fit:
res[[1]] <- fit.model.sca(data = data[dvars], inits = inits[pvars], bounds = bounds, model = model, optim = optim)#, additional = "log_q_correction")
res[[1]]$data <- data


#res[[1]] <- fit.model.sca(data = data[dvars], inits = inits[pvars], bounds = bounds, random = "recDevs", model = model)


# Fit with RV-2S model:
initnew <- update.parameters(inits[pvars], res[[1]]$fixed)
initnew$log_q_rv <- -1.6
data$S_rv_block <- 1*(data$years <= 1985) + 3*(data$years > 1985)
initnew$log_S50 <- c(1, 1, 1)
initnew$log_S95_delta <- c(0.5, 0.5, 0.5)
res[[2]] <- fit.model.sca(data = data[dvars], inits = initnew[pvars], bounds = bounds, model = model, optim = optim)#, additional = "log_q_correction")
res[[2]]$data <- data

rep  <- sdreport(res[[2]]$obj)
rep$value
rep$sd
index <- grep("SSB_year", names(rep$value))
plot(rep$value[index]+1.96*rep$sd[index], ylim = c(-500000, 1300000))
lines(1:length(index), rep$value[index])
lines(1:length(index),rep$value[index]-1.96*rep$sd[index])

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
summary <- as.data.frame(summary)

# Constant selectivities, q-correction:
#for (i in 8:14){
#   res[[i]] <- fit.model.sca(data = res[[i-7]]$data[dvars], inits = res[[i-7]]$inits[pvars], bounds = bounds, model = model, optim = optim, additional = "log_q_correction")
#   res[[i]]$data <- res[[i-7]]$data
#}

# Model summary results:
SS <- unlist(lapply(res, function(x) x$theta$value))   # Sum-of-squares.
summary$lnL <- -SS[1:7]
p  <- unlist(lapply(res, function(x) length(x$theta$par)))  # Number of parameters.
summary$k <- p[1:7]
qrv <- unlist(lapply(res, function(x) exp(x$theta$par["log_q_rv"])))  
summary$q <- qrv

#log_qrv <- unlist(lapply(res, function(x) x$theta$par["log_q_rv"]))  
#cbind(qrv[1:5], qrv[6:10])
AIC <- 2 * p + 2 * SS
summary$AIC <- AIC
#cbind(AIC[1:5], AIC[6:10])
#qcor <- unlist(lapply(res, function(x) exp(x$theta$par["log_q_correction"])))[6:10]  

clg(); output.sca(res[[1]]$obj, res[[1]]$data) # RV 1S, Fishery 1S, 2M   No
clg(); output.sca(res[[2]]$obj, res[[2]]$data) # RV 2S, Fishery 1S, 2M   Doug Yes
clg(); output.sca(res[[3]]$obj, res[[3]]$data) # RV 1S, Fishery 3S, 2M   No
clg(); output.sca(res[[4]]$obj, res[[4]]$data) # RV 2S, Fishery 3S, 2M   Maybe, RV selectivity  is a bit weird...
clg(); output.sca(res[[5]]$obj, res[[5]]$data) # RV 2S, Fishery 1S, 3M   Maybe, RV selecetivity inverted
clg(); output.sca(res[[6]]$obj, res[[6]]$data) # RV 1S, Fishery 3S, 3M   No
clg(); output.sca(res[[7]]$obj, res[[7]]$data) # RV 2S, Fishery 3S, 3M   Doug Residuals of, RV S different.


S <- rbind(unique(res[[2]]$report$S_rv), unique(res[[2]]$report$S_fishery))
cols <- c("blue", "green", "red")


#*************************************************************************************************************

output <- file.path(fp, "Figure 39.tiff")

tiff(output, width=3500, height=1500, compression="lzw", res=300)

par(mar = c(5, 5, 1, 1))    # c(bottom, left, top, right)

plot(range(data$ages), c(0, 1), xlab = "", ylab = "", xaxt="n", yaxt="n", xaxs = "i", yaxs = "i")
axis(1, at = seq(2, 12, by = 2), cex.axis = 1.2)
axis(2, at = seq(0, 1, by = 0.2), cex.axis = 1.2, las= 1)
mtext("Age", side = 1, line = 3.5, cex = 1.6, font = 2)
mtext("Proportion", side = 2, line = 3.5, cex = 1.6, font = 2)
for (i in 1:nrow(S)){
   lines(data$ages, S[i,], col = cols[i], lwd = 2)
}
legend("topleft", 
       c("Survey 1973-1985", "Survey 1986-2016", "Fishery"),
       lwd = 2, col = cols)

dev.off()

#*************************************************************************************************************

output <- file.path(fp, "Figure 45.tiff")

tiff(output, width=3500, height=2500, compression="lzw", res=300)

#par(mar = c(5, 5, 1, 1))    # c(bottom, left, top, right)

# Model comparison::
vars <- c("N_rv_2_4", "N_rv_5_7", "N_rv_8_10", "N_rv_11_12")
  
   # Plot of age-grouped RV index:
  
   age.group <- list(2:4, 5:7, 8:10, 11:12)
   m <- repvec(rep(c(1, 3), each = 4), nrow = 3)
   m <- rbind(m, m + 1)
   m <- cbind(0, m, 0)
   m <- rbind(0, m, 0)
   #windows(width = 11, height = 7)
   par(mar = c(0, 3, 0, 0))    # c(bottom, left, top, right)
   layout(m)   
   for (i in 1:4){
      index <- names(rep$value) == vars[i]
      tmp <- data.frame(mu = rep$value[index], sigma = rep$sd[index])
      tmp$lower <- tmp$mu - 1.96 * tmp$sigma
      tmp$upper <- tmp$mu + 1.96 * tmp$sigma 

      yy <- tmp$mu / 1000
      zz <- apply(data$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum)
      plot(zz, data$years, type = "n", ylim = 1.05 * c(0, max(c(yy, zz))),
                      xlim = c(1972.5, 2018), xaxs = "i",  
                      xaxt = "n", cex.axis = 1.3, las=1)
      
      polygon(c(data$years, rev(data$years)), c(tmp$lower, rev(tmp$upper)) / 1000, col = "grey90", border = "grey70")
      points(data$years, zz, pch = 21, bg = "grey", cex = 1.5)  
      lines(data$years, yy,  col = "black", lwd = 3)

      
      #yy <- apply(report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum)      
     # for (j in 2) colnames(res[[j]]$report$N_rv) <- data$ages
      #lines(data$years, apply(res[[1]]$report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum), col = "purple", lwd = 3)

      #lines(data$years, apply(res[[5]]$report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum), col = "green", lwd = 3)
      #lines(data$years, apply(res[[7]]$report$N_rv[,as.character(age.group[[i]])] / 1000, 1, sum), col = "red", lwd = 3)
                        
      #col <- c("grey", "red")
      #if (i == 3){ 
      #   legend("topright", c("Observed", "Model 1", "Model 2", "Model 5", "Model 7"), pch = c(21, NA, NA, NA, NA), pt.bg = col, 
      #          col = c("black", "purple", "blue", "green", "red"), pt.cex = c(2, 0, 0, 0, 0), cex = 1.5, bg = "white", lwd = c(0, 3, 3, 3, 3))
      #}
      if (i %in% c(2, 4)) axis(1, at = seq(1970, 2020, by = 5),  cex.axis = 1.3)
      if (i %in% c(1)) mtext("Abundance (millions)", line = 3.5, side = 2, cex = 1.3, adj = -1.2, font=2)
      if (i %in% c(2,4)) mtext("Year", line = 3, side = 1, cex = 1.3, font=2)
      text(par("usr")[1] + 0.85 * diff(par("usr")[1:2]),
           par("usr")[3] + 0.90 * diff(par("usr")[3:4]),
           paste("Age", paste(range(age.group[[i]]), collapse = "-")), cex = 1.75)
           
   }

dev.off()

#*************************************************************************************************************
   
# SSB Figures:
index <- names(rep$value) == "log_SSB_year"
tmp <- data.frame(mu = rep$value[index], sigma = rep$sd[index])
tmp$lower <- tmp$mu - 1.0 * tmp$sigma
tmp$upper <- tmp$mu + 1.0 * tmp$sigma 
tmp <- exp(tmp)
   windows(width = 10, height = 7)
   par(mar = c(5, 6, 4, 2))
plot(data$years, tmp$mu / 1000, xlab = "Year", ylab = "", cex.lab = 1.5, cex.axis = 1.3, type = "l", 
        lwd = 2, ylim = c(0, 600), xaxt = "n", las = 1, yaxs = "i")
polygon(c(data$years, rev(data$years)), c(tmp$lower, rev(tmp$upper)) / 1000, col = "grey90", border = "grey75")
grid()
lines(data$years, tmp$mu / 1000, lwd = 2, col = "black")  
box()
axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.4)
mtext("SSB (x 1000 t)", side = 2, line = 4, cex = 1.5)

age.group <- list(2:4, 5:7, 8:10, 11:12)
SSBg <- NULL
colnames(res[[2]]$report$SSB) <- data$ages
for (i in 1:length(age.group)) SSBg <- cbind(SSBg, apply(res[[2]]$report$SSB[,as.character(age.group[[i]])], 1, sum)) 
SSBg <- cbind(0, t(apply(SSBg, 1, cumsum)))
cols <- c("skyblue", "green", "yellow", "red")
windows(width = 10)
plot(data$years, tmp$mu / 1000, xlab = "Year", ylab = "", cex.lab = 1.5, cex.axis = 1.3, type = "n", 
     lwd = 2, ylim = c(0, 600), xaxt = "n", las = 1, yaxs = "i")
for (i in 1:length(age.group)){
   polygon(c(data$years, rev(data$years)), c(SSBg[,i], rev(SSBg[,i+1])) / 1000, col = cols[i], border = cols[i])
}
axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.4)
legend("topright", rev(c("Age 2-4", "Age 5-7", "Age 8-10", "Age 11+")), pch = 22, pt.bg = rev(cols), pt.cex = 4, cex = 1.5)

# Bmsy proxy figure (not finished!):
clg()
index <- names(rep$value) == "log_SSB_year"
tmp <- data.frame(mu = rep$value[index], sigma = rep$sd[index])
tmp$lower <- tmp$mu - 1.0 * tmp$sigma
tmp$upper <- tmp$mu + 1.0 * tmp$sigma 
tmp <- exp(tmp)
   windows(width = 10, height = 7)
   par(mar = c(5, 6, 4, 2))
plot(data$years, tmp$mu / 1000, xlab = "Year", ylab = "", cex.lab = 1.5, cex.axis = 1.3, type = "l", 
        lwd = 2, ylim = c(0, 600), xaxt = "n", las = 1, yaxs = "i")
polygon(c(data$years, rev(data$years)), c(tmp$lower, rev(tmp$upper)) / 1000, col = "grey90", border = "grey75")
grid()
lines(data$years, tmp$mu / 1000, lwd = 2, col = "black")  
box()
axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.4)
mtext("SSB (x 1000 t)", side = 2, line = 4, cex = 1.5)

msy.range <- c(1973, 1994)
x <- data$years
y <- tmp$mu
Bmsy <- mean(y[x %in% msy.range[1]:msy.range[2]]) / 1000
arrows(msy.range[1], Bmsy, msy.range[2], Bmsy, length = 0.15, lwd = 2, col = "red")
arrows(msy.range[2], Bmsy, msy.range[1], Bmsy, length = 0.15, lwd = 2, col = "red")
lines(par("usr")[1:2], c(0.8 * Bmsy, 0.8 * Bmsy), lty = "dashed", lwd = 2, col = "red")
lines(par("usr")[1:2], c(0.4 * Bmsy, 0.4 * Bmsy), lty = "dashed", lwd = 2, col = "red")
text(msy.range[2], Bmsy, paste0("Bmsy"), pos = 4)
text(2009, 0.8 * Bmsy, paste0("Busr"), pos = 3)
text(2010, 0.4 * Bmsy, paste0("Blim"), pos = 3)
#legend("topright", legend = c(">= 25 cm", "Smoothed 3yr avg"), col = c("black", "green"), lwd = 2, bg = "white", cex = 1.3)
box()


# F-plot:
rep <- sdreport(res[[2]]$obj)
vars <- c("log_F2_4", "log_F5_12")
xlim <- c(15e-5, 5e-3)
for (i in 1:2){
   index <- names(rep$value) == vars[i]
   tmp <- data.frame(mu = rep$value[index], sigma = rep$sd[index])
   tmp$lower <- tmp$mu - 1.0 * tmp$sigma
   tmp$upper <- tmp$mu + 1.0 * tmp$sigma 
   tmp <- exp(tmp)
   windows(width = 9, height = 7)
   par(mar = c(5, 8, 4, 2))
   plot(data$years, tmp$mu, xlab = "Year", ylab = "", cex.lab = 1.5, cex.axis = 1.3, type = "l", 
        lwd = 2, ylim = c(0, xlim[i]), xaxt = "n", las = 1)
   polygon(c(data$years, rev(data$years)), c(tmp$lower, rev(tmp$upper)), col = "grey90", border = "grey75")
   grid()
   lines(data$years, tmp$mu, lwd = 2, col = "black")  
   box()
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.4)
  # mtext(expression(paste("F (x ", 10^-3, ")", sep = "")), side = 2, line = 3.5, cex = 1.5)
   mtext("F", side = 2, line = 5, cex = 1.5) 
   if (i == 1) legend("topright", "Age 2-4", cex = 1.6, bg = "white") 
   if (i == 2) legend("topright", "Age 5+", cex = 1.6, bg = "white") 
}

# Mortality plot:
M <- exp(res[[2]]$parameters$log_M)
rownames(M) <- paste(seq(1973, 2016, by = 5), "-",  substr(c(seq(1977, 2012, by = 5), 2016), 3, 4), sep = "")
index <- names(rep$value) == "log_M"
tmp <- data.frame(mu = rep$value[index], sigma = rep$sd[index])
tmp$lower <- tmp$mu - 1.96 * tmp$sigma
tmp$upper <- tmp$mu + 1.96 * tmp$sigma 
tmp <- exp(tmp)

windows(width = 12, height = 7)
plot(1:9, tmp$mu[1:9], xlab = "Year", cex.lab = 1.6, cex.axis = 1.3, pch = 21, bg = "grey", cex = 1.8, 
     lwd = 2, ylab = "M", ylim = c(0, 1.5), xaxt = "n", las = 1, yaxs = "i")
grid()
for (i in 1:9){
   lines(c(i,i), c(tmp$lower[i], tmp$upper[i]), lwd = 2)
   lines(c(i-0.12,i+0.12), c(tmp$lower[i], tmp$lower[i]), lwd = 2)
   lines(c(i-0.12,i+0.12), c(tmp$upper[i], tmp$upper[i]), lwd = 2)
}     
points(1:9, tmp$mu[1:9], pch = 21, bg = "grey", cex = 2.5, lwd = 2)
box()
axis(1, at = 1:9, labels = rownames(M), cex.axis = 1.4)   
axis(1, at = 2*(1:4), labels = rownames(M)[2*(1:4)], cex.axis = 1.4)   
legend("topleft", "Age 2-4", cex = 1.6, bg = "white") 
      
windows(width = 12, height = 7)
plot(1:9, tmp$mu[9+(1:9)], xlab = "Year", cex.lab = 1.6, cex.axis = 1.3, pch = 21, bg = "grey", cex = 1.8, 
     lwd = 2, ylab = "M", ylim = c(0, 1.5), xaxt = "n", las = 1, yaxs = "i")
grid()
for (i in 1:9){
   lines(c(i,i), c(tmp$lower[9+i], tmp$upper[9+i]), lwd = 2)
   lines(c(i-0.12,i+0.12), c(tmp$lower[9+i], tmp$lower[9+i]), lwd = 2)
   lines(c(i-0.12,i+0.12), c(tmp$upper[9+i], tmp$upper[9+i]), lwd = 2)
}     
points(1:9, tmp$mu[10:18], pch = 21, bg = "grey", cex = 2.5, lwd = 2)
box()
axis(1, at = 1:9, labels = rownames(M), cex.axis = 1.4)   
axis(1, at = 2*(1:4), labels = rownames(M)[2*(1:4)], cex.axis = 1.4)  
legend("topleft", "Age 5+", cex = 1.6, bg = "white") 


# Joint mortality figure:
windows(width = 12, height = 7)
plot(1:9, tmp$mu[1:9], xlab = "Year", cex.lab = 1.6, cex.axis = 1.3, pch = 21, bg = "grey", cex = 1.8, 
     lwd = 2, ylab = "Natural mortality (M)", ylim = c(0, 1.5), xaxt = "n", las = 1, yaxs = "i", type = "n")
grid()
for (i in 1:9){
   lines(c(i,i)-0.12, c(tmp$lower[i], tmp$upper[i]), lwd = 2)
   lines(c(i-0.12,i+0.12)-0.12, c(tmp$lower[i], tmp$lower[i]), lwd = 2)
   lines(c(i-0.12,i+0.12)-0.12, c(tmp$upper[i], tmp$upper[i]), lwd = 2)
}     
points((1:9)-0.12, tmp$mu[1:9], pch = 21, bg = "grey", cex = 2.5, lwd = 2)
for (i in 1:9){
   lines(c(i,i)+0.12, c(tmp$lower[9+i], tmp$upper[9+i]), lwd = 2)
   lines(c(i-0.12,i+0.12)+0.12, c(tmp$lower[9+i], tmp$lower[9+i]), lwd = 2)
   lines(c(i-0.12,i+0.12)+0.12, c(tmp$upper[9+i], tmp$upper[9+i]), lwd = 2)
}     
points((1:9)+0.12, tmp$mu[10:18], pch = 24, bg = "red", cex = 2.5, lwd = 2)
box()
axis(1, at = 1:9, labels = rownames(M), cex.axis = 1.4)   
axis(1, at = 2*(1:4), labels = rownames(M)[2*(1:4)], cex.axis = 1.4)   
legend("topleft", legend = c("Age 2-4", "Age 5+"), pch = c(21, 24), pt.bg = c("grey", "red"), pt.cex = 3.0, cex = 2, bg = "white")


# Abundance by group plot:
m <- repvec(rep(c(1, 3), each = 4), nrow = 3)
m <- rbind(m, m + 1)
m <- cbind(0, m, 0)
m <- rbind(0, m, 0)
windows(width = 11, height = 7)
par(mar = c(0, 3, 0, 0))    # c(bottom, left, top, right)
layout(m)   
vars <- c("log_N2_4", "log_N5_7", "log_N8_10", "log_N11_12")
xlim <- 1.5*c(6000, 1400, 120, 11)
age.group <- list(2:4, 5:7, 8:10, 11:12)
labels <- c("Age 2-4", "Age 5-7", "Age 8-10", "Age 11+")
for (i in 1:4){
   index <- names(rep$value) == vars[i]
   tmp <- data.frame(mu = rep$value[index], sigma = rep$sd[index])
   tmp$lower <- tmp$mu - 1.0 * tmp$sigma
   tmp$upper <- tmp$mu + 1.0 * tmp$sigma 
   tmp <- exp(tmp)
   
  # windows(width = 9, height = 7)
   #par(mar = c(5, 6, 4, 2))
   plot(data$years, tmp$mu / 1000, xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.3, type = "l", 
        lwd = 2, ylim = c(0, xlim[i]), xaxt = "n", las = 1, yaxs = "i")
   polygon(c(data$years, rev(data$years)), c(tmp$lower, rev(tmp$upper)) / 1000, col = "grey90", border = "grey75")
   grid()
   lines(data$years, tmp$mu / 1000, lwd = 2, col = "black")  
   box()
   #axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.4)
   # mtext(expression(paste("F (x ", 10^-3, ")", sep = "")), side = 2, line = 3.5, cex = 1.5)

   S <- res[[2]]$report$S_rv
   colnames(S) <- data$ages
   
   #N <- (1 / exp(res[[2]]$theta$par["log_q_rv"])) * (1 / S[,as.character(age.group[[i]])]) * data$N_rv[, as.character(age.group[[i]])]  
   #N <- apply(N, 1, sum)
   #points(data$years, N / 1000, pch = 21, bg = "grey", cex = 1.5)
   
   legend("topleft", labels[i],  cex = 1.5, bg = "white")
   if (i %in% c(2, 4)) axis(1, at = seq(1970, 2020, by = 5),  cex.axis = 1.3)
   if (i %in% c(1)) mtext("Abundance (millions)", line = 4, side = 2, cex = 1.3, adj = -1.8)
   if (i %in% c(2,4)) mtext("Year", line = 3, side = 1, cex = 1.3)
}

# Trawlable biomass of 25+cm fish from RV:
TB <- c(11716.7,14407.3,9903.3,47363.9,19920.3,63442.5,5687.4,22759.7,12450.6,24859.9,24587.9,17639.1,
        25189,7043.7,8299,12959.7,4397.9,8735.6,10195.6,8439.5,5908.9,8892,9871.4,4968.8,6234.3,5952.1,
        4618.6,3068.2,7123.6,6188.2,4614,9947.8,4371.5,4990,4437.6,4035.5,4082.9,2639.7,3845.6,2597,
        1712.6,1030.5,2346.1,1400.5,1860.7,743.8)
names(TB) <- 1971:2016
SSB <- res[[2]]$report$SSB_year
names(SSB) <- data$years
mean(TB[as.character(1973:1994)])
mean(SSB[as.character(1973:1994)])

TB[as.character(1973:1994)]
   
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

#*************************************************************************************************************

output <- file.path(fp, "Figure 51.tiff")

tiff(output, width=3000, height=3500, compression="lzw", res=300)

# Retrospective plots:
# F
m <- repvec(rep(c(1, 2), each = 4), nrow = 3)
m <- rbind(m, m + 2, m+4)
m <- cbind(0, m, 0)
m <- rbind(0, m, 0)

#windows(width = 9, height = 10)
par(mar = c(0, 7, 0, 0))    # c(bottom, left, top, right)
layout(m) 
   
# Abundance 2-4:
plot(res[[k]]$data$year, apply(res[[k]]$report$N[, 1:3], 1, sum) / 1000000, lwd = 2, col = "black", type = "l", 
     xlab = "", ylab = "", cex.lab = 1.5, xaxt = "n", cex.axis = 1.2, las=1)
mtext("Abundance age 2-4 (billions)", side = 2, line = 3.5, cex = 1.2, font = 2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], apply(resp[[i]]$report$N[, 1:3], 1, sum) / 1000000, lwd = 2, col = cols[i], lty = lty[i])
}
legend("topright", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# Abundance 5-12:
plot(res[[k]]$data$year, apply(res[[k]]$report$N[, 4:11], 1, sum) / 1000, lwd = 2, col = "black", type = "l", 
     xlab = "", ylab = "", cex.lab = 1.5, cex.axis = 1.2, xaxt = "n", las=1)
mtext("Abundance age 5+ (millions)", side = 2, line = 3.5, cex = 1.2, font = 2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], apply(resp[[i]]$report$N[, 4:11], 1, sum) / 1000, lwd = 2, col = cols[i], lty = lty[i])
}
legend("topright", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# SSB
plot(res[[k]]$data$year, res[[k]]$report$SSB_year / 1000, lwd = 2, col = "black", type = "l", ylim=c(100,550), xlab = "", ylab = "",
     cex.axis = 1.2, xaxt = "n", las=1)
mtext("SSB(x1000 tonnes)", side = 2, line = 3.5, cex = 1.2, font = 2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$SSB_year / 1000, lwd = 2, col = cols[i], lty = lty[i])
}
legend("topright", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# F
plot(res[[k]]$data$year, res[[k]]$report$F, lwd = 2, col = "black", type = "l", ylim=c(0,0.054), xlab = "", ylab = "",
     cex.axis = 1.2, xaxt = "n", las=1)
mtext("F", side = 2, line = 3.5, cex = 1.2, font = 2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$F, lwd = 2, col = cols[i], lty = lty[i])
}
legend("topright", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# M young
plot(res[[k]]$data$year, res[[k]]$report$M[,1], ylim = c(0, 1.1), lwd = 2, col = "black", type = "l",
     xlab = "", ylab = "", cex.axis = 1.2, las=1)
mtext("Year", side = 1, line = 3.5, cex = 1.2, font = 2)
mtext("M ages 2-4", side = 2, line = 3.5, cex = 1.2, font = 2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$M[,1], lwd = 2, col = cols[i], lty = lty[i])
}
legend("topleft", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

# M old
plot(res[[k]]$data$year, res[[k]]$report$M[,4], ylim = c(0.5, 1.5), lwd = 2, col = "black", type = "l", 
     xlab = "", ylab = "", cex.axis = 1.2, las=1)
mtext("Year", side = 1, line = 3.5, cex = 1.2, font = 2)
mtext("M ages 5+", side = 2, line = 3.5, cex = 1.2, font = 2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("solid", "solid", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$M[,4], lwd = 2, col = cols[i], lty = lty[i])
}
legend("topleft", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = lty, bg = "white", cex = 1.15)

dev.off()


output <- file.path(fp, "WF Figure 51 jenni 1.tiff")

tiff(output, width=3000, height=3500, compression="lzw", res=300)

par(mar = c(5,5,1,1))
    
# M young
plot(res[[k]]$data$year, res[[k]]$report$M[,1], ylim = c(0, 1.1), lwd = 2, col = "black", type = "l", lty = 2,
     xlab = "", ylab = "", cex.axis = 1.2, las=1)
mtext("Year", side = 1, line = 3.5, cex = 1.2, font = 2)
mtext("M ages 2-4", side = 2, line = 3.5, cex = 1.2, font = 2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("dashed", "dashed", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$M[,1], lwd = 2, col = cols[i], lty = lty[i])
}
legend("topleft", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = c("dashed", lty), bg = "white", cex = 1.15)

dev.off()

output <- file.path(fp, "WF Figure 51 jenni 2.tiff")

tiff(output, width=3000, height=3500, compression="lzw", res=300)

par(mar = c(5,5,1,1))

# M old
plot(res[[k]]$data$year, res[[k]]$report$M[,4], ylim = c(0.5, 1.5), lwd = 2, col = "black", type = "l", lty = 2,
     xlab = "", ylab = "", cex.axis = 1.2, las=1)
mtext("Year", side = 1, line = 3.5, cex = 1.2, font = 2)
mtext("M ages 5+", side = 2, line = 3.5, cex = 1.2, font = 2)
cols <- c("blue", "green", "yellow", "orange", "red")
lty <- c("dashed", "dashed", "dashed", "dashed", "dashed")
for (i in 1:length(resp)){
   lines(resp[[i]]$data$year[1:(n_year-i)], resp[[i]]$report$M[,4], lwd = 2, col = cols[i], lty = lty[i])
}
legend("topleft", legend = 2016:2011, lwd = 2, col = c("black", cols), lty = c("dashed", lty), bg = "white", cex = 1.15)

dev.off()

#*************************************************************************************************************




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









