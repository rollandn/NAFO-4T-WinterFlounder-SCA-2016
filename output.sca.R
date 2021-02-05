output.sca <- function(obj, data){
   # OUTPUT.SCA - Plot SCA output.
   
   source("C:/R packages/gulf/R/dbarplot.R")
             
   report <- obj$report(obj$par)
   rep  <- sdreport(obj)
   fixed <- summary(rep, "fixed")
   
   # Survey index:
   windows(width = 10)
   plot(range(data$years), c(0, 80000 / 1000), type = "n", 
        xlab = "Year", ylab = "Trawlable biomass (x 1000 MT)", xaxt = "n", yaxs = "i", cex.lab = 1.3, cex.axis = 1.2)
   grid()
   cols <- c("blue", "green")
   lines(data$years, report$B_rv_pred_year / 1000, lwd = 2, col = "red", lty = "dashed")
   lines(data$years, data$B_rv / 1000, lwd = 2, col = "black", lty = "solid")
   box()
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
   legend("topright", legend = c("Observed", "Fitted"), lwd = 2, col = c("black", "red"), lty = c("solid", "dashed"), bg = "white", cex = 1.5)
   
   # Sentinel survey index:
   if ("B_sen_pred_year" %in% names(report)){
      windows(width = 10)
      plot(range(as.numeric(rownames(data$N_sen))), c(0, 5000 / 1000), type = "n", 
           xlab = "Year", ylab = "Trawlable biomass (x 1000 MT)", xaxt = "n", yaxs = "i", cex.lab = 1.3, cex.axis = 1.2)
      grid()
      cols <- c("blue", "green")
      lines(as.numeric(rownames(data$N_sen)), report$B_sen_pred_year / 1000, lwd = 2, col = "red", lty = "dashed")
      lines(as.numeric(rownames(data$N_sen)), data$B_sen / 1000, lwd = 2, col = "black", lty = "solid")
      box()
      axis(1, at = seq(2003, 2016, by = 2), cex.axis = 1.2)
      legend("topright", legend = c("Observed", "Fitted"), lwd = 2, col = c("black", "red"), lty = c("solid", "dashed"), bg = "white", cex = 1.5)
   }
   
   # RV survey age catch proportions:
   #windows(width = 10)
   #plot(range(data$years), c(-8, log(max(max(data$P_rv), max(report$P_rv_pred)))), type = "n", 
   #     xlab = "Year", ylab = "log(proportion)", xaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
   #cols <- rainbow(length(data$ages))
   #for (i in 1:length(cols)){
   #   lines(data$years, log(data$P_rv[, i]), lwd = 2, col = cols[i])
   #   lines(data$years, log(report$P_rv_pred[, i]), lwd = 2, col = cols[i], lty = "dashed")
   #}
   #axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
   #title(main = "RV proportions")
   
   # RV survey age catch proportions bubble plot:
   windows(width = 11, height = 5)
   plot(range(data$years), range(data$ages), type = "n", xlab = "Year", ylab = "Age", xaxt = "n", yaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
   scale <- 4
   for (i in 1:nrow(data$P_rv)){
      residuals <- log(data$P_rv[i, ]) - log(report$P_rv_pred[i, ])
      index <- residuals >= 0
      points(rep(data$years[i], length(residuals))[index], data$ages[index], cex = scale * sqrt(residuals[index]), pch = 21, bg = "grey90")
      points(rep(data$years[i], length(residuals))[!index], data$ages[!index], cex = scale * sqrt(-residuals[!index]), pch = 21, bg = "black")
   }
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
   axis(2, cex.axis = 1.2)
   # title(main = "Log-scale RV proportion residuals")
   
   # Sentinel survey age catch proportions:
   if ("P_sen_pred" %in% names(report)){
      windows(width = 8)
      years <- as.numeric(rownames(data$P_sen))
      plot(range(years), c(-8, log(max(max(data$P_sen), max(report$P_sen_pred)))), 
           type = "n", xaxt = "n", xlab = "Year", ylab = "log(proportion)", cex.lab = 1.3, cex.axis = 1.2)
      cols <- rainbow(length(data$ages)) 
      for (i in 1:length(cols)){
         lines(years, log(data$P_sen[, i]), lwd = 2, col = cols[i])
         lines(years, log(report$P_sen_pred[, i]), lwd = 2, col = cols[i], lty = "dashed")
      }
      axis(1, at = seq(2003, 2016, by = 2), cex.axis = 1.2)
   
      # Sentinel survey age catch proportions bubble plot:
      windows(width = 5)
      years <- as.numeric(rownames(data$P_sen))
      plot(range(years), range(data$ages), type = "n", xlab = "Year", ylab = "Age", xaxt = "n", yaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
      scale <- 4
      for (i in 1:nrow(data$P_sen)){
         residuals <- log(data$P_sen[i, ]) - log(report$P_sen_pred[i, ])
         index <- residuals >= 0
         points(rep(years[i], length(residuals))[index], data$ages[index], cex = scale * sqrt(residuals[index]), pch = 21, bg = "grey90")
         points(rep(years[i], length(residuals))[!index], data$ages[!index], cex = scale * sqrt(-residuals[!index]), pch = 21, bg = "black")
      }
      axis(1, at = seq(2003, 2016, by = 5), cex.axis = 1.2)
      axis(2, cex.axis = 1.2)
      title(main = "Log-scale Sentinel survey proportion residuals")
   }
   
   # Fishery age catch proportions:
   #windows(width = 10)
   #plot(range(data$years), c(-9, log(max(max(data$P_fishery), max(report$P_fishery_pred)))), type = "n", 
   #     xaxt = "n", xlab = "Year", ylab = "log(proportion)", cex.lab = 1.3, cex.axis = 1.2)
   #cols <- rainbow(length(data$ages))
   #for (i in 1:length(cols)){
   #   lines(data$years, log(data$P_fishery[, i]), lwd = 2, col = cols[i])
   #   lines(data$years, log(report$P_fishery_pred[, i]), lwd = 2, col = cols[i], lty = "dashed")
   #}
   #axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
    
   # Fishery age catch proportions bubble plot:
   windows(width = 11, height = 5)
   plot(range(data$years), range(data$ages), type = "n", xlab = "Year", ylab = "Age", xaxt = "n", yaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
   scale <- 4
   for (i in 1:nrow(data$P_fishery)){
      residuals <- log(data$P_fishery[i, ]) - log(report$P_fishery_pred[i, ])
      index <- residuals >= 0
      points(rep(data$years[i], length(residuals))[index], data$ages[index], cex = scale * sqrt(residuals[index]), pch = 21, bg = "grey90")
      points(rep(data$years[i], length(residuals))[!index], data$ages[!index], cex = scale * sqrt(-residuals[!index]), pch = 21, bg = "black")
   }
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
   axis(2, cex.axis = 1.2)
   # title(main = "Log-scale fishery proportion residuals")
   
   # SSB_year:
   #windows(width = 9, height = 5)
   ##plot(data$year, report$SSB_year, xlab = "Year", ylab = "Spawning Biomass (MT)")
   #dbarplot(report$SSB_year / 1000, data$year, xlab = "Year", ylab = "Spawning Biomass (x 1000 MT)",
   #        xaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
   #axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
   
   # SSB_year:
   windows(width = 9, height = 5)
   #plot(data$year, report$SSB_year, xlab = "Year", ylab = "Spawning Biomass (MT)")
   x <- cbind(apply(report$SSB[, data$ages %in% 1:5], 1, sum), apply(report$SSB[, data$ages %in% 6:15], 1, sum)) / 1000
   dbarplot(x, data$year, xlab = "Year", ylab = "Spawning Biomass (x 1000 MT)", 
            xaxt = "n", cex.lab = 1.3, cex.axis = 1.2, width = 1, col = c("grey", "white"), legend = FALSE)
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
   legend("topright", legend = c("Ages 2-5", "Ages 6-12+"), bg = "white", pch = 22, pt.cex = 2.5, pt.bg = c("grey", "white"), cex = 1.3)
   
   # Survey selectivity curves:
   windows()
   tmp <- aggregate(report$S_rv, by = data["S_rv_block"], unique)
   tmp <- tmp[, 2:ncol(tmp)]
   cols <- rainbow(length(unique(data$S_rv_block)))
   plot(range(data$ages), c(0, 1.01), type = "n", 
        xlab = "Age", ylab = "Selectivity proportion", cex.lab = 1.3, cex.axis = 1.2, yaxs = "i")
   grid()
   for (i in 1:nrow(tmp)) lines(data$ages, tmp[i, ], col = cols[i], lwd = 2)
   str <- aggregate(list(min = data$years), by = data["S_rv_block"], min)
   str <- cbind(str, aggregate(list(max = data$years), by = data["S_rv_block"], max)["max"])
   str <- paste0(str$min, "-", str$max)
   legend("topleft", legend = str, lwd = 2, col = cols, bg = "white", cex = 1.5)
   
   # Survey catchability:
   #windows()
   #tmp <- aggregate(exp(obj$par["log_q_rv"]) * report$S_rv, by = data["S_rv_block"], unique)
   #tmp <- tmp[, 2:ncol(tmp)]
   #cols <- rainbow(length(unique(data$S_rv_block)))
   #plot(range(data$ages), c(0, 1.01), type = "n", 
   #     xlab = "Age", ylab = "Catchability", cex.lab = 1.3, cex.axis = 1.2, yaxs = "i")
   #grid()
   #for (i in 1:nrow(tmp)) lines(data$ages, tmp[i, ], col = cols[i], lwd = 2)
   #str <- aggregate(list(min = data$years), by = data["S_rv_block"], min)
   #str <- cbind(str, aggregate(list(max = data$years), by = data["S_rv_block"], max)["max"])
   #str <- paste0(str$min, "-", str$max)
   #legend("topleft", legend = str, lwd = 2, col = cols, bg = "white", cex = 1.5)   
   
   # Sentinel selectivity curves:
   if ("S_sen" %in% names(report)){
      windows()
      tmp <- aggregate(report$S_sen, by = data["S_sen_block"], unique)
      tmp <- tmp[, 2:ncol(tmp)]
      cols <- rainbow(nrow(tmp))
      plot(range(data$ages), c(0, 1.01), type = "n", 
           xlab = "Age", ylab = "Selectivity proportion", cex.lab = 1.3, cex.axis = 1.2, yaxs = "i")
      grid()
      for (i in 1:nrow(tmp)) lines(data$ages, tmp[i, ], col = cols[i], lwd = 2)
      #str <- aggregate(list(min = data$years), by = data["S_sen_block"], min)
      #str <- cbind(str, aggregate(list(max = data$years), by = data["S_sen_block"], max)["max"])
      #str <- paste0(str$min, "-", str$max)
      #legend("topleft", legend = str, lwd = 2, col = cols)      
   }
   
   # Fishery selectivity curves:
   windows()
   tmp <- aggregate(report$S_fishery, by = data["S_fishery_block"], unique)
   tmp <- tmp[, 2:ncol(tmp)]
   cols <- rainbow(nrow(tmp))
   plot(range(data$ages), c(0, 1.01), type = "n", 
        xlab = "Age", ylab = "Selectivity proportion", cex.lab = 1.3, cex.axis = 1.2, yaxs = "i")
   grid()
   for (i in 1:nrow(tmp)) lines(data$ages, tmp[i, ], col = cols[i], lwd = 2)
   str <- aggregate(list(min = data$years), by = data["S_fishery_block"], min)
   str <- cbind(str, aggregate(list(max = data$years), by = data["S_fishery_block"], max)["max"])
   str <- paste0(str$min, "-", str$max)
   legend("topleft", legend = str, lwd = 2, col = cols, bg = "white", cex = 1.5)
        
   # Recruitment deviations:
   if (all(c("init_recDevs", "recDevs") %in% rownames(fixed))){
      windows(width = 9, height = 5)
      irec <- fixed[rownames(fixed) == "init_recDevs", 1] 
      rec <- fixed[rownames(fixed) == "recDevs", 1]
      years <- c((min(data$years)-length(irec)+1):(min(data$years)-1), data$years)
      plot(years, c(rev(irec), rec), pch = 21, bg = "blue", xlab = "Year", ylab = "Recruitment Deviation")
      lines(years, c(rev(irec), rec), lwd= 2, col = "blue")
   }
   
   # Recruitment:
   windows(width = 9, height = 5)
   dbarplot(report$N[, 1] / 1000000, data$year, xlab = "Year", ylab = "Age-2 Recruitment (billions)", xaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
   
   # Recruitment:
   windows(width = 9, height = 5)
   dbarplot(report$N[, 2] / 1000000, data$year, xlab = "Year", ylab = "Age-3 Recruitment (billions)", xaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
      
   # Natural mortality estimates:
   if (all(c("log_M_juvenile", "log_M_adult") %in% rownames(fixed))){
      windows()
      Mjuv <- exp(fixed[rownames(fixed) == "log_M_juvenile", 1])
      Madt <- exp(fixed[rownames(fixed) == "log_M_adult", 1])
      plot(c(1, data$n_M), c(0, max(c(Mjuv, Madt))), type = "n", xlab = "Block", ylab = "Natural Mortality")
      lines(1:data$n_M, Mjuv, col = "blue", lwd = 2)
      lines(1:data$n_M, Madt, col = "green", lwd = 2)
      legend("topleft", c("Juvenile", "adult"), col = c("blue", "green"), lwd = 2)
   }
   if ("log_M" %in% rownames(fixed)){
      windows()
      cols <- rainbow(max(data$M_block_age))
      M <- exp(fixed[rownames(fixed) == "log_M", 1])
      CIL <- exp(fixed[rownames(fixed) == "log_M", 1] - 1.96 * fixed[rownames(fixed) == "log_M", 2])
      CIU <- exp(fixed[rownames(fixed) == "log_M", 1] + 1.96 * fixed[rownames(fixed) == "log_M", 2])
      plot(c(1, max(data$M_block_year)), c(0, max(M)), type = "n", xlab = "Block", ylab = "Natural Mortality")
      dim(M) <- c(max(data$M_block_year), max(data$M_block_age))
      dim(CIL) <- c(max(data$M_block_year), max(data$M_block_age))
      dim(CIU) <- c(max(data$M_block_year), max(data$M_block_age))
      for (i in 1:dim(M)[2]){
         lines(1:dim(M)[1], M[,i], col = cols[i], lwd = 2)
         lines(1:dim(M)[1], CIL[,i], col = cols[i], lwd = 1, lty = "dashed")
         lines(1:dim(M)[1], CIU[,i], col = cols[i], lwd = 1, lty = "dashed")
      }
      tmp <- t(as.data.frame(lapply(by(data$ages, data$M_block_age, unique), range)))
      tmp <- paste(tmp[, 1], tmp[, 2], sep = "-")
      tmp <- paste("Ages", tmp)
      legend("topleft", tmp, col = cols, lwd = 2)   
   }
   

   # Fishing mortality:
   windows(width = 9, height = 5)
   dbarplot(report$F, data$year, xlab = "Year", ylab = "F", xaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
      
   # RV residuals:
   windows(width = 10, height = 5)
   plot(range(data$years), range(data$ages), type = "n", xlab = "Year", ylab = "Age", xaxt = "n", cex.lab = 1.3, cex.axis = 1.2)
   residuals <- log(data$N_rv) - log(report$N_rv)
   scale <- 4
   for (i in 1:length(data$years)){
      index <- residuals[i,] > 0
      points(rep(data$years[i], sum(index)), data$ages[index], cex = scale * sqrt(residuals[i,index]), pch = 21, bg = "grey90")
      points(rep(data$years[i], sum(!index)), data$ages[!index], cex = scale * sqrt(-residuals[i,!index]), pch = 21, bg = "black")
   }
   axis(1, at = seq(1970, 2020, by = 5), cex.axis = 1.2)
   
   # Sentinel residuals:
   if ("N_sen" %in% names(report)){
      windows(width = 10, height = 5)
      years <- as.numeric(rownames(data$N_sen))
      plot(range(years), range(data$ages), type = "n", xlab = "Year", ylab = "Age", cex.lab = 1.3)
      residuals <- log(data$N_sen) - log(report$N_sen)
      scale <- 3
      for (i in 1:length(years)){
         index <- residuals[i,] > 0
         points(rep(years[i], sum(index)), data$ages[index], cex = scale * sqrt(residuals[i,index]), pch = 21, bg = "grey90")
         points(rep(years[i], sum(!index)), data$ages[!index], cex = scale * sqrt(-residuals[i,!index]), pch = 21, bg = "black")
      }
   }
   
   # Fishery residuals:
   #windows(width = 12, height = 6)
   #plot(range(data$years), range(data$ages), type = "n", xlab = "Year", ylab = "Age", cex.lab = 1.3)
   #residuals <- log(data$N_fishery) - log(report$N_fishery)
   #scale <- 3
   #for (i in 1:length(data$years)){
   #   index <- residuals[i,] > 0
   #   points(rep(data$years[i], sum(index)), data$ages[index], cex = scale * sqrt(residuals[i,index]), pch = 21, bg = "grey90")
   #   points(rep(data$years[i], sum(!index)), data$ages[!index], cex = scale * sqrt(-residuals[i,!index]), pch = 21, bg = "black")
   #} 
   
   # Population numbers:
   windows(width = 10, height = 5)
   plot(range(data$years), range(data$ages), type = "n", xlab = "Year", ylab = "Age", cex.lab = 1.3)
   scale <- 0.15
   for (i in 1:length(data$years)){
      points(rep(data$years[i], length(data$ages)), data$ages, cex = scale * sqrt(report$N[i,]/1000), pch = 21, bg = "grey90")
   }   
   box()
   
   windows(width = 9) 
   dbarplot(apply(report$N, 1, sum) / 1000000, data$years, xlab = "Year", ylab = "Abundance (billions)", cex.lab = 1.3)
    
   # Predicted catches plot:
   #windows(); plot(data$C_obs, res[[1]]$report$C)
   
   # Stock-Recruitment plot:
   windows(width = 10, height = 10)
   xx <- res[[1]]$report$SSB_year / 1000
   yy <- res[[1]]$report$N[, 1] * data$weight_at_age_rv[, 1] / 1000
   xx <- xx[1:(length(xx)-2)]
   yy <- yy[3:length(yy)]
   cols <- colorRampPalette(c("blue", "grey", "yellow", "red"))(length(xx))
   plot(xx, yy,
        xlab = "SSB (MT)", ylab = "Recruitment (MT)", 
        xlim = c(0, 1.1*max(xx)), ylim = c(0, 1.1*max(yy)), cex.lab = 1.3, type = "n", cex.axis = 1.2, xaxs = "i", yaxs = "i")
   grid()
   for (i in 1:(length(cols)-1)) lines(xx[i:(i+1)], yy[i:(i+1)], lwd = 2, col = cols[i]) 
   points(xx, yy, pch = 21, bg = cols, cex = 2)
   text(xx, yy, data$year[-c(1,2)], pos = 1)     
   box()
   
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
      #dbarplot(zz, data$years, ylim = 1.05 * c(0, max(c(yy, zz))),
      #                xlim = c(1972.5, 2018), xaxs = "i", legend = FALSE, 
      #                xaxt = "n", width = 1, cex.axis = 1.3, border = "grey30")
      points(data$years, zz, pch = 21, bg = "grey", cex = 1.5)
      
      lines(data$years, yy, col = "red", lwd = 2)
      col <- c("grey", "red")
      if (i == 3){ 
         legend("topright", c("Observed", "Fitted"), pch = c(21, NA), pt.bg = col, 
                col = c("black", "red"), pt.cex = c(2, 0), cex = 1.5, bg = "white", lwd = c(0, 2))
      }
      if (i %in% c(2, 4)) axis(1, at = seq(1970, 2020, by = 5),  cex.axis = 1.3)
      if (i %in% c(1)) mtext("Abundance (millions)", line = 3, side = 2, cex = 1.3, adj = -1.8)
      if (i %in% c(2,4)) mtext("Year", line = 3, side = 1, cex = 1.3)
      text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]),
           par("usr")[3] + 0.90 * diff(par("usr")[3:4]),
           paste("Age", paste(range(age.group[[i]]), collapse = "-")), cex = 1.75)
           
   }

   # Plot of population age-grouped abundance:
   dimnames(report$N) <- list(year = data$years, age = data$ages)
   age.group <- list(2:4, 5:7, 8:10, 11:12)
   m <- repvec(rep(c(1, 3), each = 4), nrow = 3)
   m <- rbind(m, m + 1)
   m <- cbind(0, m, 0)
   m <- rbind(0, m, 0)
   windows(width = 11, height = 7)
   par(mar = c(0, 3, 0, 0))    # c(bottom, left, top, right)
   layout(m)   
   for (i in 1:4){
      col <- dbarplot(report$N[,as.character(age.group[[i]])] / 1000, data$years, xlim = c(1972.5, 2018), xaxs = "i",
                      legend = FALSE, xaxt = "n", width = 1, cex.axis = 1.3, border = "grey30")
      legend("topright", paste("Age", age.group[[i]]), pch = 22, pt.bg = col, pt.cex = 3.5, cex = 1.5, bg = "white")
      if (i %in% c(2, 4)) axis(1, at = seq(1970, 2020, by = 5),  cex.axis = 1.3)
      if (i %in% c(1)) mtext("Abundance (millions)", line = 3, side = 2, cex = 1.3, adj = -1.8)
      if (i %in% c(2,4)) mtext("Year", line = 3, side = 1, cex = 1.3)
   }
   
   #abline(0, 1, col = "red", lwd = 2)
   
}