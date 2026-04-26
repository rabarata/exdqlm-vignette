install.packages("exdqlm")
library(exdqlm)


#####################
##### example 1 #####
#####################
#
model = polytrendMod(order = 2, m0 = c(mean(LakeHuron),0), C0 = 10*diag(2))
model
#
M95 = exdqlmMCMC(y = LakeHuron, p0 = 0.95, model = model,
                                     df = 0.9, dim.df = 2,
                                     PriorGamma = list(m_gam = -1, s_gam = 0.1, df_gam = 1),
                                     n.burn = 2000, n.mcmc = 3000, verbose = FALSE)
M50 = exdqlmMCMC(y = LakeHuron, p0 = 0.50, model = model,
                                      df = 0.9, dim.df = 2,
                                      PriorGamma = list(m_gam = 0, s_gam = 0.1, df_gam = 1),
                                      n.burn = 2000, n.mcmc = 3000, verbose = FALSE)
M5 = exdqlmMCMC(y = LakeHuron, p0 = 0.05, model = model,
                                    df = 0.9, dim.df = 2,
                                    PriorGamma = list(m_gam = 1, s_gam = 0.1, df_gam = 1),
                                    n.burn = 2000, n.mcmc = 3000, verbose = TRUE, verbose.every = 1000)
#
# figure 1 (ex1mcmc.png)
par(mfcol=c(2,2))
par(mar = c(5.1, 4.1, 1, 1))
keep.idx = seq(1, length(M50$samp.sigma), by = 10)
sigma.trace = coda::mcmc(M50$samp.sigma[keep.idx], thin = 10)
gamma.trace = coda::mcmc(M50$samp.gamma[keep.idx], thin = 10)
coda::traceplot(sigma.trace, main = expression(sigma), ylab="trace")
coda::densplot(sigma.trace, main = "", ylab = "density")
coda::traceplot(gamma.trace, main = expression(gamma))
coda::densplot(gamma.trace, main = "")

#
M50 = exdqlmMCMC(y = LakeHuron, p0 = 0.50, model = model,
                                      df = 0.9, dim.df = 2,
                                      gam.init = 0, fix.gamma = TRUE,
                                      n.burn = 2000, n.mcmc = 3000)
#
# figure 2 (ex1quants.png)
par(mfcol=c(1,2)) 
par(mar = c(5.1, 4.1, 2.1, 1))
plot(M95); title("Dynamic quantiles")
exdqlmPlot(M50, add = TRUE, col = "blue")
exdqlmPlot(M5, add = TRUE, col = "forest green")
legend("topright", lty = 1, bty = "n", col = c("purple","blue","forest green"),
                  legend = c(expression('p'[0]*'=0.95'),expression('p'[0]*'=0.50'),
                                                 expression('p'[0]*'=0.05')))
#
fFF = model$FF
fGG =  model$GG
#
plot(LakeHuron,  xlim = c(1952,1980), ylim = c(575.5,581.5),
                   col = "dark grey", main = "Forecasted quantiles",
                    ylab = "forecast 95% CrIs")
fc95 = exdqlmForecast(start.t = length(LakeHuron), k = 8, m1 = M95,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                  return.draws = TRUE)
fc50 = exdqlmForecast(start.t = length(LakeHuron), k = 8, m1 = M50,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                 cols = c("blue","light blue"), return.draws = TRUE)
fc05 = exdqlmForecast(start.t = length(LakeHuron), k = 8, m1 = M5,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                 cols = c("forest green","green"),return.draws = TRUE)
#
syn.obs = quantileSynthesis(
            draws_list = list(M5, M50, M95),
            p = c(0.05, 0.50, 0.95),
            T_expected = length(LakeHuron))
# 
syn.fore = quantileSynthesis(
            draws_list = list(fc05, fc50, fc95),
            p = c(0.05, 0.50, 0.95),
            T_expected = 8)

#####################
##### example 2 #####
#####################
#
dlm.trend.comp = dlm::dlmModPoly(1, m0 = mean(sunspot.year), C0 = 10)
trend.comp = as.exdqlm(dlm.trend.comp)
#
seas.comp = seasMod(p = 11, h = 1:4, C0 = 10*diag(8))
#
model = trend.comp + seas.comp
#
model$GG
#
M1 = exdqlmLDVB(y = sunspot.year, p0 = 0.85, model = model,
                                     df = c(0.9,0.85), dim.df = c(1,8),
                                     dqlm.ind = TRUE, fix.sigma = FALSE,
                                     n.samp = 3000, verbose = FALSE)
M2 = exdqlmLDVB(y = sunspot.year, p0 = 0.85, model = model,
                                     df = c(0.9,0.85), dim.df = c(1,8),
                                     fix.sigma = FALSE, 
                                     n.samp = 3000, verbose = FALSE)
#
summary(M1)
summary(M2)
# figure 3 (ex2quant.png)
# not shown in article code chunk
par(mfcol=c(1,1)) 
plot(sunspot.year, col = "dark grey", xlim = c(1710,1979)) 
#
par(mfcol=c(1,2)) # not shown in article code chunk
plot(sunspot.year, xlim = c(1750,1850), col = "dark grey",
               ylab = "quantile 95% CrIs")
exdqlmPlot(M1, add = TRUE, col = "red")
exdqlmPlot(M2, add = TRUE, col = "blue")
legend("topleft", lty = 1, , bty = "n", col = c("red","blue"),legend = c("DQLM","exDQLM"))
#
hist(M2$samp.gamma, xlab=expression(gamma), main="")
abline(v = mean(M2$samp.gamma), col="blue")
#
# figure 4 (ex2checks.png)
par(mfrow=c(2,3)) # not shown in article code chunk
diagM1M2 = exdqlmDiagnostics(M1, M2, cols = c("red","blue"))
print(diagM1M2)
#
possible.dfs = cbind(0.9,seq(0.85,1,0.05))
possible.dfs
#
metrics <- matrix(NA_real_, nrow(possible.dfs), 2,
                  dimnames = list(NULL, c("CRPS","KL")))
ref.samp = rnorm(length(sunspot.year))
for(i in 1:nrow(possible.dfs)){
     temp.M2 = exdqlmLDVB(y = sunspot.year, p0 = 0.85, model = model,
                                             df = possible.dfs[i,], dim.df = c(1,8),
                                             sig.init = 2, fix.sigma = FALSE, 
                                             n.samp = 3000, verbose = FALSE)
     temp.check = exdqlmDiagnostics(temp.M2, plot = FALSE, ref = ref.samp)
     metrics[i,] = c(temp.check$m1.CRPS, temp.check$m1.KL)
   }
# optimal dfs based off CRPS
possible.dfs[which.min(metrics[, "CRPS"]),]
#
M1mcmc = exdqlmMCMC(y = sunspot.year, p0 = 0.85, model = model,
                df = c(0.9,0.85), dim.df = c(1,8),
                n.burn = 2000, n.mcmc = 3000, verbose = FALSE,
                dqlm.ind = TRUE,)
M2mcmc = exdqlmMCMC(y = sunspot.year, p0 = 0.85, model = model,
                df = c(0.9,0.85), dim.df = c(1,8),
                n.burn = 2000, n.mcmc = 3000, verbose = FALSE)


#####################
##### example 3 #####
##### OPTION 1 ######
#####################
#
# figure 5 (ex3data.png)
par(mfrow=c(2,1)) # not shown in article code chunk
plot(log(BTflow), col="dark grey",xlim = c(1988,2025),cex.lab=1.5,cex.axis=1.5,cex.main=1.8) # not shown in article code chunk
plot(nino34, col="dark grey",xlim = c(1988,2025),cex.lab=1.5,cex.axis=1.5,cex.main=1.8) # not shown in article code chunk
#
trend.comp = polytrendMod(1, m0=3, C0=0.1)
seas.comp = seasMod(p = 12, h = c(1, 2, 0.1469118636), C0 = diag(1,6))
model = trend.comp + seas.comp
#
reg.comp = regMod(nino34, m0 = 0, C0 = matrix(1,1,1))
model.w.reg = model + reg.comp
#
possible.lams = seq(0.1,0.9,0.1)
metrics <- matrix(NA_real_, length(possible.lams), 3,
                  dimnames = list(NULL, c("CRPS","KL","pplc")))
ref.samp = rnorm(length(BTflow))
for(i in 1:length(possible.lams)){
  temp.M2 = exdqlmTransferLDVB(y = log(BTflow), p0 = 0.15, model = model,
                               df = c(0.97,0.97), dim.df = c(1,6),
                               X = nino34, tf.df = c(0.95), lam = possible.lams[i],
                               tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                               sig.init = 0.1, gam.init = -0.1,
                               tol = 0.05, verbose = FALSE)
  temp.check = exdqlmDiagnostics(temp.M2, plot = FALSE, ref = ref.samp)
  metrics[i,] = c(temp.check$m1.CRPS, temp.check$m1.KL,temp.check$m1.pplc)
}
# not shown in article code chunk
par(mfrow=c(3,1)) 
plot(possible.lams,metrics[,1],type="l",main="CRPS")
plot(possible.lams,metrics[,2],type="l",main="KLs")
plot(possible.lams,metrics[,3],type="l",main="pplcs")
# optimal lam based off KL divergence
possible.lams[which.min(metrics[, "CRPS"])]
#
M1 = exdqlmLDVB(y = log(BTflow), p0 = 0.15, model = model.w.reg,
                df = c(0.97,0.97,0.97), dim.df = c(1,6,1),
                sig.init = 0.1, gam.init = -0.1,
                tol = 0.05, verbose = FALSE)
M2 = exdqlmTransferLDVB(y = log(BTflow), p0 = 0.15, model = model,
                        df = c(0.97,0.97), dim.df = c(1,6),
                        X = nino34, tf.df = c(0.97), lam = 0.9,
                        tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                        sig.init = 0.1, gam.init = -0.1,
                        tol = 0.05, verbose = FALSE)
#
# figure 6 (ex3quantcomps.png)
par(mfrow=c(3,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
plot(log(BTflow), col = "grey", ylim = c(1,8), xlim=c(1990,2010),
               ylab = "quantile 95% CrIs")
exdqlmPlot(M1, add = TRUE)
exdqlmPlot(M2, add = TRUE, col = "forest green")
#
plot(NA, col = "grey", ylim = c(-2,2), xlim=c(1990,2010), ylab = "seasonal components")
compPlot(M1, index = 2:7, add = TRUE)
compPlot(M2, index = 2:7, add = TRUE,
         col = "forest green")
#
plot(NA, col = "grey", ylim = c(-1.5,1.5), xlim=c(1990,2010), ylab = "Niño 3.4 components")
compPlot(M1, index = 8, add = TRUE)
compPlot(M2, index = 8:9, add =TRUE, col = "forest green")
abline(h = 0, col = "orange", lty = 3, lwd = 2)
#
# figure 7 (ex3zetapsi.png)
par(mfrow=c(2,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
compPlot(M2, index = 8, col= "forest green",
         add = FALSE, just.theta = TRUE)
abline(h = 0, col = "orange", lty = 3, lwd = 2)
title(expression(zeta[t]))
compPlot(M2, index = 9, col = "forest green",
         add = FALSE, just.theta = TRUE)
abline(h = 0, col = "orange", lty = 3, lwd = 2)
title(expression(psi[t]))
#
M2$median.kt
#
# figure 8 (ex3forecast.png)
par(mfrow=c(1,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
exdqlmForecast(start.t = length(BTflow) - 18, k = 18, M1)
exdqlmForecast(start.t = length(BTflow) - 18, k = 18, M2, add = TRUE, cols = c("forest green","green"))
#
checks.out = exdqlmDiagnostics(M1, M2,
                          col = c("purple","forest green"), plot = FALSE)
print(checks.out)


#####################
##### example 3 #####
##### OPTION 2 ######
#####################

# figure 5 (ex3data.png)
par(mfrow=c(2,1)) # not shown in article code chunk
plot(BTflow, col="dark grey",cex.lab=1.5,cex.axis=1.5,cex.main=1.8) # not shown in article code chunk
plot(BTprec, col="dark grey",cex.lab=1.5,cex.axis=1.5,cex.main=1.8) # not shown in article code chunk
#
trend.comp = polytrendMod(1,m0=3,C0=0.1)
seas.comp = seasMod(p = 12, h = c(1, 2, 0.1469118636), C0 = diag(1,6))
model = trend.comp + seas.comp
#
reg.comp = regMod(BTprec, m0 = 0, C0 = 1)
model.w.reg = model + reg.comp
#
possible.lams = seq(0.1,0.9,0.1)
metrics <- matrix(NA_real_, length(possible.lams), 3,
                  dimnames = list(NULL, c("CRPS","KL","pplc")))
ref.samp = rnorm(length(BTflow))
for(i in 3:length(possible.lams)){
  temp.M2 = exdqlmTransferLDVB(y = BTflow, p0 = 0.15, model = model,
                               df = c(0.97,0.97), dim.df = c(1,6),
                               X = BTprec, tf.df = c(0.95), lam = possible.lams[i],
                               tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                               sig.init = 0.1, gam.init = -0.1,
                               tol = 0.05, verbose = FALSE)
  temp.check = exdqlmDiagnostics(temp.M2, plot = FALSE, ref = ref.samp)
  metrics[i,] = c(temp.check$m1.CRPS, temp.check$m1.KL,temp.check$m1.pplc)
}
# not shown in article code chunk
par(mfrow=c(3,1)) 
plot(possible.lams,metrics[,1],type="l",main="CRPS")
plot(possible.lams,metrics[,2],type="l",main="KLs")
plot(possible.lams,metrics[,3],type="l",main="pplcs")
# optimal lam based off KL divergence
possible.lams[which.min(metrics[, "CRPS"])] # 0.5
#
M1 = exdqlmLDVB(y = BTflow, p0 = 0.15, model = model.w.reg,
                df = c(0.97,0.97,0.97), dim.df = c(1,6,1),
                sig.init = 0.1, gam.init = -0.1,
                tol = 0.05, verbose = FALSE)
M2 = exdqlmTransferLDVB(y = BTflow, p0 = 0.15, model = model,
                        df = c(0.97,0.97), dim.df = c(1,6),
                        X = BTprec, tf.df = c(0.97), lam = 0.5,
                        tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                        sig.init = 0.1, gam.init = -0.1,
                        tol = 0.05, verbose = FALSE)
#
# figure 6 (ex3quantcomps.png)
par(mfrow=c(3,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
plot(BTflow, col = "grey",ylim = c(0,500), xlim=c(1992,2008),
     ylab = "quantile 95% CrIs")
exdqlmPlot(M1, add = TRUE)
exdqlmPlot(M2, add = TRUE, col = "forest green")
#
plot(NA, col = "grey", ylim = c(-58,58), xlim=c(1992,2008), ylab = "seasonal components")
compPlot(M1, index = 2:7, add = TRUE)
compPlot(M2, index = 2:7, add = TRUE,
         col = "forest green")
#
plot(NA, col = "grey", ylim = c(0,300), xlim=c(1992,2008), ylab = "Niño 3.4 components")
compPlot(M1, index = 8, add = TRUE)
compPlot(M2, index = 8:9, add =TRUE, col = "forest green")
abline(h = 0, col = "orange", lty = 3, lwd = 2)
#
# figure 7 (ex3zetapsi.png)
par(mfrow=c(2,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
compPlot(M2, index = 8, col= "forest green",
         add = FALSE, just.theta = TRUE)
abline(h = 0, col = "orange", lty = 3, lwd = 2)
title(expression(zeta[t]))
compPlot(M2, index = 9, col = "forest green",
         add = FALSE, just.theta = TRUE)
abline(h = 0, col = "orange", lty = 3, lwd = 2)
title(expression(psi[t]))
#
M2$median.kt
#
# figure 8 (ex3forecast.png)
par(mfrow=c(1,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
exdqlmForecast(start.t = length(BTflow) - 18, k = 18, M1)
exdqlmForecast(start.t = length(BTflow) - 18, k = 18, M2, add = TRUE, cols = c("forest green","green"))
#
checks.out = exdqlmDiagnostics(M1, M2,
                               col = c("purple","forest green"), plot = FALSE)
print(checks.out)


