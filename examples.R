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
                                     fix.sigma = TRUE, sig.init = 0.07,
                                     PriorGamma = list(m_gam = -1, s_gam = 0.1, df_gam = 1),
                                     n.burn = 700, n.mcmc = 300, verbose = FALSE)
M50 = exdqlmMCMC(y = LakeHuron, p0 = 0.50, model = model,
                                      df = 0.9, dim.df = 2,
                                      fix.sigma = TRUE, sig.init = 0.4,
                                      PriorGamma = list(m_gam = 0, s_gam = 0.1, df_gam = 1),
                                      n.burn = 700, n.mcmc = 300, verbose = FALSE)
M5 = exdqlmMCMC(y = LakeHuron, p0 = 0.05, model = model,
                                    df = 0.9, dim.df = 2,
                                    fix.sigma = TRUE, sig.init = 0.07,
                                    PriorGamma = list(m_gam = 1, s_gam = 0.1, df_gam = 1),
                                    n.burn = 700, n.mcmc = 300, verbose = TRUE)
#
# figure 1 (ex1mcmc.png)
par(mfcol=c(1,2))
library(coda)
traceplot(M50$samp.gamma, main = "")
densplot(M50$samp.gamma, main = "")
#
M50 = exdqlmMCMC(y = LakeHuron, p0 = 0.50, model = model,
                                      df = 0.9, dim.df = 2,
                                      fix.sigma = TRUE, sig.init = 0.4,
                                      gam.init = 0, fix.gamma = TRUE,
                                      n.burn = 700, n.mcmc = 300)
#
# figure 2 (ex1quants.png)
par(mfcol=c(1,2)) # not shown in article code chunk
exdqlmPlot(y = LakeHuron, M95)
exdqlmPlot(y = LakeHuron, M50, add = TRUE, col = "blue")
exdqlmPlot(y = LakeHuron, M5, add = TRUE, col = "forest green")
legend("topright", lty = 1, col = c("purple","blue","forest green"),
                  legend = c(expression('p'[0]*'=0.95'),expression('p'[0]*'=0.50'),
                                                 expression('p'[0]*'=0.05')))
#
fFF = model$FF
fGG =  model$GG
#
plot(LakeHuron,  xlim = c(1952,1980), ylim = c(575,581),
                   col = "dark grey")
exdqlmForecast(y = LakeHuron, start.t = length(LakeHuron), k = 8, M95,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE)
exdqlmForecast(y = LakeHuron, start.t = length(LakeHuron), k = 8, M50,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                 cols = c("blue","light blue"))
exdqlmForecast(y = LakeHuron, start.t = length(LakeHuron), k = 8, M5,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                 cols = c("forest green","green"))

#####################
##### example 2 #####
#####################
#
library(dlm)
dlm.trend.comp = dlmModPoly(1, m0 = mean(sunspot.year), C0 = 10)
#
seas.comp = seasMod(p = 11, h = 1:4, C0 = 10*diag(8))
#
model = combineMods(dlm.trend.comp,seas.comp)
#
trend.comp = dlmMod(dlm.trend.comp)
model = combineMods(trend.comp,seas.comp)
model$GG
#
M1 = exdqlmISVB(y = sunspot.year, p0 = 0.85, model = model,
                                     df = c(0.9,0.85), dim.df = c(1,8),
                                     dqlm.ind = TRUE, fix.sigma = FALSE)
#
summary(M1$samp.sigma)
#
M1 = exdqlmISVB(y = sunspot.year, p0 = 0.85, model = model,
                                     df = c(0.9,0.85), dim.df = c(1,8),
                                     dqlm.ind = TRUE, sig.init = 2,
                                     verbose = FALSE)
M2 = exdqlmISVB(y = sunspot.year, p0 = 0.85, model = model,
                                     df = c(0.9,0.85), dim.df = c(1,8),
                                     sig.init = 2, verbose = FALSE)
#
# figure 3 (ex2quant.png)
par(mfcol=c(1,3)) # not shown in article code chunk
plot(sunspot.year, col = "dark grey") # not shown in article code chunk
plot(sunspot.year, xlim = c(1750,1850), col = "dark grey",
               ylab = "quantile 95% CrIs")
exdqlmPlot(y = sunspot.year, M1, add = TRUE, col = "red")
exdqlmPlot(y = sunspot.year, M2, add = TRUE, col = "blue")
#
hist(M2$samp.gamma,xlab=expression(gamma),main="")
#
# figure 4 (ex2checks.png)
par(mfrow=c(2,3)) # not shown in article code chunk
exdqlmChecks(y = sunspot.year, M1, M2, cols = c("red","blue"))
#
possible.dfs = cbind(0.9,seq(0.85,1,0.05))
possible.dfs
#
KLs <- vector("numeric")
ref.samp = rnorm(length(sunspot.year))
for(i in 1:nrow(possible.dfs)){
     temp.M2 = exdqlmISVB(y = sunspot.year, p0 = 0.85, model = model,
                                             df = possible.dfs[i,], dim.df = c(1,8),
                                             sig.init = 2, verbose = FALSE)
     temp.check = exdqlmChecks(y = sunspot.year, temp.M2,
                                             plot = FALSE, ref = ref.samp)
     KLs = c(KLs,temp.check$m1.KL)
   }
# optimal dfs based off KL divergence
possible.dfs[which.min(KLs),]


#####################
##### example 3 #####
#####################
#
# figure 5 (ex3data.png)
par(mfrow=c(2,1)) # not shown in article code chunk
plot(log(BTflow), col="dark grey",cex.lab=1.5,cex.axis=1.5,cex.main=1.8) # not shown in article code chunk
plot(nino34, col="dark grey",cex.lab=1.5,cex.axis=1.5,cex.main=1.8) # not shown in article code chunk
#
trend.comp = polytrendMod(1,m0=3,C0=0.1)
seas.comp = seasMod(p = 12, h = 1, C0 = diag(1,2))
model = combineMods(trend.comp,seas.comp)
model
#
reg.comp <- NULL
reg.comp$m0 = 0
reg.comp$C0 = 1
reg.comp$FF = matrix(nino34,nrow = 1)
reg.comp$GG = 1
#
model.w.reg = combineMods(model,reg.comp)
#
possible.lams = seq(0.5,0.85,0.025)
KLs <- vector("numeric")
ref.samp = rnorm(length(BTflow))
for(i in 1:length(possible.lams)){
  temp.M2 = transfn_exdqlmISVB(y = log(BTflow), p0 = 0.15, model = model,
                               df = c(1,0.9), dim.df = c(1,2),
                               X = nino34, tf.df = c(0.95), lam = possible.lams[i],
                               tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                               sig.init = 0.1, gam.init = - 0.1,
                               tol = 0.05, verbose = FALSE)
  temp.check = exdqlmChecks(y = log(BTflow), temp.M2,
                            plot = FALSE, ref = ref.samp)
  KLs = c(KLs,temp.check$m1.KL)
}
# optimal lam based off KL divergence
possible.lams[which.min(KLs)]
#
M1 = exdqlmISVB(y = log(BTflow), p0 = 0.15, model = model.w.reg,
                df = c(1,0.9,0.95), dim.df = c(1,2,1),
                sig.init = 0.1, gam.init = - 0.1,
                tol = 0.05, verbose = FALSE)
M2 = transfn_exdqlmISVB(y = log(BTflow), p0 = 0.15, model = model,
                        df = c(1,0.9), dim.df = c(1,2),
                        X = nino34, tf.df = c(0.95), lam = 0.85,
                        tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                        sig.init = 0.1, gam.init = - 0.1,
                        tol = 0.05, verbose = FALSE)
#
# figure 6 (ex3quantcomps.png)
par(mfrow=c(3,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
plot(log(BTflow), col = "grey", ylim = c(1,8), xlim=c(1970,1990),
               ylab = "quantile 95% CrIs")
exdqlmPlot(y = BTflow, M1, add = TRUE)
exdqlmPlot(y = BTflow, M2, add = TRUE, col = "forest green")
#
plot(NA, col = "grey", ylim = c(-1.5,1.5), xlim=c(1970,1990), ylab = "seasonal components")
compPlot(y = log(BTflow), M1, index = c(2,3), add = TRUE)
compPlot(y = log(BTflow), M2, index = c(2,3), add = TRUE,
         col = "forest green")
#
plot(NA, col = "grey", ylim = c(-0.5,1.5), xlim=c(1970,1990), ylab = "Niño 3.4 components")
compPlot(y = log(BTflow), M1, index = c(4), add = TRUE)
compPlot(y = log(BTflow), M2, index = c(4,5), add =TRUE, col = "forest green")
abline(h = 0, col = "orange", lty = 3, lwd = 2)
#
# figure 7 (ex3zetapsi.png)
par(mfrow=c(2,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
compPlot(y = log(BTflow), M2, index = 4, col= "forest green",
         add = FALSE, just.theta = TRUE)
abline(h = 0, col = "orange", lty = 3, lwd = 2)
title(expression(zeta[t]))
compPlot(y = log(BTflow), M2, index = 5, col = "forest green",
         add = FALSE, just.theta = TRUE)
abline(h = 0, col = "orange", lty = 3, lwd = 2)
title(expression(psi[t]))
#
M2$median.kt
#
# figure 8 (ex3forecast.png)
par(mfrow=c(1,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
exdqlmForecast(y = log(BTflow), start.t = length(BTflow) - 18, k = 18, M1)
exdqlmForecast(y = log(BTflow), start.t = length(BTflow) - 18, k = 18, M2, add = TRUE, cols = c("forest green","green"))
#
checks.out = exdqlmChecks(y = log(BTflow), M1, M2,
                          col = c("purple","orange"), plot = FALSE)
checks.out$m1.KL # not shown in article code chunk
checks.out$m2.KL # not shown in article code chunk
checks.out$m1.pplc # not shown in article code chunk
checks.out$m2.pplc # not shown in article code chunk

