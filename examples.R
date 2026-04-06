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
                                     n.burn = 700, n.mcmc = 300, verbose = FALSE)
M50 = exdqlmMCMC(y = LakeHuron, p0 = 0.50, model = model,
                                      df = 0.9, dim.df = 2,
                                      PriorGamma = list(m_gam = 0, s_gam = 0.1, df_gam = 1),
                                      n.burn = 700, n.mcmc = 300, verbose = FALSE)
M5 = exdqlmMCMC(y = LakeHuron, p0 = 0.05, model = model,
                                    df = 0.9, dim.df = 2,
                                    PriorGamma = list(m_gam = 1, s_gam = 0.1, df_gam = 1),
                                    n.burn = 700, n.mcmc = 300, verbose = TRUE)
#
# figure 1 (ex1mcmc.png)
par(mfcol=c(1,2))
coda::traceplot(M50$samp.gamma, main = "")
coda::densplot(M50$samp.gamma, main = "")
#
M50 = exdqlmMCMC(y = LakeHuron, p0 = 0.50, model = model,
                                      df = 0.9, dim.df = 2,
                                      gam.init = 0, fix.gamma = TRUE,
                                      n.burn = 700, n.mcmc = 300)
#
# figure 2 (ex1quants.png)
par(mfcol=c(1,2)) # not shown in article code chunk
plot(M95)
exdqlmPlot(M50, add = TRUE, col = "blue")
exdqlmPlot(M5, add = TRUE, col = "forest green")
legend("topright", lty = 1, col = c("purple","blue","forest green"),
                  legend = c(expression('p'[0]*'=0.95'),expression('p'[0]*'=0.50'),
                                                 expression('p'[0]*'=0.05')))
#
fFF = model$FF
fGG =  model$GG
#
plot(LakeHuron,  xlim = c(1952,1980), ylim = c(575,581),
                   col = "dark grey")
exdqlmForecast(start.t = length(LakeHuron), k = 8, M95,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE)
exdqlmForecast(start.t = length(LakeHuron), k = 8, M50,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                 cols = c("blue","light blue"))
exdqlmForecast(start.t = length(LakeHuron), k = 8, M5,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                 cols = c("forest green","green"))

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
M1 = exdqlmLDVB(y = sunspot.year, p0 = 0.95, model = model,
                                     df = c(0.9,0.85), dim.df = c(1,8),
                                     dqlm.ind = TRUE, fix.sigma = FALSE,
                                     verbose = FALSE)
M2 = exdqlmLDVB(y = sunspot.year, p0 = 0.95, model = model,
                                     df = c(0.9,0.85), dim.df = c(1,8),
                                     fix.sigma = FALSE, 
                                     verbose = FALSE)
#
summary(M1)
summary(M2)
# figure 3 (ex2quant.png)
par(mfcol=c(1,3)) # not shown in article code chunk
plot(sunspot.year, col = "dark grey") # not shown in article code chunk
plot(sunspot.year, xlim = c(1750,1850), col = "dark grey",
               ylab = "quantile 95% CrIs")
exdqlmPlot(M1, add = TRUE, col = "red")
exdqlmPlot(M2, add = TRUE, col = "blue")
hist(M2$samp.gamma,xlab=expression(gamma),main="")
#
# figure 4 (ex2checks.png)
par(mfrow=c(2,3)) # not shown in article code chunk
exdqlmDiagnostics(M1, M2, cols = c("red","blue"))
#
possible.dfs = cbind(0.9,seq(0.85,1,0.05))
possible.dfs
#
elbos <- vector("numeric")
ref.samp = rnorm(length(sunspot.year))
for(i in 1:nrow(possible.dfs)){
     temp.M2 = exdqlmLDVB(y = sunspot.year, p0 = 0.85, model = model,
                                             df = possible.dfs[i,], dim.df = c(1,8),
                                             fix.sigma = FALSE, 
                                             verbose = FALSE)
     # temp.check = exdqlmChecks(y = sunspot.year, temp.M2,
     #                                        plot = FALSE, ref = ref.samp)
     elbos = c(elbos,tail(temp.M2$diagnostics$elbo,n=1))
   }
# optimal dfs based off elbo
possible.dfs[which.max(elbos),]


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
model = trend.comp + seas.comp
#
reg.comp = regMod(nino34, m0 = 1, C0 = 1)
model.w.reg = model + reg.comp
#
possible.lams = seq(0.025,0.1,0.025)
elbos <- vector("numeric")
for(i in 1:length(possible.lams)){
  temp.M2 = transfn_exdqlmLDVB(y = log(BTflow), p0 = 0.15, model = model,
                               df = c(1,0.9), dim.df = c(1,2),
                               X = nino34, tf.df = c(0.95), lam = possible.lams[i],
                               tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                               fix.sigma = FALSE,gam.init = -0.1,
                               tol = 0.05, verbose = FALSE)
  # temp.check = exdqlmChecks(y = log(BTflow), temp.M2,
  #                           plot = FALSE, ref = ref.samp)
  elbos = c(elbos,tail(temp.M2$diagnostics$elbo,n=1))
}
# optimal lam based off KL divergence
plot(possible.lams,elbos,type="l")
possible.lams[which.max(elbos)]
#
M1 = exdqlmLDVB(y = log(BTflow), p0 = 0.15, model = model.w.reg,
                df = c(1,0.9,0.95), dim.df = c(1,2,1),
                fix.sigma = FALSE,
                tol = 0.05, verbose = FALSE)
M2 = transfn_exdqlmLDVB(y = log(BTflow), p0 = 0.15, model = model,
                        df = c(1,0.9), dim.df = c(1,2),
                        X = nino34, tf.df = c(0.95), lam = 0.1,
                        tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                        fix.sigma = FALSE,
                        tol = 0.05, verbose = FALSE)
#
# figure 6 (ex3quantcomps.png)
par(mfrow=c(3,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
plot(log(BTflow), col = "grey", ylim = c(1,8), xlim=c(1970,1990),
               ylab = "quantile 95% CrIs")
exdqlmPlot(M1, add = TRUE)
exdqlmPlot(M2, add = TRUE, col = "forest green")
#
plot(NA, col = "grey", ylim = c(-1.5,1.5), xlim=c(1970,1990), ylab = "seasonal components")
compPlot(M1, index = c(2,3), add = TRUE)
compPlot(M2, index = c(2,3), add = TRUE,
         col = "forest green")
#
plot(NA, col = "grey", ylim = c(-0.5,1.5), xlim=c(1970,1990), ylab = "Niño 3.4 components")
compPlot(M1, index = c(4), add = TRUE)
compPlot(M2, index = c(4,5), add =TRUE, col = "forest green")
abline(h = 0, col = "orange", lty = 3, lwd = 2)
#
# figure 7 (ex3zetapsi.png)
par(mfrow=c(2,1),mar=c(4.1,4.2,1.3,1.2)) # not shown in article code chunk
compPlot(M2, index = 4, col= "forest green",
         add = FALSE, just.theta = TRUE)
abline(h = 0, col = "orange", lty = 3, lwd = 2)
title(expression(zeta[t]))
compPlot(M2, index = 5, col = "forest green",
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

