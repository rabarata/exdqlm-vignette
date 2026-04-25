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
legend("topright", lty = 1, col = c("purple","blue","forest green"),
                  legend = c(expression('p'[0]*'=0.95'),expression('p'[0]*'=0.50'),
                                                 expression('p'[0]*'=0.05')))
#
fFF = model$FF
fGG =  model$GG
#
plot(LakeHuron,  xlim = c(1952,1980), ylim = c(575.5,581.5),
                   col = "dark grey", main = "Forecasted quantiles",
                    ylab = "forecast 95% CrIs")
exdqlmForecast(start.t = length(LakeHuron), k = 8, M95,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE)
exdqlmForecast(start.t = length(LakeHuron), k = 8, M50,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                 cols = c("blue","light blue"))
exdqlmForecast(start.t = length(LakeHuron), k = 8, M5,
                                 fFF = fFF, fGG = fGG, plot = TRUE, add = TRUE,
                                 cols = c("forest green","green"))
#
syn.obs = quantileSynthesis(
      draws_list = list(M5$samp.post.pred, M50$samp.post.pred,
                      + M95$samp.post.pred),
        p = c(0.05, 0.50, 0.95),
          T_expected = length(LakeHuron))
# 
## NOT WORKING
n.samp = M5$n.mcmc
q.fore = sweep(matrix(rnorm(8 * n.samp), 8, n.samp), 1,
               sqrt(fc05$fQ), "*") + fc05$ff
future.M5$ydraw = vapply(1:n.samp, function(j)
  rexal(8, p0 = M5$p0, mu = q.fore[, j],
          sigma = sigma.draws[j], gamma = gamma.draws[j]),
              numeric(8))

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
                                     verbose = FALSE)
M2 = exdqlmLDVB(y = sunspot.year, p0 = 0.85, model = model,
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
d = exdqlmDiagnostics(M1, M2, cols = c("red","blue"))
print(d)
#
possible.dfs = cbind(0.9,seq(0.85,1,0.05))
possible.dfs
#
# elbos <- vector("numeric")
crpss <- vector("numeric")
for(i in 1:nrow(possible.dfs)){
     temp.M2 = exdqlmLDVB(y = sunspot.year, p0 = 0.85, model = model,
                                             df = possible.dfs[i,], dim.df = c(1,8),
                                             fix.sigma = FALSE, 
                                             verbose = FALSE)
     temp.d = exdqlmDiagnostics(temp.M2, plot = FALSE, ref = ref.samp)
     #elbos = c(elbos,tail(temp.M2$diagnostics$elbo,n=1))
     crpss = c(crpss,temp.d$m1.CRPS)
   }
# optimal dfs based off elbo
plot(possible.dfs[,2],crpss,type="l")
possible.dfs[which.min(crpss),]


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
possible.lams = seq(0.975,0.99,0.005)
KLs <- vector("numeric")
crpss <- vector("numeric")
for(i in 1:length(possible.lams)){
  temp.M2 = transfn_exdqlmLDVB(y = log(BTflow), p0 = 0.15, model = model,
                               df = c(1,0.9), dim.df = c(1,2),
                               X = nino34, tf.df = c(0.95), lam = possible.lams[i],
                               tf.m0 = c(0,0), tf.C0 = diag(c(0.1,0.005),2),
                               fix.sigma = FALSE,gam.init = -0.1,
                               tol = 0.05, verbose = FALSE)
  temp.d = exdqlmDiagnostics(temp.M2, plot = FALSE)
  KLs = c(KLs,temp.d$m1.KL)
  crpss = c(crpss,temp.d$m1.CRPS)
}
# optimal lam based off KL divergence
plot(possible.lams,crpss,type="l")
plot(possible.lams,KLs,type="l")
possible.lams[which.min(KLs)]
#
M1 = exdqlmLDVB(y = log(BTflow), p0 = 0.15, model = model.w.reg,
                df = c(1,0.9,0.95), dim.df = c(1,2,1),
                fix.sigma = FALSE,
                tol = 0.05, verbose = FALSE)
M2 = transfn_exdqlmLDVB(y = log(BTflow), p0 = 0.15, model = model,
                        df = c(1,0.9), dim.df = c(1,2),
                        X = nino34, tf.df = c(0.95), lam = 0.985,
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

