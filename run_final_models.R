rm(list = ls())

library(sp)
library(INLA)
library(gdata)
library(ggplot2)


setwd("C:/HY-Data/SROTO/documents/fungalsex")
load(file= 'FungalSexMapAndMesh.Rdata')
load(file = 'AllDataFungalSex.Rdata')

names(d) 


  
head(d)

####################################################

# NUMBER OF NEW STRAINS:

# condition on existence previous year and presence this year ------------------------------------

 inds <- which((d$numberofstrains>0) & (d$numberofstrainsnextyear>0) & (d$nsamples <15)) 
 data_norm2<-d[inds,]
 data_norm2[, c("numberofstrains", "pathogenconnectivity", 
       "hostconnectivity", "nsamples", "proportionofcoinfected",  
       "numberofcoinfections", "absoluteabundance",  "relativeabundance", 
       "pl",  "area",  "abundance1", "abundance2" , "abundance3")] <- scale(data_norm2[, c("numberofstrains", "pathogenconnectivity", 
                                                                                  "hostconnectivity", "nsamples", "proportionofcoinfected",  
                                                                                  "numberofcoinfections", "absoluteabundance",  "relativeabundance", 
                                                                                  "pl",  "area",  "abundance1", "abundance2" , "abundance3")] )
 sp1_gps_matrix2 <- data_norm2[inds, c('x', 'y')]
 
 Nyears <- 3
 n.years = Nyears
 sigma0 = 1
 size = min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) 
 
 range0 = size / 5
 kappa0 = sqrt(8) / range0
 tau0 = 1 / (sqrt(4 * pi) * kappa0 * sigma0)
 
 size2 =(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) 
 range1 <- range(mesh$loc[, 1])
 range2 <- range(mesh$loc[, 2])
 
 spde = inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),   B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0), theta.prior.prec = c(0.1, 1) )
 
 spde.idx <- inla.spde.make.index("spatial", n.spde=mesh$n,  n.group= n.years)
 group.years <-   as.integer(data_norm2$year)-2011       
 
 Atemp2 <- inla.spde.make.A(mesh=mesh, loc=as.matrix(sp1_gps_matrix2), group = group.years)
 effects2 <- list(c(list(intercept=rep(1,mesh$n)), spde.idx),plan=data_norm2)

 hyper.prec = list(prec = list(param = c(1, 0.05)))
 h.spec <- list(theta=list(initial = 0.1, param =c(0,5)))
 sigma <- 5
 
 hyper=spde$f$hyper.default
 hyper$theta1$initial = hyper$theta1$initial - log(sigma)
 hyper$theta1$param[1] = hyper$theta1$param[1] - log(sigma)
 
 # scale = 57052.84 := 0.5125143 sca
 # - > 1 sca = 111319.5

 h.spec <- list(theta=list(initial = 0.1, param =c(0,5)))
 sigma <- 5
 hyper=spde$f$hyper.default
 hyper$theta1$initial = hyper$theta1$initial - log(sigma)
 hyper$theta1$param[1] = hyper$theta1$param[1] - log(sigma)
 inla.setOption(scale.model.default = TRUE)
 
 Ntrials <- rep(1, length(data_norm2$patch))
 
 # new strain emergence models (conditional on nonextinction) ----------------------------------------------------
 
 # First binomial, then for the number of new strains
 
 data <- inla.stack(data=list(Y=1*(data_norm2$newstrainsnextyear>0), Ntrials = Ntrials), A = list(Atemp2, 1), effects2, tag="obs")
 formula <- as.formula(Y ~ -1+Intercept+abundance2+abundance3+numberofcoinfections+pathogenconnectivity+numberofstrains+pl+ f(spatial, model = spde, group = spatial.group, control.group =list(model='ar1', hyper=h.spec)))
 resNewNonextinction <- inla(formula, family = 'binomial', Ntrials= data$Ntrials,  data = inla.stack.data(data), control.predictor =list(compute=TRUE,A=inla.stack.A(data)), control.fixed=list(expand.factor.strategy='inla'),  control.compute = list(dic=TRUE, waic = TRUE, cpo = TRUE), control.inla = list(int.strategy = "eb"), num.threads = 10)
 summary(resNewNonextinction)

 data <- inla.stack(data=list(Y=1*(data_norm2$newstrainsnextyear), Ntrials = Ntrials), A = list(Atemp2, 1), effects2, tag="obs")
 formula <- as.formula(Y ~ -1+Intercept+abundance2+abundance3+numberofcoinfections+pathogenconnectivity+numberofstrains+pl+ f(spatial, model = spde, group = spatial.group, control.group =list(model='ar1', hyper=h.spec)))
 resNewNonextinctionPoisson <- inla(formula, family = 'poisson', Ntrials= data$Ntrials,  data = inla.stack.data(data), control.predictor =list(compute=TRUE,A=inla.stack.A(data)), control.fixed=list(expand.factor.strategy='inla'),  control.compute = list(dic=TRUE, waic = TRUE, cpo = TRUE), control.inla = list(int.strategy = "eb"), num.threads = 10)
 summary(resNewNonextinctionPoisson)

 
 
 #########################
 
 # COINFECTION PRESENCE:
 
 # condition on presence this year ------------------------------------
 
 inds <- which(d$numberofstrainsC>0 & (d$nsamples <15))
 data_norm2<-d[inds,]
 data_norm2[, c("numberofstrains", "pathogenconnectivity", 
                "hostconnectivity", "nsamples", "proportionofcoinfected",  
                "numberofcoinfections", "absoluteabundance",  "relativeabundance", 
                "pl",  "area",  "abundance1", "abundance2" , "abundance3")] <- scale(data_norm2[, c("numberofstrains", "pathogenconnectivity", 
                                                                                                    "hostconnectivity", "nsamples", "proportionofcoinfected",  
                                                                                                    "numberofcoinfections", "absoluteabundance",  "relativeabundance", 
                                                                                                    "pl",  "area",  "abundance1", "abundance2" , "abundance3")] )
 
 sp1_gps_matrix2 <- d[inds,c('x', 'y')]
 Nobs2 = nrow(data_norm2)
 Ntrials <- rep(1, Nobs2)

 data_norm2[order(data_norm2$patch),] 
 
 n.years = Nyears
 sigma0 = 1
 size = min(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) # in meters?
 
 range0 = size / 5
 kappa0 = sqrt(8) / range0
 tau0 = 1 / (sqrt(4 * pi) * kappa0 * sigma0)
 
 spde = inla.spde2.matern(mesh, B.tau = cbind(log(tau0), -1, +1),   B.kappa = cbind(log(kappa0), 0, -1), theta.prior.mean = c(0, 0), theta.prior.prec = c(0.1, 1) )
 
 spde.idx <- inla.spde.make.index("spatial", n.spde=mesh$n,  n.group= n.years)
 group.years <-   as.integer(data_norm2$year)-2011       #.....#    #HERE OUT A VECTOR DESCRIBING YEAR FOR EVERY OBSERVATION
 
 Atemp2 <- inla.spde.make.A(mesh=mesh, loc=as.matrix(sp1_gps_matrix2), group = group.years)
 effects2 <- list(c(list(intercept=rep(1,mesh$n)), spde.idx), plan=data_norm2)
 
 hyper.prec = list(prec = list(param = c(1, 0.05)))
 h.spec <- list(theta=list(initial = 0.1, param =c(0,5)))
 sigma <- 5
 hyper=spde$f$hyper.default
 hyper$theta1$initial = hyper$theta1$initial - log(sigma)
 hyper$theta1$param[1] = hyper$theta1$param[1] - log(sigma)
 
 
 size2 =(c(diff(range(mesh$loc[, 1])), diff(range(mesh$loc[, 2])))) # in meters?
 range1 <- range(mesh$loc[, 1])
 range2 <- range(mesh$loc[, 2])
 
 # scale = 57052.84 := 0.5125143 sca
 # - > 1 sca = 111319.5
 
 
 h.spec <- list(theta=list(initial = 0.1, param =c(0,5)))
 sigma <- 5
 hyper=spde$f$hyper.default
 hyper$theta1$initial = hyper$theta1$initial - log(sigma)
 hyper$theta1$param[1] = hyper$theta1$param[1] - log(sigma)
 inla.setOption(scale.model.default = TRUE)
 

 # coinfection models, binary and poisson for the number (conditional on presence now) ----------------------------------------------------

 data <- inla.stack(data=list(Y=1*(data_norm2$numberofcoinfections>0), Ntrials = Ntrials), A = list(Atemp2, 1), effects2, tag="obs")
 formula <- as.formula(Y ~ -1+Intercept+abundance2+abundance3+pathogenconnectivity+numberofstrains+pl+f(spatial, model = spde, group = spatial.group, control.group =list(model='ar1', hyper=h.spec)))
 resCoPresenceNow <- inla(formula, family = 'binomial',Ntrials= data$Ntrials,  data = inla.stack.data(data), control.predictor =list(compute=TRUE,A=inla.stack.A(data)), control.fixed=list(expand.factor.strategy='inla'),  control.compute = list(dic=TRUE, waic = TRUE, cpo = TRUE), control.inla = list(int.strategy = "eb"), num.threads = 10)
 summary(resCoPresenceNow)
 
 data <- inla.stack(data=list(Y=1*(data_norm2$numberofcoinfectionsC), Ntrials = Ntrials), A = list(Atemp2, 1), effects2, tag="obs")
 formula <- as.formula(Y ~ -1+Intercept+abundance2+abundance3+pathogenconnectivity+numberofstrains+pl+f(spatial, model = spde, group = spatial.group, control.group =list(model='ar1', hyper=h.spec)))
 resCoPresenceNowPoisson <- inla(formula, family = 'poisson',Ntrials= data$Ntrials,  data = inla.stack.data(data), control.predictor =list(compute=TRUE,A=inla.stack.A(data)), control.fixed=list(expand.factor.strategy='inla'),  control.compute = list(dic=TRUE, waic = TRUE, cpo = TRUE), control.inla = list(int.strategy = "eb"), num.threads = 10)
 summary(resCoPresenceNowPoisson)
 
 
 # overwintering model, binary ----------------------------------------------------
 
 data <- inla.stack(data=list(Y=1*(data_norm2$numberofstrainsnextyear>0), Ntrials = Ntrials), A = list(Atemp2, 1), effects2, tag="obs")
 formula <- as.formula(Y ~ -1+Intercept+abundance2+abundance3+pathogenconnectivity+numberofcoinfections+numberofstrains+pl+f(spatial, model = spde, group = spatial.group, control.group =list(model='ar1', hyper=h.spec)))
 resOverwintering<- inla(formula, family = 'binomial',Ntrials= data$Ntrials,  data = inla.stack.data(data), control.predictor =list(compute=TRUE,A=inla.stack.A(data)), control.fixed=list(expand.factor.strategy='inla'),  control.compute = list(dic=TRUE, waic = TRUE, cpo = TRUE), control.inla = list(int.strategy = "eb"), num.threads = 10)
 summary(resOverwintering)
 
 # SAVE ALL RESULTS: 
 
 save(file = 'finalresultsFungalSexPaper.Rdata', resNewNonextinctionPoisson, resCoPresenceNowPoisson, resCoPresenceNow, resOverwintering)
 
 
 
 
 source('fungalsexfunctions.R')

 resNewNonextinction <-resNewNonextinctionPoisson
 resCoPresenceNow <-resCoPresenceNow
 resOverwintering <- resOverwintering
 
 png('FigCPOCoinfection.png', width = 5, height =5, units = 'in', res = 300)
 hist(resCoPresenceNow$cpo$cpo, xlab = 'CPO', main = 'Coinfection presence')
 dev.off()
 png('FigCPOOverwintering.png', width = 5, height =5, units = 'in', res = 300)
 hist(resOverwintering$cpo$cpo, xlab = 'CPO', main = 'Overwintering')
 dev.off()
 png('FigCPONewstrain.png', width = 5, height =5, units = 'in', res = 300)
 hist(resNewNonextinction$cpo$cpo, xlab = 'CPO', main = 'New strain')
 dev.off()
 
 png('FigCPOs.png', width = 10, height =3.33, units = 'in', res = 300)
 par(mfrow=c(1,3))
 hist(resCoPresenceNow$cpo$cpo, xlab = 'CPO', main = 'Coinfection presence')
 hist(resOverwintering$cpo$cpo, xlab = 'CPO', main = 'Overwintering')
 hist(resNewNonextinction$cpo$cpo, xlab = 'CPO', main = 'New strain')
 dev.off()
 
 
 
 sumcopresencenow <- summary(resCoPresenceNow) # (presence now)
 sumnewnonext <- summary(resNewNonextinction)
 sumoverwintering <-summary(resOverwintering)
 
 sumcopresencenow$fixed
 sumnewnonext$fixed
 sumoverwintering$fixed
 
 resultfixedOverwintering = inla.spde2.result(resOverwintering, "spatial", spde)
 resultfixedCoinfPresence = inla.spde2.result(resCoPresenceNow, "spatial", spde)
 resultfixedNewNonextinction = inla.spde2.result(resNewNonextinctionPoisson, "spatial", spde)
 
 kappaNewBin <- getINLAResult(resultfixedNewNonextinction$marginals.kappa[[1]], coords.scale=1/size)
 tauNewBin <- getINLAResult(resultfixedNewNonextinction$marginals.tau[[1]])
 varianceNewBin <- round(getINLAResult(resultfixedNewNonextinction$marginals.variance.nominal[[1]]), 2)
 rangeNewBin <- round(111319.5*getINLAResult(resultfixedNewNonextinction$marginals.range[[1]]),0)
 
 kappaCoBin <- getINLAResult( resultfixedCoinfPresence$marginals.kappa[[1]], coords.scale=1/size)
 tauCoBin <- getINLAResult( resultfixedCoinfPresence$marginals.tau[[1]])
 varianceCoNewBin <- round(getINLAResult( resultfixedCoinfPresence$marginals.variance.nominal[[1]]), 2)
 rangeCoBin <- round(111319.5*getINLAResult( resultfixedCoinfPresence$marginals.range[[1]]),0)
 
 kappaOverwintering <- getINLAResult(resultfixedOverwintering$marginals.kappa[[1]], coords.scale=1/size)
 tauOverwintering <- getINLAResult(resultfixedOverwintering$marginals.tau[[1]])
 varianceOverwintering <- round(getINLAResult(resultfixedOverwintering$marginals.variance.nominal[[1]]), 2)
 rangeOverwintering <- round(111319.5*getINLAResult(resultfixedOverwintering$marginals.range[[1]]),0)

 
 # plot the estimated effects:
 
 dev.off()
 setwd("C:/HY-Data/SROTO/documents/fungalsex/New Manuscript Figures")
 png('FigEstimatesNewStrainModel.png', width = 6, height =7, units = 'in', res = 300)
 plot_effects(resNewNonextinction, 'New strain emergence')
 dev.off()
 png('FigEstimatesCoinfectionPresenceNow.png', width = 6, height =7, units = 'in', res = 300)
 plot_effects(resCoPresenceNow, 'Coinfections (conditional on presence now)')
 dev.off()
 png('FigEstimatesOverwintering.png', width = 6, height =7, units = 'in', res = 300)
 plot_effects(resOverwintering, 'Overwintering')
 dev.off()
 
 #
 # model validation plots:
 
 resNewNonextinction <- inla.cpo(resNewNonextinction)
 resOverwintering <- inla.cpo(resOverwintering)
 resCoPresenceNow <- inla.cpo(resCoPresenceNow)
 resChasmoPresenceNow <- inla.cpo(resChasmoPresenceNow)
 
 require(gridExtra)
 
 setwd("C:/HY-Data/SROTO/documents/fungalsex/New Manuscript Figures")
 inds <- which((d$numberofstrains>0) & (d$numberofstrainsnextyear>0) & (d$nsamples <15)) 
 data_norm2<-d[inds,]
 sp1_gps_matrix2 <- d[inds,c('x', 'y')]
 dev.off()
 
 png('FigPITNew.png', width = 6, height =7, units = 'in', res = 300)
 PlotPITcorrected(resNewNonextinction, 1*(data_norm2$newstrainsnextyear>0), 'New Strains')
 dev.off()
 png('FigCPOTImeNew.png',  width = 8, height =3,  units = 'in', res = 300)
 PlotCPOtime(resNewNonextinction, data_norm2$year, 'New Strain')
 dev.off()
 png('FigCPOSpaceNew.png', width = 24, height =6, units = 'in', res = 300)
 PlotCPOspace2(resNewNonextinction, sp1_gps_matrix2, data_norm2$year, 'New strains')
 dev.off()
 
 dev.off()
 inds <- which(d$numberofstrainsC>0 & (d$nsamples <15))
 data_norm2<-d[inds,]
 sp1_gps_matrix2 <- d[inds,c('x', 'y')]
 
 png('FigPITOverwintering.png', width = 6, height =7, units = 'in', res = 300)
 PlotPITcorrected(resOverwintering, 1*(data_norm2$newstrainsnextyear>0), 'Overwintering')
 dev.off()
 png('FigCPOTImeOverwintering.png',  width = 8, height =3,  units = 'in', res = 300)
 PlotCPOtime(resOverwintering, data_norm2$year, 'Overwintering')
 dev.off()
 png('FigCPOSpaceOverwintering.png', width = 24, height =6, units = 'in', res = 300)
 PlotCPOspace2(resOverwintering, sp1_gps_matrix2, data_norm2$year, 'Overwintering')
 dev.off()
 
 png('FigPITCoinfPresenceNow.png', width = 6, height =7, units = 'in', res = 300)
 PlotPITcorrected(resCoPresenceNow, 1*(data_norm2$coinfpresence>0), 'Coinfection (presence now)')
 dev.off()
 png('FigCPOTImeCoinfPresenceNow.png', width = 8, height =3, units = 'in', res = 300)
 PlotCPOtime(resCoPresenceNow, data_norm2$year, 'Coinfection (presence now)')
 dev.off()
 png('FigCPOSpaceCoinfPresenceNow.png', width = 24, height =6, units = 'in', res = 300)
 PlotCPOspace2(resCoPresenceNow, sp1_gps_matrix2, data_norm2$year, 'Coinfection (presence now)')
 dev.off()
 
 data_norm2$chasmopresent2[is.na(data_norm2$chasmopresent2)] <- 0 
 data_norm2$chasmopresent[is.na(data_norm2$chasmopresent)] <- 0 
 png('FigPITChasmo.png', width = 6, height =7, units = 'in', res = 300)
 PlotPITcorrected(resChasmoPresenceNowZip, 1*(data_norm2$chasmopresent>0), 'Chasmothecia')
 dev.off()
 png('FigCPOChasmo.png', width = 6, height =7, units = 'in', res = 300)
 hist(resChasmoPresenceNow$cpo$cpo,  main= 'Chasmothecia',xlab = 'CPO-value', ylab = 'Frequency')
 dev.off()
 png('FigCPOTimeChasmo.png', width = 8, height =3, units = 'in', res = 300)
 PlotCPOtime(resChasmoPresenceNow, data_norm2$year, 'Chasmothecia')
 dev.off()
 png('FigCPOSpaceChasmo.png', width = 24, height =6, units = 'in', res = 300)
 PlotCPOspace2(resChasmoPresenceNow, sp1_gps_matrix2, data_norm2$year, 'Chasmothecia')
 dev.off()
 

 # plot spatial ranges, make summary tables:
 fixedOverwintering <- round(resOverwintering$summary.fixed, 2)
 fixedCoinfPresenceNow <- round(resCoPresenceNow$summary.fixed, 2)
 fixedChasmoPresenceNow <- round(resChasmoPresenceNow$summary.fixed, 2)
 fixedNewNonextinction<- round(resNewNonextinction$summary.fixed, 2)

 randomOverwintering <- round(resOverwintering$summary.hyperpar, 2)[3, c(1, 3,5)]
 randomCoinfPresenceNow <- round(resCoPresenceNow$summary.hyperpar, 2)[3, c(1, 3,5)]
 randomChasmoPresenceNow <- round(resChasmoPresenceNow$summary.hyperpar, 2)[3, c(1, 3,5)]
 randomNewNonextinction <- round(resNewNonextinction$summary.hyperpar, 2)[3, c(1, 3,5)]
 
 # plot spatial ranges:

 x <- seq(0,by = 0.001, 0.5)
 y_range_overwinter<- inla.dmarginal(x, resultfixedOverwintering$marginals.range.nominal[[1]])
 y_range_coinfpresencenow<- inla.dmarginal(x, resultfixedCoinfPresence$marginals.range.nominal[[1]])
 y_range_chasmo<- inla.dmarginal(x, resultfixedChasmo$marginals.range.nominal[[1]])
 y_range_new<- inla.dmarginal(x, resultfixedNewNonextinction$marginals.range.nominal[[1]])
 
 
 dev.off()
 setwd("C:/HY-Data/SROTO/documents/fungalsex/New Manuscript Figures")
 png('FigSpatialRangesModels.png', height = 4, width = 5, units = 'in', res = 300)
 plot(111319.5*x,  y_range_overwinter, type = 'l' , col = 'gold', xlab = 'Range (m)', ylab = 'Posterior marginal density', ylim = c(0, 25), bty="n", lwd=2) 
 lines(111319.5*x,  y_range_new+0.1, type = 'l', col = 'green' ,  lwd=2)
 lines(111319.5*x,  y_range_coinfpresencenow, type = 'l' , col = 'blue', lwd=2)
 lines(111319.5*x,  y_range_chasmo, type = 'l' , col = 'magenta', lwd=2) 
 legend('topright', fill =c('gold','green', 'blue', 'magenta'),  legend = c('Overwintering','New strain',  'Coinfection', 'Chasmothecia'))
 
 dev.off()
 # SUMMARY TABLE:
 
 varianceOverwintering <- round(getINLAResult(resultfixedOverwintering$marginals.variance.nominal[[1]]), 2)
 rangeOverwintering <- round(111319.5*getINLAResult(resultfixedOverwintering$marginals.range[[1]]),0)
 
 varianceCOpresence <- round(getINLAResult(resultfixedCoinfPresence$marginals.variance.nominal[[1]]), 2)
 rangeCOpresence <- round(111319.5*getINLAResult(resultfixedCoinfPresence$marginals.range[[1]]), 0)
 
 varianceChasmo <- round(getINLAResult(resultfixedChasmo$marginals.variance.nominal[[1]]), 2)
 rangeChasmo<- round(111319.5*getINLAResult(resultfixedChasmo$marginals.range[[1]]), 0)
 
 varianceNewnonextinction <- round(getINLAResult(resultfixedNewNonextinction$marginals.variance.nominal[[1]]), 2)
 rangeNewnonextinction <- round(111319.5*getINLAResult(resultfixedNewNonextinction$marginals.range[[1]]), 0)
 

 # collect results into a table:

 meansModels <- data.frame( 'covariate' = row.names(fixedOverwintering), 'overwintering' = paste(fixedOverwintering$mean, paste('(', paste0(fixedOverwintering$`0.025quant`, ' ,',  fixedOverwintering$`0.975quant`), ')')) , 

                            'New' = paste(fixedNewNonextinction$mean, paste('(', paste0(fixedNewNonextinction$`0.025quant`, ' ,',  fixedNewNonextinction$`0.975quant`), ')')))
 
 
 meansModelsCo <- data.frame('covariate' = row.names(fixedCoinfPresenceNow), 'coninfpresence'=paste( fixedCoinfPresenceNow$mean, paste('(', paste0( fixedCoinfPresenceNow$`0.025quant`, ' ,',   fixedCoinfPresenceNow$`0.975quant`), ')')), 
 'chasmo'=paste( fixedChasmoPresenceNow$mean, paste('(', paste0( fixedChasmoPresenceNow$`0.025quant`, ' ,',  fixedChasmoPresenceNow$`0.975quant`), ')')) )
 
 meansModels
 meansModelsCo 
 
 varianceOverwintering[, c(1, 3, 5)]
 rangeOverwintering[, c(1, 3, 5)]
 
 varianceCOpresence[, c(1, 3, 5)]
 rangeCOpresence[, c(1, 3, 5)]
 
 varianceChasmo[, c(1, 3, 5)]
 rangeChasmo[, c(1, 3, 5)]
 
 varianceNewnonextinction[, c(1, 3, 5)]
 rangeNewnonextinction[, c(1, 3, 5)]

 

 
 # temporal autocorrelations: 
 
 summary(resOverwintering)
 summary(resNewNonextinction)
 summary(resCoPresenceNow)
 summary(resChasmoPresenceNow)
 
 
 
 ## the plots:
 
 
 load(file = 'sexmanuscriptfinaldataCoinfNum3years2.Rdata')
 ggplot(d, aes(numberofcoinfections, numberofstrainsnextyear))+geom_point()+geom_smooth(method = 'lm')+xlim(c(0,10))
 d$numberofstrainsnextyear[is.na(d$numberofstrainsnextyear)] <-0 
 
 corr_eqn <- function(x,y, digits = 2) {
   corr_coef <- round(cor(x, y), digits = digits)
   r2<-summary(lm(x ~ y))$r.squared
   paste(paste("italic(r) == ", corr_coef))
 }
 
 corr_eqn2 <- function(x,y, digits = 2) {
   corr_coef <- round(cor(x, y), digits = digits)
   r2<-round(summary(lm(x ~ y))$r.squared, digits = digits)
   paste(paste("R^2 == ", r2))
 }
 
 dinf <- d[which(d$numberofstrains>0), ]
 
 
 labels = data.frame(x = 8, y = 3, label = corr_eqn(dinf$numberofcoinfections, dinf$numberofstrainsnextyear))
 labels2 = data.frame(x = 8, y = 2.5, label = corr_eqn2(dinf$numberofcoinfections, dinf$numberofstrainsnextyear))

 setwd("C:/HY-Data/SROTO/documents/fungalsex/New Manuscript Figures")
 
 png('CorrelationPlot2.png', width = 3.5, height = 3.5, units = 'in', res = 300)
 ggplot(dinf, aes(numberofcoinfections, numberofstrainsnextyear))+geom_smooth(method ='lm',  col = 'darkgrey', fill= 'lightblue3')+theme_classic()+xlab('Number of coinfections')+ylab('New strains to patch the following year')+geom_text(data = labels, aes(x = x, y = y,label = label), parse = TRUE)+geom_text(data = labels2, aes(x = x, y = y,label = label), parse = TRUE)+xlim(c(0, 10))
 dev.off()
 

 png('CorrelationPlot3.png', width = 3.5, height = 3.5, units = 'in', res = 300)
 ggplot(dinf, aes(numberofcoinfections, numberofstrainsnextyear))+geom_jitter()+geom_smooth(method ='lm',  col = 'darkgrey', fill= 'lightblue3')+theme_classic()+xlab('Number of coinfections')+ylab('New strains to patch the following year')+geom_text(data = labels, aes(x = x, y = y,label = label), parse = TRUE)+geom_text(data = labels2, aes(x = x, y = y,label = label), parse = TRUE)+xlim(c(0, 10))
 dev.off()
 
 