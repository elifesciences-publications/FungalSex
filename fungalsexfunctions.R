PlotCPOspace2 = function(res, data, years, title) {
  
  cpo_values2 <- res$cpo$cpo
  d <- data.frame('years' = years, 'x' = data[, 1], 'y' = data[,2],  'CPO' = cpo_values2)
  a <- fortify(wrld_RDP)
  bordercolor <- 'wheat4'
  
  inds <- which(years==2012)
  
  ylab <- paste('CPO-values, ', title)
  
  plot1 <- ggplot(d[inds,], aes(x = x, y = y, col= CPO))+geom_point(alpha = 0.7, cex = 2)+ geom_polygon(inherit.aes = FALSE, data = a, aes( x = long, y = lat, group= group), colour=bordercolor,  fill = NA)+theme_classic()+ scale_color_gradientn(colours = topo.colors(10))+ggtitle('2012')+ylab(ylab)+theme(text = element_text(size=16))
  
  
  inds <- which(years==2013)
  
  plot2 <- ggplot(d[inds,], aes(x = x, y = y, col= CPO))+geom_point(alpha = 0.7, cex = 2)+ geom_polygon(inherit.aes = FALSE, data = a, aes( x = long, y = lat, group= group), colour=bordercolor,  fill = NA)+theme_classic()+ scale_color_gradientn(colours = topo.colors(10))+ggtitle('2013')+ylab('')+theme(text = element_text(size=16))
  
  inds <- which(years==2014)
  
  plot3 <- ggplot(d[inds,], aes(x = x, y = y, col= CPO))+geom_point(alpha = 0.7, cex = 2)+ geom_polygon(inherit.aes = FALSE, data = a, aes( x = long, y = lat, group= group), colour=bordercolor,  fill = NA)+theme_classic()+ scale_color_gradientn(colours = topo.colors(10))+ggtitle('2014')+ylab('')+theme(text = element_text(size=16))
  
  grid.arrange(plot1, plot2, plot3,  ncol=3)
  
}





pit.one <- function(u,x,Px,Pxm1){
  
  F_u <- ifelse(u <= Pxm1 , 0, pmin(1,(u-Pxm1)/(Px-Pxm1) ) )
  F_u[x==0] <- pmin(1,u/Px)[x==0] # needless!?
  if(u == 1){
    F_u <- 1
  }
  if(u == 0){
    F_u <- 0
  }
  #print(F_u)
  return(mean(F_u))
}


pit <- function(J=10, x, Px, Pxm1){
  F_u.bar <- sapply((0:J)/J,pit.one, x=x, Px=Px, Pxm1=Pxm1)
  f_j <- J*diff(F_u.bar)
  erg <- list(breaks=(0:J)/J,counts=f_j, density=f_j,mids=(0:(J-1))/J+diff((0:J)/J)/2,
              xname="PIT",equidist=TRUE)
  class(erg) <- "histogram"
  return(erg)
}


PlotPITcorrected = function(res, prediction, title) {
  #result  = inla.cpo(res)
  vpit <- res$cpo$pit
  ## compute Pxm1 
  vcpo <- res$cpo$cpo
  Pxm1help <- vpit - vcpo
  ## be sure to avoid negative PITs
  Pxm1 <- ifelse(Pxm1help<0,0,Pxm1help)
  plot(pit(J=100, x=prediction, Px=vpit, Pxm1=Pxm1), ylim=c(0,1.5), ylab="Relative frequency", xlab="PIT-value", main = title)
}


PlotCPOtime = function(res, years, title) {
  
  #res <- inla.cpo(res)
  cpo_values2 <- res$cpo$cpo
  years_temp <- years
  years_temp[years_temp==2012] <- 1
  years_temp[years_temp==2013] <- 2
  years_temp[years_temp==2014] <- 3
  
  
  par(mfrow = c(1,3))  
  inds <- which(years_temp==1)
  print(sum(is.na(cpo_values2[inds])))
  print(sum(is.nan(cpo_values2[inds])))
  print(sum(cpo_values2[inds]))
  print(cpo_values2[inds])
  hist(cpo_values2[inds], main = '2012', xlab = 'CPO', breaks = 30, ylab = title, cex = 2, cex.lab = 1.5)
  
  inds <- which(years_temp==2)
  hist(cpo_values2[inds], main = '2013', xlab = 'CPO', breaks = 30, ylab = '' , cex = 2, cex.lab = 1.5)
  
  inds <- which(years_temp==3)
  hist(cpo_values2[inds], main = '2014', xlab = 'CPO', breaks = 30, ylab = '' , cex = 2, cex.lab = 1.5)
  
}


plot_effects <- function(res, title){
  # dev.off()
  library(ggplot2)
  a <- res$summary.fixed
  names <- row.names(a)
  names[names=="1"] <-'Abundance 1'
  names[names=="2"] <- 'Abundance 2'
  names[names=="3"] <- 'Abundance 3+'
  names[names=="abundance2"] <-'Abundance 2'
  names[names=="abundance3"] <- 'Abundance 3'
  names[names=="con1d3"] <- 'Pathogen connectivity'
  names[names=="numberofcoinfections"] <- 'Number of coinfections'
  names[names=="pathogenconnectivity"] <- 'Pathogen connectivity'
  names[names=="numberofstrains"] <- 'Number of strains'
  names[names=="pl"] <- 'log-host coverage (m^2)'
  names[names=="Sh"] <- 'Host connectivity'
  names[names=="road_PA"] <- 'Road presence'
  names[names=="logplm2"] <- 'log-host coverage (m^2)'
  names[names=="Number.of.strains"] <-'Number of strains'
  names[names=="prop_coinfection_prev"] <- 'Proportion coinfected'  
  a$Predictor <- names
  a <- data.frame(a)
  a$sig <- 1*(a$X0.025quant*a$X0.975quant>0)
  b <- a[a$sig==1,]
  ggplot(a, aes(Predictor, mean))+geom_abline(slope =0, intercept = 0,col = 'lightgrey', size = 1)+geom_point(cex =3)+theme_classic()+geom_errorbar(aes(ymin = X0.025quant, ymax = X0.975quant), width = 0, size = 1)+ theme(axis.text.x=element_text(angle = -90, vjust = 0.5, hjust=0, size = 16), axis.title=element_text(size=16))+ylab('The estimated effect')+ggtitle(title)+ theme(plot.title = element_text(lineheight=.8, face="bold", size = 16))
}



reid = function(id) {
  newid <- integer(length(id))
  hash <- list()
  count <- 0
  for (i in 1:length(id)) {
    if (is.null(hash[[as.character(id[i])]])) {
      count <- count + 1
      hash[[as.character(id[i])]] <- count
      newid[i] <- count        
    }
    else {
      newid[i] <- hash[[as.character(id[i])]]
    }
  }
  return(newid)
}


d$reID <- reid(d$patch)

getINLAResult = function(marginal, fun=identity, coords.scale=1) {
  m <- inla.tmarginal(function(x) fun(x) * coords.scale, marginal)
  e <- inla.emarginal(function(x) x, m)
  e2 <- inla.emarginal(function(x) x^2, m)
  sd <- sqrt(e2-e^2)
  q <- inla.qmarginal(c(0.025, 0.5, 0.975), m)
  mode <- inla.mmarginal(m)
  x <- data.frame(e=e, sd=sd, q1=q[1], q2=q[2], q3=q[3], mode=mode)
  colnames(x) <- c("mean", "sd", "0.025quant","0.5quant","0.975quant", "mode")
  return(x)
}

plot_rf_prob_pred <- function(res, N){
  xmean <- list()
  projgrid <- inla.mesh.projector(mesh, dims = c(400, 400),xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]) )
  cols= rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))
  cols <- (topo.colors(200))
  par(mar=c(1.2,1.2,1.2,3))
  par(mfrow=c(1,N))
  
  for (i in seq(1, N)){
    xmean[[i]] <- inla.mesh.project(projgrid, res$summary.fitted.values$mode[which(spde.idx$spatial.group==i)], xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]))
    xmean[[i]][!xy.in] <- NA
    image.plot(projgrid$x, projgrid$y, xmean[[i]], col = cols,  xaxt='n', yaxt='n', bty = "n", legend.shrink = 0.65, 
               legend.width = 1, legend.mar = 10)
  }
}

plot_rf <- function(res, N){
  xmean <- list()
  meanintercept <- res$summary.fixed['Intercept',1]
  projgrid <- inla.mesh.projector(mesh, dims = c(400, 400),xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]) )
  cols= rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))
  cols <- (topo.colors(200))
  par(mar=c(1.2,1.2,1.2,3))
  par(mfrow=c(1,N))
  
  for (i in seq(1, N)){
    xmean[[i]] <- inla.mesh.project(projgrid, res$summary.random$spatial$mode[which(spde.idx$spatial.group==i)], xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]))
    xmean[[i]][!xy.in] <- NA
    image.plot(projgrid$x, projgrid$y, xmean[[i]], col = cols,  xaxt='n', yaxt='n', bty = "n", legend.shrink = 0.65, 
               legend.width = 1, legend.mar = 10)
  }
}


plot_rf_prob <- function(res, N){
  xmean <- list()
  meanintercept <- res$summary.fixed['Intercept',1]
  projgrid <- inla.mesh.projector(mesh, dims = c(400, 400),xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]) )
  cols= rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))
  cols <- (topo.colors(200))
  par(mar=c(1.2,1.2,1.2,3))
  par(mfrow=c(1,N))
  
  for (i in seq(1, N)){
    xmean[[i]] <- inla.mesh.project(projgrid, res$summary.random$spatial$mode[which(spde.idx$spatial.group==i)], xlim = range(mesh$loc[,1]), ylim = range(mesh$loc[,2]))
    xmean[[i]][!xy.in] <- NA
    image.plot(projgrid$x, projgrid$y, exp(xmean[[i]]+meanintercept)/(1+exp(xmean[[i]]+meanintercept)), col = cols,  xaxt='n', yaxt='n', bty = "n", legend.shrink = 0.65, 
               legend.width = 1, legend.mar = 10)
  }
}