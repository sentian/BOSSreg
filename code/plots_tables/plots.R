### This file reproduces all the figures in the paper

library(ggplot2)
library(grid)
library(gridExtra)
library(zipfR)
library(reshape2)

## Specify the directory to save the plots
base = getwd()
base_plots = paste0(base, '/plots_tables/plots')
## Create the directory to save the plots
dir.create(base_plots, recursive = TRUE, showWarnings = FALSE)

### Functions to be used throughout this file --------
source(paste0(base, '/utils.R'))

## Extract the legend of a ggplot object
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## edf of BS when the true model is null, based on the result in Ye(1998)
integrand <- function(y, r, p){
  tmp = sqrt(y) * exp(-y/2) * Igamma(0.5, y/2)^(r-1) * (sqrt(pi)-Igamma(0.5, y/2))^(p-r)
}
moment.ordered.chisq1 <- function(p){
  k = 1:p
  const = choose(p,k) * k / (sqrt(2) * sqrt(pi)^p)
  integral = Vectorize(function(k) integrate(Vectorize(integrand), lower=0, upper=Inf, r=k, p=p )$value)
  return(const*integral(k))
}
edf.bs.null <- function(p){
  return(c(0, cumsum(rev(moment.ordered.chisq1(p)))))
}

### Results to be used to produce Figure 1-3 --------
## Calculate Cp, AICc and KL values for BS fits on replicaitons of data under various true models
df.cp.aicc.bs <- function(){
  ## parameters
  nrep = 2000
  n = 200
  p = 14
  snr = c(NA, rep(c(7, 0.2), 2))
  type = c('Null', rep('Orth-Sparse-Ex1', 2), rep('Orth-Dense', 2))
  type_snr = c('Null', 'Sparse, hsnr', 'Sparse, lsnr', 'Dense, hsnr', 'Dense, lsnr')
  
  edf = hdf = cp_edf = cp_hdf = aicc_edf = aicc_hdf = Errhat_kl= list()
  
  for(i in 1:length(type_snr)){
    ## Generate the dataset (1000 replications of y by fixing X)
    data = gen.data.orthx(n, p, snr[i], type[i], nrep, center.y=FALSE)
    mu = data$x %*% data$beta
    
    ## hdf (null true model may throw a warning message)
    hdf[[type_snr[i]]] = BOSSreg:::calc.hdf(data$x, NULL, data$sigma, mu)$hdf
    
    ## Fit BS on all replications
    betahat = lapply(1:nrep, function(rep){bs.orthx(data$x, data$y[,rep])}) 
    muhat = lapply(betahat, function(xx){data$x %*% xx})
    rss = do.call(rbind, lapply(1:nrep, function(rep){Matrix::colSums( sweep(muhat[[rep]], 1, data$y[,rep], '-')^2 )}))
    loss = do.call(rbind, lapply(1:nrep, function(rep){Matrix::colSums( sweep(muhat[[rep]], 1, mu, '-')^2 )}))
    
    ## edf
    edf[[type_snr[i]]] = calc.edf(muhat, data$y, data$sigma)
    
    ## Cp, AICc and KL
    cp_hdf[[type_snr[i]]] = sweep(rss, 2, 2*data$sigma^2*hdf[[type_snr[i]]], '+')
    cp_edf[[type_snr[i]]] = sweep(rss, 2, 2*data$sigma^2*edf[[type_snr[i]]], '+')
    aicc_hdf[[type_snr[i]]] = sweep(n*log(rss/n), 2, n*( n + hdf[[type_snr[i]]] ) / ( n - hdf[[type_snr[i]]] - 2), '+')
    aicc_edf[[type_snr[i]]] = sweep(n*log(rss/n), 2, n*( n + edf[[type_snr[i]]] ) / ( n - edf[[type_snr[i]]] - 2), '+')
    bias_kl = n*(n*data$sigma^2 + loss) / rss
    Errhat_kl[[type_snr[i]]] = sweep(n*log(rss/n), 2, Matrix::colMeans(bias_kl) - 1, '+')
  }
  
  return(list(para = list(n=n, p=p, nrep=nrep, snr=snr, type=type),
              df = list(edf=edf, hdf=hdf),
              cp = list(edf=cp_edf, hdf=cp_hdf),
              aicc = list(edf=aicc_edf, hdf=aicc_hdf, kl=Errhat_kl)))
}
result_fig123 = df.cp.aicc.bs()

### Figure 1: hdf(k) and edf(k) for BS --------
plot.hdf.edf.bs <- function(result){
  n = result$para$n
  p = result$para$p
  edf = result$df$edf
  hdf = result$df$hdf
  
  for(type_snr in names(edf)){
    df_toplot = data.frame(k = c(0:p, 0:p), 
                           df = c(edf[[type_snr]], hdf[[type_snr]]), 
                           type = rep(c('edf', 'hdf'), each = p+1))

    p1 = ggplot() + geom_point(data=df_toplot, aes(y = df, x = k, colour = type, shape = type), stat="identity",size=3, stroke=2)
    p1 = p1 + scale_colour_manual(name  = "",
                                   breaks=c("edf", "hdf"),
                                   labels=expression(edf, hdf),
                                   values=c("#F8766D", "#00BA38")) +
        scale_shape_manual(name  = "",
                           values=c(5, 3),
                           breaks=c("edf", "hdf"),
                           labels=expression(edf, hdf))
    p1 = p1 + geom_segment(aes(x = 0, y = 0, xend = p, yend = p), linetype='dashed',colour='black')
    p1 = p1 + theme(legend.text=element_text(size=30))
    p1 = p1 + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
    p1 = p1 + ggtitle(type_snr)
    pp[[type_snr]] = p1
  }

  # print to a file
  mylegend = g_legend(pp$`Null`)
  setEPS()
  postscript(file=paste(base_plots, '/hdf_edf_bs.eps',sep=""), height=12, width=18)
  pp0 = grid.arrange(arrangeGrob(pp$`Sparse, hsnr` + theme(legend.position="none"),
                                  pp$`Dense, hsnr` + theme(legend.position="none"),
                                  pp$`Null` + theme(legend.position="none"),
                                  pp$`Sparse, lsnr` + theme(legend.position="none"),
                                  pp$`Dense, lsnr` + theme(legend.position="none"),
                                  mylegend,
                                  nrow=2,ncol=3,heights=c(3,3)))
  print(pp0)
  dev.off()
}
plot.hdf.edf.bs(result_fig123)

### Figure 2: Cp-edf and Cp-hdf for BS --------
plot.ic <- function(result, ic=c('cp', 'aicc'), filename){
  # Parameters
  n = result$para$n
  p = result$para$p
  nrep = result$para$nrep
  
  ic = match.arg(ic)
  
  ic_avg = lapply(result[[ic]], function(xx){lapply(xx, Matrix::colMeans)})
  k_avg = lapply(result[[ic]], function(xx){lapply(xx, function(yy){round(mean(apply(yy, 1, which.min)-1), 0)})})

  pp = list()
  for(type_snr in names(ic_avg$edf)){
    if(ic == 'cp'){
      df_toplot <- data.frame(k = rep(0:p, 2), 
                              bias = c(ic_avg$edf[[type_snr]] / n, ic_avg$hdf[[type_snr]] / n),
                              type = rep(c('edf', 'hdf'), each = p+1))
    }else if(ic == 'aicc'){
      df_toplot <- data.frame(k = rep(0:p, 3), 
                              bias = c(ic_avg$edf[[type_snr]] / n, ic_avg$hdf[[type_snr]] / n, ic_avg$kl[[type_snr]] / n),
                              type = rep(c('edf', 'hdf', 'kl'), each = p+1))
    }
    x_op = k_avg$edf[[type_snr]]
    
    p1 <- ggplot(data=df_toplot, aes(y = bias, x = k)) + geom_point(aes(color=type, shape=type), stat="identity",size=2, stroke=2)
    p1 <- p1 + geom_path(data=data.frame(x=c(x_op, x_op), y=c(min(df_toplot$bias), max(df_toplot$bias))),aes(x=x, y=y), linetype=2, colour='black', size=2)
    if(grepl('boss', filename)){
      p1 = p1 + xlab(expression(k[Q]))
    }
    if(ic == 'cp'){
      p1 <- p1 + scale_colour_manual(name  = "",
                                     breaks=c("edf", "hdf"),
                                     labels=expression('Cp-edf'/n, 'Cp-hdf'/n),
                                     values=c("#F8766D", "#00BA38")) +
        scale_shape_manual(name  = "",
                           breaks=c("edf", "hdf"),
                           labels=expression('Cp-edf'/n, 'Cp-hdf'/n),
                           values=c(5, 3))
    }else if(ic == 'aicc'){
      p1 <- p1 + scale_colour_manual(name  = "",
                                     breaks=c("edf", "hdf", "kl"),
                                     labels=expression('AICc-edf'/n, 'AICc-hdf'/n, widehat(Err)[KL]/n),
                                     values=c("#F8766D", "#00BA38", "#619CFF")) +
        scale_shape_manual(name  = "",
                           breaks=c("edf", "hdf", "kl"),
                           labels=expression('AICc-edf'/n, 'AICc-hdf'/n, widehat(Err)[KL]/n),
                           values=c(5, 3, 1))
    }
    
    p1 <- p1 + theme(legend.text=element_text(size=30)) + theme(axis.title.y=element_blank())
    p1 <- p1 + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
    p1 <- p1 + ggtitle(type_snr)
    pp[[type_snr]] = p1
  }
  mylegend = g_legend(pp$`Null`)

  setEPS()
  postscript(file=paste(base_plots, filename, sep="/"), height=12, width=18)
  pp0 <- grid.arrange(arrangeGrob(pp$`Sparse, hsnr` + theme(legend.position="none"),
                                  pp$`Dense, hsnr` + theme(legend.position="none"),
                                  pp$Null + theme(legend.position="none"),
                                  pp$`Sparse, lsnr` + theme(legend.position="none"),
                                  pp$`Dense, lsnr` + theme(legend.position="none"),
                                  mylegend,
                                  nrow=2,ncol=3,heights=c(3,3)))
  print(pp0)
  dev.off()
}
plot.ic(result_fig123, filename = 'cp_edf_hdf_bs.eps')

### Figure 3: AICc-edf and AICc-hdf for BS --------
plot.ic(result_fig123, ic = 'aicc', filename = 'aicc_edf_hdf_kl_bs.eps')

### Figure 4: 
### Figure 4: Frequency distributions of the selected subset size for BS and LBS --------
plot.freqdist.bs.lbs <- function(){
  type = 'Orth-Sparse-Ex1'
  n = 200
  snr = 7
  nrep = 1000
  
  data_toplot = data.frame()
  for(p in c(30,180)){
    data = gen.data.orthx(n, p, snr, type)
    result_cp = eval.metrics.bs.lbs.cp.orthx(data)
    # Frequency distributions
    numvar_cp = lapply(result_cp, function(xx){xx$sparsistency + xx$extravariable})
    tmp_function <- function(xx){
      xx[xx>=9] = ">=9"
      xx
    }
    numvar_cp = lapply(numvar_cp, tmp_function)
    numvar_cp = lapply(numvar_cp, function(xx){factor(as.character(xx), levels=c('6','7','8','>=9'))})
  
    tmp = melt(table('bs' = numvar_cp$bs, 'lbs' = numvar_cp$lbs) / 10)
    tmp$type = paste0('p=',p)
    data_toplot = rbind(data_toplot, tmp)
  }
  data_toplot$type = factor(data_toplot$type, levels=c('p=30','p=180'))
  # Make the plot
  p1 = ggplot(data = data_toplot, mapping = aes(x=bs, y=lbs)) + geom_tile(aes(fill = value)) +
    geom_text(aes(label = value), size = 3)
  p1 = p1 + scale_fill_gradient(low = "white", high = "red")
  p1 = p1 + facet_wrap(.~ type, ncol=2)
  p1 = p1 + xlab(label='Selected subset size, BS') + ylab(label='Selected subset size, LBS')
  p1 = p1 + theme(strip.text = element_text(size=10)) + theme(legend.position="none")
  p1 = p1 + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
  setEPS()
  postscript(file=paste(base_plots, '/numvar_bs_lbs.eps',sep=""), height=3, width=6)
  print(p1)
  dev.off()
}
plot.freqdist.bs.lbs()

### Figure 5: Cp-edf and Cp-hdf for BOSS --------
## Calculate Cp values for BOSS fits on replicaitons of data under various true models
## May throw warnings of returning df of null model
cp.boss <- function(){
  # parameters
  nrep = 1000
  n = 200
  p = 14
  rho = 0.5
  snr = c(NA, rep(c(7, 0.2), 2))
  type = c('Null', rep('Sparse-Ex3', 2), rep('Dense', 2))
  type_snr = c('Null', 'Sparse, hsnr', 'Sparse, lsnr', 'Dense, hsnr', 'Dense, lsnr')
  
  cp_edf = cp_hdf = list()

  for(i in 1:length(type_snr)){
    data = gen.data.generalx(n, p, rho, snr[i], type[i], nrep)
    mu = data$x %*% data$beta
    # calculate the fit, bias and IC
    betahat = cp_hdf_tmp = list()
    for(rep in 1:nrep){
      boss_model = boss(data$x, data$y[,rep], intercept = FALSE, mu = mu, sigma = data$sigma)
      cp_hdf_tmp[[rep]] = as.numeric( boss_model$IC_boss$cp )
      betahat[[rep]] = boss_model$beta_boss
    }
    muhat = lapply(betahat, function(xx){data$x%*%xx})
    rss = Map(function(xx, jj){Matrix::colSums(sweep(xx,1,data$y[,jj])^2) }, muhat, 1:nrep)
    rss = do.call(rbind, rss)

    edf = calc.edf(muhat, data$y, data$sigma)
    cp_edf[[type_snr[i]]] = sweep(rss, 2, 2*data$sigma^2*edf, '+')
    cp_hdf[[type_snr[i]]] = do.call(rbind, cp_hdf_tmp)
  }
  return(list(para = list(n=n, p=p, nrep=nrep, snr=snr, rho=rho, type=type),
              cp = list(edf=cp_edf, hdf=cp_hdf)))
}
result_fig5 = cp.boss()
## Make the plot
plot.ic(result_fig5, filename = 'cp_edf_hdf_boss.eps')

### Figure 6: RMSE along the solution paths of BS, FS and BOSS --------
plot.solpath.lsmethods <- function(){
  # Parameters
  n = 200
  p = 30
  rho = 0.9
  snr = 7
  nrep = 1000
  type = c('Sparse-Ex3', 'Sparse-Ex4')

  pp = list()
  for(i in 1:length(type)){
    data = gen.data.generalx(n, p, rho, snr, type[i])
    mu = data$x %*% data$beta
    # Fit the models
    betahat = replicate(3, list(), simplify = FALSE)
    names(betahat) = c('boss', 'bs', 'fs')
    for(rep in 1:nrep){
      boss_model = boss(data$x, data$y[,rep], intercept = FALSE, hdf.ic.boss = FALSE)
      betahat[['boss']][[rep]] = boss_model$beta_boss
      betahat[['fs']][[rep]] = boss_model$beta_fs
      betahat[['bs']][[rep]] = bs.generalx(data$x, data$y[,rep], intercept = FALSE)
    }

    rmse = lapply(betahat, function(xx){lapply(xx, function(yy){ sqrt(colSums(sweep(data$x%*%yy,1,mu,'-')^2)/n) })})
    rmse_avg = lapply(rmse, function(xx){colMeans(do.call('rbind', xx))})

    # Make the plot
    k_start = 0
    k_end = p
    df_toplot <- data.frame(k=rep(k_start:k_end,3), rmse=c(rmse_avg$boss[(k_start+1):(k_end+1)],rmse_avg$fs[(k_start+1):(k_end+1)],rmse_avg$bs[(k_start+1):(k_end+1)]),
                            type=c(rep('BOSS',k_end-k_start+1), rep('FS',k_end-k_start+1), rep('BS',k_end-k_start+1)))
    df_toplot$type = factor(df_toplot$type, levels = c('BOSS', 'BS', 'FS'))

    p1 <- ggplot() + geom_point(data=df_toplot, aes(y = rmse, x = k, colour = type, shape = type), stat="identity",size=1.5)
    p1 <- p1 + xlab('Subset size') + ylab('RMSE')
    p1 <- p1 + scale_colour_manual(name  = "",
                                   breaks=c("BOSS", "BS", "FS"),
                                   labels=c("BOSS", "BS", "FS"),
                                   values=c("#F8766D", "#00BA38", "#619CFF")) +
      scale_shape_manual(name  = "",
                         breaks=c("BOSS", "BS", "FS"),
                         labels=c("BOSS", "BS", "FS"),
                         values=c(5, 3, 1))

    p1 <- p1 + theme(legend.text=element_text(size=5),legend.position="bottom") + theme(legend.title=element_blank())
    p1 <- p1 + theme(plot.title = element_text(size = 10, face = "bold", hjust=0.5),axis.text=element_text(size=10), axis.title=element_text(size=10,face="bold"))
    p1 <- p1 + ggtitle(type[i])
    pp[[type[i]]] = p1
  }

  mylegend = g_legend(pp$`Sparse-Ex3`)
  setEPS()
  postscript(file=paste(base_plots, '/rmse_solpath_lsmethods.eps',sep=""), height=3, width=6)
  pp0 <- grid.arrange(arrangeGrob(pp$`Sparse-Ex3` + theme(legend.position="none"),
                                  pp$`Sparse-Ex4` + theme(legend.position="none"),
                                  nrow=1,ncol=2,heights=3),
                      mylegend,nrow=2,heights=c(9,1))
  print(pp0)
  dev.off()
}
plot.solpath.lsmethods()
