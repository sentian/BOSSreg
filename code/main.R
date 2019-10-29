# plots/tables in the report
library(ggplot2)
library(grid)
library(gridExtra)
library(latex2exp)
library(xtable)
library(patternplot)
library(orderstats)
library(zipfR)
library(boss)
library(leaps)
library(reshape2)
#library(pBrackets)

#base = "/media/haha0542/HDD/Dropbox/Sen/Research/Model_selection"
base = "/Volumes/HDD/Dropbox/Sen/Research/Model_selection"
base = "/Users/sentian/Dropbox/Sen/Research/Model_selection"

### Functions --------
## extract the legend out of boxplot
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

## edf of BS when the true model is null
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

## hdf of BS when the true model is null
hdf.bs.null <- function(p){
  return(c(0, 1:p + 2*p*qnorm(1-1:p/(2*p))*dnorm(qnorm(1-1:p/(2*p)))) )
}

## hdf of BS in general
hdf.bs <- function(Q, y, sigma=NULL, mu=NULL){
  n = dim(Q)[1]
  p = dim(Q)[2]
  if(p>=n){
    stop('hdf is undefined when p>=n')
  }
  if(xor(is.null(sigma), is.null(mu))){
    stop('either sigma or beta is specified, need both or none')
  }

  # if mu and sigma are not specified, use the full multiple regression
  if(is.null(sigma)){
    beta_hat = xtmu = t(Q) %*% y # the multiple regression coef
    resid = y - Q %*% beta_hat
    sigma = sqrt(sum(resid^2)/(n-p))
    xtmu_matrix = matrix(rep(xtmu,each=p-1), ncol=p-1, byrow=T)
  }else{
    xtmu = t(Q)%*%mu
    xtmu_matrix = matrix(rep(xtmu,each=p-1), ncol=p-1, byrow=T)
  }
  tryCatch({
    # calculate the inverse function of E(k(lambda))=k, where k=1,...p-1
    inverse = function(f, lower, upper) {
      function(y) stats::uniroot(function(x){f(x) - y}, lower=lower, upper=upper)[1]
    }
    exp_size <- function(x){
      c = stats::pnorm((x-xtmu) / sigma)
      d = stats::pnorm((-x-xtmu) / sigma)
      return( sum(1 - c + d) )
    }
    inverse_exp_size = inverse(exp_size, 0, 100*max(abs(xtmu)))
    sqrt_2lambda = unlist(lapply(1:(p-1), inverse_exp_size))
    sqrt_2lambda_matrix = matrix(rep(sqrt_2lambda,each=p), nrow=p, byrow=F)

    # plug the sequence of lambda into the expression of df(lambda)
    a = stats::dnorm((sqrt_2lambda_matrix-xtmu_matrix) / sigma)
    b = stats::dnorm((-sqrt_2lambda_matrix-xtmu_matrix) / sigma)

    size = 1:(p-1)
    sdf = (sqrt_2lambda/sigma) * colSums(a + b)
    df = size + sdf
    names(df) = NULL
    return(list(hdf=c(0, df, p), sigma=sigma))
  }, error=function(e){
    warning('returns the df for a null model')
    return(list(hdf=c(0, 1:p + 2*p*stats::qnorm(1-1:p/(2*p))*stats::dnorm(stats::qnorm(1-1:p/(2*p)))), sigma=sigma))
  })
}

## BS on orthogonal X
bs.orthx <- function(x, y){
  p = dim(x)[2]
  z = t(x) %*% y
  tmp = order(-z^2)
  coef_bs_all = matrix(0, nrow=p, ncol=p+1)
  for(j in 1:p){
    coef_bs_all[tmp[1:j],j+1] = z[tmp[1:j]]
  }
  return(coef_bs_all)
}

## BS on general X
leaps.fullsample <- function(x,y,intercept=TRUE){
  p = dim(x)[2]
  n = dim(x)[1]
  # fit to the full data
  outs = leaps(x=x, y=y, nbest=1,strictly.compatible=FALSE, int=intercept) # leaps
  if(intercept){
    ass_coef <- matrix(0,nrow=p+1,ncol=p+1)
    ass_coef[1,1] = mean(y)
    for(i in 1:nrow(outs$which)){
      ass_coef[c(1,which(outs$which[i,])+1),i+1] = coef(lm(y~x[,outs$which[i,]]))
    }
  }else{
    ass_coef <- matrix(0,nrow=p,ncol=p+1)
    for(i in 1:nrow(outs$which)){
      ass_coef[outs$which[i,],i+1] = coef(lm(y~x[,outs$which[i,]]-1))
    }
  }

  return(ass_coef)
}

### Figure 1: hdf(k) and edf(k) for BS --------
plot.hdf.edf.bs <- function(){
  # parameters
  snr = c('hsnr','lsnr')
  n = 200
  p = 14

  edf = hdf = pp = list()

  # null true model
  hdf[['Null']] = hdf.bs.null(p)
  edf[['Null']] = edf.bs.null(p)

  # sparse and dense true models
  categ_types = c('true_model/orthogonal/sparse/ex1/p0_6', 'true_model/orthogonal/dense/beta_diminishingstrength/kappa_10')
  categ_names = c('Sparse', 'Dense')
  for(i in 1:length(snr)){
    for(j in 1:length(categ_types)){
      indir = paste(categ_types[j],'/n_',n,'/p_',p, sep='')
      x = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/x.txt',sep=""))))
      beta = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/beta.txt',sep=""))))
      sigma = c(unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/sd_',snr[i],'.txt',sep="")))))
      df = readRDS(paste(base, '/code/as_gram/results/',indir,'/',snr[i],'/df_bs.rds',sep=""))
      edf[[paste(categ_names[j], snr[i], sep=', ')]] = df$bestsub$edf
      hdf[[paste(categ_names[j], snr[i], sep=', ')]] = hdf.bs(x, NULL, sigma, x%*%beta)$hdf
    }
  }

  # make the plot
  for(case in names(edf)){
    df_toplot = data.frame(k=c(0:p,0:p), df=c(edf[[case]],hdf[[case]]), type=c(rep('edf',p+1), rep('hdf',p+1)))

    p1 = ggplot() + geom_point(data=df_toplot, aes(y = df, x = k, colour = type, shape = type), stat="identity",size=3,stroke=2)
    p1 = p1 + scale_colour_manual(name  = "",
                                   breaks=c("edf", "hdf"),
                                   labels=expression(edf, hdf),
                                   values=c("#F8766D", "#00BA38")) +
        scale_shape_manual(name  = "",
                           values=c(5, 3),
                           breaks=c("edf", "hdf"),
                           labels=expression(edf, hdf))

    #labels=expression(df[C](k),hdf(k))

    p1 = p1 + geom_segment(aes(x = 0, y = 0, xend = p, yend = p), linetype='dashed',colour='black')
    p1 = p1 + theme(legend.text=element_text(size=30),legend.position="bottom")
    p1 = p1 + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5),axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
    p1 = p1 + ggtitle(case)
    pp[[paste(case, '-df_k', sep='')]] = p1
  }

  # print to a file
  mylegend = g_legend(pp$`Null-df_k`)
  setEPS()
  postscript(file=paste(base, '/paper/figures/dfc_dflambda.eps',sep=""), height=12, width=18)
  pp0 = grid.arrange(arrangeGrob(pp$`Sparse, hsnr-df_k` + theme(legend.position="none"),
                                  pp$`Dense, hsnr-df_k` + theme(legend.position="none"),
                                  pp$`Null-df_k` + theme(legend.position="none"),
                                  pp$`Sparse, lsnr-df_k` + theme(legend.position="none"),
                                  pp$`Dense, lsnr-df_k` + theme(legend.position="none"),
                                  mylegend,
                                  nrow=2,ncol=3,heights=c(3,3)))
  print(pp0)
  dev.off()
}
plot.hdf.edf.bs()

### Figure 2: Cp-edf and Cp-hdf for BS --------
# calculate Cp and AICc values for each replication
cp.aicc.bs <- function(){
  # parameters
  n = 200
  p = 14
  nrep = 100000

  cp_edf = cp_hdf = aicc_edf = aicc_hdf = Errhat_kl= list()

  # parameters to read data
  type = c('Null', 'Sparse, hsnr', 'Sparse, lsnr', 'Dense, hsnr', 'Dense, lsnr')
  categ_type = c('true_model/orthogonal/sparse/ex1/p0_6', rep(c('true_model/orthogonal/sparse/ex1/p0_6','true_model/orthogonal/dense/beta_diminishingstrength/kappa_10'), each=2))
  snr = c('na', rep(c('hsnr', 'lsnr'), 2))


  # i = 2
  # indir = paste(categ_type[i],'/n_',n,'/p_',p, sep='')
  # x = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/x.txt',sep=""))))
  # sigma = c(unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/sd_',snr[i],'.txt',sep="")))))
  # # y = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/y_',snr[i],'.txt',sep=""))))
  # beta = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/beta.txt',sep=""))))
  # mu = x%*%beta
  #
  # tmp = c()
  # for(j in 1:1){
  #   y = matrix(rep(mu,each=nrep),ncol=nrep,byrow=TRUE) + scale(matrix(rnorm(n*nrep,mean=0,sd=sigma),nrow=n,ncol=nrep),center=TRUE,scale=FALSE)
  #
  #   betahat = t(x) %*% y
  #   muhat = x %*% betahat
  #   rss = colSums((y - muhat)^2)
  #
  #   tmp = c(tmp, mean(n^2 * sigma^2 / rss) + mean( n*(colSums(sweep(muhat, 1, mu, '-')^2)) / rss))
  # }

  # tmp = c()
  # for(rep in seq(1000,nrep,length.out=1000)){
  #   tmp = c(tmp, mean(n^2 * sigma^2 / rss[1:rep]) + mean( n*(colSums(sweep(muhat[,1:rep], 1, mu, '-')^2)) / rss[1:rep] ))
  # }


  # (tmp1 = mean(n^2 * sigma^2 / rss) + mean( n*(colSums(sweep(muhat, 1, mu, '-')^2)) / rss ))
  # (tmp2 = n^2 / (n-p-2) + n*p/(n-p-2))
  #
  #
  # (mean(n*log(rss/n)) + tmp1) / n
  # (mean(n*log(rss/n)) + tmp2) / n

  # orderchisq_rn = matrix(NA, nrow=100000, ncol=p)
  # for(j in 1:p){
  #   orderchisq_rn[,j] = order_rchisq(100000, 1, p-j+1, p)
  # }
  # chisq_rn = rchisq(100000, n-p, ncp = 0)
  #
  # tmp1 = rowSums(as.matrix(orderchisq_rn[,1:p]))
  # tmp2 = chisq_rn
  # biasAB = mean( tmp1 / tmp2 )
  # bias1B = n * mean( 1 / tmp2 )
  #
  # (bias = n*(biasAB + bias1B))



  for(i in 1:length(type)){
    # read data
    indir = paste(categ_type[i],'/n_',n,'/p_',p, sep='')
    x = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/x.txt',sep=""))))
    if(type[i] == 'Null'){
      sigma = 1
      y = scale(matrix(rnorm(n*nrep,mean=0,sd=sigma),nrow=n,ncol=nrep),center=TRUE,scale=FALSE)
      mu = rep(0, n)
      # degrees of freedom
      hdf = hdf.bs.null(p)
      edf = edf.bs.null(p)

      orderchisq_rn = matrix(NA, nrow=100000, ncol=p)
      for(j in 1:p){
        orderchisq_rn[,j] = order_rchisq(100000, 1, p-j+1, p)
      }
      chisq_rn = rchisq(100000, n-p, ncp = 0)
      bias1B = biasAB = rep(NA, p+1)
      for(j in 0:p){
        if(j %in% c(0,p)){
          biasAB[j+1] = j/(n-j-2)
          bias1B[j+1] = n/(n-j-2)
        }else{
          tmp1 = rowSums(as.matrix(orderchisq_rn[,1:j]))
          tmp2 = chisq_rn + rowSums(as.matrix(orderchisq_rn[,(j+1):p]))
          biasAB[j+1] = mean( tmp1 / tmp2 )
          bias1B[j+1] = n * mean( 1 / tmp2 )
        }
      }
      bias = n*(biasAB + bias1B)

    }else{
      sigma = c(unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/sd_',snr[i],'.txt',sep="")))))
      y = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/y_',snr[i],'.txt',sep=""))))
      beta = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/beta.txt',sep=""))))
      mu = x%*%beta
      set.seed(11)
      y = cbind(y, matrix(rep(mu,each=nrep-1000),ncol=nrep-1000,byrow=TRUE) + scale(matrix(rnorm(n*(nrep-1000),mean=0,sd=sigma),nrow=n,ncol=nrep-1000),center=TRUE,scale=FALSE))

      # degrees of freedom
      hdf = hdf.bs(x, NULL, sigma, mu)$hdf
      result_df = readRDS(paste(base, '/code/as_gram/results/',indir,'/',snr[i],'/df_bs.rds',sep=""))
      edf = result_df$bestsub$edf
    }

    # calculate the fit, bias and IC
    rss = matrix(NA, nrow=nrep, ncol=p+1)
    if(type[i] != 'Null'){
      bias = matrix(NA, nrow=nrep, ncol=p+1)
    }
    for(rep in 1:nrep){
      betahat = bs.orthx(x, y[,rep])
      muhat = x%*%betahat
      rss[rep,] = colSums( sweep(muhat, 1, y[,rep], '-')^2 )
      if(type[i] != 'Null'){
        bias[rep,] = n*(n*sigma^2 + colSums( sweep(muhat, 1, mu, '-')^2 )) / rss[rep,]
      }
    }
    cp_hdf[[type[i]]] = sweep(rss, 2, 2*sigma^2*hdf, '+')
    cp_edf[[type[i]]] = sweep(rss, 2, 2*sigma^2*edf, '+')
    aicc_hdf[[type[i]]] = sweep(n*log(rss/n), 2, n*( n + hdf ) / ( n - hdf - 2), '+')
    aicc_edf[[type[i]]] = sweep(n*log(rss/n), 2, n*( n + edf ) / ( n - edf - 2), '+')
    if(type[i] == 'Null'){
      Errhat_kl[[type[i]]] = sweep(n*log(rss/n), 2, bias, '+')
    }else{
      Errhat_kl[[type[i]]] = sweep(n*log(rss/n), 2, colMeans(bias) - 1, '+')
    }
  }

  return(list(cp=list(edf=cp_edf, hdf=cp_hdf),
              aicc=list(edf=aicc_edf, hdf=aicc_hdf, kl=Errhat_kl)))

}
result = cp.aicc.bs()

# average over replications
cp_result_avg = lapply(result$cp, function(xx){lapply(xx, colMeans)})
# average optimal subset size by minimizing Cp in each replication
cp_result_min = lapply(result$cp, function(xx){lapply(xx, function(yy){round(mean(apply(yy, 1, which.min)-1),0)})})


# plot the results
plot.cp.edf.hdf <- function(result_avg, result_min, filename){
  n = 200
  p = 14
  pp = list()
  for(case in names(result_avg$edf)){
    df_toplot <- data.frame(k=rep(0:p, 2), bias=c(result_avg$edf[[case]] / n,
                                                  result_avg$hdf[[case]] / n),
                           type=c(rep('edf',p+1), rep('hdf', p+1)))

    p1 <- ggplot(data=df_toplot, aes(y = bias, x = k)) + geom_point(aes(color=type, shape=type), stat="identity",size=2, stroke=2)

    x_op = result_min$edf[[case]]

    p1 <- p1 + geom_path(data=data.frame(x=c(x_op, x_op), y=c(min(df_toplot$bias), max(df_toplot$bias))),aes(x=x, y=y), linetype=2, colour='black', size=2)
    if(filename == 'boss_cp_edf_hdf.eps'){
      p1 = p1 + xlab(expression(k[Q]))
    }
    p1 <- p1 + scale_colour_manual(name  = "",
                                     breaks=c("edf", "hdf"),
                                     labels=expression(E('Cp-edf')/n, E('Cp-hdf')/n),
                                     values=c("#F8766D", "#00BA38")) +
      scale_shape_manual(name  = "",
                         breaks=c("edf", "hdf"),
                         labels=expression(E('Cp-edf')/n, E('Cp-hdf')/n),
                         values=c(5, 3))
    p1 <- p1 + theme(legend.text=element_text(size=30)) + theme(axis.title.y=element_blank())
    p1 <- p1 + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
    p1 <- p1 + ggtitle(case)
    pp[[case]] = p1
  }
  mylegend = g_legend(pp$`Null`)

  setEPS()
  postscript(file=paste(base, '/paper/figures/', filename, sep=""), height=12, width=18)
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
plot.cp.edf.hdf(cp_result_avg, cp_result_min, 'cp_edf_hdf.eps')

### Figure 3: AICc-edf and AICc-hdf for BS --------
# average over replications
aicc_result_avg = lapply(result$aicc, function(xx){lapply(xx, colMeans)})
# average optimal subset size by minimizing Cp in each replication
aicc_result_min = lapply(result$aicc, function(xx){lapply(xx, function(yy){round(mean(apply(yy, 1, which.min)-1),0)})})

# plot the results
plot.aicc.edf.hdf.kl <- function(result_avg, result_min){
  n = 200
  p = 14
  pp = list()
  for(case in names(result_avg$edf)){
    df_toplot <- data.frame(k=rep(0:p, 3), bias=c(result_avg$edf[[case]] / n,
                                                  result_avg$hdf[[case]] / n,
                                                  result_avg$kl[[case]] / n),
                            type=c(rep('edf',p+1), rep('hdf', p+1), rep('kl', p+1)))

    p1 <- ggplot(data=df_toplot, aes(y = bias, x = k)) + geom_point(aes(color=type, shape=type), stat="identity",size=2, stroke=2)

    x_op = result_min$edf[[case]]

    p1 <- p1 + geom_path(data=data.frame(x=c(x_op, x_op), y=c(min(df_toplot$bias), max(df_toplot$bias))),aes(x=x, y=y), linetype=2, colour='black', size=2)
    p1 <- p1 + scale_colour_manual(name  = "",
                                     breaks=c("edf", "hdf", "kl"),
                                     labels=expression(E('AICc-edf')/n, E('AICc-hdf')/n, E(widehat(Err)[KL])/n),
                                     values=c("#F8766D", "#00BA38", "#619CFF")) +
      scale_shape_manual(name  = "",
                         breaks=c("edf", "hdf", "kl"),
                         labels=expression(E('AICc-edf')/n, E('AICc-hdf')/n, E(widehat(Err)[KL])/n),
                         values=c(5, 3, 1))
    p1 <- p1 + theme(legend.text=element_text(size=30)) + theme(axis.title.y=element_blank())
    p1 <- p1 + theme(plot.title = element_text(size = 30, face = "bold", hjust=0.5), axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"))
    p1 <- p1 + ggtitle(case)
    pp[[case]] = p1
  }
  mylegend = g_legend(pp$Null)

  setEPS()
  postscript(file=paste(base, '/paper/figures/aicc_edf_eo_withrss.eps',sep=""), height=12, width=18)
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
plot.aicc.edf.hdf.kl(aicc_result_avg, aicc_result_min)


### Figure 4: Cp-edf and Cp-hdf for BOSS --------
cp.aicc.boss <- function(){
  # parameters
  n = 200
  p = 14
  nrep = 1000

  cp_edf = cp_hdf = list()

  # parameters to read data
  type = c('Null', 'Sparse, hsnr', 'Sparse, lsnr', 'Dense, hsnr', 'Dense, lsnr')
  categ_type = c('true_model/sparse/ex3/p0_6/rho_05', rep(c('true_model/sparse/ex3/p0_6/rho_05','true_model/dense/beta_diminishingstrength/kappa_10/rho_05'), each=2))
  snr = c('na', rep(c('hsnr', 'lsnr'), 2))

  for(i in 1:length(type)){
    # read data
    indir = paste(categ_type[i],'/n_',n,'/p_',p, sep='')
    x = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/x.txt',sep=""))))
    if(type[i] == 'Null'){
      sigma = 1
      y = scale(matrix(rnorm(n*nrep,mean=0,sd=sigma),nrow=n,ncol=nrep),center=TRUE,scale=FALSE)
      mu = rep(0, n)
    }else{
      sigma = c(unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/sd_',snr[i],'.txt',sep="")))))
      y = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/y_',snr[i],'.txt',sep=""))))
      beta = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',indir,'/beta.txt',sep=""))))
      mu = x%*%beta
    }

    # calculate the fit, bias and IC
    betahat = cp_hdf_tmp = list()
    for(rep in 1:nrep){
      boss_model = boss(x, y[,rep], intercept = FALSE, mu = mu, sigma = sigma)
      cp_hdf_tmp[[rep]] = as.numeric( boss_model$IC_boss$cp )
      betahat[[rep]] = boss_model$beta_boss
    }
    muhat = lapply(betahat, function(xx){x%*%xx})
    rss = Map(function(xx, jj){Matrix::colSums(sweep(xx,1,y[,jj])^2) }, muhat, 1:nrep)
    rss = do.call(rbind, rss)

    edf = c()
    for(k in 1:(p+1)){
      tmp = do.call(rbind, lapply(muhat, function(x){x[,k]}))
      edf = c(edf,
              sum(unlist(Map(function(xx,yy){cov(xx,yy)}, split(tmp, rep(1:ncol(tmp), each=nrow(tmp))), split(y, rep(1:nrow(y), ncol(y)))))) / sigma^2)
    }
    cp_edf[[type[i]]] = sweep(rss, 2, 2*sigma^2*edf, '+')
    cp_hdf[[type[i]]] = do.call(rbind, cp_hdf_tmp)

  }

  return(list(edf=cp_edf, hdf=cp_hdf))

}
result = cp.aicc.boss()

# average over replications
cp_result_avg = lapply(result, function(xx){lapply(xx, colMeans)})
# average optimal subset size by minimizing Cp in each replication
cp_result_min = lapply(result, function(xx){lapply(xx, function(yy){round(mean(apply(yy, 1, which.min)-1),0)})})

plot.cp.edf.hdf(cp_result_avg, cp_result_min, 'boss_cp_edf_hdf.eps')

### Figure 5: RMSE along the solution paths of FS and BOSS --------
plot.solpath.boss.fs <- function(){
  categ_type = c('true_model/sparse/ex3/p0_6', 'true_model/sparse/ex4/p0_6')
  categ_name = c('Sparse-Ex3', 'Sparse-Ex4')
  rho = 'rho_09'
  n = 200
  p = 30
  snr = 'hsnr'

  pp = list()
  for(i in 1:length(categ_type)){
    x = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',categ_type[i],'/',rho,'/n_',n,'/p_',p,'/x.txt',sep=""))))
    y = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',categ_type[i],'/',rho,'/n_',n,'/p_',p,'/y_',snr,'.txt',sep=""))))
    beta = unname(as.matrix(read.table(file=paste(base, '/code/as_gram/data/',categ_type[i],'/',rho,'/n_',n,'/p_',p,'/beta.txt',sep=""))))
    mu = x %*% beta

    nrep = dim(y)[2]

    betahat = list()
    for(rep in 1:nrep){
      boss_model = boss(x, y[,rep], intercept = FALSE, hdf.ic.boss = FALSE)
      betahat$boss[[rep]] = as.matrix(boss_model$beta_boss)
      betahat$fs[[rep]] = as.matrix(boss_model$beta_fs)
      betahat$bs[[rep]] = leaps.fullsample(x, y[,rep], intercept = FALSE)
    }

    rmse = lapply(betahat, function(xx){lapply(xx, function(yy){ sqrt(colSums(sweep(x%*%yy,1,mu,'-')^2)/n)  })})
    rmse_avg = lapply(rmse, function(xx){colMeans(do.call('rbind', xx))})

    # make the plot
    k_start = 0
    k_end = p
    df_toplot <- data.frame(k=rep(k_start:k_end,3), rmse=c(rmse_avg$boss[(k_start+1):(k_end+1)],rmse_avg$fs[(k_start+1):(k_end+1)],rmse_avg$bs[(k_start+1):(k_end+1)]),
                            type=c(rep('BOSS',k_end-k_start+1), rep('FS',k_end-k_start+1), rep('BS',k_end-k_start+1)))
    df_toplot$type = factor(df_toplot$type, levels = c('BOSS', 'BS', 'FS'))

    p1 <- ggplot() + geom_point(data=df_toplot, aes(y = rmse, x = k, colour = type, shape = type), stat="identity",size=1.5)
    p1 <- p1 + xlab('subset size') + ylab('RMSE')
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
    p1 <- p1 + ggtitle(categ_name[i])
    pp[[categ_name[i]]] = p1
  }

  mylegend = g_legend(pp$`Sparse-Ex3`)
  setEPS()
  postscript(file=paste(base, '/paper/figures/lossratio_fs_boss_k.eps',sep=""), height=3, width=6)
  pp0 <- grid.arrange(arrangeGrob(pp$`Sparse-Ex3` + theme(legend.position="none"),
                                  pp$`Sparse-Ex4` + theme(legend.position="none"),
                                  nrow=1,ncol=2,heights=3),
                      mylegend,nrow=2,heights=c(9,1))
  print(pp0)
  dev.off()
}
plot.solpath.boss.fs()




### Contingency table LBS vs BS --------
categ_type = 'true_model/orthogonal/sparse/ex1/p0_6'
snr = 'hsnr'
n = 200

data_toplot = data.frame()
for(p in c(30,180)){
  count = count + 1
  if(n == 200){
    result <- readRDS(paste(base, '/code/as_gram/results/',categ_type,'/n_',n,'/p_',p,'/',snr,'/result_bs_cp.rds',sep=""))
  }else if(n == 2000){
    result <- readRDS(paste(base, '/code/as_gram/results/',categ_type,'/n_',n,'/p_',p,'/',snr,'/result.rds',sep=""))
  }
  result_lbs <- readRDS(paste(base, '/code/as_gram/results/',categ_type,'/n_',n,'/p_',p,'/',snr,'/result_lbs_new.rds',sep=""))

  numvar_bs = result$extravariable$bestsub$ic$cov$edf + result$sparsistency$bestsub$ic$cov$edf
  numvar_lbs = result_lbs$extravariable$lbs$ic$cov$edf + result_lbs$sparsistency$lbs$ic$cov$edf
  numvar_bs[numvar_bs >= 9] = ">=9"
  numvar_lbs[numvar_lbs >= 9] = ">=9"
  numvar_bs = factor(as.character(numvar_bs), levels=c('6','7','8','>=9'))
  numvar_lbs = factor(as.character(numvar_lbs), levels=c('6','7','8','>=9'))
  # sum(numvar_bs > numvar_lbs)
  # sum(numvar_bs < numvar_lbs)

  # if(min(numvar_bs) < 6 | min(numvar_lbs) < 6){
  #   print(count)
  # }

  tmp = melt(table(numvar_bs, numvar_lbs) / 10)
  tmp$type = paste0('p=',p)
  data_toplot = rbind(data_toplot, tmp)
}
data_toplot$type = factor(data_toplot$type, levels=c('p=30','p=180'))
#names(data_toplot)[3] = 'percentage'

ggplot(melt(tmp), aes(x=numvar_bs, y=numvar_lbs)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme_bw()


p1 = ggplot(data = data_toplot, mapping = aes(x=numvar_bs,y=numvar_lbs)) + geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), size = 3)
p1 = p1 + scale_fill_gradient(low = "white", high = "red")
p1 = p1 + facet_wrap(.~ type, ncol=2)
p1 = p1 + xlab(label='optimal subset size, BS') + ylab(label='optimal subset size, LBS')
p1 = p1 + theme(strip.text = element_text(size=10)) + theme(legend.position="none")
p1 = p1 + theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold"))
setEPS()
postscript(file=paste(base, '/paper/figures/numvar_bs_lbs.eps',sep=""), height=3, width=6)
print(p1)
dev.off()


