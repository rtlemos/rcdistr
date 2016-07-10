#' Normal random variable
#'
#' @import rcvirtual 
#' @export Normal
#' @exportClass Normal
#'
Normal <- setRefClass(
  Class = 'Normal',
  contains = 'rcvirtual.random',
  methods = list (
    initialize = function(mean = 0, variance = 1, high.resolution = FALSE, 
                          logtol = -8, npt = 10){
      'Function that creates a new Normal random variable or vector.
      Mean and variance must be provided.
      With high.resolution = FALSE (default), standard normal variates 
      range between qnorm(1/(2^32-208)) = -6.23 and 6.23.
      With high.resolution = TRUE, the range is 
      qnorm(1/(2^53-208)) = -8.2 to 8.2.'
      
      #callSuper(...)
      .self$type <- 'Normal'
      if(is(mean, 'rcvirtual.random')){
        sz <- mean$size()
        .self$nr <- sz$nr 
      } else {
        .self$nr <- length(mean)
        sz <- list(nr = .self$nr, nc = 1)
      }
      .self$univariate <- (.self$nr == 1) 
      .self$nc <- 1
      .self$param <- list(mean = NULL, variance = NULL)
      .self$dexpr <- list(mean = NULL, variance = NULL)
      if(is(mean, 'rcvirtual.random')){
        .self$param[[1]] <- mean
        .self$dexpr[[1]] <- deparse(substitute(mean))
      } else {
        .self$param[[1]] <- Constant(mean)
        .self$dexpr[[1]] <- ".self$param$mean"
      }
      if(is(variance, 'rcvirtual.random')){
        .self$param[[2]] <- variance
        .self$dexpr[[2]] <- deparse(substitute(variance))
        szv <- variance$size()
        if (szv$nr != sz$nr | szv$nc != sz$nr) stop('Variance has wrong size')
      } else {
        if (!is.matrix(variance)) {
          if (length(variance) == 1) {
            variance <- diag(x = variance, nrow = sz$nr)
          } else if (length(variance) == sz$nr) {
            variance <- diag(x = variance)
          } else {
            stop('Variance has wrong size')
          }
        } else if (!nrow(variance) == sz$nr | !nrow(variance) == sz$nr) {
          stop('Variance has wrong size')
        }
        .self$param[[2]] <- Constant(variance)
        .self$dexpr[[2]] <- ".self$param$variance"
      }
      .self$rng <- rcrandom::rcrng(n.substreams = .self$nr * .self$nc, 
                                   high.precision = high.resolution)
      .self$operations.classes <- list(
        '%*%'=c('numeric','matrix','Constant'),
        '/'=NULL, 
        '*'=c('numeric','matrix','Constant'),
        '-'=c('numeric','matrix','Constant','Normal'),
        '+'=c('numeric','matrix','Constant','Normal'))
      
      mu  <- .self$param[[1]]
      tau <- .self$param[[2]]
      if(is(mu,'Constant')){
        par <- mu$iget.parameter(id=1,eval=TRUE)
      } else if(is(mu,'Normal')){
        par <- mu$iget.parameter(id=1,eval=TRUE)
      } else if(is(mu,'Continuous')){
        mean.evals 	<- mu$iget.evals()
        par <- mean.evals$coord[ which.max(mean.evals$logd) ]
      }
      if (is(tau,'Constant')){
        covR 		<- chol(tau$iget.parameter(id=1, eval=TRUE))
        precR 	 	<- solve(covR)
        logdetPR	<- sum(log(diag(precR)))
      } else if( is(tau,'Continuous') ){
        stop('TODO asap') #TODO
      } else {
        stop('Unknown type')
      }
      .self$mode <- list(par = par, covR = covR, precR = precR, 
                         logdetPR = logdetPR)
      
      if(is(mu,'Continuous') | is(tau,'Continuous')){
        ni <- 101
        .self$Z <- lapply(1:.self$nr,FUN=function(j){
          x.o   <- .self$iget.convert(mean.evals$coord[,j],j)
          ld.o  <- mean.evals$logd
          sd.o  <- 1
          cnv <- .Fortran(
            'convolution',
            nx = as.integer(length(ld.o)), nv = as.integer(length(sd.o)), 
            ni = as.integer(ni), npt=as.integer(npt),
            x_o=as.double(x.o), ld_o=as.double(ld.o), 
            sd_o=as.double(sd.o), coord=as.double(rep(0,ni)), 
            logdensity=as.double(rep(0,ni)) )
          mloc <- which.max(cnv$logdensity)
          ptI <- list(  
            left = which.min(abs(cnv$logdensity[1:mloc ] - logtol)),
            right= which.min(abs(cnv$logdensity[mloc:ni] - logtol)) + mloc - 1)
          st <- list(left = cnv$coord[ptI$left], 
                     right = cnv$coord[ptI$right])
          ptsL<- mapply(seq(4,1,-1), FUN=function(i){
            qtl <- st$left + (cnv$coord[mloc] - st$left) * (1-exp(-2+0.5*i))
            pt  <- which.min(abs(cnv$coord - qtl))
            return(pt)
          })
          ptsR  <- mapply(seq(1,4, 1), FUN=function(i){
            qtl <- st$right + (cnv$coord[mloc] - st$right) * (1-exp(-2+0.5*i))
            pt  <- which.min(abs(cnv$coord - qtl))
            return(pt)
          })
          cmax  <- cnv$coord[mloc]
          coord <- list(left = cnv$coord[ptsL],      
                        right= cnv$coord[ptsR] )
          logd  <- list(left = cnv$logdensity[ptsL], 
                        right= cnv$logdensity[ptsR])
          logFD <- list(left = (cnv$logdensity[ptsL+1]-cnv$logdensity[ptsL]) / 
                          (cnv$coord[ptsL+1]-cnv$coord[ptsL]),
                        right= (cnv$logdensity[ptsR]-cnv$logdensity[ptsR-1]) / 
                          (cnv$coord[ptsR]-cnv$coord[ptsR-1]))
          x	   <- c(coord$left, cnv$coord[mloc],      coord$right)
          y      <- c(logd$left,  cnv$logdensity[mloc], logd$right)
          d	   <- c(logFD$left, 0, logFD$right)
          Z.part <- Zunivariate$new(x=x, y=y, deriv=d)
          .self$lixo[[j]] <- list(coord=cnv$coord,logd=cnv$logdensity)
          return(Z.part)
        })
      }
    },
    
    rand = function(n=1, from.marginal=FALSE){
      'Generates n (multivariate) normal variates'
      
      #callSuper(...)
      if(from.marginal){
        u  <- if(.self$nr==1) matrix(nc=1,.self$rng$runif(n=n)) else .self$rng$runif(n=n)
        zz <- t(mapply(1:.self$nr,FUN=function(j) .self$Z[[j]]$invcdf(u[,j]) ))
        xx <- .self$mode$covR %*% zz + .self$mode$par
        x  <- if(.self$nr==1) as.numeric(xx) else t(xx) 
      } else {
        zz    <- qnorm(.self$rng$runif(n=n)) 
        mu 	  <- .self$param[[1]]
        Sigma <- .self$param[[2]]
        mn    <- if(is(mu,'Constant')) mu$iget.parameter(id=1, eval=TRUE) else mu$rand(n)
        if(is(Sigma,'Constant')){
          v <- Sigma$iget.parameter(id=1, eval=TRUE)
          x <- if(.self$nr == 1) {
            mn + sqrt(v) * zz
          } else {
            mn + chol(v) %*% zz
          }
        } else {
          v <- Sigma$rand(n)
          if(.self$nr == 1){
            x <- mn + sqrt(v) * zz
          } else {
            x <- mapply(mn, v, zz, FUN = function(mn, v, zz) {
              mn + chol(v) %*% zz
            })
          }
        }
      }
      return(x)
    },
    
    pdf = function(x=NULL,log=FALSE,verbose=FALSE,...){
      'Provides either the pdf or its evaluation for given x'
      
      #callSuper(...)
      if(verbose) cat(.self$show())
      mu  <- .self$param[[1]]
      tau <- .self$param[[2]]
      if(is(mu,'Constant') & is(tau,'Constant') ){
        mn  <- mu$iget.parameter(id=1)
        var <- tau$iget.parameter(id=1)
        if(.self$nr==1){
          myf <- function(x) dnorm(x,mn, sd=sqrt(var), log=log)
        } else {
          if(log){
            myf <- function(x) (-0.5*.self$nr) * log(2*pi) - 
              0.5*determinant(var)$modulus[1] - 0.5*crossprod(x-mn,solve(var, x-mn))
          } else {
            myf <- function(x) (2*pi)^(-0.5*.self$nr) * 
              exp(-0.5*crossprod(x-mn,solve(var, x-mn)))/sqrt(det(var))
          }
        }
      } else if(is(mu,'Normal') & is(tau,'Constant') ){
        if(is(mu$param$mean,'Constant') & is(mu$param$variance,'Constant') ){
          mu.mn  <- mu$iget.parameter(1)
          mu.var <- mu$iget.parameter(2)
          if(.self$nr==1){
            var <- .self$iget.parameter(2)+mu.var
            myf <- function(x) dnorm(x,mu.mn, sd=sqrt(var), log=log)
          } else {
            Var = .self$iget.parameter(2) + mu.var
            if(log){
              myf <- function(x) (-0.5*.self$nr) * log(2*pi) - 
                0.5*determinant(Var)$modulus[1] - 
                0.5*crossprod(x-mu.mn, solve(Var, x-mu.mn))
            } else {
              myf <- function(x)(2*pi)^(-0.5*.self$nr)/sqrt(det(Var)) * 
                exp(-0.5*crossprod(x-mu.mn,solve(Var, x-mu.mn)))
            }
          }
        } else {
          stop('Hierarchy is too deep.')
        }
      } else if(is(mu,'Constant') & is(tau,'InvGamma') ){
        #TODO correct this and eliminate integration
        if(is(tau$shape,'Constant') & is(tau$rate,'Constant') ){
          tau.s <- tau$iget.parameter(id=1)
          tau.r <- tau$iget.parameter(id=2)
          tR  <- tau.r/tau.s
          tn  <- 2*tau.s
          tm  <- .self$iget.parameter(1)
          f0  <- function(x)(tn+(x-tm)^2/tR)^(-0.5*(tn+1))
          nR  <- integrate(f0,tm-20*sqrt(tR),tm+20*sqrt(tR),subdivisions=1000L)$value
          myf <- function(x)(tn+(x-tm)^2/tR)^(-0.5*(tn+1))/nR 
        } else {
          stop('Hierarchy is too deep.')
        }
      } else if(is(mu,'Continuous') & is(tau,'Constant') ){
        myf <- function(x){
          z    <- as.numeric(.self$mode$precR %*% (x-.self$mode$par))
          eval <- prod(mapply(1:.self$nr, FUN=function(j) .self$Z[[j]]$pdf(z[j]) ))
          if(log){
            out <- log(eval) + .self$mode$logdetPR
          } else {
            out <- eval * exp(.self$mode$logdetPR)
          }
          return(out)
        }
      }
      if(!missing(x)){
        eval <- myf(x)
        return(as.numeric(eval))
      } else {
        return(myf)
      }
    },
    
    cdf = function(x=NULL,lower.tail=TRUE, log.p=FALSE, verbose=FALSE,...){
      'Computes the Normal cdf(x)'
      
      #callSuper(...)
      if(verbose) cat(.self$show())
      myf <- function(x){
        if(is(.self$param[[1]],'Constant') & is(.self$param[[2]],'Constant') ){
          mn  <- .self$iget.parameter(1, eval = TRUE)
          var <- .self$iget.parameter(2, eval = TRUE)
          if(.self$nr==1 & .self$nc==1){
            d <- pnorm(x,mn, sd=sqrt(var), lower.tail=lower.tail, log.p=log.p)
          } else {
            stop('Multivariate Normal CDF has not been implemented.')
          }
        } else {
          stop('Not implemented yet.')
        }
        return(d)
      }
      if(!is.null(x)){
        eval <- myf(x)
        return(as.numeric(eval))
      } else {
        return(myf)
      }
    },
    
    invcdf = function(p=NULL, lower.tail=TRUE, log.p=FALSE, verbose=FALSE){
      'Computes quantiles associated with input probability(ies)'
      
      #callSuper(...)
      if(verbose) cat(.self$show())
      if(.self$nr!=1) stop('Method only available for Univariate Normals.')
      mu  <- .self$param[[1]]
      tau <- .self$param[[2]]
      if(is(mu,'Constant') & is(tau,'Constant') ){
        mn  <- mu$iget.parameter(id=1)
        var <- tau$iget.parameter(id=1)
        myf <- function(x) qnorm(p, mean=mn, sd=sqrt(var),lower.tail=lower.tail, log.p=log.p)
      } else if(is(mu,'Constant') & is(tau,'InvGamma') ){
        stop('InvCDF for NormalInvGamma TODO')
      }
      if(!missing(p)){
        eval <- myf(p)
        return(as.numeric(eval))
      } else {
        return(myf)
      }
    },
    
    iget.optimization.parameters = function(){
      return(list(mean=.self$param$mean$iget.parameter(1),
                  sd=sqrt(.self$param$variance$iget.parameter(1))))
    },
    
    iset.operate = function(operation, operand, operand.name, operand.side, my.name){
      if(! .self$is.operation.allowed(operation, class(operand) )) stop('Normal: Invalid operation')
      mu    <- .self$param$mean$copy(shallow = FALSE)
      tau   <- .self$param$variance$copy(shallow = FALSE)
      munm  <- paste0(my.name,'$param$mean')
      taunm <- paste0(my.name,'$param$variance')
      switch(operation,
             '-'=,
             '+'= {
               switch( class(operand),
                       'numeric'=,
                       'matrix'=,
                       'Constant'= {mu$iset.operate(operation, operand, operand.name, 
                                                    operand.side, munm)},
                       'Normal' = {
                         mu.o  <- operand$param$mean
                         tau.o <- operand$param$variance
                         mnm.o <- paste0(operand.name,'$param$mean')
                         tnm.o <- paste0(operand.name,'$param$variance')
                         mu$iset.operate (operation,      mu.o, mnm.o, operand.side, munm)
                         tau$iset.operate(operation='+', tau.o, tnm.o, operand.side, taunm)
                       },
                       stop('Normal: Invalid operand.') #default
               )},
             '*'= switch( class(operand),
                          'numeric'=,
                          'matrix'=,
                          'Constant'= {
                            if(is(operand,'rcvirtual.random')){
                              operand2 <- operand
                              op.squared <- paste0(operand.name,'*',operand.name)
                            } else {
                              operand2 <- operand^2
                              op.squared <- paste0(operand.name,'^2')
                            }
                            mu$iset.operate (operation, operand, operand.name, operand.side, munm)
                            tau$iset.operate(operation, operand2,  op.squared, operand.side, taunm)
                          },
                          stop('Normal: Invalid operand.') #default
             ),
             '%*%'= switch( class(operand),
                            'numeric'=,
                            'matrix'=,
                            'Constant'= {
                              if(operand.side == 'right') stop('Operand on right side not allowed.')
                              operandt <- t(operand)
                              operandt.name <- paste0('t(',operand.name,')')
                              tau2 <- paste0(taunm, '$parameter(1) %*% t(', operand.name, ')')
                              mu$iset.operate (operation, operand,  operand.name, operand.side, munm)
                              tau$iset.operate(operation, operandt, operandt.name, 'right', taunm)
                              tau$iset.operate(operation, operand,  operand.name, 'left', tau2)
                            },
                            stop('Normal: Invalid operand.') #default
             ),    
             stop('Normal: Invalid operation, stage 2.') #default
      )
      .self$param$mean <- mu
      .self$param$variance <- tau
    }
  )
)