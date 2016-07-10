#' Allowed operations with random abstract types
#'
#' @param l 
#' @param r 
#'
#' @export '+.rcvirtual.random'
#' @export '-.rcvirtual.random'
#' @export '*.rcvirtual.random'
#' @export '/.rcvirtual.random'
#' @export '%.%'
#'
'+.rcvirtual.random' <- function (l, r) {
  if (missing(r)) { # this case occurs when the user writes " -x "
    argl <- '0'
    argr <- deparse(substitute(l))
    r <- l
    l <- 0
  } else {
    argl <- deparse(substitute(l))
    argr <- deparse(substitute(r))
  }
  if (!is(r, 'rcvirtual.random')){ 
    # right argument is a not a rcvirtual.random, 
    # therefore left argument must be a rcvirtual.random                 
    x <- l$copy(shallow = FALSE)
    szx <- x$size()[[1]]
    szr <- length(r)
    if (szx != szr) {
      if (szr == 1) {
        szr <- szx
        if (l$nc > 1) {
          r <- matrix(nr = l$nr, nc = l$nc, r)
          argr <- paste0('matrix(nr = ', l$nr, ', nc = ', l$nc, ', ', r, ')')
        } else {
          argr <- paste0('c(', paste(rep(r, szr), collapse = ', '), ')')
          r <- rep(r, szr)
        }
      } else { 
        stop('Incompatible dimensions.')
      }
    }
    x$iset.operate(operation = '+', operand = r, operand.name = argr, 
                   operand.side = 'right', my.name = argl)
  } else if (!is(l, 'rcvirtual.random')){ 
    # left argument is a not a rcvirtual.random, 
    # therefore only right argument is a rcvirtual.random 
    x <- r$copy(shallow = FALSE)
    szx <- x$size()[[1]]
    szl <- length(l)
    if (szx != szl) {
      if (szl == 1) {
        szl <- szx
        if (r$nc > 1) {
          l <- matrix(nr = r$nr, nc = r$nc, l)
          argl <- paste0('matrix(nr = ', r$nr, ', nc = ', r$nc, ', ', l, ')')
        } else {
          argr <- paste0('c(', paste(rep(l, szl), collapse = ', '), ')')
          l <- rep(l, szl)
        }
      } else { 
        stop('Incompatible dimensions.')
      }
    }
    x$iset.operate(operation = '+', operand = l, operand.name = argl, 
                   operand.side = 'left', my.name = argr)	
  } else { #both arguments are rcvirtual.randoms
    if(l$size()[[1]] != r$size()[[1]]) stop('Incompatible dimensions.')
    if(l$is.operation.allowed( operation='+', class(r) )) {
      x <- l$copy(shallow = FALSE)
      x$iset.operate(operation = '+', operand = r, operand.name = argr, 
                     operand.side = 'right', my.name = argl)
    } else if(r$is.operation.allowed( operation='+', class(l) )) {
      x <- r$copy(shallow = FALSE)
      x$iset.operate(operation = '+', operand=l, operand.name = argl, 
                     operand.side = 'left' , my.name = argr)
    } else { #converting one argument to Continuous, if possible
      if(!is(l,'Uniform') & !is(l,'Continuous')){
        x <- as.Continuous(l)
        x$iset.operate(operation = '+', operand = r, operand.name = argr, 
                       operand.side = 'right', my.name = argl)
      } else if(!is(r,'Uniform') & !is(r,'Continuous') ){
        x <- as.Continuous(r)
        x$iset.operate(operation = '+', operand = l, operand.name = argl, 
                       operand.side = 'left', my.name = argr)
      } else {
        stop('Unsucessful operation.')
      }
    }
  }
  return(x)
}

#Other operations are analogous to the function '+.rcvirtual.random'
#Simple character substitutions (e.g., '-' instead of '+') are performed
'-.rcvirtual.random' <- eval(
  parse(text = gsub('+', '-', deparse(get('+.rcvirtual.random')), fixed=TRUE)))
'*.rcvirtual.random' <- eval(
  parse(text = gsub('+', '*', deparse(get('+.rcvirtual.random')), fixed=TRUE)))
'/.rcvirtual.random' <- eval(
  parse(text = gsub('+', '/', deparse(get('+.rcvirtual.random')), fixed=TRUE)))
'%.%' <- eval(
  parse(text = gsub('+','%*%', deparse(get('+.rcvirtual.random')), fixed=TRUE)))

#' Transpose
#'
#' @param q 
#'
#' @export 't.rcvirtual.random'
#'
't.rcvirtual.random' <- function(q){
  if (class(q) == 'Constant'){
    x <- Constant(t(q$parameter(id = 1)))
  } else if(class(q) == 'Normal'){
    x <- Normal(t(q$parameter(id = 1)), q$parameter(id = 2))
  }
  return(x)
}

#' Exponential function
#'
#' @param x 
#'
#' @export 'exp.rcvirtual.random'
#'
'exp.rcvirtual.random' <- function(x){
  arg <- deparse(substitute(x))
  if (is(x,'Normal')){
    q <- LogNormal()
    q$param$mean$iset.dexpr(
      1, paste0('(', arg, ')$param$mean$parameter(1,eval=TRUE)'))
    q$param$variance$iset.dexpr(
      1, paste0('(', arg, ')$param$variance$parameter(1,eval=TRUE)'))
  } else {
    stop('TODO: Exponential transformation')
  }
  return(q)
}

#' Logarithmization of random variables
#'
#' @param x 
#'
#' @export Log
#'
Log <- function(x){
  # for some unknown reason, deparse(substitute(x)) doesn't work if the function 
  # is called "log.rcvirtual.random" instead of "Log"
  argx <- deparse(substitute(x)) #; print(argx)
  if (is(x,'Constant')){
    q <- Constant(log(x$parameter(id = 1, eval = TRUE)))
    q$iset.dexpr(1, paste0('log((',argx,')$parameter(id = 1, eval = TRUE))'))
  } else {
    stop('TODO: Logarithmic transformation')
  }
  return(q)
}

#' Redefine a random variable as a Continuous r.v.
#'
#' @param x 
#' @param npoints 
#' @param logtol 
#'
#' @export as.Continuous
#'
as.Continuous <- function(x, npoints = 9, logtol = -8) {
  if (is(x,'rcvirtual.random')) {
    lb <- if (is.null(x$param$lb)) NULL else x$param$lb$parameter(1, eval = TRUE) 
    ub <- if (is.null(x$param$ub)) NULL else x$param$ub$parameter(1, eval = TRUE)
    z  <- Continuous(f = x$pdf(), initial = x$invcdf(0.5), lower = lb, 
                     upper = ub, aux = x$iget.optimization.parameters(), 
                     high.resolution = x$is.high.resolution(),
                     is.f.logtransformed = FALSE, logtol = logtol, 
                     nside=ceiling((npoints - 1) / 2))
    return(z)
  } else {
    stop('Input is not of Random Abstract Type.')
  }
}