#' Constant 
#'
#' @import rcvirtual
#' @export Constant
#'
Constant <- setRefClass(
  Class = 'Constant',
  contains = 'rcvirtual.random',
  methods = list(
    initialize = function(value = 0){
      'Create a new Constant object (scalar, vector, matrix, or Constant)'
      
      #callSuper(...)
      .self$type <- 'Constant'
      if(is(value,'numeric')){
        .self$nr <- length(value)
        .self$nc <- 1
      } else if(is(value,'matrix')){
        .self$nr <- nrow(value)
        .self$nc <- ncol(value)
      } else if(is(value,'Constant')){
        sz <- value$size()
        .self$nr <- sz$nr
        .self$nc <- sz$nc
      } else {
        stop('Invalid value.')
      }
      if(is(value,'Constant')){
        .self$dexpr <- list(value = deparse(substitute(value)))
      } else {
        .self$dexpr <- list(value = as.character(value))
      }
      .self$param <- list(value = value)
      .self$operations.classes <- list(
        '%*%' = c('numeric', 'matrix', 'Constant'),
        '/' = c('numeric','matrix','Constant'), 
        '*' = c('numeric','matrix','Constant'),
        '-' = c('numeric','matrix','Constant'),
        '+' = c('numeric','matrix','Constant'))
    },
    
    pdf = function(x = NULL, log = FALSE, ...){
      'Computes pdf(x)'
      
      #callSuper(...)
      v <- .self$parameter(id = 1, eval = TRUE)
      myf <- function(x){
        if(log){
          d <- if(identical(x,v)) 0 else -Inf
        } else {
          d <- if(identical(x,v)) 1 else 0
        }
        return(d)
      }
      if(!missing(x)){
        eval <- if (is(x, 'list')) mapply(x, FUN = myf) else myf(x)
        return(as.numeric(eval))
      } else {
        return(myf)
      }
    },
    
    cdf = function(x = NULL, lower.tail = TRUE, log.p = FALSE,...) {
      'Computes cdf(x)'
      
      #callSuper(...)
      v <- .self$parameter(id = 1, eval = TRUE)
      myf <- function(x){
        if (lower.tail) {
          if (log.p) {
            d <- if(any(x < v)) -Inf else 0
          } else {
            d <- if (any(x < v)) 0 else 1
          } 
        } else {
          if (log.p) {
            d <- if (any(x > v)) -Inf else 0 
          } else {
            d <- if (any(x > v)) 0 else 1
          }
        }
        return(d)
      }
      if (!missing(x)){
        eval <- mapply(x, FUN = myf)
        return(as.numeric(eval))
      } else {
        return(myf)
      }
    },
    
    invcdf = function(p = NULL, lower.tail = TRUE, log.p = FALSE){
      'Computes quantiles associated with input probability(ies)'
      
      v <- .self$parameter(id = 1, eval = TRUE)
      myf <- function(p) {
        out <- mapply(p, FUN = function(p) if (p >= 0 & p <= 1) v else -Inf)
        return(out)
      }
      if (!missing(p)) {
        return(myf(p))
      } else {
        return(myf)
      }
    },
    
    rand = function(n = 1) {
      "Generates n replicates of this object's value"
      
      v <- .self$iget.parameter(id = 1, eval = TRUE)
      if (is(v, 'numeric')) {
        out <- mapply(1:n, FUN = function(i) v)
      } else if (is(v, 'matrix')) {
        out <- lapply(1:n, function(i) v)
      }
      return(out)
    },
    
    ###########################################
    # Private methods #########################
    ###########################################
    iset.operate = function(operation, operand, operand.name, 
                            operand.side, my.name){
      if (! .self$is.operation.allowed(operation, class(operand) )) {
        stop('Constant: Invalid operation')
      }
      #   print(operation)
      #   print(operand)
      #   print(operand.name)
      #   print(operand.side)
      #   print(my.name)
      if (is(operand, 'Constant')) {
        operand.value <- operand$iget.parameter(1)
        p1 <- paste0('(', operand.name,')$parameter(id = 1, eval = TRUE)')
      } else {
        operand.value <- operand
        p1 <- operand.name 
      }
      p2 <- paste0(my.name,'$parameter(id = 1, eval = TRUE)')
      current.value <- .self$iget.parameter(id = 1, eval = TRUE)
      .self$iset.parameter(
        id = 1, value = get(operation)(operand.value, current.value))
      if (operand.side == 'left') {
        .self$dexpr$value <- paste0(p1, operation, p2)
      } else{
        .self$dexpr$value <- paste0(p2, operation, p1)
      }
    }
  )
)