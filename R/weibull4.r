######Package weibull4 fits data with Weibull 4-parameter distribution
# Performs a Metropolis algorithm to Markov chain Monte Carlo method
# special usage for fitting COVID-19 epidemic data
# It can work with one or two, unimodal or bimodal Weibull distributions
# Part of the Metropolis-MCMC routine was extracted from the Florian Hartig blog:
# "https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/"
# For more details, see the original paper:
# Moreau, 2021. Model Assist. Statist. App. in press.
#####################################################################
## Basic function weibull4. This script calls weibull4.fit that perform
## all the operations
weibull4 <- function(x, y, shape=NA, scale=NA, loc=NA, area=NA,
                      shape2=NA, scale2=NA, loc2=NA, area2=NA, iter=1000,
                      xmax=0, modes=1, modes2=1, split=NA)
  if(is.na(split)) {
    fit <- weibull4.fit(x, y, shape=shape, scale=scale, loc=loc, area=area,
                        shape2=shape2, scale2=scale2, loc2=loc2, area2=area2,
                        iter=iter, xmax=xmax, modes=modes)
    return(fit)
  } else if(!is.na(split)){
    if(xmax > max(x)) {maximum <- xmax} else {maximum <- max(x)}
    x1 <- x[which(x <= split)]
    y1 <- y[which(x <= split)]
    x2 <- x[which(x > split)]
    y2 <- y[which(x > split)]

    fit1 <- weibull4.fit(x1, y1, shape=shape, scale=scale, loc=loc, area=area,
                         shape2=shape2, scale2=scale2, loc2=loc2, area2=area2, iter=iter,
                         xmax=maximum, modes=modes)
    fit2 <- weibull4.fit(x2, y2, shape=shape, scale=scale, loc=loc, area=area,
                         shape2=shape2, scale2=scale2, loc2=loc2, area2=area2, iter=iter,
                         xmax=xmax, modes=modes2)

    x1o <- fit1[[1]][,1]
    y1o <- fit1[[1]][,2]
    y2o <- fit2[[1]][,2]
    a1o <- fit1[[2]][1,4] + fit2[[2]][1,4]
    a2o <- fit1[[2]][1,8] + fit2[[2]][1,8]
    a1sdo <- sqrt(fit1[[2]][2,4]^2 + fit2[[2]][2,4]^2)
    a2sdo <- sqrt(fit1[[2]][2,8]^2 + fit2[[2]][2,8]^2)


    yo <- y1o + c(rep(0, length(y1o) - length(y2o)), y2o)
    fit1[[1]][,2] <- yo
    fit1[[2]][1,4] <- a1o
    fit1[[2]][1,8] <- a2o
    fit1[[2]][2,4] <- a1sdo
    fit1[[2]][2,8] <- a2sdo

    return(fit1)
  }

## weibull.fit
weibull4.fit <- function(x, y, shape=NA, scale=NA, loc=NA, area=NA,
                         shape2=NA, scale2=NA, loc2=NA, area2=NA, iter=1000,
                         xmax=0, modes=1) {
  xi <- seq(1,length(x))
  yi <- y

  if (is.null(loc) | is.na(loc)) {  # Sets the location parameter to the length on the initial no case time
    loc <- length(which(yi <= max(yi, na.rm=T)*0.05 & xi < mean(xi, na.rm=T)))
    if(loc <= 0) {loc <- 1}
  }
  if (is.null(shape) | is.na(shape)) {  # Sets the shape parameter to 2.5 (COVID-19 curve shape)
    shape <- 2.5
  }
  if (is.null(scale) | is.na(scale)) {    # set the initial scale parameters to about the mode of PDF
    p_scale <- seq(1, length(na.exclude(yi)))
    scale <- mean(p_scale[-c(1:loc)], na.rm=T) - loc
  }
  if (is.null(area) | is.na(area)) {   # Sets the area to the total number of cases or deaths
    area <- sum(yi, na.rm=T)
    }
  if (is.null(loc2) | is.na(loc2)) {  # Sets the location parameter to the length on the initial no case time
    loc2 <- loc + scale
  }
  if (is.null(shape2) | is.na(shape2)) {  # Sets the shape parameter to 2.5 (COVID-19 curve shape)
    shape2 <- shape
  }
  if (is.null(scale2) | is.na(scale2)) {    # set the initial scale parameters to about the mode of PDF
    scale2 <- scale
  }
  if (is.null(area2) | is.na(area2)) {   # Sets the area to the total number of cases or deaths
    area2 <- area
  }

  if(modes == 2) {
  sd <- sd(weibull4.build(xi, shape, scale, loc, area, shape2, scale2, loc2, area2,
                    modes=modes) - yi, na.rm=T)
  } else if(modes == 1) {
    sd <- sd(weibull4.build(xi, shape, scale, loc, area, modes=modes) - yi, na.rm=T)
  }

  startvalue <- c(shape, scale, loc, area, shape2, scale2, loc2, area2, sd)
  chain <- run_metropolis_MCMC(xi, yi, startvalue, iter, modes)
  burnIn <- iter/2
  acceptance <- 1-mean(duplicated(chain[-(1:burnIn),]))

  param <- rbind(apply(chain[-(1:burnIn),], 2, mean), apply(chain[-(1:burnIn),], 2, sd))
  colnames(param) <- c("shape", "scale", "location", "area",
                       "shape2", "scale2", "location2", "area2", "sd")
  rownames(param) <- c("mean", "sd")

  if (xmax != 0) {
    y.out <- c(yi, rep(NA, xmax - max(xi)))
    x.out <- c(x, seq(max(x)+1, xmax, 1))
  } else {
    y.out <- yi
    x.out <- x
  }

  fit <- weibull4.build(x.out, mean(chain[,1]), mean(chain[,2]), mean(chain[,3]), mean(chain[,4]),
                  mean(chain[,5]), mean(chain[,6]), mean(chain[,7]), mean(chain[,8]),
                  modes=modes)

  return(list(cbind(x.out, fit), param, chain))
}

#script2
likelihood <- function(x, y, param, modes){
  shape <- param[1]
  scale <- param[2]
  loc <- param[3]
  area <- param[4]
  shape2 <- param[5]
  scale2 <- param[6]
  loc2 <- param[7]
  area2 <- param[8]
  sd <- param[9]

  pred <- weibull4.build(x, shape, scale, loc, area,
                   shape2, scale2, loc2, area2, modes=modes)
  singlelikelihoods <- dnorm(y, mean=pred, sd=sd, log=T)
  sumll <- sum(singlelikelihoods, na.rm=T)
  return(sumll)
}

#script3
# Prior distribution
prior <- function(param){
  shape <- param[1]
  scale <- param[2]
  loc <- param[3]
  area <- param[4]
  shape2 <- param[5]
  scale2 <- param[6]
  loc2 <- param[7]
  area2 <- param[8]
  sd <- param[9]
  shapeprior <- dunif(shape, min=1, max=5, log=T)
  scaleprior <- dnorm(scale, sd=scale/2, log=T)
  locprior <- dunif(loc, min=1, max=loc*2, log=T)
  areaprior <- dunif(area, min=area/2, max=area*2, log=T)
  shapeprior2 <- dunif(shape2, min=1, max=5, log=T)
  scaleprior2 <- dnorm(scale2, sd=scale2/2, log=T)
  locprior2 <- dunif(loc2, min=1, max=loc2*2, log=T)
  areaprior2 <- dunif(area2, min=area2/2, max=area2*2, log=T)
  sdprior <- dunif(sd, min=1, max=sd*2, log=T)
  return(shapeprior + scaleprior + locprior + areaprior +
           shapeprior2 + scaleprior2 + locprior2 + areaprior2 + sdprior)
}

#script4
posterior <- function(x, y, param, modes){
  return (likelihood(x, y, param, modes) + prior(param))
}

#script5
######## Metropolis algorithm ################
proposalfunction <- function(param){
  return(rnorm(9, mean=param, sd=param*0.015))
}

run_metropolis_MCMC <- function(x, y, startvalue, iterations, modes){
  chain <- array(dim = c(iterations+1, 9))
  chain[1,] <- startvalue
  for (i in 1:iterations){
    proposal <- proposalfunction(chain[i,])

    probab <- exp(posterior(x, y, proposal, modes) - posterior(x, y, chain[i,], modes))
    if (runif(1) < probab){
      chain[i+1,] <- proposal
#      plot(xi,yi)     ## Uncomment these lines to plot a new curve for each iteration
#      lines(xi, weibull4.build(xi,chain[i,1],chain[i,2],chain[i,3],chain[i,4]),col="red")
                       ## This will slow down the function, so be careful
    }else{
      chain[i+1,] <- chain[i,]
    }
  }
  return(chain)
}

#######Function weibull4 builds a weibull 4-parameter PDF plot
weibull4.build <- function(x=seq(0,1,length.out=10), shape=2.5, scale=1, loc=0, area=20,
                     shape2=5, scale2=2, loc2=6, area2=1, modes=1) {
  x <- seq(1,length(x),1)
  x01 <- x[c(1:loc)]
  x11 <- x[-c(1:loc)]
  x02 <- x[c(1:loc2)]
  x12 <- x[-c(1:loc2)]

  a1 <- c(rep(0,length(x01)), area * (shape/scale) * ((x11-loc)/scale)^(shape-1) * exp(-((x11-loc)/scale)^shape))
  a2 <- c(rep(0,length(x02)), area2 * (shape2/scale2) * ((x12-loc2)/scale2)^(shape2-1) * exp(-((x12-loc2)/scale2)^shape2))

  if(modes == 2) {
  a <- a1 + a2
  } else if(modes == 1) {
    a <- a1
  }
  return(a)
}

### END #######################
