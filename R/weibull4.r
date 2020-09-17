######Function weibull4.fit fit data with Weibull 4-parameter distribution
# Performs a Metropolis algorithm to Markov chain Monte Carlo method
# special usage for fitting COVID-19 epidemic data
# Part of the Metropolis-MCMC routine was extracted from the Florian Hartig blog:
# "https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/"
#####################################################################

weibull4.fit <- function(x, y, shape=NA, scale=NA, loc=NA, area=NA, iter=1000, xmax=0) {
  xi <- seq(1,length(x))
  yi <- y

  if (is.null(loc) | is.na(loc)) {  # Sets the location parameter to the length on the initial no case time
    loc <- length(which(yi <= max(yi, na.rm=T)*0.05 & xi < mean(xi, na.rm=T)))
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

  sd <- sd(weibull4(xi, shape, scale, loc, area) - yi, na.rm=T)

  startvalue <- c(shape, scale, loc, area, sd)
  chain <- run_metropolis_MCMC(xi, yi, startvalue, iter)
  burnIn <- iter/2
  acceptance <- 1-mean(duplicated(chain[-(1:burnIn),]))

  param <- rbind(apply(chain[-(1:burnIn),], 2, mean), apply(chain[-(1:burnIn),], 2, sd))
  colnames(param) <- c("shape", "scale", "location", "area", "sd")
  rownames(param) <- c("mean", "sd")

  if (xmax != 0) {
    y.out <- c(yi, rep(NA, xmax - max(xi)))
    x.out <- c(x, seq(max(x)+1, xmax, 1))
  } else {
    y.out <- yi
    x.out <- x
  }

  fit <- weibull4(x.out, mean(chain[,1]), mean(chain[,2]), mean(chain[,3]), mean(chain[,4]))

  return(list(cbind(x.out, fit), param, chain))
}

#script2
likelihood <- function(x, y, param){
  shape <- param[1]
  scale <- param[2]
  loc <- param[3]
  area <- param[4]
  sd <- param[5]

  pred <- weibull4(x, shape, scale, loc, area)
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
  sd <- param[5]
  shapeprior <- dunif(shape, min=1, max=5, log=T)
  scaleprior <- dnorm(scale, sd=scale/2, log=T)
  locprior <- dunif(loc, min=1, max=loc*2, log=T)
  areaprior <- dunif(area, min=area/2, max=area*2, log=T)
  sdprior <- dunif(sd, min=1, max=sd*2, log=T)
  return(shapeprior + scaleprior + locprior + areaprior + sdprior)
}

#script4
posterior <- function(x, y, param){
  return (likelihood(x, y, param) + prior(param))
}

#script5
######## Metropolis algorithm ################
proposalfunction <- function(param){
  return(rnorm(5, mean=param, sd=param*0.015))
}

run_metropolis_MCMC <- function(x, y, startvalue, iterations){
  chain <- array(dim = c(iterations+1, 5))
  chain[1,] <- startvalue
  for (i in 1:iterations){
    proposal <- proposalfunction(chain[i,])

    probab <- exp(posterior(x, y, proposal) - posterior(x, y, chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] <- proposal
#      plot(xi,yi)     ## Uncomment these lines to plot a new curve for each iteration
#      lines(xi, weibull4(xi,chain[i,1],chain[i,2],chain[i,3],chain[i,4]),col="red")
                       ## This will slow down the function, so be careful
    }else{
      chain[i+1,] <- chain[i,]
    }
  }
  return(chain)
}

#######Function weibull4 builds a weibull 4-parameter PDF plot
weibull4 <- function(x=seq(0,1,length.out=10), shape=2.5, scale=1, loc=0, area=1) {
  x <- seq(1,length(x),1)
  x0 <- x[c(1:loc)]
  x1 <- x[-c(1:loc)]

  a <- c(rep(0,length(x0)), area * (shape/scale) * ((x1-loc)/scale)^(shape-1) * exp(-((x1-loc)/scale)^shape))
  return(a)
}

### END #######################
