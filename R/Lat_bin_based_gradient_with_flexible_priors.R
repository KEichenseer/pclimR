gradient <- function(x, coeff, sdy) { # sigma is labelled "sdy"
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]
  return(A + max(c(K-A,0))/(1+(exp(Q*(x-M)))) + rnorm(length(x),0,sdy))
}

loglik <- function(x, ymean, yest, sdyest, coeff, sdy) {
  A = coeff[1]
  K = coeff[2]
  M = coeff[3]
  Q = coeff[4]

  ll1 <- sum(dnorm(ymean, yest, sdyest,log=TRUE),na.rm=T)
  pred = A + max(c(K-A,0))/(1+(exp(Q*(x-M))))
  ll2 <- sum(dnorm(yest, mean = pred, sd = sdy, log = TRUE))
  return(ll1+ll2)
}

logprior <- function(coeff, yest, prior_df) {
  return(sum(c(
    dunif(coeff[1], -4, coeff[2], log = TRUE),
    dunif(coeff[2], coeff[1], 40, log = TRUE),
    dnorm(coeff[3], 45, 10, log = TRUE),
    dlnorm(coeff[4], -2, 1, log = TRUE),
    sum(sapply(1:nbin, function(b) dunif(yest[b], -4, 40, log = TRUE)+
                 extra_priors[[b]](yest[b], prior_df[b,2],prior_df[b,3],log=TRUE))))))
}
# comment: The "extra prior" could also be in the likelihood, right? Treating them as data, so to speak
# on this note, uncertainties could be incorporated into individual data points
logposterior <- function(x, ymean, yest, sdyest, coeff, sdy, prior_df){
  return (loglik(x, ymean, yest, sdyest, coeff, sdy) + logprior(coeff,yest,prior_df))
}

MH_propose_coeff <- function(coeff, prop_sd_coeff){
  return(rnorm(4,mean = coeff, sd= prop_sd_coeff))
}
MH_propose_yest <- function(yest, prop_sd_yest){
  return(mvnfast::rmvn(1,mu = yest, sigma = prop_sd_yest))
}


# function for plotting the 95 % CI shading
error_polygon <- function(x,en,ep,color) {
  polygon( c(x[1], x, x[length(x)], x[length(x)], rev(x), x[1]),
           c((ep)[1],ep, (ep)[length(ep)], (en)[length(en)], rev(en), (en)[1]),
           border = NA, col = color)
}


# Main MCMCM function
run_MCMC <- function(nIter, x, yobs, prior_df, coeff_inits, sdy_init, yest_inits, sdyest_inits,
                     prop_sd_coeff, prop_sd_yest, nbin){
  ### Initialisation
  coefficients = array(dim = c(nIter,4)) # set up array to store coefficients
  coefficients[1,] = coeff_inits # initialise coefficients
  sdy = rep(NA_real_,nIter) # set up vector to store sdy
  sdy[1] = sdy_init # intialise sdy

  yestimate = array(dim = c(nIter,nbin)) # set up array to store coefficients
  yestimate[1,] = yest_inits # initialise coefficients

  sdyest = array(dim = c(nIter,nbin)) # set up vector to store sdy
  sdyest[1,] = sdyest_inits # intialise sdy

  A_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  B_sdy = 1 # parameter for the prior on the inverse gamma distribution of sdy
  shape_sdy <- A_sdy+nbin/2 # shape parameter for the inverse gamma


  #### Investigate these: Need to be broad for single obser
  A_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
  B_sdyest = 1 # parameter for the prior on the inverse gamma distribution of sdyest
  ####
  nest = sapply(yobs,length) # number of observations per bin
  nest[which(is.na(yobs))] = NA
  shape_sdyest =  A_sdyest+(nest)/2 # shape parameter for the inverse gamma
  ### n-1?!
  ymean = sapply(yobs,mean) # y mean of each bin

  logpost = rep(NA,nIter) # store the log posterior for

  extra_priors <- list(NA,nbin)
  for(i in 1:nbin){
    if(!(is.na(prior_df[i,1]))) {
      extra_priors[[i]] = get(prior_df[i,1])
    } else  extra_priors[[i]] = function(a,b,c,log) return(0)
  }

  ### The MCMC loop
  for (i in 2:nIter){
    ## 1. Gibbs step to estimate sdy
    sdy[i] = sqrt(1/rgamma(
      1,shape_sdy,B_sdy+0.5*sum((yestimate[i-1,]-gradient(x,coefficients[i-1,],0))^2)))

    ## 2. Gibbs step to estimate sdyest
    # for(j in 1:nbin) sdyest[i,j] = sqrt(1/rgamma(
    #    1,shape_sdyest[j],B_sdyest+0.5*sum((yobs[[j]]-yestimate[i-1,j])^2)))
    for(j in 1:nbin) sdyest[i,j] = sqrt(1/rgamma(
      1,max(c(shape_sdyest[j],A_sdyest),na.rm=T),
      max(c(B_sdyest+0.5*(ymean[j]-yestimate[i-1,j])^2, B_sdyest),na.rm=T)))
    #https://stats.stackexchange.com/questions/525674/gibbs-sampler-for-normal-and-inverse-gamma-distribution-in-r
    # https://stats.stackexchange.com/questions/266665/gibbs-sampler-examples-in-r

    ## 3. Metropolis-Hastings step to estimate yest
    proposal_yest = MH_propose_yest(yestimate[i-1,],prop_sd_yest =  prop_sd_yest) # new proposed values

    ## 4. Metropolis-Hastings step to estimate the regression coefficients
    proposal_coeff = MH_propose_coeff(coefficients[i-1,],prop_sd =  prop_sd_coeff) # new proposed values
    if(any(proposal_coeff[4] <= 0)) HR = 0 else # Q needs to be >0
      # Hastings ratio of the proposal
      logpostold = logposterior(x = x, ymean = ymean, yest = yestimate[i-1,],
                                sdyest = sdyest[i,], coeff = coefficients[i-1,],
                                sdy = sdy[i], prior_df = prior_df)
    logpostnew = logposterior(x = x, ymean = ymean, yest = proposal_yest,
                              sdyest = sdyest[i,], coeff = proposal_coeff,
                              sdy = sdy[i], prior_df = prior_df)
    HR = exp(logpostnew -
               logpostold)
    # accept proposal with probability = min(HR,1)
    if (runif(1) < HR){
      coefficients[i,] = proposal_coeff
      yestimate[i,] = proposal_yest
      logpost[i] = logpostnew
      # if proposal is rejected, keep the values from the previous iteration
    }else{
      coefficients[i,] = coefficients[i-1,]
      yestimate[i,] = yestimate[i-1,]
      logpost[i] = logpostold

    }
  } # end of the MCMC loop

  ###  Function output
  output = list(data.frame(A = coefficients[,1],
                           K = coefficients[,2],
                           M = coefficients[,3],
                           Q = coefficients[,4],
                           sdy = sdy,
                           logpost = logpost),
                yestimate = yestimate,
                sdyest = sdyest)
  return(output)
}


####
####  Testing the model
####

#
# Create test data
#

nbin = 6 # number of latitudinal bins / locations

x <- seq(10,60,10) # latitudes of samples
npb <- c(0,2,4,3,1,0) # number of samples per bin/location
ym <- c(30,30,20,18,18,8) # mean temperature of bin/location
y <- lapply(1:nbin, function(x) rnorm(npb[x],ym[x],2)) # create samples as a normal distirbution
y[which(npb==0)] <- NA # set to NA if no samples are from bin/location

## Show d18O data
plot(0,0,xlim = c(10,70), ylim = c(-1,35))
for(i in (1:nbin)) if(npb[i] !=0)points(rep(x[i],npb[i]),y[[i]])

## Additional priors (geological information)
prior_df <- data.frame(matrix(NA,nrow=nbin, ncol = 3))
prior_df[1:3,1] <- "dunif"
prior_df[1:3,2] <- 21
prior_df[1:3,3] <- 35
prior_df[6,1] <- "dnorm"
prior_df[6,2:3] <- c(5,2)
class(prior_df[,3])
#extra_priors(1,1,prior_df)


## Set intitial values and proposal properties for MCMC
prop_sd_yest <- matrix(0.01,nrow = nbin, ncol = nbin) # proposal covariance matrix
diag(prop_sd_yest) <- 1 # set sd to 1
prop_sd_coeff <- c(1,1,1,0.02) # proposal sd for regression coefficients

coeff_inits = c(0,30,45,0.1) # initial values for regression coefficients
yest_inits = rep(25,6) # initial values for yest (need to be within the range allowed by priors!)
sdyest_inits = rep(2,nbin) # initial values for sdyest
sdy_init = 1 # initial values for sdy

burnin = 20001
nIter = 100000

## run the model
m1 <-  run_MCMC(nIter = nIter, x = x, yobs = y, prior_df = prior_df,
                 coeff_inits = coeff_inits,
                 sdy_init = 1, yest_inits = yest_inits,
                 sdyest_inits = sdyest_inits,
                 prop_sd_coeff=prop_sd_coeff, prop_sd_yest=prop_sd_yest,
                 nbin = nbin)



### Plot output and data

burnin = burnin
nIter = nIter
mn <- m22
plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(mn[[1]][burnin:nIter,1:4],2,median), 0),
     ylim = c(-1,35), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
### Calculate confidence intervals
latitude <- 0:90
sample_it <- sample((burnin+1):nIter,3000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))
error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:nbin) if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dunif")  lines(x[c(i,i)],c(prior_df[i,2],prior_df[i,3]), lwd = 7, col = rgb(0,0.9,0,0.5))
for(i in 1:nbin) {if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dnorm")  {beanplot::beanplot(rnorm(1000,prior_df[i,2],prior_df[i,3]), add = T,
                                                                                            at = x[i], maxwidth = 5, side = "second",
                                                                                            what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.5), border = NA)}}


points(x,apply(mn[[2]][burnin:nIter,],2,median), col = "red", pch = 19, cex = 1.25)
sapply(1:nbin, function(a) points(c(x[a],x[a]),quantile(mn[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(1,0,0,0.75), lwd =2))
#sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(mn[[2]][burnin:nIter,a])+c(-2,2)*median(mn[[3]][burnin:nIter,a])), type = "l", col= "red"))

for(i in (1:nbin)) if(npb[i] !=0) points(rep(x[i],npb[i])+rnorm(npb[i],0,0.25),y[[i]], col = rgb(0,0.3,1,0.55), pch = 17)

legend("topright",legend=c("regression line", "temperature estimate","dO18 data", "coral reef range (uniform prior)", "additional info. (normal prior)"),
       cex = 0.8, col = c("black", "red", rgb(0,0.3,1,0.55),  rgb(0,0.9,0,0.5), rgb(0,0.7,0.7,0.5)), lwd = c(2,NA,NA,4,4),
       pch = c(NA,19,17,NA,NA), pt.cex = c(NA,1.25,1,NA,NA))



mn <- output
ind <- seq(1,50000,10)
plot(ind,mn[[1]][ind,1], type = "l")
plot(ind,mn[[2]][ind,2], type = "l")

plot(ind,mn[[2]][ind,5], type = "l")

plot(mn[[2]][ind,1],mn[[2]][ind,2])
abline(a=0,b=1)
plot(m[[1]][ind,2], type = "l")

plot(m[[3]][ind,6])

ind <- 5
hist(rnorm(10000,median(m[[2]][,ind]),median(m[[3]][,ind])),100)
abline(v = y[[ind]], lwd = 2, col = "turquoise")

burnin = 20000+1
nIter = 100000
plot(seq(10,70,0.1), gradient(seq(10,70,0.1), apply(m[[1]][burnin:nIter,1:4],2,median), 0), ylim = c(-1,35))
points(x,apply(m[[2]][burnin:nIter,],2,median), col = "red", pch = 19)
sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(m[[2]][burnin:nIter,a])+c(-2,2)*median(m[[3]][burnin:nIter,a])), type = "l", col= "red"))
for(i in 1:nbin) points(rep(x[i],7),y[[i]], col = rgb(0,0.8,0.6,0.7))

points(x,apply(m[[2]][burnin:nIter,],2,median), col = rgb(0.5,1,0,0.5), pch = 19)
points(x,apply(m2[[2]][burnin:nIter,],2,median), col = rgb(0.5,0,1,0.5), pch = 19)

burnin = 20000+1
nIter = 100000
plot(seq(10,70,0.1), gradient(seq(10,70,0.1), apply(m2[[1]][burnin:nIter,1:4],2,median), 0), ylim = c(-1,35))
points(x,apply(m2[[2]][burnin:nIter,],2,median), col = "red", pch = 19)
sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(m2[[2]][burnin:nIter,a])+c(-2,2)*median(m2[[3]][burnin:nIter,a])), type = "l", col= "red"))
for(i in 1:nbin) points(rep(x[i],7),y[[i]], col = rgb(0,0.8,0.6,0.7))


plot(seq(10,70,0.1), gradient(seq(10,70,0.1), apply(m3[[1]][burnin:nIter,1:4],2,median), 0), ylim = c(-1,35))
points(x,apply(m3[[2]][burnin:nIter,],2,median), col = "red", pch = 19)
sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(m3[[2]][burnin:nIter,a])+c(-2,2)*median(m3[[3]][burnin:nIter,a])), type = "l", col= "red"))
for(i in 1:nbin) points(rep(x[i],7),y[[i]], col = rgb(0,0.8,0.6,0.7))

burnin = 20000+1
nIter = 100000
mn <- m22
plot(seq(0,90,0.1), gradient(seq(0,90,0.1), apply(mn[[1]][burnin:nIter,1:4],2,median), 0),
     ylim = c(-1,35), type = "l", lwd = 3, ylab = "Temperature", xlab = "Latitude")
### Calculate confidence intervals
latitude <- 0:90
sample_it <- sample((burnin+1):nIter,3000)
grad_025 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.025))
grad_975 <- sapply(1:length(latitude), function(f) quantile(apply(mn[[1]][sample_it,1:4],1,function(a) gradient(x=latitude[f], coeff =  a, sdy = 0)), probs = 0.975))
error_polygon(latitude,grad_025,grad_975,rgb(0,0,0,0.15))

for(i in 1:nbin) if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dunif")  lines(x[c(i,i)],c(prior_df[i,2],prior_df[i,3]), lwd = 7, col = rgb(0,0.9,0,0.5))
for(i in 1:nbin) {if(!(is.na(prior_df[i,1])) & prior_df[i,1]=="dnorm")  {beanplot::beanplot(rnorm(1000,prior_df[i,2],prior_df[i,3]), add = T,
                                                                                            at = x[i], maxwidth = 5, side = "second",
                                                                                            what = c(0,1,0,0), col = rgb(0,0.7,0.7,0.5), border = NA)}}


points(x,apply(mn[[2]][burnin:nIter,],2,median), col = "red", pch = 19, cex = 1.25)
sapply(1:nbin, function(a) points(c(x[a],x[a]),quantile(mn[[2]][burnin:nIter,a], probs = c(0.025,0.975)), type = "l", col= rgb(1,0,0,0.75), lwd =2))
#sapply(1:nbin, function(a) points(c(x[a],x[a]),c(median(mn[[2]][burnin:nIter,a])+c(-2,2)*median(mn[[3]][burnin:nIter,a])), type = "l", col= "red"))

for(i in (1:nbin)) if(npb[i] !=0) points(rep(x[i],npb[i])+rnorm(npb[i],0,0.25),y[[i]], col = rgb(0,0.3,1,0.55), pch = 17)

legend("topright",legend=c("regression line", "temperature estimate","dO18 data", "coral reef range (uniform prior)", "additional info. (normal prior)"),
       cex = 0.8, col = c("black", "red", rgb(0,0.3,1,0.55),  rgb(0,0.9,0,0.5), rgb(0,0.7,0.7,0.5)), lwd = c(2,NA,NA,4,4),
       pch = c(NA,19,17,NA,NA), pt.cex = c(NA,1.25,1,NA,NA))



###  Something is wrong with the yestimate: they move more with larger sample size. fix!
