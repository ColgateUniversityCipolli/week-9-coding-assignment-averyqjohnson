###########################################################################
# HW 9
# Avery Johnson
# MATH 240 - SPRING 2025
###########################################################################
library(tidyverse)

################################################################################
# Precipitation in Madison County
################################################################################
dat.precip <- read_csv(file = "agacis.csv")

#####################################
# Clean Data
#####################################
dat.precip.long <- dat.precip |>    
  dplyr::select(-Annual) |>                   # Remove annual column 
  pivot_longer(cols = c(Jan, Feb, Mar, Apr,   # pivot the column data into one col
                        May, Jun, Jul, Aug, 
                        Sep, Oct, Nov, Dec), 
               values_to = "Precipitation",   # store the values in Precipitation
               names_to = "Month") |>         # store the months in Month
  mutate(Precipitation = case_when(Precipitation == "M" ~ NA_character_,
                                   TRUE                 ~ Precipitation))|>
  mutate(Precipitation = as.numeric(Precipitation))

#####################################
# Summarize the data
#####################################

library(e1071)
dat.precip.long |>
  summarize(
    mean = mean(Precipitation, na.rm=T),
    sd = sd(Precipitation, na.rm=T),
    min = min(Precipitation, na.rm=T),
    max = max(Precipitation, na.rm=T),
    skew = skewness(Precipitation, na.rm=T),
    kurt = kurtosis(Precipitation, na.rm=T)
  )

###########################################################################
# Weibull
###########################################################################
llweibull <- function(par, data, neg=F){
  # a <- par[1]
  # sigma <- par[2]
  a <- exp(par[1]) # go from (-inf,inf) to (0,inf)
  sigma <- exp(par[2]) # go from (-inf,inf) to (0,inf)
  
  ll <- sum(log(dweibull(x=data, shape=a, scale=sigma)), na.rm=T)
  
  return(ifelse(neg, -ll, ll))
}
MLEs <- optim(fn = llweibull,
              par = c(1,1),
              data = dat.precip.long$Precipitation,
              neg=T)

(MLEs$par <- exp(MLEs$par)) # transform

###########################################################################
# a) Compute the MLEs for these data using a Gamma distribution
###########################################################################
llgamma <- function(par, data, neg=F){
  alpha <- exp(par[1])
  beta <- exp(par[2])
  
  ll <- sum(log(dgamma(x=data, shape=alpha, scale=beta)), na.rm=T)
  
  return(ifelse(neg, -ll, ll))
}

MLEs_gamma <- optim(fn = llgamma,
              par = c(1,1),
              data = dat.precip.long$Precipitation,
              neg=T)

(MLEs_gamma$par <- exp(MLEs_gamma$par)) # transform

###########################################################################
# b) Compute the MLEs for these data using the Log-Normal distribution
###########################################################################

lllognorm <- function(par, data, neg=F){
  mu <- par[1]
  sigma <- exp(par[2])
  
  ll <- sum(log(dlnorm(x=data, meanlog=mu, sdlog=sigma)), na.rm=T)
  
  return(ifelse(neg, -ll, ll))
}

MLEs_lognorm <- optim(fn = lllognorm,
                    par = c(1,1),
                    data = dat.precip.long$Precipitation,
                    neg=T)

(MLEs_lognorm$par[1])
(exp(MLEs_lognorm$par[2]))
###########################################################################
# c)  Compute the likelihood ratio to compare the Weibull and the Gamma dist
###########################################################################

# compute the log-likelihood values at MLEs
ll_weibull <- -2166.496
ll_gamma <- -MLEs_gamma$value

#compute the likelihood ratio
(ratio_w_g <- exp(ll_weibull-ll_gamma))

###########################################################################
# d)  Compute the likelihood ratio to compare the Weibull and the Lognormal dist
###########################################################################

# compute the log-likelihood values at MLEs
ll_lognorm <- -MLEs_lognorm$value

#compute the likelihood ratio
(ratio_w_ln <- exp(ll_weibull-ll_lognorm))

###########################################################################
# e)  Compute the likelihood ratio to compare the Gamma and the Lognormal dist
###########################################################################
(ratio_g_ln <- exp(ll_gamma-ll_lognorm))
