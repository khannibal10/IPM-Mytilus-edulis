##########################################################################
# Mytilus edulis Integral Projection Model
# Statistical models describing survivorship, growth & fecundity
#   & Functions to build IPM kernels P, F, and K
# Created by Steve Dudgeon
# Created on  04 / 09 /2021
# edited on   /  /2021
##########################################################################


# load libraries---------------------



# set working directory-------------------------
file.choose()




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the parameter vector based on statistical models of survivorship, growth and fecundity 
   # as functions of size as mussel shell length (mm). 
   # each parameter is given a name so the formulae, below, are easier to
   # read. We'll use 'model.term' scheme to name elements of the vector
m.par <- c(## survival
  surv.int  = -4.52e+0,
  surv.z    =  1.52e-1,
  ## growth 
  grow.int  =  9.00e+0,
  grow.z    =  8.89e-1,
  grow.sd   =  1.95e-0,
  ## recruitment intensity in "open" population as a function of August water temperature (°C)
  repr.int  = 1.9756e+10, # intercept from exponential model fit of recruitment per m^2 vs August Temp
  repr.z    =  -9.529e-1, # beta coefficient (slope) from exponential model fit of recruitment per m^2 vs August Temp              
  ## recruit or not
  recr.int  = -2.94e-0 ,  # linear predictor to get probability of 0.05 of successful recruitment through to year 1
  ## recruit size
  rcsz.int  =  2.19e-1,  #  from Katie's data
  #rcsz.z    =  constant,  # from Katie's data
  rcsz.sd   =  5.21e-2,  # from Katie's data
  #August seawater temperature during settlement, °C
  tempC  = 14)  # from Communications Biology 2020 paper

## functions -----------------------------------
## Define the the three demographic functions, we pass a vector of parameters "m.par" to each function

## Growth function, given you are size z now returns the pdf of size z1 next time

g_z1z <- function(z1, z, m.par)
{
  mu <- m.par["grow.int"] + m.par["grow.z"] * z           # mean size next year
  sig <- m.par["grow.sd"]                                 # sd about mean
  p.den.grow <- dnorm(z1, mean = mu, sd = sig)            # pdf that you are size z1 given you were size z
  return(p.den.grow)
}

## Survival function, logistic regression

s_z <- function(z, m.par)
{
  linear.p <- m.par["surv.int"] + m.par["surv.z"] * z       # linear predictor
  p <- (0.38/(1+exp(-1.5*linear.p)))+0.05                   # logistic transformation to probability
  return(p)
}

## Reproduction function, exponential function of decline in arriving recruits in year t with increasing temp (°C)

R_t <- function(tempC, m.par)
{
  tempC <- 14
  r_t <- m.par["repr.int"] * exp(m.par["repr.z"] * tempC)  # exponential decline in arrivals with increasing temp
  return(r_t)
}

## Recruitment probability function (from birth in spring of year t through to first summer year t+1), 
## logistic regression-intercept only to get 5% survivorship (can modify to include covariate)

pr_z <- function(m.par)
{
  linear.p <- m.par["recr.int"]                             # linear predictor
  p <- 1/(1+exp(-linear.p))                                 # logistic transformation to probability
  return(p)
}

## Recruit size function

c_z1z <- function(z1, z, m.par)
{
  mu <- m.par["rcsz.int"]                                 # mean size of recruits next year
  sig <- m.par["rcsz.sd"]                                 # sd about mean
  z <- dnorm(z1, mean = mu, sd = sig)            # pdf that offspring are size z1 given you were size z
  return(z)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P, F, and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Define the survival kernel
P_z1z <- function (z1, z, m.par) {
  
  return( s_z(z, m.par) * g_z1z(z1, z, m.par) )
  
}

## Define the reproduction kernel
F_z1z <- function (z1, z, m.par, tempC) {
  
  return( R_t(tempC, m.par) * pr_z(m.par) * c_z1z(z1, z, m.par) )
  
}

## Build the discretized kernel
mk_K <- function(m, m.par, L, U) {
  
  # mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  P <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
  return(list(P = P, meshpts = meshpts))
}

