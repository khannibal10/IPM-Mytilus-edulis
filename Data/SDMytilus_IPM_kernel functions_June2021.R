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





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Define the demographic functions and parameters
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the parameter vector based on statistical models of survivorship, growth and fecundity 
# as functions of size as mussel shell length (mm). 
# each parameter is given a name so the formulae, below, are easier to
# read. We'll use 'model.term' scheme to name elements of the vector

m.par <- c(## survival
  surv.int  = -4.52e+0 ,
  surv.z    =  1.52e-1 ,
  ## growth 
  Linf = 80 ,  # maximum size of Mytilus edulis
  k = 0.11878354 ,  # von Bertalanffy growth parameter, k
  grow.sd   =  1.63e-2 , # standard error of the growth parameter, k
 
   ## recruitment intensity in "open" population as a function of August water temperature (°C)
             
  ## recruit probability intercept
  recr.int  = 1.9756e+10 ,  
  ## recruitment slope versus SST °C
  recr.T = -9.529e-1, # beta coefficient (slope) from log10 transformed Poisson model fit of recruitment per m^2 vs August Temp 
  tempC  = 14 , #August seawater temperature during settlement, °C from Communications Biology 2020 paper
  p_r = 0.05 ,  # survisorship of recruits from year t to t+1
  ## recruit size
  rcsz.int  =  2.19e-1 ,  #  from Katie's data
  #rcsz.z    =  constant,  # from Katie's data
  rcsz.sd   =  5.21e-2 )  # from Katie's data        


## functions -----------------------------------
## Define the the three demographic functions, we pass a vector of parameters "m.par" to each function

## Growth function, given you are size z now returns the pdf of size z1 next time

g_z1z <- function(z1, z, m.par)
{
  mu <- m.par["Linf"] - (m.par["Linf"] - z ) * exp(-m.par["k"])  # mean size next year
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

## Define the reproduction function
## Reproduction function, exponential function of decline in arriving recruits in year t with increasing temp (°C)

R_t <- function(z1)
{
  linear.p <- m.par["recr.int"] * exp( m.par["recr.T"] * m.par["tempC"] )     # exponential predictor
  p <- linear.p * m.par["p_r"]        # recruitment to year `1 from those in previous year at the given August Temp`
  return(p)
}

## Recruit size function

c_z1z <- function(z1, z, m.par)
{
  mu <- m.par["rcsz.int"]                                 # mean size of recruits next year
  sig <- m.par["rcsz.sd"]                                 # sd about mean
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)            # pdf that offspring are size z1 given you were size z
  return(p.den.rcsz)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Functions to build IPM kernels P and K
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Define the survival kernel
P_z1z <- function (z1, z, m.par) {
  
  return( s_z(z, m.par) * g_z1z(z1, z, m.par) )
  
}

## Define the reproduction distribution
F_z1z <- function (z, m.par) {
  
  return( R_t(m.par) * c_z1z(z1, z, m.par) )
  
}

## Build the discretized kernel
mk_K <- function(m, m.par, L, U) {
  
  # mesh points 
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h
  K <- h * (outer(meshpts, meshpts, P_z1z, m.par = m.par))
  return(list(K = K, meshpts = meshpts))
}

