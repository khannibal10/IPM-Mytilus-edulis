##########################################################################
# Mytilus edulis Integral Projection Model - Control Script
# Calculations using Functions of IPM kernels P, F, and K
# Created by Steve Dudgeon
# Created on  04 / 14 /2021
# edited on  04 / 15 /2021
##########################################################################
## clear workspace------------------------------
rm(list=ls());gc()

# load libraries---------------------



# set working directory-------------------------
file.choose()
setwd("/Users/sd51881/Documents/Steve/NSF_LTREB_2016_2021/Mussel_patches_Maine_Scotland_Germany/integral projection model_ME/")

## run the utility functions
source("Standard Graphical Pars.R")

## run the Mytilus kernel functions
source("Mytilus_IPM_kernel functions.R")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Construct kernels and compute population growth rates  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nBigMatrix <- 100

min.size <- 0
max.size <- 100

# Run model during August SST
tempC <- 14

# Run IPM
IPM.est <- mk_K(nBigMatrix, m.par, min.size, max.size)

## calculate the population growth rate
lam.est <- Re(eigen(IPM.est$K)$values[1])
lam.est

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 2 - Construct kernels, compute useful summary quantities and plot size by age distribution
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## extract the mesh points
meshpts <- IPM.est$meshpts

## normalised stable size distribution
w.est <- Re(eigen(IPM.est$K)$vectors[, 1])
stable.z.dist.est <- w.est/sum(w.est)
stable.z.dist.est
## uncertain about this next line given that we modeled open population
stable.z.repr.dist.est <- s_z(meshpts, m.par) * R_t(meshpts, m.par) * w.est/sum(s_z(meshpts, m.par) * R_t(meshpts, m.par) * w.est)
stable.z.repr.dist.est
## mean size
mean.z.est <- sum(stable.z.dist.est * meshpts)
mean.z.est

## mean size at reproduction - again, uncertain about this given open population
mean.z.repr.est <- sum(stable.z.repr.dist.est * meshpts)
mean.z.repr.est

## variance in size
var.z.est <- sum(stable.z.dist.est * meshpts^2) - mean.z.est^2
var.z.est


## compute the size distribution for each age class with a loop
z.dist.by.age <- list()
z.dist.by.age[[1]] <- IPM.est$F %*% stable.z.dist.est/lam.est
for (i in 2:50) z.dist.by.age[[i]] <- IPM.est$P %*% z.dist.by.age[[i - 1]]/lam.est

## build a little helper function to compute the means & variances
mk_moments <- function(z.dist, meshpts) {
  z.dist <- z.dist/sum(z.dist)
  mean.z <- sum(z.dist * meshpts)
  var.z <- sum(z.dist * meshpts^2) - mean.z^2
  return(c(mean = mean.z, sd = sqrt(var.z)))
}

## 
postscript(file = "MytilusSizeAge.eps", w = 6, h = 6, horizontal = FALSE, onefile = FALSE, 
           paper = "special")
set_graph_pars(ptype = "panel1")
## age = 0 (new recruits)
z.dist <- z.dist.by.age[[1]]
plot(meshpts, z.dist, type = "n", xlab = "Length, z", ylab = "Density", xlim = c(0, 100), ylim = c(0, 0.0001))
## age = 1, 2, 3 and 4
for (A in 1:5) {
  z.dist <- z.dist.by.age[[A + 1]]
  lines(meshpts, z.dist)
  moments <- round(mk_moments(z.dist, meshpts), 2)
  text(x = moments["mean"], y = max(z.dist) + 5e-05, pos = 4, cex = 0.75, labels = paste("A = ", 
                                                                                         A, " (mean = ", moments["mean"], ", s.d. = ", moments["sd"], ")", sep = ""))
}
dev.off()


