## ########################################################################
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
#setwd("/Users/khannibal/Documents/IPM-Mytilus-edulis/IPM-Mytilus-edulis.Rproj")
## run the utility functions
source("Standard Graphical Pars.R")

## run the Mytilus kernel functions
source("Data/SDMytilus_IPM_kernel functions_June2021.R")

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Section 1 - Construct kernels and compute population growth rates  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nBigMatrix <- 100

min.size <- 0
max.size <- 100

# initial size distribution
#Original z<-numeric(100)
Z <- seq(0, 100, length.out=100) #this is creating an initial numeric vector/range of numbers dynamically, so use seq() instead.
Z 


# Run IPM
#Original  IPM.est <- mk_K(nBigMatrix, m.par, min.size, max.size) * z +F_z1z
IPM.est <- mk_K(nBigMatrix, m.par, min.size, max.size)[Z] #non-numeric argument to binary operator with *, so using []

#Adding Z and F
#OG:IPM.est + F_z1z======Error in IPM.est + F_z1z : non-numeric argument to binary operator
IPM.est<-IPM.est + (F_z1z())
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
# ERROR Error in R_t(meshpts, m.par) : unused argument (m.par)

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
plot(meshpts, z.dist, type = "n", xlab = "Length, z", ylab = "Density", xlim = c(0, 100), ylim = c(0, 0.00000002))
## age = 1, 2, 3 and 4
for (A in 0:49) {
  z.dist <- z.dist.by.age[[A + 1]]
  lines(meshpts, z.dist)
  moments <- round(mk_moments(z.dist, meshpts), 2)
  text(x = moments["mean"], y = max(z.dist) + 5e-05, pos = 4, cex = 0.75, labels = paste("A = ", 
                                                                                         A, " (mean = ", moments["mean"], ", s.d. = ", moments["sd"], ")", sep = ""))
}
dev.off()


