## =================================================================================================================================================================================================================

## R code to generate Figure S1 for manuscript titled "How public reaction to disease information across scales and the impacts of vector control methods influence disease prevalence and control efficacy"

## ==================================================================================================================================================================================================================


library(ggplot2)
library(reshape)
library(here)

## set the current directory as default to save all the following results
set_here()

N <- 9
rep <- 50

epsilon <- c(100, 100)

##find CI for each time point
confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}


##################################################################################################################################

# K = 500 local  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-local-N-9-p-0.01-t-300-K=500-strengthonly")

Ks500locallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks500localadultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks500CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()
for ( i in 1:dim(Ks500localadultsize)[1])
{
  if (Ks500localadultsize[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks500localadultsize)[1]-KK[length(KK)]+1

diff <- c(diff(KK), lastlength)

Ks500meansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks500local <- Ks500localadultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks500local <- Ks500localadultsize[KK[rep]:dim(Ks500localadultsize)[1],1]
}

for (i in 1:(rep-1))
{
  Ks500dataset1 <- Ks500locallarvaesize[KK[i]:(KK[i+1]-1),]
  Ks500dataset2 <- Ks500localadultsize[KK[i]:(KK[i+1]-1),]
  Ks500datasize1 <- Ks500dataset1[,-(which.min(Ks500CIdistance[i,])+1)] 
  Ks500datasize2 <- Ks500dataset2[,-(which.min(Ks500CIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks500meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500datasize1[,2:ncol(Ks500datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500datasize2[,2:ncol(Ks500datasize2)])/(N-1))
  Ks500meansizesum1<-approx(Ks500datasize1[,1],Ks500meansizesum, x_axis_Ks500local)$y
  
  Ks500meansizetotal[i,]<-Ks500meansizesum1
}


Ks500dataset1 <- Ks500locallarvaesize[KK[rep]:attributes(as.matrix(Ks500locallarvaesize))$dim[1],]
Ks500dataset2 <- Ks500localadultsize[KK[rep]:attributes(as.matrix(Ks500localadultsize))$dim[1],]
Ks500datasize1 <- Ks500dataset1[,-(which.min(Ks500CIdistance[rep,])+1)] 
Ks500datasize2 <- Ks500dataset2[,-(which.min(Ks500CIdistance[rep,])+1)] 
##this 1 is the time column
Ks500meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500datasize1[,2:ncol(Ks500datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500datasize2[,2:ncol(Ks500datasize2)])/(N-1))
Ks500meansizesum1<-approx(Ks500datasize1[,1],Ks500meansizesum, x_axis_Ks500local)$y

Ks500meansizetotal[rep,]<-Ks500meansizesum1

Ks500averagesize <- colSums(Ks500meansizetotal)/rep

Ks500sizeCIlow<-c()
Ks500sizeCIhigh<-c()

for(i in 1:length(Ks500averagesize))
{
  Ks500sizeCIlow[i]<-confidence_interval(Ks500meansizetotal[,i], 0.95)[1][[1]]
  Ks500sizeCIhigh[i]<-confidence_interval(Ks500meansizetotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis SF-local-N-9-p-0.01-t-1000-K=500-time+strength+control-9-23new
setwd("~\\SF-local-N-9-p-0.01-t-300-K=500-time+strength+control")

KT500locallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT500localadultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT500CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()
for ( i in 1:dim(KT500localadultsize)[1])
{
  if (KT500localadultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT500localadultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT500meansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT500local <- KT500localadultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT500local <- KT500localadultsize[KK[rep]:dim(KT500localadultsize)[1],1]
}


for (i in 1:(rep-1))
{
  KT500dataset1 <- KT500locallarvaesize[KK[i]:(KK[i+1]-1),]
  KT500dataset2 <- KT500localadultsize[KK[i]:(KK[i+1]-1),]
  KT500datasize1 <- KT500dataset1[,-(which.min(KT500CIdistance[i,])+1)] 
  KT500datasize2 <- KT500dataset2[,-(which.min(KT500CIdistance[i,])+1)] 
  ##this 1 is the time column
  KT500meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500datasize1[,2:ncol(KT500datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500datasize2[,2:ncol(KT500datasize2)])/(N-1))
  KT500meansizesum1<-approx(KT500datasize1[,1],KT500meansizesum, x_axis_KT500local)$y
  
  KT500meansizetotal[i,]<-KT500meansizesum1
}


KT500dataset1 <- KT500locallarvaesize[KK[rep]:attributes(as.matrix(KT500locallarvaesize))$dim[1],]
KT500dataset2 <- KT500localadultsize[KK[rep]:attributes(as.matrix(KT500localadultsize))$dim[1],]
KT500datasize1 <- KT500dataset1[,-(which.min(KT500CIdistance[rep,])+1)] 
KT500datasize2 <- KT500dataset2[,-(which.min(KT500CIdistance[rep,])+1)] 
##this 1 is the time column
KT500meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500datasize1[,2:ncol(KT500datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500datasize2[,2:ncol(KT500datasize2)])/(N-1))
KT500meansizesum1<-approx(KT500datasize1[,1],KT500meansizesum, x_axis_KT500local)$y


KT500meansizetotal[rep,]<-KT500meansizesum1


KT500averagesize <- colSums(KT500meansizetotal)/rep

KT500sizeCIlow<-c()
KT500sizeCIhigh<-c()

for(i in 1:length(KT500averagesize))
{
  KT500sizeCIlow[i]<-confidence_interval(KT500meansizetotal[,i], 0.95)[1][[1]]
  KT500sizeCIhigh[i]<-confidence_interval(KT500meansizetotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 500 neighbor  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=500-strengthonly")

Ks500neighborlarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks500neighboradultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks500neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks500neighboradultsize)[1])
{
  if (Ks500neighboradultsize[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks500neighboradultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks500neighbormeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks500neighbor <- Ks500neighboradultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks500neighbor <- Ks500neighboradultsize[KK[rep]:dim(Ks500neighboradultsize)[1],1]
}

for (i in 1:(rep-1))
{
  Ks500neighbordataset1 <- Ks500neighborlarvaesize[KK[i]:(KK[i+1]-1),]
  Ks500neighbordataset2 <- Ks500neighboradultsize[KK[i]:(KK[i+1]-1),]
  Ks500neighbordatasize1 <- Ks500neighbordataset1[,-(which.min(Ks500neighborCIdistance[i,])+1)] 
  Ks500neighbordatasize2 <- Ks500neighbordataset2[,-(which.min(Ks500neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks500neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500neighbordatasize1[,2:ncol(Ks500neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500neighbordatasize2[,2:ncol(Ks500neighbordatasize2)])/(N-1))
  Ks500neighbormeansizesum1<-approx(Ks500neighbordatasize1[,1],Ks500neighbormeansizesum, x_axis_Ks500neighbor)$y
  
  Ks500neighbormeansizetotal[i,]<-Ks500neighbormeansizesum1
}


Ks500neighbordataset1 <- Ks500neighborlarvaesize[KK[rep]:attributes(as.matrix(Ks500neighborlarvaesize))$dim[1],]
Ks500neighbordataset2 <- Ks500neighboradultsize[KK[rep]:attributes(as.matrix(Ks500neighboradultsize))$dim[1],]
Ks500neighbordatasize1 <- Ks500neighbordataset1[,-(which.min(Ks500neighborCIdistance[rep,])+1)] 
Ks500neighbordatasize2 <- Ks500neighbordataset2[,-(which.min(Ks500neighborCIdistance[rep,])+1)] 
##this 1 is the time column
Ks500neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500neighbordatasize1[,2:ncol(Ks500neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500neighbordatasize2[,2:ncol(Ks500neighbordatasize2)])/(N-1))
Ks500neighbormeansizesum1<-approx(Ks500neighbordatasize1[,1],Ks500neighbormeansizesum, x_axis_Ks500neighbor)$y


Ks500neighbormeansizetotal[rep,]<-Ks500neighbormeansizesum1


Ks500neighboraveragesize <- colSums(Ks500neighbormeansizetotal)/rep

Ks500neighborsizeCIlow<-c()
Ks500neighborsizeCIhigh<-c()

for(i in 1:length(Ks500neighboraveragesize))
{
  Ks500neighborsizeCIlow[i]<-confidence_interval(Ks500neighbormeansizetotal[,i], 0.95)[1][[1]]
  Ks500neighborsizeCIhigh[i]<-confidence_interval(Ks500neighbormeansizetotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=500-time+strength+control")

KT500neighborlarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT500neighboradultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT500neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT500neighboradultsize)[1])
{
  if (KT500neighboradultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT500neighboradultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT500neighbormeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT500neighbor <- KT500neighboradultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT500neighbor <- KT500neighboradultsize[KK[rep]:dim(KT500neighboradultsize)[1],1]
}

for (i in 1:(rep-1))
{
  KT500neighbordataset1 <- KT500neighborlarvaesize[KK[i]:(KK[i+1]-1),]
  KT500neighbordataset2 <- KT500neighboradultsize[KK[i]:(KK[i+1]-1),]
  KT500neighbordatasize1 <- KT500neighbordataset1[,-(which.min(KT500neighborCIdistance[i,])+1)] 
  KT500neighbordatasize2 <- KT500neighbordataset2[,-(which.min(KT500neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT500neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500neighbordatasize1[,2:ncol(KT500neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500neighbordatasize2[,2:ncol(KT500neighbordatasize2)])/(N-1))
  KT500neighbormeansizesum1<-approx(KT500neighbordatasize1[,1],KT500neighbormeansizesum, x_axis_KT500neighbor)$y
  
  KT500neighbormeansizetotal[i,]<-KT500neighbormeansizesum1
}

KT500neighbordataset1 <- KT500neighborlarvaesize[KK[rep]:attributes(as.matrix(KT500neighborlarvaesize))$dim[1],]
KT500neighbordataset2 <- KT500neighboradultsize[KK[rep]:attributes(as.matrix(KT500neighboradultsize))$dim[1],]
KT500neighbordatasize1 <- KT500neighbordataset1[,-(which.min(KT500neighborCIdistance[rep,])+1)] 
KT500neighbordatasize2 <- KT500neighbordataset2[,-(which.min(KT500neighborCIdistance[rep,])+1)] 
##this 1 is the time column
KT500neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500neighbordatasize1[,2:ncol(KT500neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500neighbordatasize2[,2:ncol(KT500neighbordatasize2)])/(N-1))
KT500neighbormeansizesum1<-approx(KT500neighbordatasize1[,1],KT500neighbormeansizesum, x_axis_KT500neighbor)$y


KT500neighbormeansizetotal[rep,]<-KT500neighbormeansizesum1


KT500neighboraveragesize <- colSums(KT500neighbormeansizetotal)/rep

KT500neighborsizeCIlow<-c()
KT500neighborsizeCIhigh<-c()

for(i in 1:length(KT500neighboraveragesize))
{
  KT500neighborsizeCIlow[i]<-confidence_interval(KT500neighbormeansizetotal[,i], 0.95)[1][[1]]
  KT500neighborsizeCIhigh[i]<-confidence_interval(KT500neighbormeansizetotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 500 global  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=500-strengthonly")

Ks500globallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks500globaladultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks500globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks500globaladultsize)[1])
{
  if (Ks500globaladultsize[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks500globaladultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks500globalmeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks500global <- Ks500globaladultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks500global <- Ks500globaladultsize[KK[rep]:dim(Ks500globaladultsize)[1],1]
}

for (i in 1:(rep-1))
{
  Ks500globaldataset1 <- Ks500globallarvaesize[KK[i]:(KK[i+1]-1),]
  Ks500globaldataset2 <- Ks500globaladultsize[KK[i]:(KK[i+1]-1),]
  Ks500globaldatasize1 <- Ks500globaldataset1[,-(which.min(Ks500globalCIdistance[i,])+1)] 
  Ks500globaldatasize2 <- Ks500globaldataset2[,-(which.min(Ks500globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks500globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500globaldatasize1[,2:ncol(Ks500globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500globaldatasize2[,2:ncol(Ks500globaldatasize2)])/(N-1))
  Ks500globalmeansizesum1<-approx(Ks500globaldatasize1[,1],Ks500globalmeansizesum, x_axis_Ks500global)$y
  
  Ks500globalmeansizetotal[i,]<-Ks500globalmeansizesum1
}

Ks500globaldataset1 <- Ks500globallarvaesize[KK[rep]:attributes(as.matrix(Ks500globallarvaesize))$dim[1],]
Ks500globaldataset2 <- Ks500globaladultsize[KK[rep]:attributes(as.matrix(Ks500globaladultsize))$dim[1],]
Ks500globaldatasize1 <- Ks500globaldataset1[,-(which.min(Ks500globalCIdistance[rep,])+1)] 
Ks500globaldatasize2 <- Ks500globaldataset2[,-(which.min(Ks500globalCIdistance[rep,])+1)] 
##this 1 is the time column
Ks500globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500globaldatasize1[,2:ncol(Ks500globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500globaldatasize2[,2:ncol(Ks500globaldatasize2)])/(N-1))
Ks500globalmeansizesum1<-approx(Ks500globaldatasize1[,1],Ks500globalmeansizesum, x_axis_Ks500global)$y


Ks500globalmeansizetotal[rep,]<-Ks500globalmeansizesum1


Ks500globalaveragesize <- colSums(Ks500globalmeansizetotal)/rep

Ks500globalsizeCIlow<-c()
Ks500globalsizeCIhigh<-c()

for(i in 1:length(Ks500globalaveragesize))
{
  Ks500globalsizeCIlow[i]<-confidence_interval(Ks500globalmeansizetotal[,i], 0.95)[1][[1]]
  Ks500globalsizeCIhigh[i]<-confidence_interval(Ks500globalmeansizetotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=500-time+strength+control")

KT500globallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT500globaladultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT500globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT500globaladultsize)[1])
{
  if (KT500globaladultsize[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT500globaladultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT500globalmeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT500global <- KT500globaladultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT500global <- KT500globaladultsize[KK[rep]:dim(KT500globaladultsize)[1],1]
}

for (i in 1:(rep-1))
{
  KT500globaldataset1 <- KT500globallarvaesize[KK[i]:(KK[i+1]-1),]
  KT500globaldataset2 <- KT500globaladultsize[KK[i]:(KK[i+1]-1),]
  KT500globaldatasize1 <- KT500globaldataset1[,-(which.min(KT500globalCIdistance[i,])+1)] 
  KT500globaldatasize2 <- KT500globaldataset2[,-(which.min(KT500globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT500globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500globaldatasize1[,2:ncol(KT500globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500globaldatasize2[,2:ncol(KT500globaldatasize2)])/(N-1))
  KT500globalmeansizesum1<-approx(KT500globaldatasize1[,1],KT500globalmeansizesum, x_axis_KT500global)$y
  
  KT500globalmeansizetotal[i,]<-KT500globalmeansizesum1
}

KT500globaldataset1 <- KT500globallarvaesize[KK[rep]:attributes(as.matrix(KT500globallarvaesize))$dim[1],]
KT500globaldataset2 <- KT500globaladultsize[KK[rep]:attributes(as.matrix(KT500globaladultsize))$dim[1],]
KT500globaldatasize1 <- KT500globaldataset1[,-(which.min(KT500globalCIdistance[rep,])+1)] 
KT500globaldatasize2 <- KT500globaldataset2[,-(which.min(KT500globalCIdistance[rep,])+1)] 
##this 1 is the time column
KT500globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500globaldatasize1[,2:ncol(KT500globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500globaldatasize2[,2:ncol(KT500globaldatasize2)])/(N-1))
KT500globalmeansizesum1<-approx(KT500globaldatasize1[,1],KT500globalmeansizesum, x_axis_KT500global)$y


KT500globalmeansizetotal[rep,]<-KT500globalmeansizesum1


KT500globalaveragesize <- colSums(KT500globalmeansizetotal)/rep

KT500globalsizeCIlow<-c()
KT500globalsizeCIhigh<-c()

for(i in 1:length(KT500globalaveragesize))
{
  KT500globalsizeCIlow[i]<-confidence_interval(KT500globalmeansizetotal[,i], 0.95)[1][[1]]
  KT500globalsizeCIhigh[i]<-confidence_interval(KT500globalmeansizetotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 local  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-local-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800locallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks800localadultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks800CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()
for ( i in 1:dim(Ks800localadultsize)[1])
{
  if (Ks800localadultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(Ks800localadultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800meansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800local <- Ks800localadultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800local <- Ks800localadultsize[KK[rep]:dim(Ks800localadultsize)[1],1]
}


for (i in 1:(rep-1))
{
  Ks800dataset1 <- Ks800locallarvaesize[KK[i]:(KK[i+1]-1),]
  Ks800dataset2 <- Ks800localadultsize[KK[i]:(KK[i+1]-1),]
  Ks800datasize1 <- Ks800dataset1[,-(which.min(Ks800CIdistance[i,])+1)] 
  Ks800datasize2 <- Ks800dataset2[,-(which.min(Ks800CIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks800meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800datasize1[,2:ncol(Ks800datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800datasize2[,2:ncol(Ks800datasize2)])/(N-1))
  Ks800meansizesum1<-approx(Ks800datasize1[,1],Ks800meansizesum, x_axis_Ks800local)$y
  
  Ks800meansizetotal[i,]<-Ks800meansizesum1
}

Ks800dataset1 <- Ks800locallarvaesize[KK[rep]:attributes(as.matrix(Ks800locallarvaesize))$dim[1],]
Ks800dataset2 <- Ks800localadultsize[KK[rep]:attributes(as.matrix(Ks800localadultsize))$dim[1],]
Ks800datasize1 <- Ks800dataset1[,-(which.min(Ks800CIdistance[rep,])+1)] 
Ks800datasize2 <- Ks800dataset2[,-(which.min(Ks800CIdistance[rep,])+1)] 
##this 1 is the time column
Ks800meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800datasize1[,2:ncol(Ks800datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800datasize2[,2:ncol(Ks800datasize2)])/(N-1))
Ks800meansizesum1<-approx(Ks800datasize1[,1],Ks800meansizesum, x_axis_Ks800local)$y


Ks800meansizetotal[rep,]<-Ks800meansizesum1


Ks800averagesize <- colSums(Ks800meansizetotal)/rep

Ks800sizeCIlow<-c()
Ks800sizeCIhigh<-c()

for(i in 1:length(Ks800averagesize))
{
  Ks800sizeCIlow[i]<-confidence_interval(Ks800meansizetotal[,i], 0.95)[1][[1]]
  Ks800sizeCIhigh[i]<-confidence_interval(Ks800meansizetotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-local-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800locallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT800localadultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT800CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()
for ( i in 1:dim(KT800localadultsize)[1])
{
  if (KT800localadultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT800localadultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800meansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800local <- KT800localadultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800local <- KT800localadultsize[KK[rep]:dim(KT800localadultsize)[1],1]
}


for (i in 1:(rep-1))
{
  KT800dataset1 <- KT800locallarvaesize[KK[i]:(KK[i+1]-1),]
  KT800dataset2 <- KT800localadultsize[KK[i]:(KK[i+1]-1),]
  KT800datasize1 <- KT800dataset1[,-(which.min(KT800CIdistance[i,])+1)] 
  KT800datasize2 <- KT800dataset2[,-(which.min(KT800CIdistance[i,])+1)] 
  ##this 1 is the time column
  KT800meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800datasize1[,2:ncol(KT800datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800datasize2[,2:ncol(KT800datasize2)])/(N-1))
  KT800meansizesum1<-approx(KT800datasize1[,1],KT800meansizesum, x_axis_KT800local)$y
  
  KT800meansizetotal[i,]<-KT800meansizesum1
}


KT800dataset1 <- KT800locallarvaesize[KK[rep]:attributes(as.matrix(KT800locallarvaesize))$dim[1],]
KT800dataset2 <- KT800localadultsize[KK[rep]:attributes(as.matrix(KT800localadultsize))$dim[1],]
KT800datasize1 <- KT800dataset1[,-(which.min(KT800CIdistance[rep,])+1)] 
KT800datasize2 <- KT800dataset2[,-(which.min(KT800CIdistance[rep,])+1)] 
##this 1 is the time column
KT800meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800datasize1[,2:ncol(KT800datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800datasize2[,2:ncol(KT800datasize2)])/(N-1))
KT800meansizesum1<-approx(KT800datasize1[,1],KT800meansizesum, x_axis_KT800local)$y


KT800meansizetotal[rep,]<-KT800meansizesum1


KT800averagesize <- colSums(KT800meansizetotal)/rep

KT800sizeCIlow<-c()
KT800sizeCIhigh<-c()

for(i in 1:length(KT800averagesize))
{
  KT800sizeCIlow[i]<-confidence_interval(KT800meansizetotal[,i], 0.95)[1][[1]]
  KT800sizeCIhigh[i]<-confidence_interval(KT800meansizetotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 neighbor  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800neighborlarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks800neighboradultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks800neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks800neighboradultsize)[1])
{
  if (Ks800neighboradultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(Ks800neighboradultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800neighbormeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800neighbor <- Ks800neighboradultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800neighbor <- Ks800neighboradultsize[KK[rep]:dim(Ks800neighboradultsize)[1],1]
}


for (i in 1:(rep-1))
{
  Ks800neighbordataset1 <- Ks800neighborlarvaesize[KK[i]:(KK[i+1]-1),]
  Ks800neighbordataset2 <- Ks800neighboradultsize[KK[i]:(KK[i+1]-1),]
  Ks800neighbordatasize1 <- Ks800neighbordataset1[,-(which.min(Ks800neighborCIdistance[i,])+1)] 
  Ks800neighbordatasize2 <- Ks800neighbordataset2[,-(which.min(Ks800neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks800neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800neighbordatasize1[,2:ncol(Ks800neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800neighbordatasize2[,2:ncol(Ks800neighbordatasize2)])/(N-1))
  Ks800neighbormeansizesum1<-approx(Ks800neighbordatasize1[,1],Ks800neighbormeansizesum, x_axis_Ks800neighbor)$y
  
  Ks800neighbormeansizetotal[i,]<-Ks800neighbormeansizesum1
}

Ks800neighbordataset1 <- Ks800neighborlarvaesize[KK[rep]:attributes(as.matrix(Ks800neighborlarvaesize))$dim[1],]
Ks800neighbordataset2 <- Ks800neighboradultsize[KK[rep]:attributes(as.matrix(Ks800neighboradultsize))$dim[1],]
Ks800neighbordatasize1 <- Ks800neighbordataset1[,-(which.min(Ks800neighborCIdistance[rep,])+1)] 
Ks800neighbordatasize2 <- Ks800neighbordataset2[,-(which.min(Ks800neighborCIdistance[rep,])+1)] 
##this 1 is the time column
Ks800neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800neighbordatasize1[,2:ncol(Ks800neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800neighbordatasize2[,2:ncol(Ks800neighbordatasize2)])/(N-1))
Ks800neighbormeansizesum1<-approx(Ks800neighbordatasize1[,1],Ks800neighbormeansizesum, x_axis_Ks800neighbor)$y


Ks800neighbormeansizetotal[rep,]<-Ks800neighbormeansizesum1


Ks800neighboraveragesize <- colSums(Ks800neighbormeansizetotal)/rep

Ks800neighborsizeCIlow<-c()
Ks800neighborsizeCIhigh<-c()

for(i in 1:length(Ks800neighboraveragesize))
{
  Ks800neighborsizeCIlow[i]<-confidence_interval(Ks800neighbormeansizetotal[,i], 0.95)[1][[1]]
  Ks800neighborsizeCIhigh[i]<-confidence_interval(Ks800neighbormeansizetotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800neighborlarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT800neighboradultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT800neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT800neighboradultsize)[1])
{
  if (KT800neighboradultsize[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT800neighboradultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800neighbormeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800neighbor <- KT800neighboradultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800neighbor <- KT800neighboradultsize[KK[rep]:dim(KT800neighboradultsize)[1],1]
}


for (i in 1:(rep-1))
{
  KT800neighbordataset1 <- KT800neighborlarvaesize[KK[i]:(KK[i+1]-1),]
  KT800neighbordataset2 <- KT800neighboradultsize[KK[i]:(KK[i+1]-1),]
  KT800neighbordatasize1 <- KT800neighbordataset1[,-(which.min(KT800neighborCIdistance[i,])+1)] 
  KT800neighbordatasize2 <- KT800neighbordataset2[,-(which.min(KT800neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT800neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800neighbordatasize1[,2:ncol(KT800neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800neighbordatasize2[,2:ncol(KT800neighbordatasize2)])/(N-1))
  KT800neighbormeansizesum1<-approx(KT800neighbordatasize1[,1],KT800neighbormeansizesum, x_axis_KT800neighbor)$y
  
  KT800neighbormeansizetotal[i,]<-KT800neighbormeansizesum1
}

KT800neighbordataset1 <- KT800neighborlarvaesize[KK[rep]:attributes(as.matrix(KT800neighborlarvaesize))$dim[1],]
KT800neighbordataset2 <- KT800neighboradultsize[KK[rep]:attributes(as.matrix(KT800neighboradultsize))$dim[1],]
KT800neighbordatasize1 <- KT800neighbordataset1[,-(which.min(KT800neighborCIdistance[rep,])+1)] 
KT800neighbordatasize2 <- KT800neighbordataset2[,-(which.min(KT800neighborCIdistance[rep,])+1)] 
##this 1 is the time column
KT800neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800neighbordatasize1[,2:ncol(KT800neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800neighbordatasize2[,2:ncol(KT800neighbordatasize2)])/(N-1))
KT800neighbormeansizesum1<-approx(KT800neighbordatasize1[,1],KT800neighbormeansizesum, x_axis_KT800neighbor)$y


KT800neighbormeansizetotal[rep,]<-KT800neighbormeansizesum1


KT800neighboraveragesize <- colSums(KT800neighbormeansizetotal)/rep

KT800neighborsizeCIlow<-c()
KT800neighborsizeCIhigh<-c()

for(i in 1:length(KT800neighboraveragesize))
{
  KT800neighborsizeCIlow[i]<-confidence_interval(KT800neighbormeansizetotal[,i], 0.95)[1][[1]]
  KT800neighborsizeCIhigh[i]<-confidence_interval(KT800neighbormeansizetotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 global  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800globallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks800globaladultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks800globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()
for ( i in 1:dim(Ks800globaladultsize)[1])
{
  if (Ks800globaladultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(Ks800globaladultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800globalmeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800global <- Ks800globaladultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800global <- Ks800globaladultsize[KK[rep]:dim(Ks800globaladultsize)[1],1]
}



for (i in 1:(rep-1))
{
  Ks800globaldataset1 <- Ks800globallarvaesize[KK[i]:(KK[i+1]-1),]
  Ks800globaldataset2 <- Ks800globaladultsize[KK[i]:(KK[i+1]-1),]
  Ks800globaldatasize1 <- Ks800globaldataset1[,-(which.min(Ks800globalCIdistance[i,])+1)] 
  Ks800globaldatasize2 <- Ks800globaldataset2[,-(which.min(Ks800globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks800globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800globaldatasize1[,2:ncol(Ks800globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800globaldatasize2[,2:ncol(Ks800globaldatasize2)])/(N-1))
  Ks800globalmeansizesum1<-approx(Ks800globaldatasize1[,1],Ks800globalmeansizesum, x_axis_Ks800global)$y
  
  Ks800globalmeansizetotal[i,]<-Ks800globalmeansizesum1
}

Ks800globaldataset1 <- Ks800globallarvaesize[KK[rep]:attributes(as.matrix(Ks800globallarvaesize))$dim[1],]
Ks800globaldataset2 <- Ks800globaladultsize[KK[rep]:attributes(as.matrix(Ks800globaladultsize))$dim[1],]
Ks800globaldatasize1 <- Ks800globaldataset1[,-(which.min(Ks800globalCIdistance[rep,])+1)] 
Ks800globaldatasize2 <- Ks800globaldataset2[,-(which.min(Ks800globalCIdistance[rep,])+1)] 
##this 1 is the time column
Ks800globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800globaldatasize1[,2:ncol(Ks800globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800globaldatasize2[,2:ncol(Ks800globaldatasize2)])/(N-1))
Ks800globalmeansizesum1<-approx(Ks800globaldatasize1[,1],Ks800globalmeansizesum, x_axis_Ks800global)$y


Ks800globalmeansizetotal[rep,]<-Ks800globalmeansizesum1


Ks800globalaveragesize <- colSums(Ks800globalmeansizetotal)/rep

Ks800globalsizeCIlow<-c()
Ks800globalsizeCIhigh<-c()

for(i in 1:length(Ks800globalaveragesize))
{
  Ks800globalsizeCIlow[i]<-confidence_interval(Ks800globalmeansizetotal[,i], 0.95)[1][[1]]
  Ks800globalsizeCIhigh[i]<-confidence_interval(Ks800globalmeansizetotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800globallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT800globaladultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT800globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT800globaladultsize)[1])
{
  if (KT800globaladultsize[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT800globaladultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800globalmeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800global <- KT800globaladultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800global <- KT800globaladultsize[KK[rep]:dim(KT800globaladultsize)[1],1]
}


for (i in 1:(rep-1))
{
  KT800globaldataset1 <- KT800globallarvaesize[KK[i]:(KK[i+1]-1),]
  KT800globaldataset2 <- KT800globaladultsize[KK[i]:(KK[i+1]-1),]
  KT800globaldatasize1 <- KT800globaldataset1[,-(which.min(KT800globalCIdistance[i,])+1)] 
  KT800globaldatasize2 <- KT800globaldataset2[,-(which.min(KT800globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT800globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800globaldatasize1[,2:ncol(KT800globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800globaldatasize2[,2:ncol(KT800globaldatasize2)])/(N-1))
  KT800globalmeansizesum1<-approx(KT800globaldatasize1[,1],KT800globalmeansizesum, x_axis_KT800global)$y
  
  KT800globalmeansizetotal[i,]<-KT800globalmeansizesum1
}

KT800globaldataset1 <- KT800globallarvaesize[KK[rep]:attributes(as.matrix(KT800globallarvaesize))$dim[1],]
KT800globaldataset2 <- KT800globaladultsize[KK[rep]:attributes(as.matrix(KT800globaladultsize))$dim[1],]
KT800globaldatasize1 <- KT800globaldataset1[,-(which.min(KT800globalCIdistance[rep,])+1)] 
KT800globaldatasize2 <- KT800globaldataset2[,-(which.min(KT800globalCIdistance[rep,])+1)] 
##this 1 is the time column
KT800globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800globaldatasize1[,2:ncol(KT800globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800globaldatasize2[,2:ncol(KT800globaldatasize2)])/(N-1))
KT800globalmeansizesum1<-approx(KT800globaldatasize1[,1],KT800globalmeansizesum, x_axis_KT800global)$y


KT800globalmeansizetotal[rep,]<-KT800globalmeansizesum1


KT800globalaveragesize <- colSums(KT800globalmeansizetotal)/rep

KT800globalsizeCIlow<-c()
KT800globalsizeCIhigh<-c()

for(i in 1:length(KT800globalaveragesize))
{
  KT800globalsizeCIlow[i]<-confidence_interval(KT800globalmeansizetotal[,i], 0.95)[1][[1]]
  KT800globalsizeCIhigh[i]<-confidence_interval(KT800globalmeansizetotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 2000 local  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-local-N-9-p-0.01-t-300-K=2000-strengthonly")

#SF-local-N-9-p-0.01-t-300-epsilon1-100-epsilon2-100-K=2000-strengthmean
Ks2000locallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks2000localadultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks2000CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


#aveKs2000localcontrol <- rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*newK2000locallarvaesize[,2:ncol(newK2000locallarvaesize)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*newK2000localadultsize[,2:ncol(newK2000locallarvaesize)])/(N-1)
#Ks2000CIread<-as.data.frame(Ks2000CIread)
####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks2000localadultsize)[1])
{
  if (Ks2000localadultsize[i,1]==0)
    KK <-c(KK, i)
}

#Ks2000meansizetotal<-matrix(, nrow = rep, ncol = diff(KK)[1])

lastlength <- dim(Ks2000localadultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000meansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000local <- Ks2000localadultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000local <- Ks2000localadultsize[KK[rep]:dim(Ks2000localadultsize)[1],1]
}



for (i in 1:(rep-1))
{
  Ks2000dataset1 <- Ks2000locallarvaesize[KK[i]:(KK[i+1]-1),]
  Ks2000dataset2 <- Ks2000localadultsize[KK[i]:(KK[i+1]-1),]
  Ks2000datasize1 <- Ks2000dataset1[,-(which.min(Ks2000CIdistance[i,])+1)] 
  Ks2000datasize2 <- Ks2000dataset2[,-(which.min(Ks2000CIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks2000meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000datasize1[,2:ncol(Ks2000datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000datasize2[,2:ncol(Ks2000datasize2)])/(N-1))
  Ks2000meansizesum1<-approx(Ks2000datasize1[,1],Ks2000meansizesum, x_axis_Ks2000local)$y
  
  Ks2000meansizetotal[i,]<-Ks2000meansizesum1
}

##last row of Ks2000meanInfectiontotal

Ks2000dataset1 <- Ks2000locallarvaesize[KK[rep]:attributes(as.matrix(Ks2000locallarvaesize))$dim[1],]
Ks2000dataset2 <- Ks2000localadultsize[KK[rep]:attributes(as.matrix(Ks2000localadultsize))$dim[1],]
Ks2000datasize1 <- Ks2000dataset1[,-(which.min(Ks2000CIdistance[rep,])+1)] 
Ks2000datasize2 <- Ks2000dataset2[,-(which.min(Ks2000CIdistance[rep,])+1)] 
##this 1 is the time column
Ks2000meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000datasize1[,2:ncol(Ks2000datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000datasize2[,2:ncol(Ks2000datasize2)])/(N-1))
Ks2000meansizesum1<-approx(Ks2000datasize1[,1],Ks2000meansizesum, x_axis_Ks2000local)$y


Ks2000meansizetotal[rep,]<-Ks2000meansizesum1


Ks2000averagesize <- colSums(Ks2000meansizetotal)/rep

Ks2000sizeCIlow<-c()
Ks2000sizeCIhigh<-c()

for(i in 1:length(Ks2000averagesize))
{
  Ks2000sizeCIlow[i]<-confidence_interval(Ks2000meansizetotal[,i], 0.95)[1][[1]]
  Ks2000sizeCIhigh[i]<-confidence_interval(Ks2000meansizetotal[,i],0.95)[2][[1]]
}


#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-local-N-9-p-0.01-t-300-K=2000-time+strength+control")

KT2000locallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT2000localadultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT2000CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT2000localadultsize)[1])
{
  if (KT2000localadultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT2000localadultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000meansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000local <- KT2000localadultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000local <- KT2000localadultsize[KK[rep]:dim(KT2000localadultsize)[1],1]
}


for (i in 1:(rep-1))
{
  KT2000dataset1 <- KT2000locallarvaesize[KK[i]:(KK[i+1]-1),]
  KT2000dataset2 <- KT2000localadultsize[KK[i]:(KK[i+1]-1),]
  KT2000datasize1 <- KT2000dataset1[,-(which.min(KT2000CIdistance[i,])+1)] 
  KT2000datasize2 <- KT2000dataset2[,-(which.min(KT2000CIdistance[i,])+1)] 
  ##this 1 is the time column
  KT2000meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000datasize1[,2:ncol(KT2000datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000datasize2[,2:ncol(KT2000datasize2)])/(N-1))
  KT2000meansizesum1<-approx(KT2000datasize1[,1],KT2000meansizesum, x_axis_KT2000local)$y
  
  KT2000meansizetotal[i,]<-KT2000meansizesum1
}


KT2000dataset1 <- KT2000locallarvaesize[KK[rep]:attributes(as.matrix(KT2000locallarvaesize))$dim[1],]
KT2000dataset2 <- KT2000localadultsize[KK[rep]:attributes(as.matrix(KT2000localadultsize))$dim[1],]
KT2000datasize1 <- KT2000dataset1[,-(which.min(KT2000CIdistance[rep,])+1)] 
KT2000datasize2 <- KT2000dataset2[,-(which.min(KT2000CIdistance[rep,])+1)] 
##this 1 is the time column
KT2000meansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000datasize1[,2:ncol(KT2000datasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000datasize2[,2:ncol(KT2000datasize2)])/(N-1))
KT2000meansizesum1<-approx(KT2000datasize1[,1],KT2000meansizesum, x_axis_KT2000local)$y


KT2000meansizetotal[rep,]<-KT2000meansizesum1


KT2000averagesize <- colSums(KT2000meansizetotal)/rep

KT2000sizeCIlow<-c()
KT2000sizeCIhigh<-c()

for(i in 1:length(KT2000averagesize))
{
  KT2000sizeCIlow[i]<-confidence_interval(KT2000meansizetotal[,i], 0.95)[1][[1]]
  KT2000sizeCIhigh[i]<-confidence_interval(KT2000meansizetotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 2000 neighbor  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=2000-strengthonly-10-17")

Ks2000neighborlarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks2000neighboradultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks2000neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks2000neighboradultsize)[1])
{
  if (Ks2000neighboradultsize[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks2000neighboradultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000neighbormeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000neighbor <- Ks2000neighboradultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000neighbor <- Ks2000neighboradultsize[KK[rep]:dim(Ks2000neighboradultsize)[1],1]
}

for (i in 1:(rep-1))
{
  Ks2000neighbordataset1 <- Ks2000neighborlarvaesize[KK[i]:(KK[i+1]-1),]
  Ks2000neighbordataset2 <- Ks2000neighboradultsize[KK[i]:(KK[i+1]-1),]
  Ks2000neighbordatasize1 <- Ks2000neighbordataset1[,-(which.min(Ks2000neighborCIdistance[i,])+1)] 
  Ks2000neighbordatasize2 <- Ks2000neighbordataset2[,-(which.min(Ks2000neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks2000neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000neighbordatasize1[,2:ncol(Ks2000neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000neighbordatasize2[,2:ncol(Ks2000neighbordatasize2)])/(N-1))
  Ks2000neighbormeansizesum1<-approx(Ks2000neighbordatasize1[,1],Ks2000neighbormeansizesum, x_axis_Ks2000neighbor)$y
  
  Ks2000neighbormeansizetotal[i,]<-Ks2000neighbormeansizesum1
}

Ks2000neighbordataset1 <- Ks2000neighborlarvaesize[KK[rep]:attributes(as.matrix(Ks2000neighborlarvaesize))$dim[1],]
Ks2000neighbordataset2 <- Ks2000neighboradultsize[KK[rep]:attributes(as.matrix(Ks2000neighboradultsize))$dim[1],]
Ks2000neighbordatasize1 <- Ks2000neighbordataset1[,-(which.min(Ks2000neighborCIdistance[rep,])+1)] 
Ks2000neighbordatasize2 <- Ks2000neighbordataset2[,-(which.min(Ks2000neighborCIdistance[rep,])+1)] 
##this 1 is the time column
Ks2000neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000neighbordatasize1[,2:ncol(Ks2000neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000neighbordatasize2[,2:ncol(Ks2000neighbordatasize2)])/(N-1))
Ks2000neighbormeansizesum1<-approx(Ks2000neighbordatasize1[,1],Ks2000neighbormeansizesum, x_axis_Ks2000neighbor)$y


Ks2000neighbormeansizetotal[rep,]<-Ks2000neighbormeansizesum1


Ks2000neighboraveragesize <- colSums(Ks2000neighbormeansizetotal)/rep

Ks2000neighborsizeCIlow<-c()
Ks2000neighborsizeCIhigh<-c()

for(i in 1:length(Ks2000neighboraveragesize))
{
  Ks2000neighborsizeCIlow[i]<-confidence_interval(Ks2000neighbormeansizetotal[,i], 0.95)[1][[1]]
  Ks2000neighborsizeCIhigh[i]<-confidence_interval(Ks2000neighbormeansizetotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=2000-time+strength+control")

KT2000neighborlarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT2000neighboradultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT2000neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT2000neighboradultsize)[1])
{
  if (KT2000neighboradultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT2000neighboradultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000neighbormeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000neighbor <- KT2000neighboradultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000neighbor <- KT2000neighboradultsize[KK[rep]:dim(KT2000neighboradultsize)[1],1]
}


for (i in 1:(rep-1))
{
  KT2000neighbordataset1 <- KT2000neighborlarvaesize[KK[i]:(KK[i+1]-1),]
  KT2000neighbordataset2 <- KT2000neighboradultsize[KK[i]:(KK[i+1]-1),]
  KT2000neighbordatasize1 <- KT2000neighbordataset1[,-(which.min(KT2000neighborCIdistance[i,])+1)] 
  KT2000neighbordatasize2 <- KT2000neighbordataset2[,-(which.min(KT2000neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT2000neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000neighbordatasize1[,2:ncol(KT2000neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000neighbordatasize2[,2:ncol(KT2000neighbordatasize2)])/(N-1))
  KT2000neighbormeansizesum1<-approx(KT2000neighbordatasize1[,1],KT2000neighbormeansizesum, x_axis_KT2000neighbor)$y
  
  KT2000neighbormeansizetotal[i,]<-KT2000neighbormeansizesum1
}

KT2000neighbordataset1 <- KT2000neighborlarvaesize[KK[rep]:attributes(as.matrix(KT2000neighborlarvaesize))$dim[1],]
KT2000neighbordataset2 <- KT2000neighboradultsize[KK[rep]:attributes(as.matrix(KT2000neighboradultsize))$dim[1],]
KT2000neighbordatasize1 <- KT2000neighbordataset1[,-(which.min(KT2000neighborCIdistance[rep,])+1)] 
KT2000neighbordatasize2 <- KT2000neighbordataset2[,-(which.min(KT2000neighborCIdistance[rep,])+1)] 
##this 1 is the time column
KT2000neighbormeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000neighbordatasize1[,2:ncol(KT2000neighbordatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000neighbordatasize2[,2:ncol(KT2000neighbordatasize2)])/(N-1))
KT2000neighbormeansizesum1<-approx(KT2000neighbordatasize1[,1],KT2000neighbormeansizesum, x_axis_KT2000neighbor)$y


KT2000neighbormeansizetotal[rep,]<-KT2000neighbormeansizesum1


KT2000neighboraveragesize <- colSums(KT2000neighbormeansizetotal)/rep

KT2000neighborsizeCIlow<-c()
KT2000neighborsizeCIhigh<-c()

for(i in 1:length(KT2000neighboraveragesize))
{
  KT2000neighborsizeCIlow[i]<-confidence_interval(KT2000neighbormeansizetotal[,i], 0.95)[1][[1]]
  KT2000neighborsizeCIhigh[i]<-confidence_interval(KT2000neighbormeansizetotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 2000 global  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=2000-strengthonly-10-17")

Ks2000globallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

Ks2000globaladultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

Ks2000globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks2000globaladultsize)[1])
{
  if (Ks2000globaladultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(Ks2000globaladultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000globalmeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000global <- Ks2000globaladultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000global <- Ks2000globaladultsize[KK[rep]:dim(Ks2000globaladultsize)[1],1]
}


for (i in 1:(rep-1))
{
  Ks2000globaldataset1 <- Ks2000globallarvaesize[KK[i]:(KK[i+1]-1),]
  Ks2000globaldataset2 <- Ks2000globaladultsize[KK[i]:(KK[i+1]-1),]
  Ks2000globaldatasize1 <- Ks2000globaldataset1[,-(which.min(Ks2000globalCIdistance[i,])+1)] 
  Ks2000globaldatasize2 <- Ks2000globaldataset2[,-(which.min(Ks2000globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks2000globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000globaldatasize1[,2:ncol(Ks2000globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000globaldatasize2[,2:ncol(Ks2000globaldatasize2)])/(N-1))
  Ks2000globalmeansizesum1<-approx(Ks2000globaldatasize1[,1],Ks2000globalmeansizesum, x_axis_Ks2000global)$y
  
  Ks2000globalmeansizetotal[i,]<-Ks2000globalmeansizesum1
}


Ks2000globaldataset1 <- Ks2000globallarvaesize[KK[rep]:attributes(as.matrix(Ks2000globallarvaesize))$dim[1],]
Ks2000globaldataset2 <- Ks2000globaladultsize[KK[rep]:attributes(as.matrix(Ks2000globaladultsize))$dim[1],]
Ks2000globaldatasize1 <- Ks2000globaldataset1[,-(which.min(Ks2000globalCIdistance[rep,])+1)] 
Ks2000globaldatasize2 <- Ks2000globaldataset2[,-(which.min(Ks2000globalCIdistance[rep,])+1)] 
##this 1 is the time column
Ks2000globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000globaldatasize1[,2:ncol(Ks2000globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000globaldatasize2[,2:ncol(Ks2000globaldatasize2)])/(N-1))
Ks2000globalmeansizesum1<-approx(Ks2000globaldatasize1[,1],Ks2000globalmeansizesum, x_axis_Ks2000global)$y


Ks2000globalmeansizetotal[rep,]<-Ks2000globalmeansizesum1


Ks2000globalaveragesize <- colSums(Ks2000globalmeansizetotal)/rep

Ks2000globalsizeCIlow<-c()
Ks2000globalsizeCIhigh<-c()

for(i in 1:length(Ks2000globalaveragesize))
{
  Ks2000globalsizeCIlow[i]<-confidence_interval(Ks2000globalmeansizetotal[,i], 0.95)[1][[1]]
  Ks2000globalsizeCIhigh[i]<-confidence_interval(Ks2000globalmeansizetotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=2000-time+strength+control")


KT2000globallarvaesize <- read.csv("2-a1--a2--g1--g2--e1--e2--MOSQUITOES-LARVAE.csv",header=FALSE)

KT2000globaladultsize <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-MOSQUITOES.csv", header=FALSE)

KT2000globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()
for ( i in 1:dim(KT2000globaladultsize)[1])
{
  if (KT2000globaladultsize[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT2000globaladultsize)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000globalmeansizetotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000global <- KT2000globaladultsize[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000global <- KT2000globaladultsize[KK[rep]:dim(KT2000globaladultsize)[1],1]
}


for (i in 1:(rep-1))
{
  KT2000globaldataset1 <- KT2000globallarvaesize[KK[i]:(KK[i+1]-1),]
  KT2000globaldataset2 <- KT2000globaladultsize[KK[i]:(KK[i+1]-1),]
  KT2000globaldatasize1 <- KT2000globaldataset1[,-(which.min(KT2000globalCIdistance[i,])+1)] 
  KT2000globaldatasize2 <- KT2000globaldataset2[,-(which.min(KT2000globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT2000globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000globaldatasize1[,2:ncol(KT2000globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000globaldatasize2[,2:ncol(KT2000globaldatasize2)])/(N-1))
  KT2000globalmeansizesum1<-approx(KT2000globaldatasize1[,1],KT2000globalmeansizesum, x_axis_KT2000global)$y
  
  KT2000globalmeansizetotal[i,]<-KT2000globalmeansizesum1
}


KT2000globaldataset1 <- KT2000globallarvaesize[KK[rep]:attributes(as.matrix(KT2000globallarvaesize))$dim[1],]
KT2000globaldataset2 <- KT2000globaladultsize[KK[rep]:attributes(as.matrix(KT2000globaladultsize))$dim[1],]
KT2000globaldatasize1 <- KT2000globaldataset1[,-(which.min(KT2000globalCIdistance[rep,])+1)] 
KT2000globaldatasize2 <- KT2000globaldataset2[,-(which.min(KT2000globalCIdistance[rep,])+1)] 
##this 1 is the time column
KT2000globalmeansizesum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000globaldatasize1[,2:ncol(KT2000globaldatasize1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000globaldatasize2[,2:ncol(KT2000globaldatasize2)])/(N-1))
KT2000globalmeansizesum1<-approx(KT2000globaldatasize1[,1],KT2000globalmeansizesum, x_axis_KT2000global)$y


KT2000globalmeansizetotal[rep,]<-KT2000globalmeansizesum1


KT2000globalaveragesize <- colSums(KT2000globalmeansizetotal)/rep

KT2000globalsizeCIlow<-c()
KT2000globalsizeCIhigh<-c()

for(i in 1:length(KT2000globalaveragesize))
{
  KT2000globalsizeCIlow[i]<-confidence_interval(KT2000globalmeansizetotal[,i], 0.95)[1][[1]]
  KT2000globalsizeCIhigh[i]<-confidence_interval(KT2000globalmeansizetotal[,i],0.95)[2][[1]]
}

###################################################################################################################################

#PLOT the above all figures (K = 500 and K = 2000)

####################################################################################################################################

########################################################################################
##########################################################################################

tiff("Figure_S1.tiff", width=10,height=10, units='in',res=600)

par(mfrow=c(3,3))
par(mar=c(2,3,3,3))

### K = 500 local
plot(x_axis_KT500local,KT500averagesize , type = "l", col = "red", lwd = 2, ylim = c(0,0.1))
points(x_axis_KT500local,KT500sizeCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT500local,KT500sizeCIhigh, type = "l", lty = 2, col = "red")

points(x_axis_Ks500local,Ks500averagesize, type = "l", col = "blue", lwd = 2)

### K = 500 neighbor
plot(x_axis_KT500neighbor,KT500neighboraveragesize, type = "l", col = "red", lwd = 2, ylim = c(0,0.1))
points(x_axis_KT500neighbor,KT500neighborsizeCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT500neighbor,KT500neighborsizeCIhigh, type = "l", lty =  2, col = "red")

points(x_axis_Ks500neighbor,Ks500neighboraveragesize, type = "l", col = "blue", lwd = 2)


### K = 500 global

plot(x_axis_KT500global,KT500globalaveragesize, type = "l", col = "red", lwd = 2, ylim = c(0,0.1))
points(x_axis_KT500global,KT500globalsizeCIlow, type = "l", lty=2, col = "red")
points(x_axis_KT500global,KT500globalsizeCIhigh, type = "l", lty=2, col = "red")

points(x_axis_Ks500global,Ks500globalaveragesize, type = "l", col = "blue", lwd = 2)


### K = 800 local
plot(x_axis_KT800local,KT800averagesize, type = "l", col = "red", lwd = 2, ylim = c(0,0.1))
points(x_axis_KT800local,KT800sizeCIlow, type = "l", lty=2, col = "red")
points(x_axis_KT800local,KT800sizeCIhigh, type = "l", lty=2, col = "red")

points(x_axis_Ks800local,Ks800averagesize, type = "l", col = "blue", lwd = 2)

### K = 800 neighbor
plot(x_axis_KT800neighbor,KT800neighboraveragesize, type = "l", col = "red", lwd = 2, ylim = c(0,0.1))
points(x_axis_KT800neighbor,KT800neighborsizeCIlow, type = "l", lty=2, col = "red")
points(x_axis_KT800neighbor,KT800neighborsizeCIhigh, type = "l", lty=2, col = "red")

points(x_axis_Ks800neighbor,Ks800neighboraveragesize, type = "l", col = "blue", lwd = 2)


### K = 800 global

plot(x_axis_KT800global,KT800globalaveragesize, type = "l", col = "red", lwd = 2, ylim = c(0,0.1))
points(x_axis_KT800global,KT800globalsizeCIlow, type = "l", lty=2, col = "red")
points(x_axis_KT800global,KT800globalsizeCIhigh, type = "l", lty=2, col = "red")

points(x_axis_Ks800global,Ks800globalaveragesize, type = "l", col = "blue", lwd = 2)


### K = 2000 local
plot(x_axis_KT2000local,KT2000averagesize, type = "l", col = "red", lwd = 2, ylim = c(0,0.15))
points(x_axis_KT2000local,KT2000sizeCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT2000local,KT2000sizeCIhigh, type = "l", lty = 2, col = "red")

points(x_axis_Ks2000local,Ks2000averagesize, type = "l", col = "blue", lwd = 2)

### K = 2000 neighbor
plot(x_axis_KT2000neighbor,KT2000neighboraveragesize, type = "l", col = "red", lwd = 2, ylim = c(0,0.15))
points(x_axis_KT2000neighbor,KT2000neighborsizeCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT2000neighbor,KT2000neighborsizeCIhigh, type = "l", lty = 2, col = "red")

points(x_axis_Ks2000neighbor,Ks2000neighboraveragesize, type = "l", col = "blue", lwd = 2)


### K = 2000 global

plot(x_axis_KT2000global,KT2000globalaveragesize, type = "l", col = "red", lwd = 2, ylim = c(0,0.15))
points(x_axis_KT2000global,KT2000globalsizeCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT2000global,KT2000globalsizeCIhigh, type = "l", lty = 2, col = "red")

points(x_axis_Ks2000global,Ks2000globalaveragesize, type = "l", col = "blue", lwd = 2)

dev.off()


