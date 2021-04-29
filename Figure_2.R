## =================================================================================================================================================================================================================

## R code to generate Figure 2 for manuscript titled "How public reaction to disease information across scales and the impacts of vector control methods influence disease prevalence and control efficacy"

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
Ks500CIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks500CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks500CIread)[1])
{
  if (Ks500CIread[i,1]==0)
    KK <-c(KK, i)
}

Ks500meanInfectiontotal<-matrix(, nrow = rep, ncol = diff(KK)[1])

for (i in 1:(rep-1))
{
  Ks500dataset <- Ks500CIread[KK[i]:(KK[i+1]-1),]
  Ks500dataInfection <- Ks500dataset[,-(which.min(Ks500CIdistance[i,])+1)]  ##this 1 is the time column
  Ks500meanInfection1 <- as.vector(rowSums(Ks500dataInfection[,2:N])/(N-1))
  Ks500meanInfectiontotal[i,]<-Ks500meanInfection1
}

##last row of Ks500meanInfectiontotal
Ks500dataset <- Ks500CIread[KK[rep]:attributes(as.matrix(Ks500CIread))$dim[1],]
Ks500dataInfection <- Ks500dataset[,-(which.min(Ks500CIdistance[rep,])+1)]  ##this 1 is the time column
Ks500meanInfection1 <- as.vector(rowSums(Ks500dataInfection[,2:N])/(N-1))
Ks500meanInfectiontotal[rep,]<-Ks500meanInfection1


Ks500averageInfection <- colSums(Ks500meanInfectiontotal)/rep

Ks500InfectionCIlow<-c()
Ks500InfectionCIhigh<-c()

for(i in 1:length(Ks500averageInfection))
{
  Ks500InfectionCIlow[i]<-confidence_interval(Ks500meanInfectiontotal[,i], 0.95)[1][[1]]
  Ks500InfectionCIhigh[i]<-confidence_interval(Ks500meanInfectiontotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-local-N-9-p-0.01-t-300-K=500-time+strength+control")

KT500CIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT500CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT500CIread)[1])
{
  if (KT500CIread[i,1]==0)
    KK <-c(KK, i)
}

KT500meanInfectiontotal<-matrix(, nrow = rep, ncol = diff(KK)[1])

for (i in 1:(rep-1))
{
  KT500dataset <- KT500CIread[KK[i]:(KK[i+1]-1),]
  KT500dataInfection <- KT500dataset[,-(which.min(KT500CIdistance[i,])+1)]  ##this 1 is the time column
  KT500meanInfection1 <- as.vector(rowSums(KT500dataInfection[,2:N])/(N-1))
  KT500meanInfectiontotal[i,]<-KT500meanInfection1
}

##last row of KT500meanInfectiontotal
KT500dataset <- KT500CIread[KK[rep]:attributes(as.matrix(KT500CIread))$dim[1],]
KT500dataInfection <- KT500dataset[,-(which.min(KT500CIdistance[rep,])+1)]  ##this 1 is the time column
KT500meanInfection1 <- as.vector(rowSums(KT500dataInfection[,2:N])/(N-1))
KT500meanInfectiontotal[rep,]<-KT500meanInfection1


KT500averageInfection <- colSums(KT500meanInfectiontotal)/rep

KT500InfectionCIlow<-c()
KT500InfectionCIhigh<-c()

for(i in 1:length(KT500averageInfection))
{
  KT500InfectionCIlow[i]<-confidence_interval(KT500meanInfectiontotal[,i], 0.95)[1][[1]]
  KT500InfectionCIhigh[i]<-confidence_interval(KT500meanInfectiontotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 500 neighbor  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=500-strengthonly")

Ks500neighborCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks500neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks500neighborCIread)[1])
{
  if (Ks500neighborCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks500neighborCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks500neighbormeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks500neighbor <- Ks500neighborCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks500neighbor <- Ks500neighborCIread[KK[rep]:dim(Ks500neighborCIread)[1],1]
}

for (i in 1:(rep-1))
{
  Ks500neighbordataset <- Ks500neighborCIread[KK[i]:(KK[i+1]-1),]
  Ks500neighbordataInfection <- Ks500neighbordataset[,-(which.min(Ks500neighborCIdistance[i,])+1)]  ##this 1 is the time column
  Ks500neighbormeanInfection1 <- as.vector(rowSums(Ks500neighbordataInfection[,2:N])/(N-1))
  Ks500neighbormeanInfection2<-approx(Ks500neighbordataset[,1],Ks500neighbormeanInfection1, x_axis_Ks500neighbor)$y
  
  Ks500neighbormeanInfectiontotal[i,]<-Ks500neighbormeanInfection2
}

##last row of Ks500neighbormeanInfectiontotal
Ks500neighbordataset <- Ks500neighborCIread[KK[rep]:attributes(as.matrix(Ks500neighborCIread))$dim[1],]
Ks500neighbordataInfection <- Ks500neighbordataset[,-(which.min(Ks500neighborCIdistance[rep,])+1)]  ##this 1 is the time column
Ks500neighbormeanInfection1 <- as.vector(rowSums(Ks500neighbordataInfection[,2:N])/(N-1))
Ks500neighbormeanInfection2<-approx(Ks500neighbordataset[,1],Ks500neighbormeanInfection1, x_axis_Ks500neighbor)$y

Ks500neighbormeanInfectiontotal[rep,]<-Ks500neighbormeanInfection2

Ks500neighboraverageInfection <- colSums(Ks500neighbormeanInfectiontotal)/rep

Ks500neighborInfectionCIlow<-c()
Ks500neighborInfectionCIhigh<-c()

for(i in 1:length(Ks500neighboraverageInfection))
{
  Ks500neighborInfectionCIlow[i]<-confidence_interval(Ks500neighbormeanInfectiontotal[,i], 0.95)[1][[1]]
  Ks500neighborInfectionCIhigh[i]<-confidence_interval(Ks500neighbormeanInfectiontotal[,i],0.95)[2][[1]]
}


#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=500-time+strength+control")

KT500neighborCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT500neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT500neighborCIread)[1])
{
  if (KT500neighborCIread[i,1]==0)
    KK <-c(KK, i)
}

KT500neighbormeanInfectiontotal<-matrix(, nrow = rep, ncol = diff(KK)[1])

for (i in 1:(rep-1))
{
  KT500neighbordataset <- KT500neighborCIread[KK[i]:(KK[i+1]-1),]
  KT500neighbordataInfection <- KT500neighbordataset[,-(which.min(KT500neighborCIdistance[i,])+1)]  ##this 1 is the time column
  KT500neighbormeanInfection1 <- as.vector(rowSums(KT500neighbordataInfection[,2:N])/(N-1))
  KT500neighbormeanInfectiontotal[i,]<-KT500neighbormeanInfection1
}

KT500neighbordataset <- KT500neighborCIread[KK[rep]:attributes(as.matrix(KT500neighborCIread))$dim[1],]
KT500neighbordataInfection <- KT500neighbordataset[,-(which.min(KT500neighborCIdistance[rep,])+1)]  ##this 1 is the time column
KT500neighbormeanInfection1 <- as.vector(rowSums(KT500neighbordataInfection[,2:N])/(N-1))
KT500neighbormeanInfectiontotal[rep,]<-KT500neighbormeanInfection1


KT500neighboraverageInfection <- colSums(KT500neighbormeanInfectiontotal)/rep

KT500neighborInfectionCIlow<-c()
KT500neighborInfectionCIhigh<-c()

for(i in 1:length(KT500neighboraverageInfection))
{
  KT500neighborInfectionCIlow[i]<-confidence_interval(KT500neighbormeanInfectiontotal[,i], 0.95)[1][[1]]
  KT500neighborInfectionCIhigh[i]<-confidence_interval(KT500neighbormeanInfectiontotal[,i],0.95)[2][[1]]
}


##################################################################################################################################

# K = 500 global  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=500-strengthonly")

Ks500globalCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks500globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks500globalCIread)[1])
{
  if (Ks500globalCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks500globalCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks500globalmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks500global <- Ks500globalCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks500global <- Ks500globalCIread[KK[rep]:dim(Ks500globalCIread)[1],1]
}


for (i in 1:(rep-1))
{
  Ks500globaldataset <- Ks500globalCIread[KK[i]:(KK[i+1]-1),]
  Ks500globaldataInfection <- Ks500globaldataset[,-(which.min(Ks500globalCIdistance[i,])+1)]  ##this 1 is the time column
  Ks500globalmeanInfection1 <- as.vector(rowSums(Ks500globaldataInfection[,2:N])/(N-1))
  Ks500globalmeanInfection2<-approx(Ks500globaldataset[,1],Ks500globalmeanInfection1, x_axis_Ks500global)$y
  
  Ks500globalmeanInfectiontotal[i,]<-Ks500globalmeanInfection2
}

Ks500globaldataset <- Ks500globalCIread[KK[rep]:attributes(as.matrix(Ks500globalCIread))$dim[1],]
Ks500globaldataInfection <- Ks500globaldataset[,-(which.min(Ks500globalCIdistance[rep,])+1)]  ##this 1 is the time column
Ks500globalmeanInfection1 <- as.vector(rowSums(Ks500globaldataInfection[,2:N])/(N-1))
Ks500globalmeanInfection2<-approx(Ks500globaldataset[,1],Ks500globalmeanInfection1, x_axis_Ks500global)$y

Ks500globalmeanInfectiontotal[rep,]<-Ks500globalmeanInfection2

Ks500globalaverageInfection <- colSums(Ks500globalmeanInfectiontotal)/rep

Ks500globalInfectionCIlow<-c()
Ks500globalInfectionCIhigh<-c()

for(i in 1:length(Ks500globalaverageInfection))
{
  Ks500globalInfectionCIlow[i]<-confidence_interval(Ks500globalmeanInfectiontotal[,i], 0.95)[1][[1]]
  Ks500globalInfectionCIhigh[i]<-confidence_interval(Ks500globalmeanInfectiontotal[,i],0.95)[2][[1]]
}

# #####strength +time
#######################################################
###define the tspan in x axis
setwd("~\\SF-global-N-9-p-0.01-t-300-K=500-time+strength+control")

KT500globalCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT500globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT500globalCIread)[1])
{
  if (KT500globalCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT500globalCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT500globalmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT500global <- KT500globalCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT500global <- KT500globalCIread[KK[rep]:dim(KT500globalCIread)[1],1]
}


for (i in 1:(rep-1))
{
  KT500globaldataset <- KT500globalCIread[KK[i]:(KK[i+1]-1),]
  KT500globaldataInfection <- KT500globaldataset[,-(which.min(KT500globalCIdistance[i,])+1)]  ##this 1 is the time column
  KT500globalmeanInfection1 <- as.vector(rowSums(KT500globaldataInfection[,2:N])/(N-1))
  KT500globalmeanInfection2<-approx(KT500globaldataset[,1],KT500globalmeanInfection1, x_axis_KT500global)$y
  
  KT500globalmeanInfectiontotal[i,]<-KT500globalmeanInfection2
}

KT500globaldataset <- KT500globalCIread[KK[rep]:attributes(as.matrix(KT500globalCIread))$dim[1],]
KT500globaldataInfection <- KT500globaldataset[,-(which.min(KT500globalCIdistance[rep,])+1)]  ##this 1 is the time column
KT500globalmeanInfection1 <- as.vector(rowSums(KT500globaldataInfection[,2:N])/(N-1))
KT500globalmeanInfection2<-approx(KT500globaldataset[,1],KT500globalmeanInfection1, x_axis_KT500global)$y

KT500globalmeanInfectiontotal[rep,]<-KT500globalmeanInfection2

KT500globalaverageInfection <- colSums(KT500globalmeanInfectiontotal)/rep

KT500globalInfectionCIlow<-c()
KT500globalInfectionCIhigh<-c()

for(i in 1:length(KT500globalaverageInfection))
{
  KT500globalInfectionCIlow[i]<-confidence_interval(KT500globalmeanInfectiontotal[,i], 0.95)[1][[1]]
  KT500globalInfectionCIhigh[i]<-confidence_interval(KT500globalmeanInfectiontotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 local  (strength first, added timing later)

##################################################################################################################################

#####only strength
############################################################
###define the tspan in x axis
setwd("~\\SF-local-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800localCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks800localCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks800localCIread)[1])
{
  if (Ks800localCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks800localCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800localmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800local <- Ks800localCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800local <- Ks800localCIread[KK[rep]:dim(Ks800localCIread)[1],1]
}

for (i in 1:(rep-1))
{
  Ks800localdataset <- Ks800localCIread[KK[i]:(KK[i+1]-1),]
  Ks800localdataInfection <- Ks800localdataset[,-(which.min(Ks800localCIdistance[i,])+1)]  ##this 1 is the time column
  Ks800localmeanInfection1 <- as.vector(rowSums(Ks800localdataInfection[,2:N])/(N-1))
  Ks800localmeanInfection2<-approx(Ks800localdataset[,1],Ks800localmeanInfection1, x_axis_Ks800local)$y
  
  Ks800localmeanInfectiontotal[i,]<-Ks800localmeanInfection2
}

Ks800localdataset <- Ks800localCIread[KK[rep]:attributes(as.matrix(Ks800localCIread))$dim[1],]
Ks800localdataInfection <- Ks800localdataset[,-(which.min(Ks800localCIdistance[rep,])+1)]  ##this 1 is the time column
Ks800localmeanInfection1 <- as.vector(rowSums(Ks800localdataInfection[,2:N])/(N-1))
Ks800localmeanInfection2<-approx(Ks800localdataset[,1],Ks800localmeanInfection1, x_axis_Ks800local)$y

Ks800localmeanInfectiontotal[rep,]<-Ks800localmeanInfection2

Ks800localaverageInfection <- colSums(Ks800localmeanInfectiontotal)/rep

Ks800localInfectionCIlow<-c()
Ks800localInfectionCIhigh<-c()

for(i in 1:length(Ks800localaverageInfection))
{
  Ks800localInfectionCIlow[i]<-confidence_interval(Ks800localmeanInfectiontotal[,i], 0.95)[1][[1]]
  Ks800localInfectionCIhigh[i]<-confidence_interval(Ks800localmeanInfectiontotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-local-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800localCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT800localCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT800localCIread)[1])
{
  if (KT800localCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT800localCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800localmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800local <- KT800localCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800local <- KT800localCIread[KK[rep]:dim(KT800localCIread)[1],1]
}


for (i in 1:(rep-1))
{
  KT800localdataset <- KT800localCIread[KK[i]:(KK[i+1]-1),]
  KT800localdataInfection <- KT800localdataset[,-(which.min(KT800localCIdistance[i,])+1)]  ##this 1 is the time column
  KT800localmeanInfection1 <- as.vector(rowSums(KT800localdataInfection[,2:N])/(N-1))
  KT800localmeanInfection2<-approx(KT800localdataset[,1],KT800localmeanInfection1, x_axis_KT800local)$y
  
  KT800localmeanInfectiontotal[i,]<-KT800localmeanInfection2
}

KT800localdataset <- KT800localCIread[KK[rep]:attributes(as.matrix(KT800localCIread))$dim[1],]
KT800localdataInfection <- KT800localdataset[,-(which.min(KT800localCIdistance[rep,])+1)]  ##this 1 is the time column
KT800localmeanInfection1 <- as.vector(rowSums(KT800localdataInfection[,2:N])/(N-1))
KT800localmeanInfection2<-approx(KT800localdataset[,1],KT800localmeanInfection1, x_axis_KT800local)$y

KT800localmeanInfectiontotal[rep,]<-KT800localmeanInfection2

KT800localaverageInfection <- colSums(KT800localmeanInfectiontotal)/rep

KT800localInfectionCIlow<-c()
KT800localInfectionCIhigh<-c()

for(i in 1:length(KT800localaverageInfection))
{
  KT800localInfectionCIlow[i]<-confidence_interval(KT800localmeanInfectiontotal[,i], 0.95)[1][[1]]
  KT800localInfectionCIhigh[i]<-confidence_interval(KT800localmeanInfectiontotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 neighbor  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800neighborCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks800neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks800neighborCIread)[1])
{
  if (Ks800neighborCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks800neighborCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800neighbormeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800neighbor <- Ks800neighborCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800neighbor <- Ks800neighborCIread[KK[rep]:dim(Ks800neighborCIread)[1],1]
}


for (i in 1:(rep-1))
{
  Ks800neighbordataset <- Ks800neighborCIread[KK[i]:(KK[i+1]-1),]
  Ks800neighbordataInfection <- Ks800neighbordataset[,-(which.min(Ks800neighborCIdistance[i,])+1)]  ##this 1 is the time column
  Ks800neighbormeanInfection1 <- as.vector(rowSums(Ks800neighbordataInfection[,2:N])/(N-1))
  Ks800neighbormeanInfection2<-approx(Ks800neighbordataset[,1],Ks800neighbormeanInfection1, x_axis_Ks800neighbor)$y
  
  Ks800neighbormeanInfectiontotal[i,]<-Ks800neighbormeanInfection2
}

Ks800neighbordataset <- Ks800neighborCIread[KK[rep]:attributes(as.matrix(Ks800neighborCIread))$dim[1],]
Ks800neighbordataInfection <- Ks800neighbordataset[,-(which.min(Ks800neighborCIdistance[rep,])+1)]  ##this 1 is the time column
Ks800neighbormeanInfection1 <- as.vector(rowSums(Ks800neighbordataInfection[,2:N])/(N-1))
Ks800neighbormeanInfection2<-approx(Ks800neighbordataset[,1],Ks800neighbormeanInfection1, x_axis_Ks800neighbor)$y

Ks800neighbormeanInfectiontotal[rep,]<-Ks800neighbormeanInfection2

Ks800neighboraverageInfection <- colSums(Ks800neighbormeanInfectiontotal)/rep

Ks800neighborInfectionCIlow<-c()
Ks800neighborInfectionCIhigh<-c()

for(i in 1:length(Ks800neighboraverageInfection))
{
  Ks800neighborInfectionCIlow[i]<-confidence_interval(Ks800neighbormeanInfectiontotal[,i], 0.95)[1][[1]]
  Ks800neighborInfectionCIhigh[i]<-confidence_interval(Ks800neighbormeanInfectiontotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800neighborCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT800neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT800neighborCIread)[1])
{
  if (KT800neighborCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT800neighborCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800neighbormeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800neighbor <- KT800neighborCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800neighbor <- KT800neighborCIread[KK[rep]:dim(KT800neighborCIread)[1],1]
}


for (i in 1:(rep-1))
{
  KT800neighbordataset <- KT800neighborCIread[KK[i]:(KK[i+1]-1),]
  KT800neighbordataInfection <- KT800neighbordataset[,-(which.min(KT800neighborCIdistance[i,])+1)]  ##this 1 is the time column
  KT800neighbormeanInfection1 <- as.vector(rowSums(KT800neighbordataInfection[,2:N])/(N-1))
  KT800neighbormeanInfection2<-approx(KT800neighbordataset[,1],KT800neighbormeanInfection1, x_axis_KT800neighbor)$y
  
  KT800neighbormeanInfectiontotal[i,]<-KT800neighbormeanInfection2
}

KT800neighbordataset <- KT800neighborCIread[KK[rep]:attributes(as.matrix(KT800neighborCIread))$dim[1],]
KT800neighbordataInfection <- KT800neighbordataset[,-(which.min(KT800neighborCIdistance[rep,])+1)]  ##this 1 is the time column
KT800neighbormeanInfection1 <- as.vector(rowSums(KT800neighbordataInfection[,2:N])/(N-1))
KT800neighbormeanInfection2<-approx(KT800neighbordataset[,1],KT800neighbormeanInfection1, x_axis_KT800neighbor)$y

KT800neighbormeanInfectiontotal[rep,]<-KT800neighbormeanInfection2

KT800neighboraverageInfection <- colSums(KT800neighbormeanInfectiontotal)/rep

KT800neighborInfectionCIlow<-c()
KT800neighborInfectionCIhigh<-c()

for(i in 1:length(KT800neighboraverageInfection))
{
  KT800neighborInfectionCIlow[i]<-confidence_interval(KT800neighbormeanInfectiontotal[,i], 0.95)[1][[1]]
  KT800neighborInfectionCIhigh[i]<-confidence_interval(KT800neighbormeanInfectiontotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 global  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800globalCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks800globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks800globalCIread)[1])
{
  if (Ks800globalCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks800globalCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800globalmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800global <- Ks800globalCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800global <- Ks800globalCIread[KK[rep]:dim(Ks800globalCIread)[1],1]
}


for (i in 1:(rep-1))
{
  Ks800globaldataset <- Ks800globalCIread[KK[i]:(KK[i+1]-1),]
  Ks800globaldataInfection <- Ks800globaldataset[,-(which.min(Ks800globalCIdistance[i,])+1)]  ##this 1 is the time column
  Ks800globalmeanInfection1 <- as.vector(rowSums(Ks800globaldataInfection[,2:N])/(N-1))
  Ks800globalmeanInfection2<-approx(Ks800globaldataset[,1],Ks800globalmeanInfection1, x_axis_Ks800global)$y
  
  Ks800globalmeanInfectiontotal[i,]<-Ks800globalmeanInfection2
}

Ks800globaldataset <- Ks800globalCIread[KK[rep]:attributes(as.matrix(Ks800globalCIread))$dim[1],]
Ks800globaldataInfection <- Ks800globaldataset[,-(which.min(Ks800globalCIdistance[rep,])+1)]  ##this 1 is the time column
Ks800globalmeanInfection1 <- as.vector(rowSums(Ks800globaldataInfection[,2:N])/(N-1))
Ks800globalmeanInfection2<-approx(Ks800globaldataset[,1],Ks800globalmeanInfection1, x_axis_Ks800global)$y

Ks800globalmeanInfectiontotal[rep,]<-Ks800globalmeanInfection2

Ks800globalaverageInfection <- colSums(Ks800globalmeanInfectiontotal)/rep

Ks800globalInfectionCIlow<-c()
Ks800globalInfectionCIhigh<-c()

for(i in 1:length(Ks800globalaverageInfection))
{
  Ks800globalInfectionCIlow[i]<-confidence_interval(Ks800globalmeanInfectiontotal[,i], 0.95)[1][[1]]
  Ks800globalInfectionCIhigh[i]<-confidence_interval(Ks800globalmeanInfectiontotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-global-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800globalCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT800globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT800globalCIread)[1])
{
  if (KT800globalCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT800globalCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800globalmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800global <- KT800globalCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800global <- KT800globalCIread[KK[rep]:dim(KT800globalCIread)[1],1]
}


for (i in 1:(rep-1))
{
  KT800globaldataset <- KT800globalCIread[KK[i]:(KK[i+1]-1),]
  KT800globaldataInfection <- KT800globaldataset[,-(which.min(KT800globalCIdistance[i,])+1)]  ##this 1 is the time column
  KT800globalmeanInfection1 <- as.vector(rowSums(KT800globaldataInfection[,2:N])/(N-1))
  KT800globalmeanInfection2<-approx(KT800globaldataset[,1],KT800globalmeanInfection1, x_axis_KT800global)$y
  
  KT800globalmeanInfectiontotal[i,]<-KT800globalmeanInfection2
}

KT800globaldataset <- KT800globalCIread[KK[rep]:attributes(as.matrix(KT800globalCIread))$dim[1],]
KT800globaldataInfection <- KT800globaldataset[,-(which.min(KT800globalCIdistance[rep,])+1)]  ##this 1 is the time column
KT800globalmeanInfection1 <- as.vector(rowSums(KT800globaldataInfection[,2:N])/(N-1))
KT800globalmeanInfection2<-approx(KT800globaldataset[,1],KT800globalmeanInfection1, x_axis_KT800global)$y

KT800globalmeanInfectiontotal[rep,]<-KT800globalmeanInfection2

KT800globalaverageInfection <- colSums(KT800globalmeanInfectiontotal)/rep

KT800globalInfectionCIlow<-c()
KT800globalInfectionCIhigh<-c()

for(i in 1:length(KT800globalaverageInfection))
{
  KT800globalInfectionCIlow[i]<-confidence_interval(KT800globalmeanInfectiontotal[,i], 0.95)[1][[1]]
  KT800globalInfectionCIhigh[i]<-confidence_interval(KT800globalmeanInfectiontotal[,i],0.95)[2][[1]]
}


##################################################################################################################################

# K = 2000 local  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-local-N-9-p-0.01-t-300-K=2000-strengthonly")

Ks2000localCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks2000localCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks2000localCIread)[1])
{
  if (Ks2000localCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks2000localCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000localmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000local <- Ks2000localCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000local <-Ks2000localCIread[KK[rep]:dim(Ks2000localCIread)[1],1]
}


for (i in 1:(rep-1))
{
  Ks2000localdataset <- Ks2000localCIread[KK[i]:(KK[i+1]-1),]
  Ks2000localdataInfection <- Ks2000localdataset[,-(which.min(Ks2000localCIdistance[i,])+1)]  ##this 1 is the time column
  Ks2000localmeanInfection1 <- as.vector(rowSums(Ks2000localdataInfection[,2:N])/(N-1))
  Ks2000localmeanInfection2<-approx(Ks2000localdataset[,1],Ks2000localmeanInfection1, x_axis_Ks2000local)$y
  
  Ks2000localmeanInfectiontotal[i,]<-Ks2000localmeanInfection2
}

Ks2000localdataset <- Ks2000localCIread[KK[rep]:attributes(as.matrix(Ks2000localCIread))$dim[1],]
Ks2000localdataInfection <- Ks2000localdataset[,-(which.min(Ks2000localCIdistance[rep,])+1)]  ##this 1 is the time column
Ks2000localmeanInfection1 <- as.vector(rowSums(Ks2000localdataInfection[,2:N])/(N-1))
Ks2000localmeanInfection2<-approx(Ks2000localdataset[,1],Ks2000localmeanInfection1, x_axis_Ks2000local)$y

Ks2000localmeanInfectiontotal[rep,]<-Ks2000localmeanInfection2

Ks2000localaverageInfection <- colSums(Ks2000localmeanInfectiontotal)/rep

Ks2000localInfectionCIlow<-c()
Ks2000localInfectionCIhigh<-c()

for(i in 1:length(Ks2000localaverageInfection))
{
  Ks2000localInfectionCIlow[i]<-confidence_interval(Ks2000localmeanInfectiontotal[,i], 0.95)[1][[1]]
  Ks2000localInfectionCIhigh[i]<-confidence_interval(Ks2000localmeanInfectiontotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\Zika Project\\SF-local-N-9-p-0.01-t-300-K=2000-time+strength+control")

KT2000localCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT2000localCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT2000localCIread)[1])
{
  if (KT2000localCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT2000localCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000localmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000local <- KT2000localCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000local <-KT2000localCIread[KK[rep]:dim(KT2000localCIread)[1],1]
}


for (i in 1:(rep-1))
{
  KT2000localdataset <- KT2000localCIread[KK[i]:(KK[i+1]-1),]
  KT2000localdataInfection <- KT2000localdataset[,-(which.min(KT2000localCIdistance[i,])+1)]  ##this 1 is the time column
  KT2000localmeanInfection1 <- as.vector(rowSums(KT2000localdataInfection[,2:N])/(N-1))
  KT2000localmeanInfection2<-approx(KT2000localdataset[,1],KT2000localmeanInfection1, x_axis_KT2000local)$y
  
  KT2000localmeanInfectiontotal[i,]<-KT2000localmeanInfection2
}

KT2000localdataset <- KT2000localCIread[KK[rep]:attributes(as.matrix(KT2000localCIread))$dim[1],]
KT2000localdataInfection <- KT2000localdataset[,-(which.min(KT2000localCIdistance[rep,])+1)]  ##this 1 is the time column
KT2000localmeanInfection1 <- as.vector(rowSums(KT2000localdataInfection[,2:N])/(N-1))
KT2000localmeanInfection2<-approx(KT2000localdataset[,1],KT2000localmeanInfection1, x_axis_KT2000local)$y

KT2000localmeanInfectiontotal[rep,]<-KT2000localmeanInfection2

KT2000localaverageInfection <- colSums(KT2000localmeanInfectiontotal)/rep

KT2000localInfectionCIlow<-c()
KT2000localInfectionCIhigh<-c()

for(i in 1:length(KT2000localaverageInfection))
{
  KT2000localInfectionCIlow[i]<-confidence_interval(KT2000localmeanInfectiontotal[,i], 0.95)[1][[1]]
  KT2000localInfectionCIhigh[i]<-confidence_interval(KT2000localmeanInfectiontotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 2000 neighbor  (strength first, added timing later)

##################################################################################################################################

setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=2000-strengthonly")

Ks2000neighborCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks2000neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks2000neighborCIread)[1])
{
  if (Ks2000neighborCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks2000neighborCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000neighbormeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000neighbor <- Ks2000neighborCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000neighbor <-Ks2000neighborCIread[KK[rep]:dim(Ks2000neighborCIread)[1],1]
}


for (i in 1:(rep-1))
{
  Ks2000neighbordataset <- Ks2000neighborCIread[KK[i]:(KK[i+1]-1),]
  Ks2000neighbordataInfection <- Ks2000neighbordataset[,-(which.min(Ks2000neighborCIdistance[i,])+1)]  ##this 1 is the time column
  Ks2000neighbormeanInfection1 <- as.vector(rowSums(Ks2000neighbordataInfection[,2:N])/(N-1))
  Ks2000neighbormeanInfection2<-approx(Ks2000neighbordataset[,1],Ks2000neighbormeanInfection1, x_axis_Ks2000neighbor)$y
  
  Ks2000neighbormeanInfectiontotal[i,]<-Ks2000neighbormeanInfection2
}

Ks2000neighbordataset <- Ks2000neighborCIread[KK[rep]:attributes(as.matrix(Ks2000neighborCIread))$dim[1],]
Ks2000neighbordataInfection <- Ks2000neighbordataset[,-(which.min(Ks2000neighborCIdistance[rep,])+1)]  ##this 1 is the time column
Ks2000neighbormeanInfection1 <- as.vector(rowSums(Ks2000neighbordataInfection[,2:N])/(N-1))
Ks2000neighbormeanInfection2<-approx(Ks2000neighbordataset[,1],Ks2000neighbormeanInfection1, x_axis_Ks2000neighbor)$y

Ks2000neighbormeanInfectiontotal[rep,]<-Ks2000neighbormeanInfection2

Ks2000neighboraverageInfection <- colSums(Ks2000neighbormeanInfectiontotal)/rep

Ks2000neighborInfectionCIlow<-c()
Ks2000neighborInfectionCIhigh<-c()

for(i in 1:length(Ks2000neighboraverageInfection))
{
  Ks2000neighborInfectionCIlow[i]<-confidence_interval(Ks2000neighbormeanInfectiontotal[,i], 0.95)[1][[1]]
  Ks2000neighborInfectionCIhigh[i]<-confidence_interval(Ks2000neighbormeanInfectiontotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=2000-time+strength+control")

KT2000neighborCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT2000neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT2000neighborCIread)[1])
{
  if (KT2000neighborCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT2000neighborCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000neighbormeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000neighbor <- KT2000neighborCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000neighbor <-KT2000neighborCIread[KK[rep]:dim(KT2000neighborCIread)[1],1]
}


for (i in 1:(rep-1))
{
  KT2000neighbordataset <- KT2000neighborCIread[KK[i]:(KK[i+1]-1),]
  KT2000neighbordataInfection <- KT2000neighbordataset[,-(which.min(KT2000neighborCIdistance[i,])+1)]  ##this 1 is the time column
  KT2000neighbormeanInfection1 <- as.vector(rowSums(KT2000neighbordataInfection[,2:N])/(N-1))
  KT2000neighbormeanInfection2<-approx(KT2000neighbordataset[,1],KT2000neighbormeanInfection1, x_axis_KT2000neighbor)$y
  
  KT2000neighbormeanInfectiontotal[i,]<-KT2000neighbormeanInfection2
}

KT2000neighbordataset <- KT2000neighborCIread[KK[rep]:attributes(as.matrix(KT2000neighborCIread))$dim[1],]
KT2000neighbordataInfection <- KT2000neighbordataset[,-(which.min(KT2000neighborCIdistance[rep,])+1)]  ##this 1 is the time column
KT2000neighbormeanInfection1 <- as.vector(rowSums(KT2000neighbordataInfection[,2:N])/(N-1))
KT2000neighbormeanInfection2<-approx(KT2000neighbordataset[,1],KT2000neighbormeanInfection1, x_axis_KT2000neighbor)$y

KT2000neighbormeanInfectiontotal[rep,]<-KT2000neighbormeanInfection2

KT2000neighboraverageInfection <- colSums(KT2000neighbormeanInfectiontotal)/rep

KT2000neighborInfectionCIlow<-c()
KT2000neighborInfectionCIhigh<-c()

for(i in 1:length(KT2000neighboraverageInfection))
{
  KT2000neighborInfectionCIlow[i]<-confidence_interval(KT2000neighbormeanInfectiontotal[,i], 0.95)[1][[1]]
  KT2000neighborInfectionCIhigh[i]<-confidence_interval(KT2000neighbormeanInfectiontotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 2000 global  (strength first, added timing later)

##################################################################################################################################

setwd("~\\SF-global-N-9-p-0.01-t-300-K=2000-strengthonly")

Ks2000globalCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks2000globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(Ks2000globalCIread)[1])
{
  if (Ks2000globalCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks2000globalCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000globalmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000global <- Ks2000globalCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000global <-Ks2000globalCIread[KK[rep]:dim(Ks2000globalCIread)[1],1]
}


for (i in 1:(rep-1))
{
  Ks2000globaldataset <- Ks2000globalCIread[KK[i]:(KK[i+1]-1),]
  Ks2000globaldataInfection <- Ks2000globaldataset[,-(which.min(Ks2000globalCIdistance[i,])+1)]  ##this 1 is the time column
  Ks2000globalmeanInfection1 <- as.vector(rowSums(Ks2000globaldataInfection[,2:N])/(N-1))
  Ks2000globalmeanInfection2<-approx(Ks2000globaldataset[,1],Ks2000globalmeanInfection1, x_axis_Ks2000global)$y
  
  Ks2000globalmeanInfectiontotal[i,]<-Ks2000globalmeanInfection2
}

Ks2000globaldataset <- Ks2000globalCIread[KK[rep]:attributes(as.matrix(Ks2000globalCIread))$dim[1],]
Ks2000globaldataInfection <- Ks2000globaldataset[,-(which.min(Ks2000globalCIdistance[rep,])+1)]  ##this 1 is the time column
Ks2000globalmeanInfection1 <- as.vector(rowSums(Ks2000globaldataInfection[,2:N])/(N-1))
Ks2000globalmeanInfection2<-approx(Ks2000globaldataset[,1],Ks2000globalmeanInfection1, x_axis_Ks2000global)$y

Ks2000globalmeanInfectiontotal[rep,]<-Ks2000globalmeanInfection2

Ks2000globalaverageInfection <- colSums(Ks2000globalmeanInfectiontotal)/rep

Ks2000globalInfectionCIlow <- c()
Ks2000globalInfectionCIhigh <- c()

for(i in 1:length(Ks2000globalaverageInfection))
{
  Ks2000globalInfectionCIlow[i]<-confidence_interval(Ks2000globalmeanInfectiontotal[,i], 0.95)[1][[1]]
  Ks2000globalInfectionCIhigh[i]<-confidence_interval(Ks2000globalmeanInfectiontotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-global-N-9-p-0.01-t-300-K=2000-time+strength+control")

KT2000globalCIread<-read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
KT2000globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

KK<-c()
for ( i in 1:dim(KT2000globalCIread)[1])
{
  if (KT2000globalCIread[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT2000globalCIread)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000globalmeanInfectiontotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000global <- KT2000globalCIread[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000global <-KT2000globalCIread[KK[rep]:dim(KT2000globalCIread)[1],1]
}


for (i in 1:(rep-1))
{
  KT2000globaldataset <- KT2000globalCIread[KK[i]:(KK[i+1]-1),]
  KT2000globaldataInfection <- KT2000globaldataset[,-(which.min(KT2000globalCIdistance[i,])+1)]  ##this 1 is the time column
  KT2000globalmeanInfection1 <- as.vector(rowSums(KT2000globaldataInfection[,2:N])/(N-1))
  KT2000globalmeanInfection2<-approx(KT2000globaldataset[,1],KT2000globalmeanInfection1, x_axis_KT2000global)$y
  
  KT2000globalmeanInfectiontotal[i,]<-KT2000globalmeanInfection2
}

KT2000globaldataset <- KT2000globalCIread[KK[rep]:attributes(as.matrix(KT2000globalCIread))$dim[1],]
KT2000globaldataInfection <- KT2000globaldataset[,-(which.min(KT2000globalCIdistance[rep,])+1)]  ##this 1 is the time column
KT2000globalmeanInfection1 <- as.vector(rowSums(KT2000globaldataInfection[,2:N])/(N-1))
KT2000globalmeanInfection2<-approx(KT2000globaldataset[,1],KT2000globalmeanInfection1, x_axis_KT2000global)$y

KT2000globalmeanInfectiontotal[rep,]<-KT2000globalmeanInfection2

KT2000globalaverageInfection <- colSums(KT2000globalmeanInfectiontotal)/rep

KT2000globalInfectionCIlow <- c()
KT2000globalInfectionCIhigh <- c()

for(i in 1:length(KT2000globalaverageInfection))
{
  KT2000globalInfectionCIlow[i]<-confidence_interval(KT2000globalmeanInfectiontotal[,i], 0.95)[1][[1]]
  KT2000globalInfectionCIhigh[i]<-confidence_interval(KT2000globalmeanInfectiontotal[,i],0.95)[2][[1]]
}

###################################################################################################################################

#### No Control

####################################################################################################################################


#####################
### K = 500

setwd("~\\SF-local-N-9-p-0.01-t-300-K=500-nocontrol")

Ks500nocontrollocalInfection <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks500nocontrollocaldistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)
Ks500nocontrollocaldistance_new <- apply(Ks500nocontrollocaldistance,2,function(x){rep(x,each = 601)})

pp <- apply(Ks500nocontrollocaldistance_new,1,function(x){which.min(x)+1})

c<-data.frame(matrix(NA,nrow(Ks500nocontrollocalInfection),(ncol(Ks500nocontrollocalInfection)-1)))

for(i in 1:nrow(Ks500nocontrollocalInfection))
    {
    c[i,] = Ks500nocontrollocalInfection[i,-pp[i]]  
      
    }

c1 <- c[,-1]

newKs500nocontrolglobalmeaninfection <- apply(c1,1,mean)
newKs500nocontrollocalmeaninfection <- newKs500nocontrolglobalmeaninfection
newKs500nocontrolneighbormeaninfection <- newKs500nocontrolglobalmeaninfection

K500globalnocontrol <- c()

for (i in 1:601)
{
  q = 0
  
  for(j in 0:(rep-1))
  {
    
    q <- sum(q, newKs500nocontrolglobalmeaninfection[i+601*j])
  } 
  
  K500globalnocontrol[i]<-q/rep

}

K500localnocontrol <- K500globalnocontrol
K500neighbornocontrol <- K500globalnocontrol

#################################################################
## K = 800

setwd("~\\SF-local-N-9-p-0.01-t-300-K=800-nocontrol-9-23")

Ks800nocontrollocalInfection <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks800nocontrollocaldistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)
Ks800nocontrollocaldistance_new <- apply(Ks800nocontrollocaldistance,2,function(x){rep(x,each = 601)})

pp <- apply(Ks800nocontrollocaldistance_new,1,function(x){which.min(x)+1})

c<-data.frame(matrix(NA,nrow(Ks800nocontrollocalInfection),(ncol(Ks800nocontrollocalInfection)-1)))

for(i in 1:nrow(Ks800nocontrollocalInfection))
{
  c[i,] = Ks800nocontrollocalInfection[i,-pp[i]]  
  
}


c1 <- c[,-1]

newKs800nocontrolglobalmeaninfection <- apply(c1,1,mean)

newKs800nocontrollocalmeaninfection <- newKs800nocontrolglobalmeaninfection
newKs800nocontrolneighbormeaninfection <- newKs800nocontrolglobalmeaninfection

K800globalnocontrol <- c()

for (i in 1:601)
{
  q = 0
  
  for(j in 0:(rep-1))
  {
    
    q <- sum(q, newKs800nocontrolglobalmeaninfection[i+601*j])
  } 
  
  K800globalnocontrol[i]<-q/rep
  
}

K800localnocontrol <- K800globalnocontrol
K800neighbornocontrol <- K800globalnocontrol

#####################################################################
### K = 2000

setwd("~\\SF-local-N-9-p-0.01-t-300-K=2000-nocontrol")

Ks2000nocontrollocalInfection <- read.csv("2-a1--a2--g1--g2--e1--e2--INFECTED-HUMANS.csv",header=FALSE)
Ks2000nocontrollocaldistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)
Ks2000nocontrollocaldistance_new <- apply(Ks2000nocontrollocaldistance,2,function(x){rep(x,each = 601)})

pp <- apply(Ks2000nocontrollocaldistance_new,1,function(x){which.min(x)+1})

c<-data.frame(matrix(NA,nrow(Ks2000nocontrollocalInfection),(ncol(Ks2000nocontrollocalInfection)-1)))

for(i in 1:nrow(Ks2000nocontrollocalInfection))
{
  c[i,] = Ks2000nocontrollocalInfection[i,-pp[i]]  
  
}

c1 <- c[,-1]

newKs2000nocontrolglobalmeaninfection <- apply(c1,1,mean)
newKs2000nocontrollocalmeaninfection <- newKs2000nocontrolglobalmeaninfection
newKs2000nocontrolneighbormeaninfection <- newKs2000nocontrolglobalmeaninfection

K2000globalnocontrol <- c()

for (i in 1:601)
{
  q = 0
  
  for(j in 0:(rep-1))
  {
    
    q <- sum(q, newKs2000nocontrolglobalmeaninfection[i+601*j])
  } 
  
  K2000globalnocontrol[i]<-q/rep
  
}

K2000localnocontrol <- K2000globalnocontrol
K2000neighbornocontrol <- K2000globalnocontrol


########################################################################################
##########################################################################################

tiff("Figure_2.tiff", width=10,height=10, units='in',res=600)

par(mfrow=c(3,3))
par(mar=c(2,3,3,3))

### K = 500 local
plot(KT500dataset[,1],K500localnocontrol, type = "l", col = "black")

points(KT500dataset[,1], KT500averageInfection, type = "l", col = "red", lwd = 2)
points(KT500dataset[,1], Ks500averageInfection, type = "l", col = "blue", lwd = 2)

### K = 500 neighbor

plot(KT500neighbordataset[,1],K500neighbornocontrol, type = "l", col = "black")

points(KT500neighbordataset[,1], KT500neighboraverageInfection, type = "l", col = "red", lwd = 2) 
points(KT500neighbordataset[,1], KT500neighborInfectionCIlow, type = "l", lty = 2, col = "red")
points(KT500neighbordataset[,1], KT500neighborInfectionCIhigh, type = "l", lty = 2, col = "red")
points(KT500neighbordataset[,1], Ks500neighboraverageInfection, type = "l", col = "blue", lwd = 2)

### K = 500 global


plot(x_axis_KT500global,K500globalnocontrol, type = "l", col = "black")

points(x_axis_KT500global, KT500globalaverageInfection, type = "l", col = "red", lwd = 2)
points(x_axis_KT500global, KT500globalInfectionCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT500global, KT500globalInfectionCIhigh, type = "l", lty = 2, col = "red")
points(x_axis_KT500global, Ks500globalaverageInfection, type = "l", col = "blue", lwd = 2)


### K = 800 local
plot(x_axis_KT800local,K800localnocontrol, type = "l", col = "black")

points(x_axis_KT800local, KT800localaverageInfection, type = "l", col = "red", lwd = 2) 
points(x_axis_KT800local, Ks800localaverageInfection, type = "l", col = "blue", lwd = 2)


### K = 800 neighbor
plot(x_axis_KT800neighbor,K800neighbornocontrol, type = "l", col = "black")

points(x_axis_KT800neighbor, KT800neighboraverageInfection, type = "l", col="red", lwd = 2) 
points(x_axis_KT800neighbor, KT800neighborInfectionCIlow, type = "l", lty = 2, col="red")
points(x_axis_KT800neighbor, KT800neighborInfectionCIhigh, type = "l", lty = 2, col="red")
points(x_axis_KT800neighbor, Ks800neighboraverageInfection, type = "l", col = "blue", lwd = 2)

### K = 800 global

plot(x_axis_KT800global,K800globalnocontrol, type = "l", col = "black")

points(x_axis_KT800global, KT800globalaverageInfection, type = "l",col="red", lwd = 2)
points(x_axis_KT800global, KT800globalInfectionCIlow, type = "l", lty = 2,col="red")
points(x_axis_KT800global, KT800globalInfectionCIhigh, type = "l", lty = 2, col="red")
points(x_axis_KT800global, Ks800globalaverageInfection, type = "l", col = "blue", lwd = 2)

### K = 2000 local
plot(x_axis_KT2000local,K2000globalnocontrol, type = "l", col = "black", ylim = c(0, 0.01))

points(x_axis_KT2000local, KT2000localaverageInfection, type = "l", col = "red", lwd = 2) 
points(x_axis_KT2000local, Ks2000localaverageInfection, type = "l", col = "blue", lwd = 2)

### K = 2000 neighbor
plot(x_axis_KT2000neighbor,K2000neighbornocontrol, type = "l", col = "black", ylim = c(0, 0.01))

points(x_axis_KT2000neighbor, KT2000neighboraverageInfection, type = "l", col="red", lwd = 2) 
points(x_axis_KT2000neighbor, KT2000neighborInfectionCIlow, type = "l", lty = 2, col="red")
points(x_axis_KT2000neighbor, KT2000neighborInfectionCIhigh, type = "l", lty = 2, col="red")
points(x_axis_KT2000neighbor, Ks2000neighboraverageInfection, type = "l", col = "blue", lwd = 2)

### K = 2000 global

plot(x_axis_KT2000global,K2000globalnocontrol, type = "l", col = "black", ylim = c(0, 0.01) )

points(x_axis_KT2000global, KT2000globalaverageInfection, type = "l",col="red", lwd = 2)
points(x_axis_KT2000global, KT2000globalInfectionCIlow, type = "l", lty = 2,col="red")
points(x_axis_KT2000global, KT2000globalInfectionCIhigh, type = "l", lty = 2, col="red")
points(x_axis_KT2000global, Ks2000globalaverageInfection, type = "l", col = "blue", lwd = 2)

dev.off()



