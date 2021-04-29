## ==========================================================================================================================================================================================================================
##
## R code to generate Figure 1 and 3 for manuscript titled "How public reaction to disease information across scales and the impacts of vector control methods influence disease prevalence and control efficacy"

## ==================================================================================================================================================================================================================


library(ggplot2)
library(reshape)
library(here)

## set the current directory as default to save all the following results
set_here()

##all csv. files are from the matlab code named "Figure 1.m"

N <- 9
rep <-50

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

#################################################################################################################################
#################################################################################################################################

## Infection Section (used for calculation of control efficacy in Figure 3)

#################################################################################################################################
##################################################################################################################################

# K = 500 local  (strength first with Ks notation, added timing later with KT notation) 

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


##################################################################################################################################

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

##################################################################################################################################
##################################################################################################################################

## Control effort: this section is to calculate the control effort at diff breeding capacities and information 

##################################################################################################################################
##################################################################################################################################

# K = 500 local  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-local-N-9-p-0.01-t-300-K=500-strengthonly")

Ks500localcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks500localcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks500CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)



KK<-c()
for ( i in 1:dim(Ks500localcontroladult)[1])
{
  if (Ks500localcontroladult[i,1]==0)
    KK <-c(KK, i)
}

#Ks500meanControltotal<-matrix(, nrow = rep, ncol = diff(KK)[1])

lastlength <- dim(Ks500localcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks500meanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks500local <- Ks500localcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks500local <- Ks500localcontroladult[KK[rep]:dim(Ks500localcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks500dataset1 <- Ks500localcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks500dataset2 <- Ks500localcontroladult[KK[i]:(KK[i+1]-1),]
  Ks500dataControl1 <- Ks500dataset1[,-(which.min(Ks500CIdistance[i,])+1)] 
  Ks500dataControl2 <- Ks500dataset2[,-(which.min(Ks500CIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks500meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500dataControl1[,2:ncol(Ks500dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500dataControl2[,2:ncol(Ks500dataControl2)])/(N-1))
  Ks500meanControlsum1<-approx(Ks500dataControl1[,1],Ks500meanControlsum, x_axis_Ks500local)$y
  
  Ks500meanControltotal[i,]<-Ks500meanControlsum1
}

##last row of Ks500meanInfectiontotal

Ks500dataset1 <- Ks500localcontrollarvae[KK[rep]:attributes(as.matrix(Ks500localcontrollarvae))$dim[1],]
Ks500dataset2 <- Ks500localcontroladult[KK[rep]:attributes(as.matrix(Ks500localcontroladult))$dim[1],]
Ks500dataControl1 <- Ks500dataset1[,-(which.min(Ks500CIdistance[rep,])+1)] 
Ks500dataControl2 <- Ks500dataset2[,-(which.min(Ks500CIdistance[rep,])+1)] 
##this 1 is the time column
Ks500meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500dataControl1[,2:ncol(Ks500dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500dataControl2[,2:ncol(Ks500dataControl2)])/(N-1))
Ks500meanControlsum1<-approx(Ks500dataControl1[,1],Ks500meanControlsum, x_axis_Ks500local)$y


Ks500meanControltotal[rep,]<-Ks500meanControlsum1


Ks500averageControl <- colSums(Ks500meanControltotal)/rep

Ks500ControlCIlow<-c()
Ks500ControlCIhigh<-c()

for(i in 1:length(Ks500averageControl))
{
  Ks500ControlCIlow[i]<-confidence_interval(Ks500meanControltotal[,i], 0.95)[1][[1]]
  Ks500ControlCIhigh[i]<-confidence_interval(Ks500meanControltotal[,i],0.95)[2][[1]]
}


#####strength +time
#######################################################

###define the tspan in x axis
setwd("~//SF-local-N-9-p-0.01-t-300-K=500-strengthonly")


KT500localcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT500localcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT500CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()

for ( i in 1:dim(KT500localcontroladult)[1])
{
  if (KT500localcontroladult[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT500localcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT500meanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT500local <- KT500localcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT500local <- KT500localcontroladult[KK[rep]:dim(KT500localcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  KT500dataset1 <- KT500localcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT500dataset2 <- KT500localcontroladult[KK[i]:(KK[i+1]-1),]
  KT500dataControl1 <- KT500dataset1[,-(which.min(KT500CIdistance[i,])+1)] 
  KT500dataControl2 <- KT500dataset2[,-(which.min(KT500CIdistance[i,])+1)] 
  ##this 1 is the time column
  KT500meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500dataControl1[,2:ncol(KT500dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500dataControl2[,2:ncol(KT500dataControl2)])/(N-1))
  KT500meanControlsum1<-approx(KT500dataControl1[,1],KT500meanControlsum, x_axis_KT500local)$y
  
  KT500meanControltotal[i,]<-KT500meanControlsum1
}

##last row of KT500meanInfectiontotal

KT500dataset1 <- KT500localcontrollarvae[KK[rep]:attributes(as.matrix(KT500localcontrollarvae))$dim[1],]
KT500dataset2 <- KT500localcontroladult[KK[rep]:attributes(as.matrix(KT500localcontroladult))$dim[1],]
KT500dataControl1 <- KT500dataset1[,-(which.min(KT500CIdistance[rep,])+1)] 
KT500dataControl2 <- KT500dataset2[,-(which.min(KT500CIdistance[rep,])+1)] 
##this 1 is the time column
KT500meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500dataControl1[,2:ncol(KT500dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500dataControl2[,2:ncol(KT500dataControl2)])/(N-1))
KT500meanControlsum1<-approx(KT500dataControl1[,1],KT500meanControlsum, x_axis_KT500local)$y


KT500meanControltotal[rep,]<-KT500meanControlsum1


KT500averageControl <- colSums(KT500meanControltotal)/rep

KT500ControlCIlow<-c()
KT500ControlCIhigh<-c()

for(i in 1:length(KT500averageControl))
{
  KT500ControlCIlow[i]<-confidence_interval(KT500meanControltotal[,i], 0.95)[1][[1]]
  KT500ControlCIhigh[i]<-confidence_interval(KT500meanControltotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 500 neighbor  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################

setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=500-strengthonly")

Ks500neighborcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks500neighborcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks500neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()
for ( i in 1:dim(Ks500neighborcontroladult)[1])
{
  if (Ks500neighborcontroladult[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks500neighborcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks500neighbormeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks500neighbor <- Ks500neighborcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks500neighbor <- Ks500neighborcontroladult[KK[rep]:dim(Ks500neighborcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks500neighbordataset1 <- Ks500neighborcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks500neighbordataset2 <- Ks500neighborcontroladult[KK[i]:(KK[i+1]-1),]
  Ks500neighbordataControl1 <- Ks500neighbordataset1[,-(which.min(Ks500neighborCIdistance[i,])+1)] 
  Ks500neighbordataControl2 <- Ks500neighbordataset2[,-(which.min(Ks500neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks500neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500neighbordataControl1[,2:ncol(Ks500neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500neighbordataControl2[,2:ncol(Ks500neighbordataControl2)])/(N-1))
  Ks500neighbormeanControlsum1<-approx(Ks500neighbordataControl1[,1],Ks500neighbormeanControlsum, x_axis_Ks500neighbor)$y
  
  Ks500neighbormeanControltotal[i,]<-Ks500neighbormeanControlsum1
}

##last row of Ks500meanInfectiontotal

Ks500neighbordataset1 <- Ks500neighborcontrollarvae[KK[rep]:attributes(as.matrix(Ks500neighborcontrollarvae))$dim[1],]
Ks500neighbordataset2 <- Ks500neighborcontroladult[KK[rep]:attributes(as.matrix(Ks500neighborcontroladult))$dim[1],]
Ks500neighbordataControl1 <- Ks500neighbordataset1[,-(which.min(Ks500neighborCIdistance[rep,])+1)] 
Ks500neighbordataControl2 <- Ks500neighbordataset2[,-(which.min(Ks500neighborCIdistance[rep,])+1)] 
##this 1 is the time column
Ks500neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500neighbordataControl1[,2:ncol(Ks500neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500neighbordataControl2[,2:ncol(Ks500neighbordataControl2)])/(N-1))
Ks500neighbormeanControlsum1<-approx(Ks500neighbordataControl1[,1],Ks500neighbormeanControlsum, x_axis_Ks500neighbor)$y


Ks500neighbormeanControltotal[rep,]<-Ks500neighbormeanControlsum1


Ks500neighboraverageControl <- colSums(Ks500neighbormeanControltotal)/rep

Ks500neighborControlCIlow<-c()
Ks500neighborControlCIhigh<-c()

for(i in 1:length(Ks500neighboraverageControl))
{
  Ks500neighborControlCIlow[i]<-confidence_interval(Ks500neighbormeanControltotal[,i], 0.95)[1][[1]]
  Ks500neighborControlCIhigh[i]<-confidence_interval(Ks500neighbormeanControltotal[,i],0.95)[2][[1]]
}

#####only strength
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=500-time+strength+control")

KT500neighborcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT500neighborcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT500neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


KK<-c()
for ( i in 1:dim(KT500neighborcontroladult)[1])
{
  if (KT500neighborcontroladult[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT500neighborcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT500neighbormeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT500neighbor <- KT500neighborcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT500neighbor <- KT500neighborcontroladult[KK[rep]:dim(KT500neighborcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  KT500neighbordataset1 <- KT500neighborcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT500neighbordataset2 <- KT500neighborcontroladult[KK[i]:(KK[i+1]-1),]
  KT500neighbordataControl1 <- KT500neighbordataset1[,-(which.min(KT500neighborCIdistance[i,])+1)] 
  KT500neighbordataControl2 <- KT500neighbordataset2[,-(which.min(KT500neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT500neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500neighbordataControl1[,2:ncol(KT500neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500neighbordataControl2[,2:ncol(KT500neighbordataControl2)])/(N-1))
  KT500neighbormeanControlsum1<-approx(KT500neighbordataControl1[,1],KT500neighbormeanControlsum, x_axis_KT500neighbor)$y
  
  KT500neighbormeanControltotal[i,]<-KT500neighbormeanControlsum1
}


KT500neighbordataset1 <- KT500neighborcontrollarvae[KK[rep]:attributes(as.matrix(KT500neighborcontrollarvae))$dim[1],]
KT500neighbordataset2 <- KT500neighborcontroladult[KK[rep]:attributes(as.matrix(KT500neighborcontroladult))$dim[1],]
KT500neighbordataControl1 <- KT500neighbordataset1[,-(which.min(KT500neighborCIdistance[rep,])+1)] 
KT500neighbordataControl2 <- KT500neighbordataset2[,-(which.min(KT500neighborCIdistance[rep,])+1)] 
##this 1 is the time column
KT500neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500neighbordataControl1[,2:ncol(KT500neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500neighbordataControl2[,2:ncol(KT500neighbordataControl2)])/(N-1))
KT500neighbormeanControlsum1<-approx(KT500neighbordataControl1[,1],KT500neighbormeanControlsum, x_axis_KT500neighbor)$y


KT500neighbormeanControltotal[rep,]<-KT500neighbormeanControlsum1


KT500neighboraverageControl <- colSums(KT500neighbormeanControltotal)/rep

KT500neighborControlCIlow<-c()
KT500neighborControlCIhigh<-c()

for(i in 1:length(KT500neighboraverageControl))
{
  KT500neighborControlCIlow[i]<-confidence_interval(KT500neighbormeanControltotal[,i], 0.95)[1][[1]]
  KT500neighborControlCIhigh[i]<-confidence_interval(KT500neighbormeanControltotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 500 global  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=500-strengthonly")

Ks500globalcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks500globalcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks500globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks500globalcontroladult)[1])
{
  if (Ks500globalcontroladult[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(Ks500globalcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks500globalmeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks500global <- Ks500globalcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks500global <- Ks500globalcontroladult[KK[rep]:dim(Ks500globalcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks500globaldataset1 <- Ks500globalcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks500globaldataset2 <- Ks500globalcontroladult[KK[i]:(KK[i+1]-1),]
  Ks500globaldataControl1 <- Ks500globaldataset1[,-(which.min(Ks500globalCIdistance[i,])+1)] 
  Ks500globaldataControl2 <- Ks500globaldataset2[,-(which.min(Ks500globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks500globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500globaldataControl1[,2:ncol(Ks500globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500globaldataControl2[,2:ncol(Ks500globaldataControl2)])/(N-1))
  Ks500globalmeanControlsum1<-approx(Ks500globaldataControl1[,1],Ks500globalmeanControlsum, x_axis_Ks500global)$y
  
  Ks500globalmeanControltotal[i,]<-Ks500globalmeanControlsum1
}

##last row of Ks500meanInfectiontotal

Ks500globaldataset1 <- Ks500globalcontrollarvae[KK[rep]:attributes(as.matrix(Ks500globalcontrollarvae))$dim[1],]
Ks500globaldataset2 <- Ks500globalcontroladult[KK[rep]:attributes(as.matrix(Ks500globalcontroladult))$dim[1],]
Ks500globaldataControl1 <- Ks500globaldataset1[,-(which.min(Ks500globalCIdistance[rep,])+1)] 
Ks500globaldataControl2 <- Ks500globaldataset2[,-(which.min(Ks500globalCIdistance[rep,])+1)] 
##this 1 is the time column
Ks500globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500globaldataControl1[,2:ncol(Ks500globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks500globaldataControl2[,2:ncol(Ks500globaldataControl2)])/(N-1))
Ks500globalmeanControlsum1<-approx(Ks500globaldataControl1[,1],Ks500globalmeanControlsum, x_axis_Ks500global)$y


Ks500globalmeanControltotal[rep,]<-Ks500globalmeanControlsum1


Ks500globalaverageControl <- colSums(Ks500globalmeanControltotal)/rep

Ks500globalControlCIlow<-c()
Ks500globalControlCIhigh<-c()

for(i in 1:length(Ks500globalaverageControl))
{
  Ks500globalControlCIlow[i]<-confidence_interval(Ks500globalmeanControltotal[,i], 0.95)[1][[1]]
  Ks500globalControlCIhigh[i]<-confidence_interval(Ks500globalmeanControltotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=500-time+strength+control")

KT500globalcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT500globalcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT500globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT500globalcontroladult)[1])
{
  if (KT500globalcontroladult[i,1]==0)
    KK <-c(KK, i)
}

#KT500meanControltotal<-matrix(, nrow = rep, ncol = diff(KK)[1])

lastlength <- dim(KT500globalcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT500globalmeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT500global <- KT500globalcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT500global <- KT500globalcontroladult[KK[rep]:dim(KT500globalcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  KT500globaldataset1 <- KT500globalcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT500globaldataset2 <- KT500globalcontroladult[KK[i]:(KK[i+1]-1),]
  KT500globaldataControl1 <- KT500globaldataset1[,-(which.min(KT500globalCIdistance[i,])+1)] 
  KT500globaldataControl2 <- KT500globaldataset2[,-(which.min(KT500globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT500globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500globaldataControl1[,2:ncol(KT500globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500globaldataControl2[,2:ncol(KT500globaldataControl2)])/(N-1))
  KT500globalmeanControlsum1<-approx(KT500globaldataControl1[,1],KT500globalmeanControlsum, x_axis_KT500global)$y
  
  KT500globalmeanControltotal[i,]<-KT500globalmeanControlsum1
}


KT500globaldataset1 <- KT500globalcontrollarvae[KK[rep]:attributes(as.matrix(KT500globalcontrollarvae))$dim[1],]
KT500globaldataset2 <- KT500globalcontroladult[KK[rep]:attributes(as.matrix(KT500globalcontroladult))$dim[1],]
KT500globaldataControl1 <- KT500globaldataset1[,-(which.min(KT500globalCIdistance[rep,])+1)] 
KT500globaldataControl2 <- KT500globaldataset2[,-(which.min(KT500globalCIdistance[rep,])+1)] 
##this 1 is the time column
KT500globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500globaldataControl1[,2:ncol(KT500globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT500globaldataControl2[,2:ncol(KT500globaldataControl2)])/(N-1))
KT500globalmeanControlsum1<-approx(KT500globaldataControl1[,1],KT500globalmeanControlsum, x_axis_KT500global)$y


KT500globalmeanControltotal[rep,]<-KT500globalmeanControlsum1


KT500globalaverageControl <- colSums(KT500globalmeanControltotal)/rep

KT500globalControlCIlow<-c()
KT500globalControlCIhigh<-c()

for(i in 1:length(KT500globalaverageControl))
{
  KT500globalControlCIlow[i]<-confidence_interval(KT500globalmeanControltotal[,i], 0.95)[1][[1]]
  KT500globalControlCIhigh[i]<-confidence_interval(KT500globalmeanControltotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 local  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-local-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800localcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks800localcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks800CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks800localcontroladult)[1])
{
  if (Ks800localcontroladult[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(Ks800localcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800meanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800local <- Ks800localcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800local <- Ks800localcontroladult[KK[rep]:dim(Ks800localcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks800dataset1 <- Ks800localcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks800dataset2 <- Ks800localcontroladult[KK[i]:(KK[i+1]-1),]
  Ks800dataControl1 <- Ks800dataset1[,-(which.min(Ks800CIdistance[i,])+1)] 
  Ks800dataControl2 <- Ks800dataset2[,-(which.min(Ks800CIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks800meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800dataControl1[,2:ncol(Ks800dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800dataControl2[,2:ncol(Ks800dataControl2)])/(N-1))
  Ks800meanControlsum1<-approx(Ks800dataControl1[,1],Ks800meanControlsum, x_axis_Ks800local)$y
  
  Ks800meanControltotal[i,]<-Ks800meanControlsum1
}


Ks800dataset1 <- Ks800localcontrollarvae[KK[rep]:attributes(as.matrix(Ks800localcontrollarvae))$dim[1],]
Ks800dataset2 <- Ks800localcontroladult[KK[rep]:attributes(as.matrix(Ks800localcontroladult))$dim[1],]
Ks800dataControl1 <- Ks800dataset1[,-(which.min(Ks800CIdistance[rep,])+1)] 
Ks800dataControl2 <- Ks800dataset2[,-(which.min(Ks800CIdistance[rep,])+1)] 
##this 1 is the time column
Ks800meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800dataControl1[,2:ncol(Ks800dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800dataControl2[,2:ncol(Ks800dataControl2)])/(N-1))
Ks800meanControlsum1<-approx(Ks800dataControl1[,1],Ks800meanControlsum, x_axis_Ks800local)$y


Ks800meanControltotal[rep,]<-Ks800meanControlsum1


Ks800averageControl <- colSums(Ks800meanControltotal)/rep

Ks800ControlCIlow<-c()
Ks800ControlCIhigh<-c()

for(i in 1:length(Ks800averageControl))
{
  Ks800ControlCIlow[i]<-confidence_interval(Ks800meanControltotal[,i], 0.95)[1][[1]]
  Ks800ControlCIhigh[i]<-confidence_interval(Ks800meanControltotal[,i],0.95)[2][[1]]
}


#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-local-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800localcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT800localcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT800CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)


####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT800localcontroladult)[1])
{
  if (KT800localcontroladult[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT800localcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800meanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800local <- KT800localcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800local <- KT800localcontroladult[KK[rep]:dim(KT800localcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  KT800dataset1 <- KT800localcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT800dataset2 <- KT800localcontroladult[KK[i]:(KK[i+1]-1),]
  KT800dataControl1 <- KT800dataset1[,-(which.min(KT800CIdistance[i,])+1)] 
  KT800dataControl2 <- KT800dataset2[,-(which.min(KT800CIdistance[i,])+1)] 
  ##this 1 is the time column
  KT800meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800dataControl1[,2:ncol(KT800dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800dataControl2[,2:ncol(KT800dataControl2)])/(N-1))
  KT800meanControlsum1<-approx(KT800dataControl1[,1],KT800meanControlsum, x_axis_KT800local)$y
  
  KT800meanControltotal[i,]<-KT800meanControlsum1
}


KT800dataset1 <- KT800localcontrollarvae[KK[rep]:attributes(as.matrix(KT800localcontrollarvae))$dim[1],]
KT800dataset2 <- KT800localcontroladult[KK[rep]:attributes(as.matrix(KT800localcontroladult))$dim[1],]
KT800dataControl1 <- KT800dataset1[,-(which.min(KT800CIdistance[rep,])+1)] 
KT800dataControl2 <- KT800dataset2[,-(which.min(KT800CIdistance[rep,])+1)] 
##this 1 is the time column
KT800meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800dataControl1[,2:ncol(KT800dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800dataControl2[,2:ncol(KT800dataControl2)])/(N-1))
KT800meanControlsum1<-approx(KT800dataControl1[,1],KT800meanControlsum, x_axis_KT800local)$y


KT800meanControltotal[rep,]<-KT800meanControlsum1


KT800averageControl <- colSums(KT800meanControltotal)/rep

KT800ControlCIlow<-c()
KT800ControlCIhigh<-c()

for(i in 1:length(KT800averageControl))
{
  KT800ControlCIlow[i]<-confidence_interval(KT800meanControltotal[,i], 0.95)[1][[1]]
  KT800ControlCIhigh[i]<-confidence_interval(KT800meanControltotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 neighbor  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800neighborcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks800neighborcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks800neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks800neighborcontroladult)[1])
{
  if (Ks800neighborcontroladult[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks800neighborcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800neighbormeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800neighbor <- Ks800neighborcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800neighbor <- Ks800neighborcontroladult[KK[rep]:dim(Ks800neighborcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks800neighbordataset1 <- Ks800neighborcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks800neighbordataset2 <- Ks800neighborcontroladult[KK[i]:(KK[i+1]-1),]
  Ks800neighbordataControl1 <- Ks800neighbordataset1[,-(which.min(Ks800neighborCIdistance[i,])+1)] 
  Ks800neighbordataControl2 <- Ks800neighbordataset2[,-(which.min(Ks800neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks800neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800neighbordataControl1[,2:ncol(Ks800neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800neighbordataControl2[,2:ncol(Ks800neighbordataControl2)])/(N-1))
  Ks800neighbormeanControlsum1<-approx(Ks800neighbordataControl1[,1],Ks800neighbormeanControlsum, x_axis_Ks800neighbor)$y
  
  Ks800neighbormeanControltotal[i,]<-Ks800neighbormeanControlsum1
}


Ks800neighbordataset1 <- Ks800neighborcontrollarvae[KK[rep]:attributes(as.matrix(Ks800neighborcontrollarvae))$dim[1],]
Ks800neighbordataset2 <- Ks800neighborcontroladult[KK[rep]:attributes(as.matrix(Ks800neighborcontroladult))$dim[1],]
Ks800neighbordataControl1 <- Ks800neighbordataset1[,-(which.min(Ks800neighborCIdistance[rep,])+1)] 
Ks800neighbordataControl2 <- Ks800neighbordataset2[,-(which.min(Ks800neighborCIdistance[rep,])+1)] 
##this 1 is the time column
Ks800neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800neighbordataControl1[,2:ncol(Ks800neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800neighbordataControl2[,2:ncol(Ks800neighbordataControl2)])/(N-1))
Ks800neighbormeanControlsum1<-approx(Ks800neighbordataControl1[,1],Ks800neighbormeanControlsum, x_axis_Ks800neighbor)$y


Ks800neighbormeanControltotal[rep,]<-Ks800neighbormeanControlsum1


Ks800neighboraverageControl <- colSums(Ks800neighbormeanControltotal)/rep

Ks800neighborControlCIlow<-c()
Ks800neighborControlCIhigh<-c()

for(i in 1:length(Ks800neighboraverageControl))
{
  Ks800neighborControlCIlow[i]<-confidence_interval(Ks800neighbormeanControltotal[,i], 0.95)[1][[1]]
  Ks800neighborControlCIhigh[i]<-confidence_interval(Ks800neighbormeanControltotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800neighborcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT800neighborcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT800neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT800neighborcontroladult)[1])
{
  if (KT800neighborcontroladult[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(KT800neighborcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800neighbormeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800neighbor <- KT800neighborcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800neighbor <- KT800neighborcontroladult[KK[rep]:dim(KT800neighborcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  KT800neighbordataset1 <- KT800neighborcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT800neighbordataset2 <- KT800neighborcontroladult[KK[i]:(KK[i+1]-1),]
  KT800neighbordataControl1 <- KT800neighbordataset1[,-(which.min(KT800neighborCIdistance[i,])+1)] 
  KT800neighbordataControl2 <- KT800neighbordataset2[,-(which.min(KT800neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT800neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800neighbordataControl1[,2:ncol(KT800neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800neighbordataControl2[,2:ncol(KT800neighbordataControl2)])/(N-1))
  KT800neighbormeanControlsum1<-approx(KT800neighbordataControl1[,1],KT800neighbormeanControlsum, x_axis_KT800neighbor)$y
  
  KT800neighbormeanControltotal[i,]<-KT800neighbormeanControlsum1
}

KT800neighbordataset1 <- KT800neighborcontrollarvae[KK[rep]:attributes(as.matrix(KT800neighborcontrollarvae))$dim[1],]
KT800neighbordataset2 <- KT800neighborcontroladult[KK[rep]:attributes(as.matrix(KT800neighborcontroladult))$dim[1],]
KT800neighbordataControl1 <- KT800neighbordataset1[,-(which.min(KT800neighborCIdistance[rep,])+1)] 
KT800neighbordataControl2 <- KT800neighbordataset2[,-(which.min(KT800neighborCIdistance[rep,])+1)] 
##this 1 is the time column
KT800neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800neighbordataControl1[,2:ncol(KT800neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800neighbordataControl2[,2:ncol(KT800neighbordataControl2)])/(N-1))
KT800neighbormeanControlsum1<-approx(KT800neighbordataControl1[,1],KT800neighbormeanControlsum, x_axis_KT800neighbor)$y

KT800neighbormeanControltotal[rep,]<-KT800neighbormeanControlsum1

KT800neighboraverageControl <- colSums(KT800neighbormeanControltotal)/rep

KT800neighborControlCIlow<-c()
KT800neighborControlCIhigh<-c()

for(i in 1:length(KT800neighboraverageControl))
{
  KT800neighborControlCIlow[i]<-confidence_interval(KT800neighbormeanControltotal[,i], 0.95)[1][[1]]
  KT800neighborControlCIhigh[i]<-confidence_interval(KT800neighbormeanControltotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 800 global  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=800-strengthonly")

Ks800globalcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks800globalcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks800globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks800globalcontroladult)[1])
{
  if (Ks800globalcontroladult[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(Ks800globalcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks800globalmeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks800global <- Ks800globalcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks800global <- Ks800globalcontroladult[KK[rep]:dim(Ks800globalcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks800globaldataset1 <- Ks800globalcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks800globaldataset2 <- Ks800globalcontroladult[KK[i]:(KK[i+1]-1),]
  Ks800globaldataControl1 <- Ks800globaldataset1[,-(which.min(Ks800globalCIdistance[i,])+1)] 
  Ks800globaldataControl2 <- Ks800globaldataset2[,-(which.min(Ks800globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks800globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800globaldataControl1[,2:ncol(Ks800globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800globaldataControl2[,2:ncol(Ks800globaldataControl2)])/(N-1))
  Ks800globalmeanControlsum1<-approx(Ks800globaldataControl1[,1],Ks800globalmeanControlsum, x_axis_Ks800global)$y
  
  Ks800globalmeanControltotal[i,]<-Ks800globalmeanControlsum1
}


Ks800globaldataset1 <- Ks800globalcontrollarvae[KK[rep]:attributes(as.matrix(Ks800globalcontrollarvae))$dim[1],]
Ks800globaldataset2 <- Ks800globalcontroladult[KK[rep]:attributes(as.matrix(Ks800globalcontroladult))$dim[1],]
Ks800globaldataControl1 <- Ks800globaldataset1[,-(which.min(Ks800globalCIdistance[rep,])+1)] 
Ks800globaldataControl2 <- Ks800globaldataset2[,-(which.min(Ks800globalCIdistance[rep,])+1)] 
##this 1 is the time column
Ks800globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800globaldataControl1[,2:ncol(Ks800globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks800globaldataControl2[,2:ncol(Ks800globaldataControl2)])/(N-1))
Ks800globalmeanControlsum1<-approx(Ks800globaldataControl1[,1],Ks800globalmeanControlsum, x_axis_Ks800global)$y


Ks800globalmeanControltotal[rep,]<-Ks800globalmeanControlsum1


Ks800globalaverageControl <- colSums(Ks800globalmeanControltotal)/rep

Ks800globalControlCIlow<-c()
Ks800globalControlCIhigh<-c()

for(i in 1:length(Ks800globalaverageControl))
{
  Ks800globalControlCIlow[i]<-confidence_interval(Ks800globalmeanControltotal[,i], 0.95)[1][[1]]
  Ks800globalControlCIhigh[i]<-confidence_interval(Ks800globalmeanControltotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=800-time+strength+control")

KT800globalcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT800globalcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT800globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT800globalcontroladult)[1])
{
  if (KT800globalcontroladult[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT800globalcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT800globalmeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT800global <- KT800globalcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT800global <- KT800globalcontroladult[KK[rep]:dim(KT800globalcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  KT800globaldataset1 <- KT800globalcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT800globaldataset2 <- KT800globalcontroladult[KK[i]:(KK[i+1]-1),]
  KT800globaldataControl1 <- KT800globaldataset1[,-(which.min(KT800globalCIdistance[i,])+1)] 
  KT800globaldataControl2 <- KT800globaldataset2[,-(which.min(KT800globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT800globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800globaldataControl1[,2:ncol(KT800globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800globaldataControl2[,2:ncol(KT800globaldataControl2)])/(N-1))
  KT800globalmeanControlsum1<-approx(KT800globaldataControl1[,1],KT800globalmeanControlsum, x_axis_KT800global)$y
  
  KT800globalmeanControltotal[i,]<-KT800globalmeanControlsum1
}

##last row of KT800meanInfectiontotal

KT800globaldataset1 <- KT800globalcontrollarvae[KK[rep]:attributes(as.matrix(KT800globalcontrollarvae))$dim[1],]
KT800globaldataset2 <- KT800globalcontroladult[KK[rep]:attributes(as.matrix(KT800globalcontroladult))$dim[1],]
KT800globaldataControl1 <- KT800globaldataset1[,-(which.min(KT800globalCIdistance[rep,])+1)] 
KT800globaldataControl2 <- KT800globaldataset2[,-(which.min(KT800globalCIdistance[rep,])+1)] 
##this 1 is the time column
KT800globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800globaldataControl1[,2:ncol(KT800globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT800globaldataControl2[,2:ncol(KT800globaldataControl2)])/(N-1))
KT800globalmeanControlsum1<-approx(KT800globaldataControl1[,1],KT800globalmeanControlsum, x_axis_KT800global)$y


KT800globalmeanControltotal[rep,]<-KT800globalmeanControlsum1


KT800globalaverageControl <- colSums(KT800globalmeanControltotal)/rep

KT800globalControlCIlow<-c()
KT800globalControlCIhigh<-c()

for(i in 1:length(KT800globalaverageControl))
{
  KT800globalControlCIlow[i]<-confidence_interval(KT800globalmeanControltotal[,i], 0.95)[1][[1]]
  KT800globalControlCIhigh[i]<-confidence_interval(KT800globalmeanControltotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 2000 local  (strength first, added timing later)

##################################################################################################################################
setwd("~\\SF-local-N-9-p-0.01-t-300-K=2000-strengthonly")

Ks2000localcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks2000localcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks2000CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks2000localcontroladult)[1])
{
  if (Ks2000localcontroladult[i,1]==0)
    KK <-c(KK, i)
}

#Ks2000meanControltotal<-matrix(, nrow = rep, ncol = diff(KK)[1])

lastlength <- dim(Ks2000localcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000meanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000local <- Ks2000localcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000local <- Ks2000localcontroladult[KK[rep]:dim(Ks2000localcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks2000dataset1 <- Ks2000localcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks2000dataset2 <- Ks2000localcontroladult[KK[i]:(KK[i+1]-1),]
  Ks2000dataControl1 <- Ks2000dataset1[,-(which.min(Ks2000CIdistance[i,])+1)] 
  Ks2000dataControl2 <- Ks2000dataset2[,-(which.min(Ks2000CIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks2000meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000dataControl1[,2:ncol(Ks2000dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000dataControl2[,2:ncol(Ks2000dataControl2)])/(N-1))
  Ks2000meanControlsum1<-approx(Ks2000dataControl1[,1],Ks2000meanControlsum, x_axis_Ks2000local)$y
  
  Ks2000meanControltotal[i,]<-Ks2000meanControlsum1
}


Ks2000dataset1 <- Ks2000localcontrollarvae[KK[rep]:attributes(as.matrix(Ks2000localcontrollarvae))$dim[1],]
Ks2000dataset2 <- Ks2000localcontroladult[KK[rep]:attributes(as.matrix(Ks2000localcontroladult))$dim[1],]
Ks2000dataControl1 <- Ks2000dataset1[,-(which.min(Ks2000CIdistance[rep,])+1)] 
Ks2000dataControl2 <- Ks2000dataset2[,-(which.min(Ks2000CIdistance[rep,])+1)] 
##this 1 is the time column
Ks2000meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000dataControl1[,2:ncol(Ks2000dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000dataControl2[,2:ncol(Ks2000dataControl2)])/(N-1))
Ks2000meanControlsum1<-approx(Ks2000dataControl1[,1],Ks2000meanControlsum, x_axis_Ks2000local)$y


Ks2000meanControltotal[rep,]<-Ks2000meanControlsum1


Ks2000averageControl <- colSums(Ks2000meanControltotal)/rep

Ks2000ControlCIlow<-c()
Ks2000ControlCIhigh<-c()

for(i in 1:length(Ks2000averageControl))
{
  Ks2000ControlCIlow[i]<-confidence_interval(Ks2000meanControltotal[,i], 0.95)[1][[1]]
  Ks2000ControlCIhigh[i]<-confidence_interval(Ks2000meanControltotal[,i],0.95)[2][[1]]
}

#####strength +time
#######################################################

###define the tspan in x axis
setwd("~\\SF-local-N-9-p-0.01-t-300-K=2000-time+strength+control")

KT2000localcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT2000localcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT2000CIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT2000localcontroladult)[1])
{
  if (KT2000localcontroladult[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT2000localcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000meanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000local <- KT2000localcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000local <- KT2000localcontroladult[KK[rep]:dim(KT2000localcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  KT2000dataset1 <- KT2000localcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT2000dataset2 <- KT2000localcontroladult[KK[i]:(KK[i+1]-1),]
  KT2000dataControl1 <- KT2000dataset1[,-(which.min(KT2000CIdistance[i,])+1)] 
  KT2000dataControl2 <- KT2000dataset2[,-(which.min(KT2000CIdistance[i,])+1)] 
  ##this 1 is the time column
  KT2000meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000dataControl1[,2:ncol(KT2000dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000dataControl2[,2:ncol(KT2000dataControl2)])/(N-1))
  KT2000meanControlsum1<-approx(KT2000dataControl1[,1],KT2000meanControlsum, x_axis_KT2000local)$y
  
  KT2000meanControltotal[i,]<-KT2000meanControlsum1
}


KT2000dataset1 <- KT2000localcontrollarvae[KK[rep]:attributes(as.matrix(KT2000localcontrollarvae))$dim[1],]
KT2000dataset2 <- KT2000localcontroladult[KK[rep]:attributes(as.matrix(KT2000localcontroladult))$dim[1],]
KT2000dataControl1 <- KT2000dataset1[,-(which.min(KT2000CIdistance[rep,])+1)] 
KT2000dataControl2 <- KT2000dataset2[,-(which.min(KT2000CIdistance[rep,])+1)] 
##this 1 is the time column
KT2000meanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000dataControl1[,2:ncol(KT2000dataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000dataControl2[,2:ncol(KT2000dataControl2)])/(N-1))
KT2000meanControlsum1<-approx(KT2000dataControl1[,1],KT2000meanControlsum, x_axis_KT2000local)$y


KT2000meanControltotal[rep,]<-KT2000meanControlsum1


KT2000averageControl <- colSums(KT2000meanControltotal)/rep

KT2000ControlCIlow<-c()
KT2000ControlCIhigh<-c()

for(i in 1:length(KT2000averageControl))
{
  KT2000ControlCIlow[i]<-confidence_interval(KT2000meanControltotal[,i], 0.95)[1][[1]]
  KT2000ControlCIhigh[i]<-confidence_interval(KT2000meanControltotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 2000 neighbor  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=2000-strengthonly")

Ks2000neighborcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks2000neighborcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks2000neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks2000neighborcontroladult)[1])
{
  if (Ks2000neighborcontroladult[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(Ks2000neighborcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000neighbormeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000neighbor <- Ks2000neighborcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000neighbor <- Ks2000neighborcontroladult[KK[rep]:dim(Ks2000neighborcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks2000neighbordataset1 <- Ks2000neighborcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks2000neighbordataset2 <- Ks2000neighborcontroladult[KK[i]:(KK[i+1]-1),]
  Ks2000neighbordataControl1 <- Ks2000neighbordataset1[,-(which.min(Ks2000neighborCIdistance[i,])+1)] 
  Ks2000neighbordataControl2 <- Ks2000neighbordataset2[,-(which.min(Ks2000neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks2000neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000neighbordataControl1[,2:ncol(Ks2000neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000neighbordataControl2[,2:ncol(Ks2000neighbordataControl2)])/(N-1))
  Ks2000neighbormeanControlsum1<-approx(Ks2000neighbordataControl1[,1],Ks2000neighbormeanControlsum, x_axis_Ks2000neighbor)$y
  
  Ks2000neighbormeanControltotal[i,]<-Ks2000neighbormeanControlsum1
}


Ks2000neighbordataset1 <- Ks2000neighborcontrollarvae[KK[rep]:attributes(as.matrix(Ks2000neighborcontrollarvae))$dim[1],]
Ks2000neighbordataset2 <- Ks2000neighborcontroladult[KK[rep]:attributes(as.matrix(Ks2000neighborcontroladult))$dim[1],]
Ks2000neighbordataControl1 <- Ks2000neighbordataset1[,-(which.min(Ks2000neighborCIdistance[rep,])+1)] 
Ks2000neighbordataControl2 <- Ks2000neighbordataset2[,-(which.min(Ks2000neighborCIdistance[rep,])+1)] 
##this 1 is the time column
Ks2000neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000neighbordataControl1[,2:ncol(Ks2000neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000neighbordataControl2[,2:ncol(Ks2000neighbordataControl2)])/(N-1))
Ks2000neighbormeanControlsum1<-approx(Ks2000neighbordataControl1[,1],Ks2000neighbormeanControlsum, x_axis_Ks2000neighbor)$y


Ks2000neighbormeanControltotal[rep,]<-Ks2000neighbormeanControlsum1


Ks2000neighboraverageControl <- colSums(Ks2000neighbormeanControltotal)/rep

Ks2000neighborControlCIlow<-c()
Ks2000neighborControlCIhigh<-c()

for(i in 1:length(Ks2000neighboraverageControl))
{
  Ks2000neighborControlCIlow[i]<-confidence_interval(Ks2000neighbormeanControltotal[,i], 0.95)[1][[1]]
  Ks2000neighborControlCIhigh[i]<-confidence_interval(Ks2000neighbormeanControltotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-neighbor-N-9-p-0.01-t-300-K=2000-time+strength+control-9-23new")

KT2000neighborcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT2000neighborcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT2000neighborCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT2000neighborcontroladult)[1])
{
  if (KT2000neighborcontroladult[i,1]==0)
    KK <-c(KK, i)
}

lastlength <- dim(KT2000neighborcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000neighbormeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000neighbor <- KT2000neighborcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000neighbor <- KT2000neighborcontroladult[KK[rep]:dim(KT2000neighborcontroladult)[1],1]
}



for (i in 1:(rep-1))
{
  KT2000neighbordataset1 <- KT2000neighborcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT2000neighbordataset2 <- KT2000neighborcontroladult[KK[i]:(KK[i+1]-1),]
  KT2000neighbordataControl1 <- KT2000neighbordataset1[,-(which.min(KT2000neighborCIdistance[i,])+1)] 
  KT2000neighbordataControl2 <- KT2000neighbordataset2[,-(which.min(KT2000neighborCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT2000neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000neighbordataControl1[,2:ncol(KT2000neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000neighbordataControl2[,2:ncol(KT2000neighbordataControl2)])/(N-1))
  KT2000neighbormeanControlsum1<-approx(KT2000neighbordataControl1[,1],KT2000neighbormeanControlsum, x_axis_KT2000neighbor)$y
  
  KT2000neighbormeanControltotal[i,]<-KT2000neighbormeanControlsum1
}

KT2000neighbordataset1 <- KT2000neighborcontrollarvae[KK[rep]:attributes(as.matrix(KT2000neighborcontrollarvae))$dim[1],]
KT2000neighbordataset2 <- KT2000neighborcontroladult[KK[rep]:attributes(as.matrix(KT2000neighborcontroladult))$dim[1],]
KT2000neighbordataControl1 <- KT2000neighbordataset1[,-(which.min(KT2000neighborCIdistance[rep,])+1)] 
KT2000neighbordataControl2 <- KT2000neighbordataset2[,-(which.min(KT2000neighborCIdistance[rep,])+1)] 
##this 1 is the time column
KT2000neighbormeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000neighbordataControl1[,2:ncol(KT2000neighbordataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000neighbordataControl2[,2:ncol(KT2000neighbordataControl2)])/(N-1))
KT2000neighbormeanControlsum1<-approx(KT2000neighbordataControl1[,1],KT2000neighbormeanControlsum, x_axis_KT2000neighbor)$y


KT2000neighbormeanControltotal[rep,]<-KT2000neighbormeanControlsum1


KT2000neighboraverageControl <- colSums(KT2000neighbormeanControltotal)/rep

KT2000neighborControlCIlow<-c()
KT2000neighborControlCIhigh<-c()

for(i in 1:length(KT2000neighboraverageControl))
{
  KT2000neighborControlCIlow[i]<-confidence_interval(KT2000neighbormeanControltotal[,i], 0.95)[1][[1]]
  KT2000neighborControlCIhigh[i]<-confidence_interval(KT2000neighbormeanControltotal[,i],0.95)[2][[1]]
}

##################################################################################################################################

# K = 2000 global  (strength first, added timing later)

##################################################################################################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=2000-strengthonly")

Ks2000globalcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

Ks2000globalcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

Ks2000globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(Ks2000globalcontroladult)[1])
{
  if (Ks2000globalcontroladult[i,1]==0)
    KK <-c(KK, i)
}


lastlength <- dim(Ks2000globalcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

Ks2000globalmeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_Ks2000global <- Ks2000globalcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_Ks2000global <- Ks2000globalcontroladult[KK[rep]:dim(Ks2000globalcontroladult)[1],1]
}


for (i in 1:(rep-1))
{
  Ks2000globaldataset1 <- Ks2000globalcontrollarvae[KK[i]:(KK[i+1]-1),]
  Ks2000globaldataset2 <- Ks2000globalcontroladult[KK[i]:(KK[i+1]-1),]
  Ks2000globaldataControl1 <- Ks2000globaldataset1[,-(which.min(Ks2000globalCIdistance[i,])+1)] 
  Ks2000globaldataControl2 <- Ks2000globaldataset2[,-(which.min(Ks2000globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  Ks2000globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000globaldataControl1[,2:ncol(Ks2000globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000globaldataControl2[,2:ncol(Ks2000globaldataControl2)])/(N-1))
  Ks2000globalmeanControlsum1<-approx(Ks2000globaldataControl1[,1],Ks2000globalmeanControlsum, x_axis_Ks2000global)$y
  
  Ks2000globalmeanControltotal[i,]<-Ks2000globalmeanControlsum1
}

Ks2000globaldataset1 <- Ks2000globalcontrollarvae[KK[rep]:attributes(as.matrix(Ks2000globalcontrollarvae))$dim[1],]
Ks2000globaldataset2 <- Ks2000globalcontroladult[KK[rep]:attributes(as.matrix(Ks2000globalcontroladult))$dim[1],]
Ks2000globaldataControl1 <- Ks2000globaldataset1[,-(which.min(Ks2000globalCIdistance[rep,])+1)] 
Ks2000globaldataControl2 <- Ks2000globaldataset2[,-(which.min(Ks2000globalCIdistance[rep,])+1)] 
##this 1 is the time column
Ks2000globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000globaldataControl1[,2:ncol(Ks2000globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*Ks2000globaldataControl2[,2:ncol(Ks2000globaldataControl2)])/(N-1))
Ks2000globalmeanControlsum1<-approx(Ks2000globaldataControl1[,1],Ks2000globalmeanControlsum, x_axis_Ks2000global)$y


Ks2000globalmeanControltotal[rep,]<-Ks2000globalmeanControlsum1


Ks2000globalaverageControl <- colSums(Ks2000globalmeanControltotal)/rep

Ks2000globalControlCIlow<-c()
Ks2000globalControlCIhigh<-c()

for(i in 1:length(Ks2000globalaverageControl))
{
  Ks2000globalControlCIlow[i]<-confidence_interval(Ks2000globalmeanControltotal[,i], 0.95)[1][[1]]
  Ks2000globalControlCIhigh[i]<-confidence_interval(Ks2000globalmeanControltotal[,i],0.95)[2][[1]]
}

#####only strength
############################################################
##################################################################################################################################
setwd("~\\SF-global-N-9-p-0.01-t-300-K=2000-time+strength+control")

KT2000globalcontrollarvae <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICICDES-LARVAE.csv",header=FALSE)

KT2000globalcontroladult <- read.csv("2-a1--a2--g1--g2--e1--e2--PESTICIDES-ADULTS.csv", header=FALSE)

KT2000globalCIdistance <- read.csv("2-a1--a2--g1--g2--e1--e2--DISTANCE-PATCH.csv",header = FALSE)

####divde rows whenever all 0s met
KK<-c()
for ( i in 1:dim(KT2000globalcontroladult)[1])
{
  if (KT2000globalcontroladult[i,1]==0)
    KK <-c(KK, i)
}

#KT2000meanControltotal<-matrix(, nrow = rep, ncol = diff(KK)[1])

lastlength <- dim(KT2000globalcontroladult)[1]-KK[length(KK)]+1
diff <- c(diff(KK), lastlength)

KT2000globalmeanControltotal<-matrix(, nrow = rep, ncol = max(diff))

index <- which.max(diff)

if(index < rep)
{
  x_axis_KT2000global <- KT2000globalcontroladult[KK[index]:(KK[index+1]-1),1]
} else
{
  x_axis_KT2000global <- KT2000globalcontroladult[KK[rep]:dim(KT2000globalcontroladult)[1],1]
}

for (i in 1:(rep-1))
{
  KT2000globaldataset1 <- KT2000globalcontrollarvae[KK[i]:(KK[i+1]-1),]
  KT2000globaldataset2 <- KT2000globalcontroladult[KK[i]:(KK[i+1]-1),]
  KT2000globaldataControl1 <- KT2000globaldataset1[,-(which.min(KT2000globalCIdistance[i,])+1)] 
  KT2000globaldataControl2 <- KT2000globaldataset2[,-(which.min(KT2000globalCIdistance[i,])+1)] 
  ##this 1 is the time column
  KT2000globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000globaldataControl1[,2:ncol(KT2000globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000globaldataControl2[,2:ncol(KT2000globaldataControl2)])/(N-1))
  KT2000globalmeanControlsum1<-approx(KT2000globaldataControl1[,1],KT2000globalmeanControlsum, x_axis_KT2000global)$y
  
  KT2000globalmeanControltotal[i,]<-KT2000globalmeanControlsum1
}


KT2000globaldataset1 <- KT2000globalcontrollarvae[KK[rep]:attributes(as.matrix(KT2000globalcontrollarvae))$dim[1],]
KT2000globaldataset2 <- KT2000globalcontroladult[KK[rep]:attributes(as.matrix(KT2000globalcontroladult))$dim[1],]
KT2000globaldataControl1 <- KT2000globaldataset1[,-(which.min(KT2000globalCIdistance[rep,])+1)] 
KT2000globaldataControl2 <- KT2000globaldataset2[,-(which.min(KT2000globalCIdistance[rep,])+1)] 
##this 1 is the time column
KT2000globalmeanControlsum <- as.vector(rowSums((epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000globaldataControl1[,2:ncol(KT2000globaldataControl1)]+(epsilon[1]+epsilon[2])/(epsilon[1]+epsilon[2])*KT2000globaldataControl2[,2:ncol(KT2000globaldataControl2)])/(N-1))
KT2000globalmeanControlsum1<-approx(KT2000globaldataControl1[,1],KT2000globalmeanControlsum, x_axis_KT2000global)$y


KT2000globalmeanControltotal[rep,]<-KT2000globalmeanControlsum1


KT2000globalaverageControl <- colSums(KT2000globalmeanControltotal)/rep

KT2000globalControlCIlow<-c()
KT2000globalControlCIhigh<-c()

for(i in 1:length(KT2000globalaverageControl))
{
  KT2000globalControlCIlow[i]<-confidence_interval(KT2000globalmeanControltotal[,i], 0.95)[1][[1]]
  KT2000globalControlCIhigh[i]<-confidence_interval(KT2000globalmeanControltotal[,i],0.95)[2][[1]]
}

###################################################################################################################################

#PLOT the above all figures (K = 500, 800 and K = 2000)

####################################################################################################################################

########################################################################################
##########################################################################################

tiff("Figure_1.tiff", width=10,height=10, units='in',res=600)

par(mfrow=c(3,3))
par(mar=c(2,3,3,3))

### K = 500 local
plot(x_axis_KT500local,KT500averageControl , type = "l", col = "red", lwd = 2, ylim = c(0,0.042))
points(x_axis_KT500local,KT500ControlCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT500local,KT500ControlCIhigh, type = "l", lty = 2, col = "red")

points(x_axis_Ks500local,Ks500averageControl, type = "l", col = "blue", lwd = 2)

### K = 500 neighbor
plot(x_axis_KT500neighbor,KT500neighboraverageControl, type = "l", col = "red", lwd = 2, ylim = c(0,0.042))
points(x_axis_KT500neighbor,KT500neighborControlCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT500neighbor,KT500neighborControlCIhigh, type = "l", lty =  2, col = "red")

points(x_axis_Ks500neighbor,Ks500neighboraverageControl, type = "l", col = "blue", lwd = 2)


### K = 500 global

plot(x_axis_KT500global,KT500globalaverageControl, type = "l", col = "red", lwd = 2, ylim = c(0,0.042))
points(x_axis_KT500global,KT500globalControlCIlow, type = "l", lty=2, col = "red")
points(x_axis_KT500global,KT500globalControlCIhigh, type = "l", lty=2, col = "red")

points(x_axis_Ks500global,Ks500globalaverageControl, type = "l", col = "blue", lwd = 2)


### K = 800 local
plot(x_axis_KT800local,KT800averageControl, type = "l", col = "red", lwd = 2, ylim = c(0,0.04))
points(x_axis_KT800local,KT800ControlCIlow, type = "l", lty=2, col = "red")
points(x_axis_KT800local,KT800ControlCIhigh, type = "l", lty=2, col = "red")

points(x_axis_Ks800local,Ks800averageControl, type = "l", col = "blue", lwd = 2)

### K = 800 neighbor
plot(x_axis_KT800neighbor,KT800neighboraverageControl, type = "l", col = "red", lwd = 2, ylim = c(0,0.04))
points(x_axis_KT800neighbor,KT800neighborControlCIlow, type = "l", lty=2, col = "red")
points(x_axis_KT800neighbor,KT800neighborControlCIhigh, type = "l", lty=2, col = "red")

points(x_axis_Ks800neighbor,Ks800neighboraverageControl, type = "l", col = "blue", lwd = 2)


### K = 800 global

plot(x_axis_KT800global,KT800globalaverageControl, type = "l", col = "red", lwd = 2, ylim = c(0,0.04))
points(x_axis_KT800global,KT800globalControlCIlow, type = "l", lty=2, col = "red")
points(x_axis_KT800global,KT800globalControlCIhigh, type = "l", lty=2, col = "red")

points(x_axis_Ks800global,Ks800globalaverageControl, type = "l", col = "blue", lwd = 2)


### K = 2000 local
plot(x_axis_KT2000local,KT2000averageControl, type = "l", col = "red", lwd = 2, ylim = c(0,0.25))
points(x_axis_KT2000local,KT2000ControlCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT2000local,KT2000ControlCIhigh, type = "l", lty = 2, col = "red")

points(x_axis_Ks2000local,Ks2000averageControl, type = "l", col = "blue", lwd = 2)

### K = 2000 neighbor
plot(x_axis_KT2000neighbor,KT2000neighboraverageControl, type = "l", col = "red", lwd = 2, ylim = c(0,0.25))
points(x_axis_KT2000neighbor,KT2000neighborControlCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT2000neighbor,KT2000neighborControlCIhigh, type = "l", lty = 2, col = "red")

points(x_axis_Ks2000neighbor,Ks2000neighboraverageControl, type = "l", col = "blue", lwd = 2)


### K = 2000 global

plot(x_axis_KT2000global,KT2000globalaverageControl, type = "l", col = "red", lwd = 2, ylim = c(0,0.25))
points(x_axis_KT2000global,KT2000globalControlCIlow, type = "l", lty = 2, col = "red")
points(x_axis_KT2000global,KT2000globalControlCIhigh, type = "l", lty = 2, col = "red")

points(x_axis_Ks2000global,Ks2000globalaverageControl, type = "l", col = "blue", lwd = 2)

dev.off()


####################################################################################################################################

## Results from NO Control at three levels of breeding capacity. This section is used for calculating control efficacy in Figure 3.

####################################################################################################################################
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

setwd("~\\SF-local-N-9-p-0.01-t-300-K=800-nocontrol")

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


K500nocontrol = rep(K500localnocontrol[601],length(cont))
K800nocontrol = rep(K800localnocontrol[601],length(cont))
K2000nocontrol = rep(K2000localnocontrol[601],length(cont))



############################################################################################
############################################################################################

## Index for control efficacy = (Infectionwithoutcontrol - InfectionwithControl)/MosquitoControl

############################################################################################
############################################################################################

###### K = 500 Local

IndexKs500Local = -(Ks500averageInfection - K500localnocontrol)/Ks500averageControl


IndexKT500Local = -(KT500averageInfection - K500localnocontrol)/KT500averageControl

###### K = 500 Neighbor

IndexKs500neighbor = -(Ks500neighboraverageInfection - K500neighbornocontrol)/Ks500neighboraverageControl


IndexKT500neighbor = -(KT500neighboraverageInfection - K500neighbornocontrol)/KT500neighboraverageControl

###### K = 500 Global

IndexKs500global = -(Ks500globalaverageInfection - K500globalnocontrol)/Ks500globalaverageControl


IndexKT500global = -(KT500globalaverageInfection - K500globalnocontrol)/KT500globalaverageControl


######## K = 800 Local

IndexKs800local = -(Ks800localaverageInfection - K800localnocontrol)/Ks800averageControl


IndexKT800local = -(KT800localaverageInfection - K800localnocontrol)/KT800averageControl

######## K = 800 neighbor

IndexKs800neighbor = -(Ks800neighboraverageInfection - K800neighbornocontrol)/Ks800neighboraverageControl


IndexKT800neighbor = -(KT800neighboraverageInfection - K800neighbornocontrol)/KT800neighboraverageControl


######## K = 800 global

IndexKs800global = -(Ks800globalaverageInfection - K800globalnocontrol)/Ks800globalaverageControl


IndexKT800global = -(KT800globalaverageInfection - K800globalnocontrol)/KT800globalaverageControl


######## K = 2000 Local


IndexKs2000local = -(KT2000localaverageInfection - K2000localnocontrol)/Ks2000averageControl


IndexKT2000local = -(KT2000localaverageInfection - K2000localnocontrol)/KT2000averageControl

######## K = 2000 neighbor


IndexKs2000neighbor = -(Ks2000neighboraverageInfection - K2000neighbornocontrol)/Ks2000neighboraverageControl


IndexKT2000neighbor = -(KT2000neighboraverageInfection - K2000neighbornocontrol)/KT2000neighboraverageControl


######## K = 2000 global


#newKs2000global <- approx(x_axis_KT2000global, K2000globalnocontrol, x_axis_Ks2000global)$y


IndexKs2000global = -(Ks2000globalaverageInfection - K2000globalnocontrol)/Ks2000globalaverageControl


IndexKT2000global = -(KT2000globalaverageInfection - K2000globalnocontrol)/KT2000globalaverageControl


################################################################################################################

tiff("Figure_3.tiff", width=10,height=10, units='in',res=600)

par(mfrow=c(3,3))
par(mar=c(2,3,3,3))

plot(x_axis_Ks500local,IndexKs500Local,type = "l", col = "blue", lwd = 2, ylim = c(0,0.007))
points(x_axis_KT500local,IndexKT500Local,type = "l", col = "red", lwd = 2)

plot(x_axis_Ks500neighbor,IndexKs500neighbor,type = "l", col = "blue", lwd = 2, ylim = c(0,0.007))
points(x_axis_KT500neighbor,IndexKT500neighbor,type = "l", col = "red", lwd = 2)

plot(x_axis_Ks500global,IndexKs500global,type = "l", col = "blue", lwd = 2, ylim = c(0,0.007))
points(x_axis_KT500global,IndexKT500global,type = "l", col ="red", lwd = 2)

plot(x_axis_Ks800local,IndexKs800local,type = "l", col = "blue", lwd = 2, ylim = c(0,0.0165))
points(x_axis_KT800local,IndexKT800local,type = "l", col = "red", lwd = 2)

plot(x_axis_Ks800neighbor,IndexKs800neighbor,type = "l", col = "blue", lwd = 2, ylim = c(0,0.0165))
points(x_axis_KT800neighbor,IndexKT800neighbor,type = "l", col = "red", lwd = 2)

plot(x_axis_Ks800global,IndexKs800global,type = "l", col = "blue", lwd = 2, ylim = c(0,0.0165))
points(x_axis_KT800global,IndexKT800global,type = "l", col = "red", lwd = 2)

plot(x_axis_Ks2000local,IndexKs2000local,type = "l", col = "blue", lwd = 2, ylim = c(0,0.6))
points(x_axis_KT2000local,IndexKT2000local,type = "l", col = "red", lwd = 2)

plot(x_axis_Ks2000neighbor,IndexKs2000neighbor,type = "l", col = "blue", lwd = 2, ylim = c(0,0.6))
points(x_axis_KT2000neighbor,IndexKT2000neighbor,type = "l", col = "red", lwd = 2)


plot(x_axis_Ks2000global,IndexKs2000global,type = "l", col = "blue", lwd = 2, ylim = c(0,0.6))
points(x_axis_KT2000global,IndexKT2000global,type = "l", col = "red", lwd = 2)


dev.off()