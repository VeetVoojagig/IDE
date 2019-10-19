##########################################################################################################
## Simulations of adult weight management programmes - truncated and non-truncated - evaluated for IDEs ##
##########################################################################################################

rm(list=ls())
require(MASS); require(rpsychi); require(psych); require(R2MLwiN); require(foreign); library(lavaan); require(ggplot2)

##########################
## Setup all parameters ##
##########################

Nrep    <- 10 # Number or repeated simulations - ultimately set this to 10k
Nsmp    <- 500   # Sample size of typical weight-management study
Npop    <- 5000  # Population sample size that is simulated before subsampling to obtain typical study size
Thold   <- 30    # Selection threshold adopted (obese adults)
seed    <- 13; set.seed(seed)
version <- paste0("-Adult-Nrep-",Nrep,"-01.csv")
CovCon  <- t(matrix(c(as.character(c(rep(2,4))),"Grp0","Grp1","Grp0","GrpT0","Grp1",
                      "GrpT0","GrpT1","GrpT1"),ncol=3))                       # Cross-group covariance constraints
rand.c  <- as.matrix(rep(0,9),ncol=,1); rand.c[3,] <- 1; rand.c[5,] <- -1     # IDE covariance constraint

## Population & Subsample Control (no truncation & no intervention) ##
C0.mean1 <- 26.4   # Mean = Median BMI (males 25-34, health survey England)
C0.mean2 <- 26.4   # No change in mean between periods as no RTM or intervention
C0.sd    <- 5.0    # derived from 0.25*sqrt(8011) from health survey for England; same between periods
C0.r     <- 0.85   # Correlation between t1 and t2 
C0.Corr  <- matrix(c(1, C0.r, C0.r, 1), nrow=2)    # Correlation matrix
C0.Sigma <- r2cov(sd=c(C0.sd,C0.sd), R=C0.Corr)    # Covariance  matrix

## Population Null (no truncation & null IDE) ##
T0.mean1 <- 26.4  # Initial mean = MedianBMI (males 25-34, health survey England)
T0.mean2 <- 22.0  # Reduction of 4.4 BMI due to intervention
T0.sd    <- 5.0   # Same SD for both measurements as no RTM or IDE
T0.r     <- 0.85  # Correlation between t1 and t2
T0.Corr  <- matrix(c(1, T0.r, T0.r, 1), nrow=2)    # Correlation matrix
T0.Sigma <- r2cov(sd=c(T0.sd,T0.sd), R=T0.Corr)    # Covariance  matrix

## Population IDE (truncation & non-zero IDE) ##
T1.mean1 <- 26.4  # Initial mean = MedianBMI (males 25-34, health survey England)
T1.mean2 <- 22.0  # Reduction of 4.4 BMI due to intervention
T1.sd1   <- 5.0   # Initial SD
T1.sd2   <- 4.0   # Lower SD by 1 BMI at t2 due to IDE
T1.r     <- 0.85  # Correlation between t1 and t2
T1.Corr  <- matrix(c(1, T1.r, T1.r, 1), nrow=2)    # Correlation matrix
T1.Sigma <- r2cov(sd=c(T1.sd1,T1.sd2), R=T1.Corr)  # Covariance  matrix

#########################
## Setup all Functions ##
#########################

# Function to generate correlated t1 and t2 measures for control or treatment group
SimDat <- function(N,Mu,Sig,Emp=FALSE) {
  data <- data.frame(mvrnorm(n=N, mu=Mu, Sigma=Sig, empirical=Emp)); names(data) <- c("T1","T2")
  data <- data.frame(data, Diff = data$T2 - data$T1, Summ = data$T2 + data$T1)
  return(data) }

# Function to generate long data for single group MLwiN analysis
MyDataLong <- function(data) {
  MyData   <- data.frame(data[,1:2]); names(MyData) <- c("T1","T2")
  MyData.long <- reshape(MyData, varying=c("T1","T2"), v.names="BMI", timevar="Time", 
                         idvar="Id", times=as.numeric(c("-0.5","0.5")), direction="long")
  MyData.long <- MyData.long[order(MyData.long$Id,MyData.long$Time),] 
  return (MyData.long) }

# Function to generate long data for combined group MLwiN analysis
MyDataLong2 <- function(data1,data2) {
  n1 <- dim(data1)[1]; n2 <- dim(data2)[1]
  MyData <- rbind(data.frame(Grp0=rep(0,n1),data1[,1:2]),data.frame(Grp0=rep(1,n2),data2[,1:2]))
  names(MyData)[2:3] <- c("T1","T2")
  MyData.long <- reshape(MyData, varying=c("T1","T2"), v.names="BMI", timevar="Time", idvar="Id", 
                         times=as.numeric(c("-0.5","0.5")), direction="long")
  MyData.long$Grp1  <- 1 - MyData.long$Grp0                                   # Create Grp1  variable
  MyData.long$GrpT0 <- MyData.long$Time*MyData.long$Grp0                      # Create GrpT0 variable
  MyData.long$GrpT1 <- MyData.long$Time*MyData.long$Grp1                      # Create GrpT1 variable
  MyData.long <- MyData.long[order(MyData.long$Id,MyData.long$Time),] 
  return (MyData.long) }

# Function to undertake MC on single group
EstMCg  <- function(MyData) {
  Rho   <- cor.test(MyData[,3],MyData[,1])
  MCg   <- data.frame(MCr=Rho$estimate,MCp=Rho$p.value)
  return (MCg) }

# Function to undertake Oldham on single group
EstOLg  <- function(MyData) {
  Rho   <- cor.test(MyData[,3],MyData[,4])
  OLg   <- data.frame(OLr=Rho$estimate,OLp=Rho$p.value)
  return (OLg) }

# Function to undertake MLM on single group
EstMLg  <- function(MyData.long) {
  n     <- dim(MyData.long)[1]/2
  ML    <- runMLwiN(BMI ~ 1 + Time + (1 + Time| Id) + (1 | Time), data=MyData.long, 
                    estoptions = list(EstM = 0, Meth=0, resi.store = TRUE))
  Rho   <- ML@RP[2]/sqrt(ML@RP[1]*ML@RP[3])
  MLg   <- data.frame(MLr=Rho,MLp=r.test(n,Rho)$p)
  return (MLg) }

# Function to undertake MLM on two groups
EstML2  <- function(MyData.long) {
  n     <- dim(MyData.long)[1]/2
  # Unconstrained model
  MLu   <- runMLwiN(BMI ~ Grp0 + Grp1 + GrpT0 + GrpT1 + (Grp0 + Grp1 + GrpT0 + GrpT1 | Id) + (Grp0 + Grp1 | Time), 
                    data=MyData.long, estoptions=list(clre=CovCon))
  # Constrained model
  MLc   <- runMLwiN(BMI ~ Grp0 + Grp1 + GrpT0 + GrpT1 + (Grp0 + Grp1 + GrpT0 + GrpT1 | Id) + (Grp0 + Grp1 | Time), 
                    data=MyData.long,estoptions=list(clre=CovCon, constraints=list(random.ui=rand.c,random.ci=c(0))))
  Rho1  <- MLu@RP[3]/sqrt(MLu@RP[1]*MLu@RP[4])
  Rho2  <- MLu@RP[5]/sqrt(MLu@RP[2]*MLu@RP[6])
  ML2   <- data.frame(MLr=Rho1-Rho2,MLp=pchisq(abs(abs(2*as.numeric(logLik(MLu))-2*as.numeric(logLik(MLc)))),1,lower.tail=FALSE))
  return (ML2) }

# Function to contrast two independent correlations
EstCor  <- function(Rho1,Rho2,n) {
  t1    <- SIN::fisherz(Rho1)
  t2    <- SIN::fisherz(Rho2)
  Evalr <- Rho2 - Rho1
  tDiff <- (t1-t2)/(sqrt((1/(n-3))+(1/(n-3))))
  pDiff <- 2*(1-pnorm(abs(tDiff),0,1))
  Eval  <- data.frame(Evalr,pDiff)
  return (Eval) }

# Function to expedite Mills Ratio
Mills <- function(x) {
  LogMill <- pnorm(x, lower.tail=FALSE, log.p=TRUE) - dnorm(x, log=TRUE) 
  return(exp(LogMill)) } 

# Function to estimate baseline SD for truncated sample
tSDx <- function(t,SDx) {
  R <- Mills(t)
  xVar <- 1 + (t/R) - (1/(R*R))
  xSD  <- sqrt(xVar)
  return(xSD) }

# Function to estimate follow-up SD for truncated sample
tSDy <- function(t,SDx,Rho) {
  R <- Mills(t)
  yVar <- 1 + Rho*Rho*((t/R) - (1/(R*R)))
  ySD  <- sqrt(yVar)
  return(ySD) }

##############################################################
##############################################################
##                                                          ##
## SIMULATION - generate population & sample truncated data ## 
##                                                          ##
##############################################################
##############################################################

## Establish arrays for datasets to be Simulated ##
PopCH0 <- PopTH0 <- PopTH1 <- array(dim = c(Nrep,Npop,4)) # Population Control / Intervention No IDE / Intervention Non-Zero IDE
SmpCH0 <- SmpTH0 <- SmpTH1 <- array(dim = c(Nrep,Nsmp,4)) # Subsample  Control / Intervention No IDE / Intervention Non-Zero IDE
## Populate the arrays
for (itn in 1:Nrep) {
  PopCH0[itn,,] <- as.matrix(SimDat(Npop, c(C0.mean1,C0.mean2), C0.Sigma))
  PopTH0[itn,,] <- as.matrix(SimDat(Npop, c(T0.mean1,T0.mean2), T0.Sigma))
  PopTH1[itn,,] <- as.matrix(SimDat(Npop, c(T1.mean1,T1.mean2), T1.Sigma))
  # Subsample selection according to threshold & subsample size (Nsmp)
  idCH0 <- sample(sort(unique((1:Npop)*(PopCH0[itn,,1]>Thold)))[-1],Nsmp); SmpCH0[itn,,] <- PopCH0[itn,idCH0,]
  idTH0 <- sample(sort(unique((1:Npop)*(PopTH0[itn,,1]>Thold)))[-1],Nsmp); SmpTH0[itn,,] <- PopTH0[itn,idTH0,]
  idTH1 <- sample(sort(unique((1:Npop)*(PopTH1[itn,,1]>Thold)))[-1],Nsmp); SmpTH1[itn,,] <- PopTH1[itn,idTH1,] }

################################
## Examine the simulated data ##
################################

RepMean  <- NULL; for (itn in 1:Nrep) RepMean <- rbind(RepMean,apply(PopCH0[itn,,1:2],2,mean))
pCH0mean <- round(apply(RepMean,2,mean),1)
RepMean  <- NULL; for (itn in 1:Nrep) RepMean <- rbind(RepMean,apply(PopTH0[itn,,1:2],2,mean))
pTH0mean <- round(apply(RepMean,2,mean),1)
RepMean  <- NULL; for (itn in 1:Nrep) RepMean <- rbind(RepMean,apply(PopTH1[itn,,1:2],2,mean))
pTH1mean <- round(apply(RepMean,2,mean),1)
RepMean  <- NULL; for (itn in 1:Nrep) RepMean <- rbind(RepMean,apply(SmpCH0[itn,,1:2],2,mean))
sCH0mean <- round(apply(RepMean,2,mean),1)
RepMean  <- NULL; for (itn in 1:Nrep) RepMean <- rbind(RepMean,apply(SmpTH0[itn,,1:2],2,mean))
sTH0mean <- round(apply(RepMean,2,mean),1)
RepMean  <- NULL; for (itn in 1:Nrep) RepMean <- rbind(RepMean,apply(SmpTH1[itn,,1:2],2,mean))
sTH1mean <- round(apply(RepMean,2,mean),1)
RepSD    <- NULL; for (itn in 1:Nrep) RepSD   <- rbind(RepSD,  apply(PopCH0[itn,,1:2],2,sd))
pCH0sd   <- round(apply(RepSD,  2,mean),1)
RepSD    <- NULL; for (itn in 1:Nrep) RepSD   <- rbind(RepSD,  apply(PopTH0[itn,,1:2],2,sd))
pTH0sd   <- round(apply(RepSD,  2,mean),1)
RepSD    <- NULL; for (itn in 1:Nrep) RepSD   <- rbind(RepSD,  apply(PopTH1[itn,,1:2],2,sd))
pTH1sd   <- round(apply(RepSD,  2,mean),1)
RepSD    <- NULL; for (itn in 1:Nrep) RepSD   <- rbind(RepSD,  apply(SmpCH0[itn,,1:2],2,sd))
sCH0sd   <- round(apply(RepSD,  2,mean),1)
RepSD    <- NULL; for (itn in 1:Nrep) RepSD   <- rbind(RepSD,  apply(SmpTH0[itn,,1:2],2,sd))
sTH0sd   <- round(apply(RepSD,  2,mean),1)
RepSD    <- NULL; for (itn in 1:Nrep) RepSD   <- rbind(RepSD,  apply(SmpTH1[itn,,1:2],2,sd))
sTH1sd   <- round(apply(RepSD,  2,mean),1)
RepRho   <- NULL; for (itn in 1:Nrep) RepRho  <- rbind(RepRho, cor(PopCH0[itn,,1],PopCH0[itn,,2]))
pCH0r    <- round(apply(RepRho, 2,mean),2)
RepRho   <- NULL; for (itn in 1:Nrep) RepRho  <- rbind(RepRho, cor(PopTH0[itn,,1],PopTH0[itn,,2]))
pTH0r    <- round(apply(RepRho, 2,mean),2)
RepRho   <- NULL; for (itn in 1:Nrep) RepRho  <- rbind(RepRho, cor(PopTH1[itn,,1],PopTH1[itn,,2]))
pTH1r    <- round(apply(RepRho, 2,mean),2)
RepRho   <- NULL; for (itn in 1:Nrep) RepRho  <- rbind(RepRho, cor(SmpCH0[itn,,1],SmpCH0[itn,,2]))
sCH0r    <- round(apply(RepRho, 2,mean),2)
RepRho   <- NULL; for (itn in 1:Nrep) RepRho  <- rbind(RepRho, cor(SmpTH0[itn,,1],SmpTH0[itn,,2]))
sTH0r    <- round(apply(RepRho, 2,mean),2)
RepRho   <- NULL; for (itn in 1:Nrep) RepRho  <- rbind(RepRho, cor(SmpTH1[itn,,1],SmpTH1[itn,,2]))
sTH1r    <- round(apply(RepRho, 2,mean),2)
RepMin   <- NULL; for (itn in 1:Nrep) RepMin  <- rbind(RepMin, apply(PopCH0[itn,,1:2],2,min))
pCH0min  <- round(apply(RepMin, 2,mean),1)
RepMin   <- NULL; for (itn in 1:Nrep) RepMin  <- rbind(RepMin, apply(PopTH0[itn,,1:2],2,min))
pTH0min  <- round(apply(RepMin, 2,mean),1)
RepMin   <- NULL; for (itn in 1:Nrep) RepMin  <- rbind(RepMin, apply(PopTH1[itn,,1:2],2,min))
pTH1min  <- round(apply(RepMin, 2,mean),1)
RepMin   <- NULL; for (itn in 1:Nrep) RepMin  <- rbind(RepMin, apply(SmpCH0[itn,,1:2],2,min))
sCH0min  <- round(apply(RepMin, 2,mean),1)
RepMin   <- NULL; for (itn in 1:Nrep) RepMin  <- rbind(RepMin, apply(SmpTH0[itn,,1:2],2,min))
sTH0min  <- round(apply(RepMin, 2,mean),1)
RepMin   <- NULL; for (itn in 1:Nrep) RepMin  <- rbind(RepMin, apply(SmpTH1[itn,,1:2],2,min))
sTH1min  <- round(apply(RepMin, 2,mean),1)
filename <- paste0("Simulation",version)
write.table("Baseline mean BMI (population / sample)",              filename,row.names=FALSE,col.names=FALSE)
write.table(t(c(pCH0mean[1],sCH0mean[1])),sep=",",                  filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table("Follow-up mean BMI (population / sample) - no IDE",    filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(t(c(pTH0mean[2],sTH0mean[2])),sep=",",                  filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table("Follow-up mean BMI (population / sample) - true IDE",  filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(t(c(pTH1mean[2],sTH1mean[2])),sep=",",                  filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table("Baseline BMI SD (population / sample)",                filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(t(c(pCH0sd[1],sCH0sd[1])),sep=",",                      filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table("Follow-up BMI SD (population / sample) - control",     filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(t(c(pCH0sd[2],sCH0sd[2])),sep=",",                      filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table("Follow-up BMI SD (population / sample) - no IDE",      filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(t(c(pTH0sd[2],sTH0sd[2])),sep=",",                      filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table("Follow-up BMI SD (population / sample) - true IDE",    filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(t(c(pTH1sd[2],sTH1sd[2])),sep=",",                      filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table("Correlation (population / sample)",                    filename,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(t(c(pTH0r,sTH0r)),sep=",",                              filename,row.names=FALSE,col.names=FALSE,append=TRUE)


#################################################
#################################################
##                                             ##
## Execute analyses & generate summary outputs ##
##                                             ##
#################################################
#################################################

######################################
## Naive MC analyses - single group ##
######################################

PC0MCg <- PT0MCg <- PT1MCg <- SC0MCg <- ST0MCg <- ST1MCg <- NULL
for (itn in 1:Nrep) {
  PC0MCg <- rbind(PC0MCg,EstMCg(PopCH0[itn,,]))
  PT0MCg <- rbind(PT0MCg,EstMCg(PopTH0[itn,,]))
  PT1MCg <- rbind(PT1MCg,EstMCg(PopTH1[itn,,]))
  SC0MCg <- rbind(SC0MCg,EstMCg(SmpCH0[itn,,]))
  ST0MCg <- rbind(ST0MCg,EstMCg(SmpTH0[itn,,]))
  ST1MCg <- rbind(ST1MCg,EstMCg(SmpTH1[itn,,])) }
EvalMC   <- data.frame(PC0MCg,PT0MCg,PT1MCg,SC0MCg,ST0MCg,ST1MCg)
names(EvalMC) <- c("PC0MCr","PC0MCp","PT0MCr","PT0MCp","PT1MCr","PT1MCp","SC0MCr","SC0MCp","ST0MCr","ST0MCp","ST1MCr","ST1MCp")
filename <- paste0("EvalMC",version); write.csv(EvalMC,filename,row.names=FALSE)

####################################
## Oldham analyses - single group ##
####################################

PC0OLg <- PT0OLg <- PT1OLg <- SC0OLg <- ST0OLg <- ST1OLg <- NULL
for (itn in 1:Nrep) {
  PC0OLg <- rbind(PC0OLg,EstOLg(PopCH0[itn,,]))
  PT0OLg <- rbind(PT0OLg,EstOLg(PopTH0[itn,,]))
  PT1OLg <- rbind(PT1OLg,EstOLg(PopTH1[itn,,]))
  SC0OLg <- rbind(SC0OLg,EstOLg(SmpCH0[itn,,]))
  ST0OLg <- rbind(ST0OLg,EstOLg(SmpTH0[itn,,]))
  ST1OLg <- rbind(ST1OLg,EstOLg(SmpTH1[itn,,])) }
EvalOL   <- data.frame(PC0OLg,PT0OLg,PT1OLg,SC0OLg,ST0OLg,ST1OLg)
names(EvalOL) <- c("PC0OLr","PC0OLp","PT0OLr","PT0OLp","PT1OLr","PT1OLp","SC0OLr","SC0OLp","ST0OLr","ST0OLp","ST1OLr","ST1OLp")
filename <- paste0("EvalOL",version); write.csv(EvalOL,filename,row.names=FALSE)

########################################
## Multilevel analyses - single group ##
########################################

PC0MLg <- NULL
for (itn in 1:Nrep) { PC0MLg <- rbind(PC0MLg,EstMLg(MyDataLong(PopCH0[itn,,]))) }
filename <- paste0("PC0MLg",version); write.csv(PC0MLg,filename,row.names=FALSE)
PT0MLg   <- NULL
for (itn in 1:Nrep) { PT0MLg <- rbind(PT0MLg,EstMLg(MyDataLong(PopTH0[itn,,]))) }
filename <- paste0("PT0MLg",version); write.csv(PT0MLg,filename,row.names=FALSE)
PT1MLg   <- NULL
for (itn in 1:Nrep) { PT1MLg <- rbind(PT1MLg,EstMLg(MyDataLong(PopTH1[itn,,]))) }
filename <- paste0("PT1MLg",version); write.csv(PT1MLg,filename,row.names=FALSE)
SC0MLg   <- NULL
for (itn in 1:Nrep) { SC0MLg <- rbind(SC0MLg,EstMLg(MyDataLong(SmpCH0[itn,,]))) }
filename <- paste0("SC0MLg",version); write.csv(SC0MLg,filename,row.names=FALSE)
ST0MLg   <- NULL
for (itn in 1:Nrep) { ST0MLg <- rbind(ST0MLg,EstMLg(MyDataLong(SmpTH0[itn,,]))) }
filename <- paste0("ST0MLg",version); write.csv(ST0MLg,filename,row.names=FALSE)
ST1MLg   <- NULL
for (itn in 1:Nrep) { ST1MLg <- rbind(ST1MLg,EstMLg(MyDataLong(SmpTH1[itn,,]))) }
filename <- paste0("ST1MLg",version); write.csv(ST1MLg,filename,row.names=FALSE)

##########################################################################
## Multilevel analyses - combined control & intervention group analysis ##
##########################################################################

PCT0ML2  <- NULL
for (itn in 1:Nrep) { PCT0ML2 <- rbind(PCT0ML2,EstML2(MyDataLong2(PopCH0[itn,,],PopTH0[itn,,]))) }
filename <- paste0("PCT0ML2",version); write.csv(PCT0ML2,filename,row.names=FALSE)
PCT1ML2  <- NULL
for (itn in 1:Nrep) { PCT1ML2 <- rbind(PCT1ML2,EstML2(MyDataLong2(PopCH0[itn,,],PopTH1[itn,,]))) }
filename <- paste0("PCT1ML2",version); write.csv(PCT1ML2,filename,row.names=FALSE)
SCT0ML2  <- NULL
for (itn in 1:Nrep) { SCT0ML2 <- rbind(SCT0ML2,EstML2(MyDataLong2(SmpCH0[itn,,],SmpTH0[itn,,]))) }
filename <- paste0("SCT0ML2",version); write.csv(SCT0ML2,filename,row.names=FALSE)
SCT1ML2  <- NULL
for (itn in 1:Nrep) { SCT1ML2 <- rbind(SCT1ML2,EstML2(MyDataLong2(SmpCH0[itn,,],SmpTH1[itn,,]))) }
filename <- paste0("SCT1ML2",version); write.csv(SCT1ML2,filename,row.names=FALSE)

#########################################################################
## Collate previously evaluated single group multilevel analyses - ML1 ##
#########################################################################

filename <- paste0("PC0MLg",version); PC0MLg <- read.csv(filename)
filename <- paste0("PT0MLg",version); PT0MLg <- read.csv(filename)
filename <- paste0("PT1MLg",version); PT1MLg <- read.csv(filename)
filename <- paste0("SC0MLg",version); SC0MLg <- read.csv(filename)
filename <- paste0("ST0MLg",version); ST0MLg <- read.csv(filename)
filename <- paste0("ST1MLg",version); ST1MLg <- read.csv(filename)
EvalML1  <- data.frame(PC0MLg,PT0MLg,PT1MLg,SC0MLg,ST0MLg,ST1MLg)
names(EvalML1) <- c("PC0MLr","PC0MLp","PT0MLr","PT0MLp","PT1MLr","PT1MLp","SC0MLr","SC0MLp","ST0MLr","ST0MLp","ST1MLr","ST1MLp")
filename       <- paste0("EvalML1",version); write.csv(EvalML1,filename,row.names=FALSE)

###########################################################################
## Collate previously evaluated combined group multilevel analyses - ML2 ##
###########################################################################

filename <- paste0("PCT0ML2",version); PCT0ML2 <- read.csv(filename)
filename <- paste0("PCT1ML2",version); PCT1ML2 <- read.csv(filename)
filename <- paste0("SCT0ML2",version); SCT0ML2 <- read.csv(filename)
filename <- paste0("SCT1ML2",version); SCT1ML2 <- read.csv(filename)
EvalML2  <- data.frame(PCT0ML2,PCT1ML2,SCT0ML2,SCT1ML2)
names(EvalML2) <- c("PCT0MLr","PCT0MLp","PCT1MLr","PCT1MLp","SCT0MLr","SCT0MLp","SCT1MLr","SCT1MLp")
filename       <- paste0("EvalML2",version); write.csv(EvalML2,filename,row.names=FALSE)

#################################
# Pairwise contrasts of groups ##
#################################

filename <- paste0("EvalMC",version);  EvalMC  <- read.csv(filename)
filename <- paste0("EvalOL",version);  EvalOL  <- read.csv(filename)
filename <- paste0("EvalML1",version); EvalML  <- read.csv(filename)
PCT0MCc  <- PCT1MCc <- SCT0MCc <- SCT1MCc <- NULL
PCT0OLc  <- PCT1OLc <- SCT0OLc <- SCT1OLc <- NULL
PCT0MLc  <- PCT1MLc <- SCT0MLc <- SCT1MLc <- NULL 
for (itn in 1:Nrep) {
  PCT0MCc <- rbind(PCT0MCc,EstCor(EvalMC[itn,1],EvalMC[itn,3], Npop))
  PCT1MCc <- rbind(PCT1MCc,EstCor(EvalMC[itn,1],EvalMC[itn,5], Npop))
  SCT0MCc <- rbind(SCT0MCc,EstCor(EvalMC[itn,7],EvalMC[itn,9], Nsmp))
  SCT1MCc <- rbind(SCT1MCc,EstCor(EvalMC[itn,7],EvalMC[itn,11],Nsmp))
  PCT0OLc <- rbind(PCT0OLc,EstCor(EvalOL[itn,1],EvalOL[itn,3], Npop))
  PCT1OLc <- rbind(PCT1OLc,EstCor(EvalOL[itn,1],EvalOL[itn,5], Npop))
  SCT0OLc <- rbind(SCT0OLc,EstCor(EvalOL[itn,7],EvalOL[itn,9], Nsmp))
  SCT1OLc <- rbind(SCT1OLc,EstCor(EvalOL[itn,7],EvalOL[itn,11],Nsmp))
  PCT0MLc <- rbind(PCT0MLc,EstCor(EvalML[itn,1],EvalML[itn,3], Npop))
  PCT1MLc <- rbind(PCT1MLc,EstCor(EvalML[itn,1],EvalML[itn,5], Npop))
  SCT0MLc <- rbind(SCT0MLc,EstCor(EvalML[itn,7],EvalML[itn,9], Nsmp))
  SCT1MLc <- rbind(SCT1MLc,EstCor(EvalML[itn,7],EvalML[itn,11],Nsmp)) }
EvalPZ    <- data.frame(PCT0MCc,PCT1MCc,SCT0MCc,SCT1MCc,PCT0OLc,PCT1OLc,SCT0OLc,SCT1OLc,PCT0MLc,PCT1MLc,SCT0MLc,SCT1MLc)
names(EvalPZ) <- c("PCT0MCr","PCT0MCp","PCT1MCr","PCT1MCp","SCT0MCr","SCT0MCp","SCT1MCr","SCT1MCp","PCT0OLr","PCT0OLp",
                   "PCT1OLr","PCT1OLp","SCT0OLr","SCT0OLp","SCT1OLr","SCT1OLp","PCT0MLr","PCT0MLp","PCT1MLr","PCT1MLp",
                   "SCT0MLr","SCT0MLp","SCT1MLr","SCT1MLp")
filename <- paste0("EvalPZ",version); write.csv(EvalPZ,filename,row.names=FALSE)

#########################
## Calculate summaries ##
#########################

filename <- paste0("EvalMC",version);  EvalMC <- read.csv(filename)
filename <- paste0("EvalOL",version);  EvalOL <- read.csv(filename)
filename <- paste0("EvalML1",version); EvalML <- read.csv(filename)
filename <- paste0("EvalPZ",version);  EvalPZ <- read.csv(filename)
filename <- paste0("EvalML2",version); EvalMM <- read.csv(filename)

# Summarise single group evaluations of correlation & its p-value
Cval   <- grep("r",names(EvalMC)) # Match correlations
MCr    <- apply(EvalMC[,Cval],2,function(x){quantile(x,c(0.025,0.5,0.975))})
OLr    <- apply(EvalOL[,Cval],2,function(x){quantile(x,c(0.025,0.5,0.975))})
MLr    <- apply(EvalML[,Cval],2,function(x){quantile(x,c(0.025,0.5,0.975))})
Cor.g  <- rbind(MCr,OLr,MLr)       
rownames(Cor.g) <- paste(c(rep("MCr",3),rep("OLr",3),rep("MLr",3)),rownames(MCr))
colnames(Cor.g) <- substr(colnames(Cor.g),1,3)
Pval   <- grep("p",names(EvalMC)) # Match p-values
MCp    <- 100*apply(EvalMC[,Pval],2,function(x){sum(x<=0.05)})/Nrep
OLp    <- 100*apply(EvalOL[,Pval],2,function(x){sum(x<=0.05)})/Nrep
MLp    <- 100*apply(EvalML[,Pval],2,function(x){sum(x<=0.05)})/Nrep
Pval.g <- rbind(MCp,OLp,MLp)       
colnames(Pval.g) <- colnames(Cor.g)

# Summarise pairwise group difference in correlation & corresponding p-value or MLM LRT p-value
Cval   <- grep("r",names(EvalPZ)) # Match correlations
PZr    <- apply(EvalPZ[,Cval],2,function(x){quantile(x,c(0.025,0.5,0.975))})
Cval   <- grep("r",names(EvalMM)) # Match correlations
MMr    <- apply(EvalMM[,Cval],2,function(x){quantile(x,c(0.025,0.5,0.975))})
Cor.p  <- rbind(PZr[,1:4],PZr[,5:8],PZr[,9:12],MMr)
rownames(Cor.p) <- paste(c(rep("MCr",3),rep("OLr",3),rep("MLr",3),rep("MMr",3)),rownames(PZr))
colnames(Cor.p) <- substr(colnames(Cor.p),1,4)
Pval   <- grep("p",names(EvalPZ)) # Match p-values
PZp    <- 100*apply(EvalPZ[,Pval],2,function(x){sum(x<=0.05)})/Nrep
Pval   <- grep("p",names(EvalMM)) # Match p-values
MMp    <- 100*apply(EvalMM[,Pval],2,function(x){sum(x<=0.05)})/Nrep
Pval.p <- rbind(PZp[1:4],PZp[5:8],PZp[9:12],MMp)
rownames(Pval.p) <- c(rownames(Pval.g),"MMp")
colnames(Pval.p) <- substr(colnames(Pval.p),1,4)

filename <- paste0("Summary",version)
write.table("Single group correlations",filename,row.names=FALSE,col.names=FALSE)
write.table(t(c("Est",colnames(Cor.g))),filename,sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(Cor.g,filename,sep=",",append=TRUE,row.names=TRUE,col.names=FALSE)
write.table("",filename,append=TRUE,row.names=FALSE,col.names=FALSE)
write.table("Single group p-values",filename,append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(t(c("Est",colnames(Pval.g))),filename,sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(Pval.g,filename,sep=",",append=TRUE,row.names=TRUE,col.names=FALSE)
write.table("",filename,append=TRUE,row.names=FALSE,col.names=FALSE)

write.table("Pairwise group correlation differences for each method",filename,append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(t(c("Est",colnames(Cor.p))),filename,sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(Cor.p,filename,sep=",",append=TRUE,row.names=TRUE,col.names=FALSE)
write.table("",filename,append=TRUE,row.names=FALSE,col.names=FALSE)
write.table("Pairwise p-values",filename,append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(t(c("Est",colnames(Pval.p))),filename,sep=",",append=TRUE,row.names=FALSE,col.names=FALSE)
write.table(Pval.p,filename,sep=",",append=TRUE,row.names=TRUE,col.names=FALSE)
write.table("",filename,append=TRUE,row.names=FALSE,col.names=FALSE) 

#########################################
## Plot the consequences of truncation ##
#########################################

Summary <- NULL
for (Rho in seq(0.05,0.95,0.1)) {
  Corr    <- matrix(c(1, Rho, Rho, 1), nrow=2)  # Correlation matrix
  Sigma   <- r2cov(sd=c(1,1), R=Corr)         # Covariance  matrix
  Pop     <- data.frame(SimDat(Npop,c(1,1),Sigma,Emp=TRUE))[,1:2]; names(Pop) <- c("PopX","PopY")
  PopRho  <- cor(Pop$PopX,Pop$PopY)
  for (t in seq(-2,3,0.05)) {
    Smp     <- Pop[Pop$PopX>t,]; names(Smp) <- c("SmpX","SmpY")
    EstSDx  <- tSDx(t,SD)
    EstSDy  <- tSDy(t,SD,Rho)
    EstInd  <- EstSDy/EstSDx
    Summary <- rbind(Summary,data.frame(t=t,Rho=round(Rho,2),EstSDx=EstSDx,EstSDy=EstSDy,EstInd=EstInd)) } }
Summary$r <- as.factor(Summary$Rho)

ggplot(Summary) + geom_smooth(aes(x=t,y=EstInd,color=r),method=loess,se=FALSE) + 
  geom_hline(yintercept=1, linetype="solid", color = "black", size=1) +
  scale_x_continuous(name="Z-score threshold of population truncation at baseline") +
  scale_y_continuous(name="Ratio of follow-up to baseline standard deviation") +
  theme(axis.title.x = element_text(face="bold", colour="black", size=18),
        axis.text.x  = element_text(vjust=0.5, size=16),
        axis.title.y = element_text(face="bold", colour="black", size=18),
        axis.text.y  = element_text(vjust=0.5, size=16),
        legend.title = element_text(colour="black", size=18, face="italic"),
        legend.text = element_text(colour="black", size = 16))

round(100*(Summary[(Summary$t==0)&(Summary$Rho==0.95),]$EstInd-1),1)
Thold <- 30     # Selection threshold adopted (obese adults)
Mu    <- 26.4   # Mean = Median BMI (males 25-34, health survey England)
SD    <- 5.0    # derived from 0.25*sqrt(8011) from health survey for England; same between periods
(T <- (Thold-Mu)/SD)
Rho <- 0.95
EstSDx  <- tSDx(T,SD)
EstSDy  <- tSDy(T,SD,Rho)
EstInd  <- EstSDy/EstSDx
round(100*(EstInd-1),1)