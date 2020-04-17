#######################################################################################
##                                                                                   ##
## Simulation code - data from the National Center for Health Statistics. Percentile ##
## Data Files with LMS Values. Centres for Disease Control & Prevention; 2009.       ##
## From https://www.cdc.gov/growthcharts/percentile_data_files.htm & from 1st author ##
##                                                                                   ##
#######################################################################################

# Require packages
require(foreign); require(MASS); require(rpsychi); require(reshape2); require(ggplot2); require(nlme)
require(R2MLwiN); require(lmtest); require(psych); require(cocor)

# Set directory to find MLwiN software - ensure this is compatible on your computer
options(MLwiN_path = "C:/Program Files/MLwiN v3.03/")

######################
## Ascertain values ##
######################

wtage1      <- read.csv("wtage0-36.csv")            # Read in 0 to 36 months data
wtage2      <- wtage1[ c(14,32,52,70), ]            # Restrict to required ages only
wtage2$mean <- wtage2$M; print(wtage2$mean)         # Create mean = median variable
wtage2$SD   <- wtage2$M*wtage2$S; print(wtage2$SD)  # Create sd variable where sd = M*S
wtage3      <- read.csv("wtage2-20.csv")            # Read in 2-20 years data
wtage4      <- wtage3[ c(8,26,226,244), ]           # Restrict to required ages only
wtage4$mean <- wtage4$M; print(wtage4$mean)         # Create mean = median variable
wtage4$SD   <- wtage4$M*wtage4$S; print(wtage4$SD)  # Create sd variable where sd = M*S

#################
## Simulations ##
#################

# Set seed
set.seed(56079)

# Set number of iterations
iterations <- 100000

# Create vector for p values for ide present and absent for MLM & Oldham's
p.value.ide.mlm          <- vector(mode="numeric", length=iterations)
p.value.MLM.z.test.ide   <- vector(mode="numeric", length=iterations)
p.value.noide.mlm        <- vector(mode="numeric", length=iterations)
p.value.MLM.z.test.noide <- vector(mode="numeric", length=iterations)
p.value.ide.old.c        <- vector(mode="numeric", length=iterations)
p.value.ide.old.i        <- vector(mode="numeric", length=iterations)
p.value.z.test.ide       <- vector(mode="numeric", length=iterations)
p.value.noide.old.c      <- vector(mode="numeric", length=iterations)
p.value.noide.old.i      <- vector(mode="numeric", length=iterations)
p.value.z.test.noide     <- vector(mode="numeric", length=iterations)

# Create vectors for MLM & Oldham's correlations for both scenarios
MLM.cor.ide.c   <- vector(mode="numeric", length=iterations)
MLM.cor.ide.i   <- vector(mode="numeric", length=iterations)
old.cor.ide.c   <- vector(mode="numeric", length=iterations)
old.cor.ide.i   <- vector(mode="numeric", length=iterations)
MLM.cor.noide.c <- vector(mode="numeric", length=iterations)
MLM.cor.noide.i <- vector(mode="numeric", length=iterations)
old.cor.noide.c <- vector(mode="numeric", length=iterations)
old.cor.noide.i <- vector(mode="numeric", length=iterations)

# Create for loop for ide present
for (i in 1:iterations) {
  # Create id for 1000 children
  id.i <- seq(1, 1000, by=1)
  
  # Create group variable
  group.i <- rep(seq(0,1), each=500)
  
  # Set autocorrelation
  rho.i <- 0.9
  rho12.i <- rho.i^1.5
  rho23.i <- rho.i^1.5
  rho13.i <- rho.i^3
  
  # Set multipliers
  MultC.i <- (15.88-9.67)/3    # Girls' data - approx increase 2.07kg per year
  MultI.i <- (14.17-9.67)/3    # Intervention causing only increase 1.5kg per year
  
  # Create mu vectors
  MuC.i <- c(9.67,9.67+MultC.i*1.5,9.67+MultC.i*3)
  MuI.i <- c(9.67,9.67+MultI.i*1.5,9.67+MultI.i*3)
  
  # Create covariance matrix from correlations and sds
  SigC.i <- r2cov(sd=c(1.0,1.5,2.0),R=lower2R(c(rho12.i,rho13.i,rho23.i)))    # SDs from girls' data
  SigI.i <- r2cov(sd=c(1.0,1.375,1.75),R=lower2R(c(rho12.i,rho13.i,rho23.i))) # SDs decreased proportionally to mean
  
  # Create multivariate normal weight data
  WtC.i <- mvrnorm(500,MuC.i,SigC.i)
  WtI.i <- mvrnorm(500,MuI.i,SigI.i)
  
  # Create control and intervention dataframes with id and group
  c.data.i <- data.frame(id.i[1:500],    group.i[1:500],    WtC.i, stringsAsFactors = FALSE)
  i.data.i <- data.frame(id.i[501:1000], group.i[501:1000], WtI.i, stringsAsFactors = FALSE)
  
  # Rename id variable in both dataframes
  names(c.data.i)[names(c.data.i)=="id.i.1.500."]    <- "id.i"
  names(i.data.i)[names(i.data.i)=="id.i.501.1000."] <- "id.i"
  
  # Rename group variable in both dataframes
  names(c.data.i)[names(c.data.i)=="group.i.1.500."]    <- "group.i"
  names(i.data.i)[names(i.data.i)=="group.i.501.1000."] <- "group.i"
  
  # Create stacked dataframe
  full.data.i <- rbind(c.data.i, i.data.i)
  
  # Rename weight variables in stacked dataframe
  names(full.data.i)[names(full.data.i)=="X1"] <- "wt1.i"
  names(full.data.i)[names(full.data.i)=="X2"] <- "wt2.5.i"
  names(full.data.i)[names(full.data.i)=="X3"] <- "wt4.i"
  
  # Reshape into long form with time centered
  long.data.i <- reshape(full.data.i, varying = c("wt1.i","wt2.5.i","wt4.i"), v.names = "weight.i", 
                         timevar = "year.i", idvar="id.i", times = c("-1.5","0","1.5"), direction = "long")
  
  # Convert character year variable to numeric year variable
  long.data.i$year.i <- as.numeric(long.data.i$year.i)
  
  # Create group0 and group1 variables
  long.data.i$group0.i <- 1-long.data.i$group.i
  long.data.i$group1.i <- long.data.i$group.i
  
  # Create yeargroup0 and yeargroup1 variables
  long.data.i$yeargroup0.i <- long.data.i$year.i*long.data.i$group0.i
  long.data.i$yeargroup1.i <- long.data.i$year.i*long.data.i$group1.i
  
  # Sort dataset by id then year
  ide.data.i <- long.data.i[order(long.data.i$id.i, long.data.i$year.i), ]
  
  ######################
  ## Multilevel Model ##
  ######################
  
  # Create covariance matrix for constraining covariances
  covmatrix.i      <- matrix(NA, nrow=3, ncol=5)
  covmatrix.i[1,1] <- 2
  covmatrix.i[2,1] <- "group0.i"
  covmatrix.i[3,1] <- "group1.i"
  covmatrix.i[1,2] <- 2
  covmatrix.i[2,2] <- "group1.i"
  covmatrix.i[3,2] <- "yeargroup0.i"
  covmatrix.i[1,3] <- 2
  covmatrix.i[2,3] <- "group0.i"
  covmatrix.i[3,3] <- "yeargroup1.i"
  covmatrix.i[1,4] <- 2
  covmatrix.i[2,4] <- "yeargroup0.i"
  covmatrix.i[3,4] <- "yeargroup1.i"
  covmatrix.i[1,5] <- 1
  covmatrix.i[2,5] <- "group0.i"
  covmatrix.i[3,5] <- "group1.i"
  
  # Run model where covariances allowed to differ (unconstrained)
  moduncon.i <- runMLwiN(weight.i ~ group0.i + group1.i + yeargroup0.i + yeargroup1.i + 
                           (group0.i + group1.i + yeargroup0.i + yeargroup1.i | id.i) + 
                           (group0.i + group1.i | year.i), data=ide.data.i, 
                         estoptions=list(EstM=0, Meth=0, reset=c(2,2), clre=covmatrix.i))
  
  # Calculate correlations for each group separately from covariances
  MLM.ide.c <- as.vector((moduncon.i@RP["RP2_cov_group0_i_yeargroup0_i"])/
                           ((sqrt(moduncon.i@RP["RP2_var_group0_i"])*(sqrt(moduncon.i@RP["RP2_var_yeargroup0_i"])))))
  MLM.ide.i <- as.vector((moduncon.i@RP["RP2_cov_group1_i_yeargroup1_i"])/
                           ((sqrt(moduncon.i@RP["RP2_var_group1_i"])*(sqrt(moduncon.i@RP["RP2_var_yeargroup1_i"])))))
  
  # Run Fisher's z test on MLM correlations
  MLM.z.test.ide            <- cocor.indep.groups(MLM.ide.c, MLM.ide.i, 500, 500, 
                                                  alternative="two.sided", test="fisher1925", return.htest=TRUE)
  p.value.MLM.z.test.ide[i] <- signif(MLM.z.test.ide$fisher1925$p.value,4)
  
  # Store MLM correlations
  MLM.cor.ide.c[i] <- MLM.ide.c
  MLM.cor.ide.i[i] <- MLM.ide.i
  
  # Run nested model with constrained covariances
  modcon.i <- runMLwiN(weight.i ~ group0.i + group1.i + yeargroup0.i + yeargroup1.i + 
                         (group0.i + group1.i + yeargroup0.i + yeargroup1.i | id.i) + 
                         (group0.i + group1.i | year.i), data=ide.data.i, 
                       estoptions=list(EstM=0, Meth=0, reset=c(2,2), clre=covmatrix.i, 
                                       constraints=list(random.ui=matrix(c(0,0,1,0,-1,0,0,0), nrow=8), random.ci=c(0))))
  
  # Extract LR test value and calculate p value
  lr.value.i         <- as.numeric(-2*logLik(modcon.i) - -2*logLik(moduncon.i))
  p.value.ide.mlm[i] <- signif((1-pchisq(abs(lr.value.i),1)),2)
  
  #####################
  ## Oldham's Method ##
  #####################
  
  # Calculate correlation between difference and average for control group
  diff.ide.c     <- as.vector(c.data.i$X1 - c.data.i$X3)
  ave.ide.c      <- as.vector((c.data.i$X1 + c.data.i$X3)/2)
  old.ide.c.test <- cor.test(diff.ide.c, ave.ide.c, alternative = "two.sided", method = "pearson")
  old.ide.c      <- cor(diff.ide.c, ave.ide.c)
  old.ide.c.z    <- fisherz(old.ide.c)
  
  # Calculate correlation between difference and average for intervention group
  diff.ide.i     <- as.vector(i.data.i$X1 - i.data.i$X3)
  ave.ide.i      <- as.vector((i.data.i$X1 + i.data.i$X3)/2)
  old.ide.i.test <- cor.test(diff.ide.i, ave.ide.i, alternative = "two.sided", method = "pearson")
  old.ide.i      <- cor(diff.ide.i, ave.ide.i)
  old.ide.i.z    <- fisherz(old.ide.i)
  
  # Store Oldham's correlations
  old.cor.ide.c[i] <- old.ide.c
  old.cor.ide.i[i] <- old.ide.i
  
  # Extract p values from each correlation
  p.value.ide.old.c[i] <- signif(old.ide.c.test$p.value,4)
  p.value.ide.old.i[i] <- signif(old.ide.i.test$p.value,4)
  
  # Use Fisher's z test to assess difference between 2 Oldham's correlations
  z.test.ide            <- cocor.indep.groups(old.ide.c, old.ide.i, 500, 500, 
                                              alternative="two.sided", test="fisher1925", return.htest=TRUE)
  p.value.z.test.ide[i] <- signif(z.test.ide$fisher1925$p.value,4)
}

# Create for loop for ide absent
for (i in 1:iterations) {
  # Create id for 1000 children
  id.n <- seq(1, 1000, by=1)
  
  # Create group variable
  group.n <- rep(seq(0,1), each=500)
  
  # Set autocorrelation
  rho.n   <- 0.9
  rho12.n <- rho.n^1.5
  rho23.n <- rho.n^1.5
  rho13.n <- rho.n^3
  
  # Set multipliers
  MultC.n <- (15.88-9.67)/3   # Girls' data - approx increase 2.07kg per year
  MultI.n <- (14.17-9.67)/3   # no IDE causing increase 1.5kg per year for everyone
  
  # Create mu vectors
  MuC.n <- c(9.67,9.67+MultC.n*1.5,9.67+MultC.n*3)
  MuI.n <- c(9.67,9.67+MultI.n*1.5,9.67+MultI.n*3)
  
  #Create covariance matrix from correlations and sds
  SigC.n <- r2cov(sd=c(1.0,1.5,2.0),R=lower2R(c(rho12.n,rho13.n,rho23.n)))  # SDs from girls' data
  SigI.n <- r2cov(sd=c(1.0,1.5,2.0),R=lower2R(c(rho12.n,rho13.n,rho23.n)))  # SDs same as no IDE
  
  # Create multivariate normal weight data
  WtC.n <- mvrnorm(500,MuC.n,SigC.n)
  WtI.n <- mvrnorm(500,MuI.n,SigI.n)
  
  # Create control and intervention dataframes with id and group
  c.data.n <- data.frame(id.n[1:500],    group.n[1:500],    WtC.n, stringsAsFactors = FALSE)
  i.data.n <- data.frame(id.n[501:1000], group.n[501:1000], WtI.n, stringsAsFactors = FALSE)
  
  # Rename id variable in both dataframes
  names(c.data.n)[names(c.data.n)=="id.n.1.500."]    <- "id.n"
  names(i.data.n)[names(i.data.n)=="id.n.501.1000."] <- "id.n"
  
  # Rename group variable in both dataframes
  names(c.data.n)[names(c.data.n)=="group.n.1.500."]    <- "group.n"
  names(i.data.n)[names(i.data.n)=="group.n.501.1000."] <- "group.n"
  
  # Create stacked dataframe
  full.data.n <- rbind(c.data.n, i.data.n)
  
  # Rename weight variables in stacked dataframe
  names(full.data.n)[names(full.data.n)=="X1"] <- "wt1.n"
  names(full.data.n)[names(full.data.n)=="X2"] <- "wt2.5.n"
  names(full.data.n)[names(full.data.n)=="X3"] <- "wt4.n"
  
  # Reshape into long form with time centered
  long.data.n <- reshape(full.data.n, varying = c("wt1.n","wt2.5.n","wt4.n"), v.names = "weight.n", 
                         timevar = "year.n", idvar="id.n", times = c("-1.5","0","1.5"), direction = "long")
  
  # Convert character year variable to numeric year variable
  long.data.n$year.n <- as.numeric(long.data.n$year.n)
  
  # Create group0 and group1 variables
  long.data.n$group0.n <- 1-long.data.n$group.n
  long.data.n$group1.n <- long.data.n$group.n
  
  # Create yeargroup0 and yeargroup1 variables
  long.data.n$yeargroup0.n <- long.data.n$year.n*long.data.n$group0.n
  long.data.n$yeargroup1.n <- long.data.n$year.n*long.data.n$group1.n
  
  # Sort dataset by id then year
  noide.data.n <- long.data.n[order(long.data.n$id.n, long.data.n$year.n), ]
  
  ######################
  ## Multilevel Model ##
  ######################
  
  # Create covariance matrix for constraining covariances
  covmatrix.n      <- matrix(NA, nrow=3, ncol=5)
  covmatrix.n[1,1] <- 2
  covmatrix.n[2,1] <- "group0.n"
  covmatrix.n[3,1] <- "group1.n"
  covmatrix.n[1,2] <- 2
  covmatrix.n[2,2] <- "group1.n"
  covmatrix.n[3,2] <- "yeargroup0.n"
  covmatrix.n[1,3] <- 2
  covmatrix.n[2,3] <- "group0.n"
  covmatrix.n[3,3] <- "yeargroup1.n"
  covmatrix.n[1,4] <- 2
  covmatrix.n[2,4] <- "yeargroup0.n"
  covmatrix.n[3,4] <- "yeargroup1.n"
  covmatrix.n[1,5] <- 1
  covmatrix.n[2,5] <- "group0.n"
  covmatrix.n[3,5] <- "group1.n"
  
  # Run model where covariances allowed to differ (unconstrained)
  moduncon.n <- runMLwiN(weight.n ~ group0.n + group1.n + yeargroup0.n + yeargroup1.n + 
                           (group0.n + group1.n + yeargroup0.n + yeargroup1.n | id.n) + 
                           (group0.n + group1.n | year.n), data=noide.data.n, 
                         estoptions=list(EstM=0, Meth=0, reset=c(2,2), clre=covmatrix.n))
  
  # Calculate correlations for each group separately from covariances
  MLM.noide.c <- as.vector((moduncon.n@RP["RP2_cov_group0_n_yeargroup0_n"])/
                             ((sqrt(moduncon.n@RP["RP2_var_group0_n"])*(sqrt(moduncon.n@RP["RP2_var_yeargroup0_n"])))))
  MLM.noide.i <- as.vector((moduncon.n@RP["RP2_cov_group1_n_yeargroup1_n"])/
                             ((sqrt(moduncon.n@RP["RP2_var_group1_n"])*(sqrt(moduncon.n@RP["RP2_var_yeargroup1_n"])))))
  
  # Run Fisher's z test on MLM correlations
  MLM.z.test.noide            <- cocor.indep.groups(MLM.noide.c, MLM.noide.i, 500, 500, 
                                                    alternative="two.sided", test="fisher1925", return.htest=TRUE)
  p.value.MLM.z.test.noide[i] <- signif(MLM.z.test.noide$fisher1925$p.value,4)
  
  # Store MLM correlations
  MLM.cor.noide.c[i] <- MLM.noide.c
  MLM.cor.noide.i[i] <- MLM.noide.i
  
  # Run nested model with constrained covariances
  modcon.n <- runMLwiN(weight.n ~ group0.n + group1.n + yeargroup0.n + yeargroup1.n + 
                         (group0.n + group1.n + yeargroup0.n + yeargroup1.n | id.n) 
                       + (group0.n + group1.n | year.n), data=noide.data.n, 
                       estoptions=list(EstM=0, Meth=0, reset=c(2,2), clre=covmatrix.n, 
                                       constraints=list(random.ui=matrix(c(0,0,1,0,-1,0,0,0), nrow=8), random.ci=c(0))))
  
  # Extract LR test value and calculate p value
  lr.value.n           <- as.numeric(-2*logLik(modcon.n) - -2*logLik(moduncon.n))
  p.value.noide.mlm[i] <- signif((1-pchisq(abs(lr.value.n),1)),2)
  
  
  #####################
  ## Oldham's Method ##
  #####################
  
  # Calculate correlation between difference and average for control group
  diff.noide.c     <- as.vector(c.data.n$X1 - c.data.n$X3)
  ave.noide.c      <- as.vector((c.data.n$X1 + c.data.n$X3)/2)
  old.noide.c.test <- cor.test(diff.noide.c, ave.noide.c, alternative = "two.sided", method = "pearson")
  old.noide.c      <- cor(diff.noide.c, ave.noide.c)
  old.noide.c.z    <- fisherz(old.noide.c)
  
  # Calculate correlation between difference and average for intervention group
  diff.noide.i     <- as.vector(i.data.n$X1 - i.data.n$X3)
  ave.noide.i      <- as.vector((i.data.n$X1 + i.data.n$X3)/2)
  old.noide.i.test <- cor.test(diff.noide.i, ave.noide.i, alternative = "two.sided", method = "pearson")
  old.noide.i      <- cor(diff.noide.i, ave.noide.i)
  old.noide.i.z    <- fisherz(old.noide.i)
  
  # Store Oldham's correlations
  old.cor.noide.c[i] <- old.noide.c
  old.cor.noide.i[i] <- old.noide.i
  
  # Extract p values
  p.value.noide.old.c[i] <- signif(old.noide.c.test$p.value,4)
  p.value.noide.old.i[i] <- signif(old.noide.i.test$p.value,4)
  
  # Use Fisher's z test to assess difference between 2 Oldham's correlations
  z.test.noide            <- cocor.indep.groups(old.noide.c, old.noide.i, 500, 500, 
                                                alternative="two.sided", test="fisher1925", return.htest=TRUE)
  p.value.z.test.noide[i] <- signif(z.test.noide$fisher1925$p.value,4)
}

# Calculate power and error rates
power.mlm.lr             <- ((sum(p.value.ide.mlm          <= 0.05))/iterations)*100
power.mlm.z              <- ((sum(p.value.MLM.z.test.ide   <= 0.05))/iterations)*100
t1error.mlm.lr           <- ((sum(p.value.noide.mlm        <= 0.05))/iterations)*100
t1error.mlm.z            <- ((sum(p.value.MLM.z.test.noide <= 0.05))/iterations)*100
t1error.old.ide.c        <- ((sum(p.value.ide.old.c        <= 0.05))/iterations)*100
power.old.ide.i          <- ((sum(p.value.ide.old.i        <= 0.05))/iterations)*100
t1error.old.noide.c      <- ((sum(p.value.noide.old.c      <= 0.05))/iterations)*100
t1error.old.noide.i      <- ((sum(p.value.noide.old.i      <= 0.05))/iterations)*100
power.old.ide.2group     <- ((sum(p.value.z.test.ide       <= 0.05))/iterations)*100
t1error.old.noide.2group <- ((sum(p.value.z.test.noide     <= 0.05))/iterations)*100

# Display results: power & error rates
power.mlm.lr 
power.mlm.z
t1error.mlm.lr 
t1error.mlm.z
t1error.old.ide.c 
power.old.ide.i 
t1error.old.noide.c 
t1error.old.noide.i 
power.old.ide.2group 
t1error.old.noide.2group 

# Display results: correlations
MLM.cor.ide.c
MLM.cor.ide.i
old.cor.ide.c
old.cor.ide.i
MLM.cor.noide.c
MLM.cor.noide.i
old.cor.noide.c
old.cor.noide.i