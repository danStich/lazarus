###########################################################################################################
#   LAZARUS CODE FOR JAGS (INFORMED MODEL)
#
#   GENERALIZED FRAMEWORK FOR ESTIMATING SURVIVAL OF FISH
#   WHEN OBSERVATION OF MORTALITY CONTAINS FALSE-POSITIVE
#   IDENTIFICATION OF DEAD FISH
# 
#   CODE MAKES EXTENSIVE USE OF CODE AND FORMATTING
#   CONVENTIONS FROM KERY AND SCHAUB (2012) BAYESIAN POPULATION ANALYSIS
#
#   CITATION FOR THIS WORK:
#   Stich, D. S. , Y. Jiao, and B. R. Murphy. 2015. 
#     Life, death, and resurrection: Accounting for state
#     uncertainty in survival estimation from radio-tagged
#     Grass Carp. North American Journal of Fisheries Management 35:321-330.
#
################################################################################
# Install and load necessary packages and software
# If you do not have jags installed on your PC, the following link will open
# a download window to install.  Select yes, then follow instructions.  You
# must uncomment to run it
# browseURL("http://sourceforge.net/projects/mcmc-jags/files/latest/download?source=files")

# Install the JAGS packages and load them
  #install.packages("R2jags")
  require(R2jags)
	
#psiABt <- c()
#psiR <- c()
#psiD <- c()
#Delta <- c()
bias <- c()

for(i in 1:1000){
  
# Define mean survival, transitions, recapture, as well as number of occasions,
# states, observations and released individuals
phiA <- 1.00   # Survival of live fish
phiB <- 1.00   # Survival of dead fish
psiAB <- runif(1, 0.05, 0.50)  # Mortality rate
psiBA <- 0.00  # Resurrection rate
pA <- runif(1, 0.50, 0.90)     # Detection in live state  
pA2 <- runif(1, 0.50, 0.90)     # Detection live fish in dead state  
pB <- runif(1, 0.50, 0.90)     # Detection in dead state 
delt <- runif(1, 0.50, 0.95)   # Correct classification rate
n.occasions <- 13
n.states <- 3
n.obs <- 3

# Create matrix of marking
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- c(100,rep(0, (n.occasions-1)))  
marked[,2] <- rep(0, n.occasions)
marked[,3] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
   # Dimension 1: state of departure
   # Dimension 2: state of arrival
   # Dimension 3: individual
   # Dimension 4: time

# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.STATE[,,i,t] <- matrix(c(
      phiA*(1-psiAB), phiA*psiAB,     1-phiA,
      phiB*psiBA,     phiB*(1-psiBA), 1-phiB,
      0,              0,              1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
   for (t in 1:(n.occasions-1)){
      PSI.OBS[,,i,t] <- matrix(c(
      pA*(delt), pA*(1-delt),  1-pA,
      0,  pB, 1-pB,
      0,  0,  1       ), nrow = n.states, byrow = TRUE)
      } #t
   } #i

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
   # Unobservable: number of state that is unobservable
   n.occasions <- dim(PSI.STATE)[4] + 1
   CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
   # Define a vector with the occasion of marking
   mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
   g <- colSums(marked)
   for (s in 1:dim(PSI.STATE)[1]){
      if (g[s]==0) next  # To avoid error message if nothing to replace
      mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
      } #s
   for (i in 1:sum(marked)){
      for (s in 1:dim(PSI.STATE)[1]){
         if (mark.occ[i,s]==0) next
         first <- mark.occ[i,s]
         CH[i,first] <- s
         CH.TRUE[i,first] <- s
         } #s
      for (t in (first+1):n.occasions){
         # Multinomial trials for state transitions
         if (first==n.occasions) next
         state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
         CH.TRUE[i,t] <- state
         # Multinomial trials for observation process
         event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
         CH[i,t] <- event
         } #t
      } #i
   # Replace the NA and the highest state number (dead) in the file by 0
   CH[is.na(CH)] <- 0
   CH[CH==dim(PSI.STATE)[1]] <- 0
   CH[CH==0] <- 3
   # id <- numeric(0)
   # for (i in 1:dim(CH)[1]){
      # z <- min(which(CH[i,]!=0))
      # ifelse(z==dim(CH)[2], id <- c(id,i), id <- c(id))
      # }
   return(CH)
   # CH: capture histories to be used
   # CH.TRUE: capture histories with perfect observation
   }

# Execute function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim
f <- c(rep(1,nrow(CH)))

# Define number of occasions in CH as number of columns
  n.occasions <- ncol(CH)

# Define number of states, including one for dead (absorbing state: 3)
  n.states <- max(CH)

# Define number of observable states (one for unobserved: 3)
  n.obs <- max(CH)

# Define number of individuals (number of rows in encounter history)
  nind <-nrow(CH)

# Set a new working directory to separate the files for different models
  if(!file.exists('C:/Laz0')) dir.create('C:/Laz0')
  setwd('C:/Laz0')	
	
# Specify the model in JAGS language
  sink("ms.jags") # save as a model file
  cat(" # concatenate the model string
    model{ # Start model specification

      # PARAMETERS DEFINITIONS:
        # phiAlive: survival of individuals in live state
        # phiDead: survival of individuals with tags transmitting mortality
        # psiDeath: probability of moving from live to possibly dead
        # psiResurrection: probability of resurrection
        # pAlive: probability of detection in live state
        # pDead: probability of detection in possibly dead state
        
      # STATES:
        # 1 Live
        # 2 Possibly dead
        # 3 Really dead
   
      # OBSERVABLE STATES:
        # 1 Seen alive
        # 2 Seen with tag transmitting mortality
        # 3 Not seen      

      # PRIORS AND CONSTRAINTS:
        # Parameter definitions for state-transition probabilities
          psiDeath <- mean.psi[1]
          psiResurrection <- mean.psi[2]
        # Parameter definitions for probability of detection
          pAlive <- mean.p[1]
          pAlive2 ~ dunif(0,1)
          pDead <- mean.p[2]
        # Priors for state-transition and detection probabilities  
          for (u in 1:2){
            mean.psi[u] ~ dunif(0,1)
            mean.p[u] ~ dunif(0,1)
          } 

        # Prior for misclassification rate
          delta ~ dunif(0,1)

          mu <- psiDeath*pAlive*delta-((1-pAlive2*delta)-(1-pAlive2))

      # Four-dimensional state-transition and detection probability matrices
      # conditional on individual i being in state h at time t
        for(i in 1:nind){ # For each individual
          for(t in 1:(n.occasions-1)){ # At during each interval
            # Define probability of state h(t+1) given h(t) [OMEGA MATRIX]
              ps[1,i,t,1] <- (1-psiDeath)
              ps[1,i,t,2] <- psiDeath
              ps[1,i,t,3] <- 0
              ps[2,i,t,1] <- 0.0000000000000000000000000000000000001
              ps[2,i,t,2] <- 1-0.0000000000000000000000000000000000001
              ps[2,i,t,3] <- 0
              ps[3,i,t,1] <- 0
              ps[3,i,t,2] <- 0 
              ps[3,i,t,3] <- 1
            # Define probability of detection given h(t) [THETA MATRIX]
              po[1,i,t,1] <- pAlive*delta
              po[1,i,t,2] <- pAlive2*(1-delta)
              po[1,i,t,3] <- 1-(pAlive)
              po[2,i,t,1] <- 0
              po[2,i,t,2] <- pDead
              po[2,i,t,3] <- 1-pDead
              po[3,i,t,1] <- 0
              po[3,i,t,2] <- 0 
              po[3,i,t,3] <- 1
          } # t
        } # i    
           
      # STATE-SPACE LIKELIHOOD:
        for(i in 1:nind){ # For each individual
          # Define state at first capture
            z[i,f[i]] <- y[i,f[i]]
            for (t in (f[i]+1):n.occasions){ # At each occasion
              # STATE PROCESS: draw S(t) given S(t-1)
                z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, ]) 
              # OBSERVATION PROCESS: draw O(t) given S(t)
                y[i,t] ~ dcat(po[z[i,t], i, t-1, ])
            } # t
        } # i
    }
  ",fill = TRUE)
  sink()
# End model specification

# Function to create known latent states z to speed computation
  known.state.ms <- function(ms, notseen){
   # notseen: label for not seen
   state <- ms
   state[state==notseen] <- NA # Change unobserved to NA, need to estimate
   for (i in 1:dim(ms)[1]){ 
      m <- min(which(!is.na(state[i,]))) # Fill in known states
      state[i,m] <- NA 
   }
     return(state)
   }

# Function to create initial values for unknown z
  ms.init.z <- function(ch, f){
    for (i in 1:dim(ch)[1]){
      ch[i,1:f[i]] <- NA
    }
    states <- max(ch, na.rm = TRUE)
    known.states <- 1:(states-1)
    v <- which(ch==states)
    ch[-v] <- NA
    ch[v] <- sample(known.states, length(v), replace = TRUE)
    return(ch)
  }

# Bundle data
  jags.data <- list(y = CH, f = f, n.occasions = dim(CH)[2], 
    nind = dim(CH)[1], z = known.state.ms(CH, 3))

# Initial values
  inits <- function(){list(
  	delta=runif(1,0,1),
    mean.psi=runif(2,0,1),
    mean.p=runif(2,0,1), 
  	pAlive2=runif(1,0,1),
    z=ms.init.z(CH,f))
  }  

# Parameters monitored
  parameters <- c("psiDeath", "delta", "mu",
   "pAlive", "pAlive2", "pDead")

# MCMC settings
  ni <- 600
  nt <- 2
  nb <- 200
  nc <- 3

# Call JAGS from R and run the model
   msJAGS <- jags(jags.data, inits, parameters, model.file="ms.jags",
     n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

  bias <- append(bias, mean(msJAGS$BUGSoutput$sims.list$mu)-psiAB)
  hist(bias, main='', xlab='Bias', col='gray87')
  mtext(side=3, paste(length(bias)), adj=1, cex=2)
  mtext(side=3, sprintf('%.6f', mean(bias)), adj=0, cex=2)
}