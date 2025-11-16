### Cartel Recruitment
### Age Cohort analysis
### 2025/09/22
### Code created by Rafael Prieto-Curiel

#######################################################################
#### Cartel recruitment
#######################################################################
# The code is divided into three sections
# S1 - Time series analysis of cartels
# S2 - Age assignment of agents
# S3 - Analysis of the model

# S1 - Time series
# The code creates a time series that corresponds to weekly steps of the model
# Five time series are created: Number of active members, recruited, incapacitated, killed and retired
# each time series indicates the number of people of that weekly step

# S2 - Age assignment of agents
# Using the time series produced in S1, an age is assigned to each person
# the age is updated each week for all active members of the cartel

# S3 - Analysis of the model
# 1.4 million males born in 1990
# get first the whole cartels with all population
# then distribute those according to some peak 

### considerations
# the code is based on functions that can be reused for other cases.
# relevant parameters
# 19300 recruits in a year
# 6500 deaths
# 5700 incapacitations
# on date 44562 (01/01/2022) we get 175000
# with the current parameters we get 2012 and 2022 data from Science

#######################################################################
#### parameters and packages used
#######################################################################
{
#require(scales)
#require(dplyr)
#require(reshape)
P0 <- 23000 #### initial population of males in 30 years
AgeCohort <- 1400000
rho = 370/175000 #net gain of 90 members each week
#### from 2005 to 2022 so 175000 in the end
kappaP = (5700/52)/175000 #### probability of arrest
kappaK = (6500/52)/175000 #### probability of killed 
sat <- (rho-kappaP-kappaK)/600000
pi_rec = 1 ### prob of successful recruitment
### here, all exposures to recruitment are considered so pi_rec = 1
### however, it is possible to alter the parameter
Equilibrium <- (rho-kappaP-kappaK)/sat
}


#######################################################################
#### functions
#######################################################################
{
  #### function that creates five time series (S1)
  #### the function runs for a determined number of steps (weeks)
  CartelTS <- function(P0,     # initial cartel size
                       kappaP, # prob of being incapacitated
                       kappaK, # prob of being killed
                       pi_rec, # prob of recruitment
                       rho,    # recruitment rate
                       t0 = 23743, #01/01/1965 and with 4227 steps we reach end of 2045
                       steps = 4227){ 
    #### Initial members of a cartel  
    Active <- P0
    #### not recruited but exposed
    NotRecruited <- 0*Active
    Recruited <- 0*Active
    Killed <- 0*Active #### killed
    Incapacitated <- 0*Active
    Exposed <- 0
    Saturated <- 0
    Date <- t0
    
    #### time vars
    At <- Active
    Rt <- Recruited
    Kt <- Killed
    It <- Incapacitated
    Nt <- NotRecruited
    Et <- Exposed
    St <- Saturated
    Dt <- Date
    
    for(k in 1:(steps+1)){
      Date <- Date + 7 # one week later
      
      # Exposed <- round(Active* rho * ((P0 - sum(Kt) - sum(It) -sum(Nt) - Active)/P0)) for finite pop
      Exposed <- round(Active* rho) # for infinite pop
      Recruited <- IsRecruited(N = Exposed,
                               pi_rec = pi_rec)
      NotRecruited <- Exposed - Recruited ##### they were exposed and were not recruited
      Active <- Active + Recruited ##### new members of a cartel
      
      #### incapacitations
      Incapacitated <- sum(runif(Active) < kappaP)
      Active <- Active - Incapacitated
      
      #### kills
      Killed <- sum(runif(Active) < kappaK)
      Active <- Active - Killed
      
      #### saturation or retirement
      Saturated <- round(sat*Active^2)
      Active <- Active - Saturated
      
      At <- c(At, Active)
      Kt <- c(Kt, Killed)
      It <- c(It, Incapacitated)
      Rt <- c(Rt, Recruited)
      Nt <- c(Nt, NotRecruited)
      Et <- c(Et, Exposed)
      St <- c(St, Saturated)
      Dt <- c(Dt, Date)
    }
    DF <- data.frame(Active = At, 
                     Killed = Kt, 
                     Incap = It, 
                     Recruit = Rt, 
                     NotRect = Nt, 
                     Exposed = Et, 
                     Sat = St,
                     Dt = Dt)
    return(DF)
  }
  
  #### function that takes N, the size and pi_rec, prob of being recruited
  ### returns the number of people who are recruited
  IsRecruited <- function(N,    # number of people who are exposed to recruitment 
                          pi_rec){
    #### a function that takes all values and returns a 1 if the person is recruited
    return(sum(runif(N) < pi_rec))
  }
  
  #### function that takes the time series and assigns a specific age for each agent
  DemogImpact <- function(Efs){
    n <- dim(Efs)[1]
    Gen <- paste("G", 1930:2035, sep = "-")
    
    ### create initial age
    ActAge <- CreateAgeRecruit(Efs$Active[1], date = Efs$Dt[1])
    
    #### to get the share by the week
    {
      counts <- data.frame(matrix(
        ncol = length(Gen) + 1,  # +1 for Dt
        nrow = 0                 # start empty
      ))
      names(counts) <- c("Dt", Gen)
    }
    
    #### create initial DFs
    KilledDF <- data.frame(Age = c(), Birth = c(), Gen = c(), dateK = c())
    IncapDF <- data.frame(Age = c(), Birth = c(), Gen = c(), dateK = c())
    RecruitDF <- data.frame(Age = c(), Birth = c(), Gen = c(), dateK = c())
    SatDF <- data.frame(Age = c(), Birth = c(), Gen = c(), dateK = c())
    
    #### each step
    for(k in 2:n){
      #cat(k, "\t -- size -- ", dim(ActAge)[1], "\n")
      ActAge$Age <- ActAge$Age + 7/365 # age them one week 
      AgeRecruits <- CreateAgeRecruit(Efs$Recruit[k], date = Efs$Dt[k])
      RecruitDF <- rbind(RecruitDF, data.frame(Age = AgeRecruits$Age,
                                               Birth = AgeRecruits$Birth,
                                               Gen = AgeRecruits$Gen, 
                                               dateK = Efs$Dt[k]))
      ActAge <- rbind(ActAge, AgeRecruits)
      
      ### update shares
      #### to get the share by the week
      {
        # frequency counts for this step
        tab <- table(ActAge$Gen)
        
        # start with a zero row
        row <- as.data.frame(t(c(Dt = Efs$Dt[k], setNames(rep(0, length(Gen)), Gen))),
                             stringsAsFactors = FALSE)
        
        # fill in the counts
        row[1, names(tab)] <- as.integer(tab)
        
        # append to accumulator
        counts <- rbind(counts, row)
      }
      
      ### sample Killed
      u <- sample.int(dim(ActAge)[1], size = Efs$Killed[k])
      KilledDF <- rbind(KilledDF, data.frame(Age = ActAge$Age[u],
                                             Birth = ActAge$Birth[u],
                                             Gen = ActAge$Gen[u], 
                                             dateK = Efs$Dt[k]))
      ActAge <- ActAge[-u,]
      
      ### sample Incapacitated
      u <- sample.int(dim(ActAge)[1], size = Efs$Incap[k])
      IncapDF <- rbind(IncapDF, data.frame(Age = ActAge$Age[u],
                                           Birth = ActAge$Birth[u],
                                           Gen = ActAge$Gen[u], 
                                           dateK = Efs$Dt[k]))
      ActAge <- ActAge[-u,]
      
      ### sample Saturated
      u <- sample.int(dim(ActAge)[1], size = Efs$Sat[k], prob = ActAge$Age^15)
      SatDF <- rbind(SatDF, data.frame(Age = ActAge$Age[u],
                                       Birth = ActAge$Birth[u],
                                       Gen = ActAge$Gen[u], 
                                       dateK = Efs$Dt[k]))
      ActAge <- ActAge[-u,]
    }
    
    return(list(ActAge = ActAge, 
                KilledDF = KilledDF, 
                IncapDF = IncapDF, 
                RecruitDF = RecruitDF, 
                SatDF = SatDF,
                AgeShare = counts))
  }
  
  #### for the sensitivity analysis, vary the MaxAge of recruitment by deltaM
  DemogImpactVariations <- function(Efs, deltaM){
    n <- dim(Efs)[1]
    Gen <- paste("G", 1930:2035, sep = "-")
    
    ### create initial age
    ActAge <- CreateAgeRecruitVariations(Efs$Active[1], 
                                         date = Efs$Dt[1],
                                         deltaM = deltaM)
    
    #### create initial DFs
    RecruitDF <- data.frame(Age = c(), Birth = c(), Gen = c(), dateK = c())
    
    #### each step
    for(k in 2:n){
      #cat(k, "\t -- size -- ", dim(ActAge)[1], "\n")
      ActAge$Age <- ActAge$Age + 7/365 # age them one week 
      AgeRecruits <- CreateAgeRecruitVariations(Efs$Recruit[k], 
                                                date = Efs$Dt[k],
                                                deltaM = deltaM)
      RecruitDF <- rbind(RecruitDF, data.frame(Age = AgeRecruits$Age,
                                               Birth = AgeRecruits$Birth,
                                               Gen = AgeRecruits$Gen, 
                                               dateK = Efs$Dt[k]))
    }
    return(sum(RecruitDF$Gen == "G-1990"))
  }
  
  #### function that creates a distribution for the age of people
  CreateAgeRecruit <- function(N,           # number of people
                               MinAge = 15, # min age for recruitment
                               MaxAge = 35, # max age for recruitment
                               date){
    x <- 0:2500
    y <- (x-MinAge*52)*(MaxAge*52-x)
    y <- pmax(y, 0)
    Age <- sample(x, size = N, replace = T, prob = y)/52
    Birth <- date - Age*365
    Gen <- paste("G", 1965+floor((Birth-23743)/365), sep = "-")
    return(data.frame(Age = Age, Birth = Birth, Gen = Gen))
  }
  
  #### for the sensitivity analysis, vary MaxAge = MaxAge + deltaM
  CreateAgeRecruitVariations <- function(N, 
                                         MinAge = 15,
                                         MaxAge = 35,
                                         date, deltaM){
    ### simulate the same process but varying the max age pm 10 years
    MaxAge <- 35 + deltaM
    x <- 0:2500
    y <- (x-MinAge*52)*(MaxAge*52-x)
    y <- pmax(y, 0)
    Age <- sample(x, size = N, replace = T, prob = y)/52
    Birth <- date - Age*365
    Gen <- paste("G", 1965+floor((Birth-23743)/365), sep = "-")
    return(data.frame(Age = Age, Birth = Birth, Gen = Gen))
  }

  #### takes the sentences and the number of people incapacitated
  #### for each case, assigns the length of a sentence
  #### then exposes people to risks and returns the outcome
  Recidiv <- function(IncapDF,  ### data frame with incapacitations
                      Efs,      ### time series
                      Sentence, ### sentences assigned to people
                      kappaP,   ### arrest
                      kappaK,   ### killed
                      sat){ ### sat * P0 = p of retirement
    f <- IncapDF$Gen == "G-1990"
    IncapDF <- IncapDF[f,]
    IncapDF$SENTENCE <- sample(Sentence$SENTENCE,
                               size = sum(f), 
                               prob = 1/Sentence$SENTENCE, ### the weight for shorter sentences
                               replace = T)
    IncapDF$ReleasedAge <- IncapDF$Age + IncapDF$SENTENCE
    IncapDF$AfterPrison <- NA
    IncapDF$AfterPrison[IncapDF$ReleasedAge > 69] <- "DEATH IN PRISON"
    u <- which(is.na(IncapDF$AfterPrison ))
    for(k in 1:length(u)){
      weeks <- (69-IncapDF$ReleasedAge[u[k]])*52+1 ### weeksLeft
      counter <- 0
      startI <- max(floor((IncapDF$dateK[u[k]]-38681)/7),0)
      while(weeks > 0){
        ToK <- runif(1) < kappaK*Efs$Active[2000]/Efs$Active[2000+counter+startI] ### person got killed
        ToI <- runif(1) < kappaP*Efs$Active[2000]/Efs$Active[2000+counter+startI] ### person got arrested
        ToS <- runif(1) < sat*Efs$Active[2000+counter+startI]
        
        if(ToS){
          IncapDF$AfterPrison[u[k]] <- "RETIRED"
          weeks <- 0
        }
        
        if(ToI){
          IncapDF$AfterPrison[u[k]] <- "REINCIDENCE"
          weeks <- 0
        }
        
        if(ToK){
          IncapDF$AfterPrison[u[k]] <- "MURDERED"
          weeks <- 0
        }
        
        weeks <- weeks - 1
        #counter <- counter + 1
      }
    }
    
    IncapDF$AfterPrison[is.na(IncapDF$AfterPrison)] <- "NATURAL DEATH"
    return(IncapDF)
  }
  
  #### the analysis of incarceration by multiplying the length of senteces
  #### by a factor Alpha
  AlphaRecidiv <- function(IncapDF, 
                           Alpha = 1, # multiplier of sentences
                           Efs,       # time series 
                           Sentence,  # observed sentences
                           kappaP,    # arrest
                           kappaK,    # killed
                           sat){ ### sat * P0 = p of retirement
    f <- IncapDF$Gen == "G-1990"
    IncapDF <- IncapDF[f,]
    IncapDF$SENTENCE <- Alpha * sample(Sentence$SENTENCE,
                                       size = sum(f), 
                                       prob = 1/Sentence$SENTENCE, ### the weight for shorter sentences
                                       replace = T)
    IncapDF$ReleasedAge <- IncapDF$Age + IncapDF$SENTENCE
    IncapDF$AfterPrison <- NA
    IncapDF$AfterPrison[IncapDF$ReleasedAge > 69] <- "DEATH IN PRISON"
    u <- which(is.na(IncapDF$AfterPrison ))
    for(k in 1:length(u)){
      weeks <- (69-IncapDF$ReleasedAge[u[k]])*52+1 ### weeksLeft
      counter <- 0
      startI <- max(floor((IncapDF$dateK[u[k]]-38681)/7),0)
      while(weeks > 0){
        ToK <- runif(1) < kappaK*Efs$Active[2000]/Efs$Active[2000+counter+startI] ### person got killed
        ToI <- runif(1) < kappaP*Efs$Active[2000]/Efs$Active[2000+counter+startI] ### person got arrested
        ToS <- runif(1) < sat*Efs$Active[2000+counter+startI]
        
        if(ToS){
          IncapDF$AfterPrison[u[k]] <- "RETIRED"
          weeks <- 0
        }
        
        if(ToI){
          IncapDF$AfterPrison[u[k]] <- "REINCIDENCE"
          weeks <- 0
        }
        
        if(ToK){
          IncapDF$AfterPrison[u[k]] <- "MURDERED"
          weeks <- 0
        }
        
        weeks <- weeks - 1
        #counter <- counter + 1
      }
    }
    
    IncapDF$AfterPrison[is.na(IncapDF$AfterPrison)] <- "NATURAL DEATH"
    return(IncapDF)
  }
  
  #### function that converts dates to a consecutive number
  YearF <- function(serial, origin = "1899-12-30") {
    # Default origin "1899-12-30" works for Windows Excel (1900 system)
    as.integer(format(as.Date(serial, origin = origin), "%Y"))
  }  

  #### function that converts a number to a date
  DateF <- function(year, origin = "1899-12-30") {
    # Convert January 1st of the given year to Excel serial number
    as.integer(as.Date(paste0(year, "-01-01")) - as.Date(origin))
  }
}

#######################################################################
#### Execute analysis (S1 - Time series)
#######################################################################
{
  Efs <- CartelTS(P0 = P0,
                  kappaP = kappaP, 
                  kappaK = kappaK, 
                  pi_rec = pi_rec,
                  rho = rho,
                  t0 = 23743,
                  steps = 4227)
  save(Efs, file = "CartelDynamics.RData")
}

#######################################################################
#### assign members to demographics (S2 - demography)
#######################################################################
{
  #### uses the time series (Efs) from S1 and assigns a specific date to each event
AgeDynamics <- DemogImpact(Efs)
save(AgeDynamics, file = "AgeDynamics.RData")
}

#######################################################################
#### analysis (S3 - analysis of the model)
#######################################################################
#### sensitivity to the age of recruitment
{
  Res <- data.frame(G1990 = c(), deltaM = c())
  for(k in 1:100){
    u <- sample(-5:5, size = 1)
    G1990 <- DemogImpactVariations(Efs, deltaM = u)
    Res <- rbind(Res, data.frame(G1990 = G1990, deltaM = u))
    cat(G1990, "\t", u, "\n")
  }
  save(Res, file = "AgeDynamicsVariations.RData")
}

#### recidivism
{
   Sentence <- read.csv("ENPOL sentences.csv")
   for(k in 1:100){
   RecidAnalysis <- Recidiv(IncapDF, 
                            Efs,### the Active members
                            Sentence,
                            kappaP, ### arrest
                            kappaK, ### killed
                            sat)
   cat(k, "\n")
   
   if(k == 1){
     dtab <- data.frame(table(RecidAnalysis$AfterPrison)/dim(RecidAnalysis)[1])
     df <- data.frame(DeathPrison = dtab$Freq[1],
                      Killed = dtab$Freq[2],
                      Reincidence = dtab$Freq[4],
                      Retire = dtab$Freq[3]+dtab$Freq[5],
                      MeanSentence = mean(RecidAnalysis$SENTENCE))
   } else{
     dtab <- data.frame(table(RecidAnalysis$AfterPrison)/dim(RecidAnalysis)[1])
     df <- rbind(df, data.frame(DeathPrison = dtab$Freq[1],
                      Killed = dtab$Freq[2],
                      Reincidence = dtab$Freq[4],
                      Retire = dtab$Freq[3]+dtab$Freq[5],
                      MeanSentence = mean(RecidAnalysis$SENTENCE)))
   }
   }
   save(df, file = "RecidivismSimulationsWithTime.RData")
}

#### impacts of incarceration
{
  CrimeSavedIncarc <- data.frame(MaxAct = c(),
                                 KilledS = c(),
                                 SaturS = c(),
                                 IncapS = c())
  
  f <- RecruitDF$Gen == "G-1990"
  R1990 <- RecruitDF[f,]
  
  f <- KilledDF$Gen == "G-1990"
  K1990 <- KilledDF[f, ]
  
  f <- IncapDF$Gen == "G-1990"
  I1990 <- IncapDF[f, ]
  
  f <- SatDF$Gen == "G-1990"
  S1990 <- SatDF[f, ]
  
  ## days of activity (no lost)
  R1990$NatDeath <- R1990$Birth+69*365
  MaxAct <- sum(R1990$NatDeath - R1990$dateK)
  
  ## days saved to homicides
  K1990$NatDeath <- K1990$Birth+69*365
  KilledS <- sum(K1990$NatDeath-K1990$dateK)
  
  ## days saved to saturation
  S1990$NatDeath <- S1990$Birth+69*365
  SaturS <- sum(S1990$NatDeath-S1990$dateK)
  
  ### saved to incarceration
  Sentence <- read.csv("ENPOL sentences.csv")
  
  RecidAnalysis <- Recidiv(IncapDF, 
                            Efs,### the Active members
                       Sentence,
                       kappaP, ### arrest
                       kappaK, ### killed
                       sat)
  IncapS <- sum(RecidAnalysis$SENTENCE*365)
  CrimeSavedIncarc <- rbind(CrimeSavedIncarc, 
                            data.frame(MaxAct = MaxAct,
                                 KilledS = KilledS,
                                 SaturS = SaturS,
                                 IncapS = IncapS))

  save(CrimeSavedIncarc, 
       file = "IncarcerationSavedCrime.RData")
  CrimeSavedIncarc/CrimeSavedIncarc$MaxAct
  CrimeSavedIncarc$IncapS/CrimeSavedIncarc$MaxAct*(1 + 0.26) ### the 0.26 is reincapacitation
}

#### alpha recidivism
{
CrimeSavedIncarcAlpha <- data.frame(Alpha = c(),
                                    MaxAct = c(),
                                    KilledS = c(),
                                    SaturS = c(),
                                    IncapS = c())

f <- RecruitDF$Gen == "G-1990"
R1990 <- RecruitDF[f,]

f <- KilledDF$Gen == "G-1990"
K1990 <- KilledDF[f, ]

f <- IncapDF$Gen == "G-1990"
I1990 <- IncapDF[f, ]

f <- SatDF$Gen == "G-1990"
S1990 <- SatDF[f, ]

## days of activity (no lost)
R1990$NatDeath <- R1990$Birth+69*365
MaxAct <- sum(R1990$NatDeath - R1990$dateK)

## days saved to homicides
K1990$NatDeath <- K1990$Birth+69*365
KilledS <- sum(K1990$NatDeath-K1990$dateK)

## days saved to saturation
S1990$NatDeath <- S1990$Birth+69*365
SaturS <- sum(S1990$NatDeath-S1990$dateK)

### saved to incarceration
Sentence <- read.csv("ENPOL sentences.csv")

for(k in 1:250){
  Alpha = 2*runif(1)
  RecidAnalysis <- AlphaRecidiv(IncapDF, 
                                Alpha = Alpha,
                                Efs,### the Active members
                                Sentence,
                                kappaP, ### arrest
                                kappaK, ### killed
                                sat)
  IncapS <- sum(RecidAnalysis$SENTENCE*365)
  CrimeSavedIncarcAlpha <- rbind(CrimeSavedIncarcAlpha, 
                            data.frame(Alpha = Alpha,
                                       MaxAct = MaxAct,
                                       KilledS = KilledS,
                                       SaturS = SaturS,
                                       IncapS = IncapS))
  cat(k, Alpha, "\t", round(IncapS/MaxAct,4), "\n")
}

save(CrimeSavedIncarcAlpha, 
     file = "IncarcerationSavedCrimeAlpha.RData")
CrimeSavedIncarcAlpha/CrimeSavedIncarcAlpha$MaxAct
CrimeSavedIncarcAlpha$IncapS/CrimeSavedIncarcAlpha$MaxAct*(1 + 0.26) ### the 0.26 is reincapacitation
}

#### analysis of time spent in prison
{
timeSpentPrison <- sum(RecidAnalysis$SENTENCE)*(1 + 0.26)
timeSpentPrison/16000

### total number of recruits in 2000-2040:
sum(RecruitDF$dateK >= 36526 & RecruitDF$dateK <= 51501 )

f <- KilledDF$Gen == "G-1990"
G1990K <- KilledDF[f,]
6500*(75-mean(KilledDF$Age))/1400000
RecruitDF$year <- YearF(RecruitDF$dateK)
RecsTable <- table(RecruitDF$Gen, RecruitDF$year)
RecruitsShare1990 <- RecsTable[41, ]/apply(RecsTable, 2, sum)

KilledDF$year <- YearF(KilledDF$dateK)
KillsTable <- table(KilledDF$Gen, KilledDF$year)
KillsShare1990 <- KillsTable[61, ]/apply(KillsTable, 2, sum)

### years in prison and crimes averted
f <- which (CrimeSavedIncarcAlpha$Alpha == 1)
CrimeSavedIncarcAlpha[f,]
CrimeSavedIncarcAlpha$MaxAct[f]/(16000*365)
CrimeSavedIncarcAlpha$MaxAct[f]/(365)

#### years of prison life
69-mean(IncapDF$Age)
sum(IncapDF$Gen == "G-1990") * (69-mean(IncapDF$Age))

mean(IncapDF$dateK - IncapDF$Birth)/365
}
