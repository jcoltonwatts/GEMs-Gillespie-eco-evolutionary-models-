# Load in data for simulations initialized with both competitors below
# the unstable equilibrium.
# Make sure that these relative paths point to the location of the data files.
noVarBelow = 
  read.csv("data/Extinctions_j1_below.csv",
                      header = FALSE)
noEvolBelow = 
  read.csv("data/Extinctions_j2_below.csv",
                      header = FALSE)
evolBelow = 
  read.csv("data/Extinctions_j3_below.csv",
                      header = FALSE)

# rename columns so that times come out as numeric in long format
colnames(noVarBelow) <- c(1:ncol(noVarBelow))
colnames(noEvolBelow) <- c(1:ncol(noEvolBelow))
colnames(evolBelow) <- c(1:ncol(evolBelow))

# add a column identifying each row as a replicate
noVarBelow$replicate = 1:nrow(noVarBelow)
noEvolBelow$replicate = 1:nrow(noEvolBelow)
evolBelow$replicate = 1:nrow(evolBelow)

# transform to long format
library(tidyr)
# gather(data, new key column, new value column, names of source columns with values)
noVarBelowLong = gather(noVarBelow, time, extinct, 1:ncol(noVarBelow)-1)
noEvolBelowLong = gather(noEvolBelow, time, extinct, 1:ncol(noEvolBelow)-1)
evolBelowLong = gather(evolBelow, time, extinct, 1:ncol(evolBelow)-1)

# add a variable to each encoding which scenario type
noVarBelowLong$scenario = "NoVar"
noEvolBelowLong$scenario = "NoEvol"
evolBelowLong$scenario = "Evol"

# now combine them
extsBelow = rbind(noVarBelowLong, noEvolBelowLong, evolBelowLong)

# make sure time variable is numeric
str(extsBelow$time) # it's not
as.numeric(extsBelow$time) # but this works
extsBelow$ntime = as.numeric(extsBelow$time) # so save it as "ntime"
str(extsBelow$ntime) # great.

# make sure that scenario is a factor, and that noVar is the reference level.
str(extsBelow$scenario)
extsBelow$scenario = factor(extsBelow$scenario, ordered = FALSE)
str(extsBelow$scenario)
extsBelow$scenario = relevel(extsBelow$scenario, ref = "NoVar")
str(extsBelow$scenario)

# fit a survival model with scenario type
library(survival)
mod1 = coxph(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsBelow)
anova(mod1)
plot(survfit(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsBelow), lty = c(1,2,3))
summary(mod1)
# check proportional hazards assumption
plot(cox.zph(mod1, terms = TRUE, transform = "identity"))

# repeat, but only comparing no evolution and evolution scenarios, to more
# easily interpret the proportional hazards change over time
mod2 = coxph(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsBelow[which(extsBelow$scenario != "NoVar"),])
anova(mod2)
plot(survfit(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsBelow[which(extsBelow$scenario != "NoVar"),]), lty = c(1,2))
summary(mod2)
plot(cox.zph(mod2, terms = TRUE), xlim = c(0, 1))


# fit an explicitly time-varying effect of scenario using tt() functionality
mod3 = coxph(Surv(ntime, extinct, type = "right") ~ scenario + tt(scenario),
             data = extsBelow[which(extsBelow$scenario != "NoVar"),] ,
             tt = function(x, t, ...) {
               mtrx <- model.matrix(~x)[,-1]
               #mtrx * log(t)
               #x * nsk(x, knots = seq(1,1200,299), Boundary.knots = FALSE)
               pspline(x*log(t+1))
             }
             )
anova(mod3)
summary(mod3)
plot(cox.zph(coxph(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsBelow[which(extsBelow$scenario != "NoVar"),]), terms = TRUE, transform = "identity")) 
abline(coef(mod3)[3], lwd=2, lty=3, col=2)
coef(mod3)[3]


# spline?
mod4 = coxph(Surv(ntime, extinct, type = "right") ~ scenario + tt(scenario),
             data = extsBelow[which(extsBelow$scenario != "NoVar"),] ,
             tt = function(x, t, ...) {
               mtrx <- model.matrix(~x)[,-1]
               #mtrx * log(t)
               #x * nsk(x, knots = seq(1,1200,299), Boundary.knots = FALSE)
               pspline(x*log(t+1))
             }
)
anova(mod3)
summary(mod3)
plot(cox.zph(coxph(Surv(ntime, extinct, type = "right") ~ scenario,
                   data = extsBelow[which(extsBelow$scenario != "NoVar"),]), terms = TRUE, transform = "identity")) 
abline(coef(mod3)[3], lwd=2, lty=3, col=2)
coef(mod3)[3]


#### Repeat for scenarios that start at the unstable ecological equilibrium ####

noVarAt = 
  read.csv("/Volumes/GoogleDrive-108856664039640963521/My Drive/PDocGEMs/gemfiles/data_202308/Extinctions_j1_at.csv",
           header = FALSE)
noEvolAt = 
  read.csv("/Volumes/GoogleDrive-108856664039640963521/My Drive/PDocGEMs/gemfiles/data_202308/Extinctions_j2_at.csv",
           header = FALSE)
evolAt = 
  read.csv("/Volumes/GoogleDrive-108856664039640963521/My Drive/PDocGEMs/gemfiles/data_202308/Extinctions_j3_at.csv",
           header = FALSE)

# rename columns so that times come out as numeric in long format
colnames(noVarAt) <- c(1:ncol(noVarAt))
colnames(noEvolAt) <- c(1:ncol(noEvolAt))
colnames(evolAt) <- c(1:ncol(evolAt))

# add a column identifying each row as a replicate
noVarAt$replicate = 1:nrow(noVarAt)
noEvolAt$replicate = 1:nrow(noEvolAt)
evolAt$replicate = 1:nrow(evolAt)

# transform to long format
# gather(data, new key column, new value column, names of source columns with values)
noVarAtLong = gather(noVarAt, time, extinct, 1:ncol(noVarAt)-1)
noEvolAtLong = gather(noEvolAt, time, extinct, 1:ncol(noEvolAt)-1)
evolAtLong = gather(evolAt, time, extinct, 1:ncol(evolAt)-1)

# add a variable to each encoding which scenario type
noVarAtLong$scenario = "NoVar"
noEvolAtLong$scenario = "NoEvol"
evolAtLong$scenario = "Evol"

# now combine them
extsAt = rbind(noVarAtLong, noEvolAtLong, evolAtLong)

# make sure time variable is numeric
str(extsAt$time) # it's not
as.numeric(extsAt$time) # but this works
extsAt$ntime = as.numeric(extsAt$time) # so save it as "ntime"
str(extsAt$ntime) # great.

# make sure that scenario is a factor, and that noVar is the reference level.
str(extsAt$scenario)
extsAt$scenario = factor(extsAt$scenario, ordered = FALSE)
str(extsAt$scenario)
extsAt$scenario = relevel(extsAt$scenario, ref = "NoVar")
str(extsAt$scenario)

# fit a survival model with scenario type
mod1At = coxph(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsAt)
anova(mod1At)
plot(survfit(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsAt), lty = c(1,2,3))
summary(mod1At)
# check proportional hazards assumption
plot(cox.zph(mod1At, terms = TRUE, transform = "identity"))

# repeat, but only comparing no evolution and evolution scenarios, to more
# easily interpret the proportional hazards change over time
mod2At = coxph(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsAt[which(extsAt$scenario != "NoVar"),])
anova(mod2At)
plot(survfit(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsAt[which(extsAt$scenario != "NoVar"),]), lty = c(1,2))
summary(mod2At)
plot(cox.zph(mod2At, terms = TRUE, transform = "identity"))



#### Repeat for scenarios that start ABOVE the unstable ecological equilibrium ####

# forgot to save the extinction csv for these starting conditions, need to calculate again from abundances
# (or traits, which become NAs when population is extinct)

noVarAbove = 
  read.csv("/Volumes/GoogleDrive-108856664039640963521/My Drive/PDocGEMs/gemfiles/data_202308/Extinctions_j1_above.csv",
           header = FALSE)
noEvolAbove = 
  read.csv("/Volumes/GoogleDrive-108856664039640963521/My Drive/PDocGEMs/gemfiles/data_202308/Extinctions_j2_above.csv",
           header = FALSE)
evolAbove = 
  read.csv("/Volumes/GoogleDrive-108856664039640963521/My Drive/PDocGEMs/gemfiles/data_202308/Extinctions_j3_above.csv",
           header = FALSE)

# rename columns so that times come out as numeric in long format
colnames(noVarAbove) <- c(1:ncol(noVarAbove))
colnames(noEvolAbove) <- c(1:ncol(noEvolAbove))
colnames(evolAbove) <- c(1:ncol(evolAbove))

# add a column identifying each row as a replicate
noVarAbove$replicate = 1:nrow(noVarAbove)
noEvolAbove$replicate = 1:nrow(noEvolAbove)
evolAbove$replicate = 1:nrow(evolAbove)

# transform to long format
# gather(data, new key column, new value column, names of source columns with values)
noVarAboveLong = gather(noVarAbove, time, extinct, 1:ncol(noVarAbove)-1)
noEvolAboveLong = gather(noEvolAbove, time, extinct, 1:ncol(noEvolAbove)-1)
evolAboveLong = gather(evolAbove, time, extinct, 1:ncol(evolAbove)-1)

# add a variable to each encoding which scenario type
noVarAboveLong$scenario = "NoVar"
noEvolAboveLong$scenario = "NoEvol"
evolAboveLong$scenario = "Evol"

# now combine them
extsAbove = rbind(noVarAboveLong, noEvolAboveLong, evolAboveLong)

# make sure time variable is numeric
str(extsAbove$time) # it's not
as.numeric(extsAbove$time) # but this works
extsAbove$ntime = as.numeric(extsAbove$time) # so save it as "ntime"
str(extsAbove$ntime) # great.

# make sure that scenario is a factor, and that noVar is the reference level.
str(extsAbove$scenario)
extsAbove$scenario = factor(extsAbove$scenario, ordered = FALSE)
str(extsAbove$scenario)
extsAbove$scenario = relevel(extsAbove$scenario, ref = "NoVar")
str(extsAbove$scenario)

# fit a survival model with scenario type
mod1Above = coxph(Surv(ntime, extinct, type = "right") ~ scenario,
               data = extsAbove)
anova(mod1Above)
plot(survfit(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsAbove), lty = c(1,2,3))
summary(mod1Above)
# check proportional hazards assumption
plot(cox.zph(mod1Above, terms = TRUE, transform = "identity"))

# repeat, but only comparing no evolution and evolution scenarios, to more
# easily interpret the proportional hazards change over time
extsAboveNoNovar = extsAbove[which(extsAbove$scenario != "NoVar"),]
extsAboveNoNovar$scenario =factor(extsAboveNoNovar$scenario )
extsAboveNoNovar$scenario = relevel(extsAboveNoNovar$scenario, ref = "NoEvol")
str(extsAboveNoNovar$scenario)
mod2Above = coxph(Surv(ntime, extinct, type = "right") ~ scenario,
               data = extsAboveNoNovar)
anova(mod2Above)
plot(survfit(Surv(ntime, extinct, type = "right") ~ scenario,
             data = extsAbove[which(extsAbove$scenario != "NoVar"),]), lty = c(1,2))
summary(mod2Above)
plot(cox.zph(mod2Above, terms = TRUE, transform = "identity"))

str(cox.zph(mod2Above, terms = TRUE, transform = "identity"))




#### Nice plot ####
# panel plot showing time-varying hazard ratio between Evol/No Evol scenarios
library(ggplot2)
library(gridExtra)

# label each df based on starting ecological conditions
extsBelow$InitialCond = "Below"
extsAt$InitialCond = "At"
extsAbove$InitialCond = "Above"

# bind them together into a new df
exts = rbind(extsBelow, extsAt, extsAbove)
# create version subsetting to only NoEvol and Evol scenarios
extsNoNoVar = exts[-which(exts$scenario == "NoVar"),]

# create a blank df to hold the HR for each time and initial condition
betaOfts = data.frame(expand.grid("time" = 1:1200,
                                  "InitialCond" = c("Below", "At", "Above"),
                                  "HazardRatio" = NA,
                                  "OddsRatio" = NA))
# create a vector of the unique values of initial conditions
InitialConds = unique(betaOfts$InitialCond)
# loop over them
for (i in 1:length(InitialConds)){
  # for ecah initial ecological condition...
  # fit the model for this initial condition and store the result of cox.zph()
  this.zph = cox.zph(coxph(Surv(ntime, extinct, type = "right") ~ scenario,
                           data = extsNoNoVar[which(extsNoNoVar$InitialCond == 
                                                      InitialConds[i]),]),
                     terms = TRUE,
                     transform = "identity")
  # loop over times and fill in the corresponding hazard, if applicable
  # (can't simply fill in, b/c lots of NAs)
  for (t in 1:max(betaOfts$time)){
    if(length(this.zph$y[which(this.zph$time == t)])>0){
      betaOfts$HazardRatio[which(betaOfts$InitialCond == InitialConds[i] &
                                   betaOfts$time == t)] = 
        this.zph$y[which(this.zph$time == t)]
    }
  }
}

resids = coxph(Surv(ntime, extinct, type = "right") ~ scenario,
      data = extsNoNoVar[which(extsNoNoVar$InitialCond == 
                                 InitialConds[i]),])$resid
resids = as.data.frame(resids)
str(resids)

ggplot(betaOfts, aes(x = time, y = HazardRatio, group = InitialCond)) + 
  geom_point() +
  geom_smooth(data = betaOfts, aes(group = InitialCond))


# use survminer package to make nice plots of cox.zph output
library("survminer")
library("ggplot2")

ggcoxzph(cox.zph(mod2, terms = T, transform = "identity"), resid = T, se = F, df = 3, ylim = c(-2,2))
ggcoxzph(cox.zph(mod2At, terms = T, transform = "identity"), resid = F, se = F, df = 3, ylim = c(-2,2))
ggcoxzph(cox.zph(mod2Above, terms = T, transform = "identity"), resid = F, se = F, df = 3, ylim = c(-2,2))

quartz()
par(mfrow = c(1,3),
    mar = c(5,4,1,1))
plot(cox.zph(mod2, terms = F, transform = "identity"), xlim = c(1,1200))#, df = 4)
plot(cox.zph(mod2At, terms = F, transform = "identity"), xlim = c(1,1200))#, df = 4)
plot(cox.zph(mod2Above, terms = F, transform = "identity"), xlim = c(1,1200))#, df = 4)

# is there any way to convert the y-axis to an odds ratio before plotting?
# maybe one step in this direction is to specify hr = T to plot
# the HR itself

quartz()
par(mfrow = c(1,3),
    mar = c(5,4,1,1))
plot(cox.zph(mod2, terms = F, transform = "identity"), xlim = c(1,1200), ylim = c(-3,3))
plot(cox.zph(mod2At, terms = F, transform = "identity"), xlim = c(1,1200), ylim = c(-3,3))
plot(cox.zph(mod2Above, terms = F, transform = "identity"), xlim = c(1,1200), ylim = c(-3,3))

summary(mod2)
summary(mod2At)     
summary(mod2Above)     
