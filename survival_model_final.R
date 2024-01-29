#####################
# 1. Initialisation #
#####################
library(emmeans)
library(logistf)
library(ggplot2)
library(lme4)
library(readxl)
library(scales)
library(survival)
library(survminer)
library(emmeans)
library(ggpubr)

dat <- read.csv("Death.csv")

# Remove the microcolony where all bees starved to death.
dat = dplyr::filter(dat, Microcolony!="5clean45")

# Fix some small inconsistencies.
# For 7bclean45 Week 3, (Dead,Alive) = (10,34), which doesn't add to 45.
# Let's assume that the missing bee escaped and is alive.
dat$Alive[dat$Microcolony=="7bclean45" & dat$Week==3] = 35
# For 7cClean15, (Dead,Alive) = (1,14) for Weeks 1 and 2,
# and (0,15) for Week 3.
# Let's ignore the extra alive bee in Week 3.
dat$Dead[dat$Microcolony=="7cClean15" & dat$Week==3] = 1
dat$Alive[dat$Microcolony=="7cClean15" & dat$Week==3] = 14

# Match "Colony" of each control to "Colony" of the paired treatments,
# since we already have the variable "Crithidia" to specify inoculation status.
dat$Colony[dat$Colony=="1clean"] = "1a"
dat$Colony[dat$Colony=="3clean"] = "3b"
dat$Colony[dat$Colony=="5clean"] = "5b"
dat$Colony[dat$Colony=="6clean"] = "6b"

# Let's work only with "complete" colonies that have all four treatments.
complete_colonies = c("1a", "3b", "6b", "7a", "7b", "7c", "8a", "8b", "8c","9a","10a")
dat = dplyr::filter(dat, Colony %in% complete_colonies)

# Convert factor variables into factors.
# But before that, keep a copy of the numerical version of "Size" as "Size_num".
dat$Size_num = dat$Size
dat$Size = as.factor(dat$Size)
dat$Crithidia = as.factor(dat$Crithidia)
dat$Trial = as.factor(dat$Trial)
dat$Biobest = as.factor(dat$Biobest)
dat$Microcolony = as.factor(dat$Microcolony)
dat$Colony = as.factor(dat$Colony)

# Rename "Fraction Dead" to "Dead_Prop" since the latter name is easier to work with.
#dat$Dead_Prop = dat$"Fraction.Dead"

# Create a subset corresponding to Week 3 data.
# This will be used when assessing mortality at the end.
dat_3 = dplyr::filter(dat, Week==3)

#################################
# 2. Preparing "long form" data #
#################################
# This is required for Cox and parametric survival models.

long_dead_dat = c()
long_alive_dat = c()

for(Mic in unique(dat$Microcolony)){
  # Create a subset corresponding to the current microcolony.
  dat_Mic = dplyr::filter(dat, Microcolony==Mic)
  dat_Mic_3 = dplyr::filter(dat_Mic, Week==3)
  
  # Right now, "Dead" is cumulative.
  # Here, we create a new column that gives the number that died each week.
  # First, make sure the rows in dat_Mic are ordered by weeks.
  dat_Mic = dat_Mic[order(dat_Mic$Week),]
  # This line should then give the number that died each week.
  dat_Mic$Dead_week = dat_Mic$Dead - c(0, dat_Mic$Dead[-nrow(dat_Mic)])

  # Now create the long form data.
  long_dead_dat =  rbind(long_dead_dat,
                         dat_Mic[rep(seq_len(nrow(dat_Mic)), dat_Mic$Dead_week), c("Trial", "Week", "Biobest", "Colony", "Microcolony", "Size", "Crithidia")])
  long_alive_dat = rbind(long_alive_dat,
                         dat_Mic_3[rep(1, dat_Mic_3$Alive), c("Trial", "Week", "Biobest", "Colony", "Microcolony", "Size", "Crithidia")])
}
long_dead_dat$Status_dead = 1
long_alive_dat$Status_dead = 0
long_dat = rbind(long_dead_dat, long_alive_dat)

#############
# 4. Models #
#############


# Cox survival models.
model_cox = coxph(Surv(Week,Status_dead) ~ Size*Crithidia + Colony, dat=long_dat)
summary(model_cox)
# Interpreting the output.
# - Note that for all the fitted coefficients below, a positive estimate means
#   higher mortality, since it means that the hazard rate is higher.
#
# - "ColonyXX" is the difference between Colony XX and Colony 1a for small
#   control colonies.
#
# - "Size45" is the difference between large and small control colonies.
#
# - "CrithidiaInoculated" is the difference between small inoculated and 
#   small control colonies.
#
# - "Size45:CrithidiaInoculated" is the ADDITIONAL difference between large
#   inoculated and large control colonies, on top of "CrithidiaInoculated".
#   This is the one telling us that colony size affects the difference in
#   mortality between control and inoculated colonies.

# Contrast tests 
# Wee Hao: This directly supports our claim that disease affects survival for large colonies. 
# A significant interaction does tell us that disease effects are different between small and large colonies, 
# but it does not directly support our claim.

pairs(emmeans(model_cox, ~ Crithidia|Size))

# with bonferroni correction 
emm <- emmeans(model_cox, ~ Crithidia|Size)
summary(pairs(emm), by=NULL, adjust="bonferroni")

# GLM for mortality at the end (Week 3).
model_glm = glm(cbind(Dead,Alive) ~ Size*Crithidia + Colony, dat=dat_3, family="binomial")
summary(model_glm)


# contrast tests, only looked at effect of inoculation within each size
emmeans(model_glm, specs=pairwise ~ Crithidia|Size, type="response")

# contrast tests with bonferroni correction 
emm1 <- emmeans(model_glm, ~ Crithidia|Size)
summary(pairs(emm1), by=NULL, adjust = "bonferroni")


