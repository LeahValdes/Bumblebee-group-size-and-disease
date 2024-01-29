#####################
# 1. Initialisation #
#####################
library(DHARMa)
library(glmmTMB)
library(readxl)
library(ggplot2)
library(ggpubr)
library(lme4)

dat <- read.csv("Dissections.csv")
# Filter out male bees.
dat = dat[dat$BeeID!="male",]
# Filter out uninfected bees.
dat = dat[dat$Crithidia>0,]


# Convert factor variables into factors.
# But before that, keep a copy of the numerical version of "Size" as "Size_num".
dat$Size_num = dat$Size
dat$Size = as.factor(dat$Size)
dat$Trial = as.factor(dat$Trial)
dat$Biobest = as.factor(dat$Biobest)
dat$Micro = as.factor(dat$Micro)
dat$Colony = as.factor(dat$Colony)

# dat_figures: Do not filter out bees without marginal cell measurements
dat_figures <- dat
# dat: Filter out bees without marginal cell measurements
dat <- na.omit(dat)


#############
# 3. Models #
#############

# Comparing infection intensity between small and large colonies
# Main effects: Colony size and body size (marginal cell length)
# Truncated negative binomial distribution to account for filtering out uninfected bees.
model_nb_glmmTMB1 = glmmTMB(Crithidia ~ Size + Marginal_Cell + (1|Colony), data=dat, family=truncated_nbinom2)
summary(model_nb_glmmTMB1)
drop1(model_nb_glmmTMB1, test="Chisq")

# Is a truncated negative binomial distribution appropriate? Are the residuals well-behaved? Yes
plot( simulateResiduals(model_nb_glmmTMB1) )

##################

# Comparing body size (marginal cell) between small and large colonies
# Using lmm because data are ~normally distributed
bodysize_model = lmer(Marginal_Cell ~ Size + (1|Colony), data = dat)
summary(bodysize_model)
Anova(bodysize_model)

# testing model assumptions: normality and homogeneity of residuals
plot(simulateResiduals(bodysize_model))
ggqqplot(residuals(bodysize_model))
