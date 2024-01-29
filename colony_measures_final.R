library(glmmTMB)
library(DHARMa)
library(emmeans)
library(dplyr)
library(ggpubr)
library(nlme)
library(lme4)
library(car)
library("cowplot")
library(gridExtra)
library(ggpubr)


# loading data only for which there are matched pairs
dat <- read.csv("Nests_matched.csv")
dat$Size <- as.factor(dat$Size)
dat$Inoculated = as.factor(dat$Inoculated)
dat$Biobest = as.factor(dat$Biobest)
dat$Colony = as.factor(dat$Colony)

# color palette for figures
caPalette <- c("#0072B2", "#D55E00", "#E69F00", "#CC79A7")



###################################
# TOTAL NEST WEIGHT MODEL
###################################

weight_model <- lmer(log(Nest_Weight) ~ Size*Inoculated + (1 | Colony), data = dat, REML = FALSE)
summary(weight_model)
Anova(weight_model)

# without interaction, because interaction was not significant
weight_model2 <- lmer(log(Nest_Weight) ~ Size + Inoculated + (1 | Colony), data = dat, REML = FALSE)
summary(weight_model2)
Anova(weight_model2)


# CHECKING ASSUMPTIONS: Normality and Homogenous Residuals
# check for normality of residuals (not response variable)
# normality for linear models can be eyeballed - shapiro test is often too precise (-Erika)
ggqqplot(residuals(weight_model))

# test if residuals are homogenous; you want to see no particular pattern 
plot(predict(weight_model), residuals(weight_model))


# contrast tests, only looked at effect of inoculation within each size
emmeans(weight_model, specs=pairwise ~ Inoculated|Size, type="response")

# with bonferonni correction 
emm1 <- emmeans(weight_model, ~ Inoculated|Size, type="response")
summary(pairs(emm1), by = NULL, adjust = "bonferroni")



###################################
# TOTAL PUPAE + MALE MODEL
###################################


fitness_model <- lmer(log(EggMale) ~ Size*Inoculated + (1 | Colony), data = dat, REML = FALSE)
summary(fitness_model)
Anova(fitness_model)

# without interaction because it was not significant 
fitness_model2 <- lmer(log(EggMale) ~ Size + Inoculated + (1 | Colony), data = dat, REML = FALSE)
summary(fitness_model2)
Anova(fitness_model2)

# testing for normality and homogenous residuals
ggqqplot(residuals(fitness_model))
plot(predict(fitness_model), residuals(fitness_model))


# contrast tests, only look at effect of inoculation within each size
emmeans(fitness_model, specs=pairwise ~ Inoculated|Size, type="response")


# with bonferonni correction 
emm2 <- emmeans(fitness_model, ~ Inoculated|Size, type="response")
summary(pairs(emm2), by = NULL, adjust = "bonferroni")


##########################
# PER CAPITA NEST WEIGHT MODEL
##########################


percap_weight_model <- lmer(log(Per_Capita_Weight) ~ Size*Inoculated + (1 | Colony), data = dat, REML = FALSE)
summary(percap_weight_model)
Anova(percap_weight_model)

# without interaction to report fixed effects because interaction was not significant
percap_weight_model2 <- lmer(log(Per_Capita_Weight) ~ Size+Inoculated + (1 | Colony), data = dat, REML = FALSE)
summary(percap_weight_model2)
Anova(percap_weight_model2)

# CHECKING ASSUMPTIONS: Normality and Homogenous Residuals
# check for normality of residuals (not response variable)
# normality for linear models can be eyeballed - shapiro test is often too precise (-Erika)
ggqqplot(residuals(percap_weight_model))

# test if residuals are homogenous; you want to see no particular pattern 
plot(predict(percap_weight_model), residuals(percap_weight_model))


# contrast tests, only looked at effect of inoculation within each size
emmeans(percap_weight_model, specs=pairwise ~ Inoculated|Size, type="response")


# with bonferonni correction 
emm3 <- emmeans(percap_weight_model, ~ Inoculated|Size, type="response")
summary(pairs(emm3), by = NULL, adjust = "bonferroni")


##########################
# PER CAPITA PUPAE+MALE MODEL
##########################


percap_eggmale_model <- lmer(log(Per_Capita_EggMale) ~ Size*Inoculated + (1 | Colony), data = dat, REML = FALSE)
summary(percap_eggmale_model)
Anova(percap_eggmale_model)


percap_eggmale_model2 <- lmer(log(Per_Capita_EggMale) ~ Size+Inoculated + (1 | Colony), data = dat, REML = FALSE)
summary(percap_eggmale_model2)
Anova(percap_eggmale_model2)


# CHECKING ASSUMPTIONS: Normality and Homogenous Residuals
ggqqplot(residuals(percap_eggmale_model))
# test if residuals are homogenous; you want to see no particular pattern 
plot(predict(percap_eggmale_model), residuals(percap_eggmale_model))


# contrast tests, only looked at effect of inoculation within each size
emmeans(percap_eggmale_model, specs=pairwise ~ Inoculated|Size, type="response")


# with bonferonni correction 
emm4 <- emmeans(percap_eggmale_model, ~ Inoculated|Size, type="response")
summary(pairs(emm4), by = NULL, adjust = "bonferroni")


# main results figures
ggarrange(weight,fitness,ncol=2,nrow=1)

