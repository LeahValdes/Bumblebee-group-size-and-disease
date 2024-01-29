#####################
# 1. Initialisation #
#####################
library(emmeans)
library(ggplot2)
library(glmmTMB)
library(lme4)
library(readxl)
library(scales)


dat <- read.csv("Colony_Proportion.csv")

# Convert factor variables into factors.
# But before that, keep a copy of the numerical version of "Size" as "Size_num".
dat$Size_num = dat$Size
dat$Size = as.factor(dat$Size)
dat$Trial = as.factor(dat$Trial)
dat$Biobest = as.factor(dat$Biobest)
dat$Micro = as.factor(dat$Micro)
dat$Colony = as.factor(dat$Colony)

# Create columns relevant to infection proportion, used for the plots.
dat$Screened = dat$Infected + dat$Healthy
dat$Infected_Prop = dat$Infected / dat$Screened

# Create a subset corresponding to Week 3 data.
# This will be used when assessing infection at the end.
dat_3 = dplyr::filter(dat, Week==3)





############
# 2. Plots #
############

caPalette <- c("#0072B2", "#D55E00", "#E69F00", "#CC79A7")


# Infection at the end, separated by colony.
ggplot(data=dat_3, 
       aes(x=Size, y=Infected_Prop, color=Colony)) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_line(aes(group=Colony), position=position_dodge(width=0.2)) +
  labs(y = "Proportion Infected", title = "Infection prevalence at the end (Week 3)") + 
  theme(plot.title = element_text(face = "bold", size = 12)) + cowplot::theme_cowplot()



# Weekly infection, separated by colony.
ggplot(data=dat, 
       aes(x=Week, y=Infected_Prop, color=Size)) +
  geom_point(size=3) +
  geom_line(aes(group=Micro)) +
  facet_wrap("Colony", ncol=5) + 
  labs(y = "Proportion Infected", title = "Infection over time") +
  scale_colour_manual(values=caPalette) + 
  theme(plot.title = element_text(face = "bold", size = 12))



# Weekly infection, pooled across colony.
p1 <- ggplot(data=dat,
       aes(x=Week, y=Infected_Prop, color=Size)) +
  geom_point(size=3, position=position_dodge(0.12)) +
  geom_smooth(method="glm", method.args=list(family="binomial"),
              aes(x=Week, y=Infected_Prop, group=Size, weight=Screened), se=F) +
  labs(y = "Proportion Infected", tag = "(a)") + 
  scale_colour_manual(values=caPalette) + 
  theme(plot.title = element_text(face = "bold", size = 12)) +
  cowplot::theme_cowplot()

p2 <- ggplot(data=dat_3, aes(x=Size, y=Infected_Prop, fill=Size)) + 
  geom_boxplot(show.legend=FALSE) + scale_fill_manual(values=caPalette) +
  geom_point(position=position_dodge(width = .75), show.legend = FALSE) + 
  scale_x_discrete(labels=c("Small\n (15 bees)", "Large\n (45 bees)")) +
  theme_light() + labs(y="Proportion Infected", tag = "(b)") + labs(fill="Size") +
  cowplot::theme_cowplot()

Figure3 <- ggarrange(p1,p2,nrow = 1,ncol=2)


ggsave(Figure3, file = "/Users/leah/Library/CloudStorage/OneDrive-CornellUniversity/Cornell/Density_Project/FIGURES/Figure3.jpg", 
       width=10, height=5.5)

###############################################
# 3. Models for infection at the end (Week 3) #
###############################################

# Mixed model 
glmer_week3 = glmer(cbind(Infected,Healthy) ~ Size + (1|Colony), dat=dat_3, family="binomial")
summary(glmer_week3)


ggqqplot(residuals(glmer_week3))
# test if residuals are homogenous; you want to see no particular pattern 
plot(predict(glmer_week3), residuals(glmer_week3))


###################################
# 4. Models for weekly infection. #
###################################

# The goal is to see if the slope is steeper for a larger colony.
# In other words, we want to see if Week:Size is significant. (it's not)

# Mixed model version.
glmer_weekly = glmer(cbind(Infected,Healthy) ~ Week*Size + (1|Colony), data=dat, family="binomial")
summary(glmer_weekly)

ggqqplot(residuals(glmer_weekly))
# test if residuals are homogenous; you want to see no particular pattern 
plot(predict(glmer_weekly), residuals(glmer_weekly))

