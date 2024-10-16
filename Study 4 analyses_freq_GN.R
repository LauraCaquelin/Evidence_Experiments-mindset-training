# Script from Yeager et al, adapted by Gustav Nilsonne

# Set working directory to source file directory
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)

# Require libraries
library(tidyr)
library(dplyr)
library(readr)
library(lme4)
#library(multibart) # Not able to install
library(tidyverse)
library(effects)

# Read data
dat_long <- read_csv("study4long.csv")

# Clean data and create variables for analyses
dat_long$tpr[dat_long$tpr<650] = NA
dat_long$tpr[dat_long$tpr>5000] = NA

dat_long$pid <- factor(dat_long$pid)

dat_long = dat_long %>% dplyr::group_by(pid) %>% 
  mutate(
         mean_base = base::mean(tpr[time<6], na.rm=TRUE),
         mean_base_sv = mean(sv[time<6], na.rm=TRUE), # toggle these depending on the outcome
         mean_base_pep = mean(pep[time<6], na.rm=TRUE))
         
dat_long <- dat_long %>%
  drop_na(
    tpr,
    #sv, # toggle these depending on the outcome
    #pep,
    cond.factor)

# Drop two subjects with <= 2 observations
goodid = dat_long %>% dplyr::summarize(n=n()) %>% filter(n>2) %>% select(pid)

# Create tsst reactivity dv:s
dat_long$tpr_react <- dat_long$tpr - dat_long$mean_base
#dat_long$sv_react <- dat_long$sv - dat_long$mean_base_sv
#dat_long$pep_react <- dat_long$pep - dat_long$mean_base_pep

# Make factors of continuous predictors that should be factors
dat_long$stressconds <- as_factor(dat_long$stressconds)
dat_long$mindsetconds <- as_factor(dat_long$mindsetconds)
dat_long$stressonly <- as_factor(dat_long$stressonly)
dat_long$msonly <- as_factor(dat_long$msonly)
dat_long$synergistic <- as_factor(dat_long$synergistic)
dat_long$racefactor <- as_factor(dat_long$s1.race)

# Set response variable
dat_long <- dat_long %>% filter(is.na(manipcheck)==F) # Toggle this depending on whether using the appraisal/well-being DVs
#dat_long <- dat_long %>% filter(is.na(wellbeing)==F) # Toggle this depending on whether using the appraisal/well-being DVs

# Mean center outcome variable
dat_long$tpr_react_mc_z <- as.vector(scale(dat_long$tpr_react))

# Modelling
# TPR is the main outcome if interest
mod_tpr_factors_unadjusted <- lmer(tpr_react_mc_z ~ stressonly + msonly + synergistic + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_factors_unadjusted)
confint(mod_tpr_factors_unadjusted)
confint(mod_tpr_factors_unadjusted, level = 0.98)
plot(mod_tpr_factors_unadjusted)
plot(Effect(c("stressonly"), mod_tpr_factors_unadjusted), ylim = c(-0.2, 1))
plot(Effect(c("msonly"), mod_tpr_factors_unadjusted), ylim = c(-0.2, 1))
plot(Effect(c("synergistic"), mod_tpr_factors_unadjusted), ylim = c(-0.2, 1))

mod_tpr_interaction_unadjusted <- lmer(tpr_react_mc_z ~ stressconds * mindsetconds + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_interaction_unadjusted)
plot(mod_tpr_interaction_unadjusted)
plot(Effect(c("stressconds", "mindsetconds"), mod_tpr_interaction_unadjusted))

mod_tpr_factors_fullyadjusted <- lmer(tpr_react_mc_z ~ stressonly + msonly + synergistic +
                                       s1.sex + s1.age + racefactor + pss.baseline + 
                                       stressmindset.baseline + fixedmindset.baseline +
                                       bothmindsets.baseline + selfesteem.baseline +
                                       testanxiety.baseline + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_factors_fullyadjusted)
plot(mod_tpr_factors_fullyadjusted)
confint(mod_tpr_factors_fullyadjusted)
confint(mod_tpr_factors_fullyadjusted, level = 0.98)
plot(Effect(c("stressonly"), mod_tpr_factors_fullyadjusted), ylim = c(-0.2, 1))
plot(Effect(c("msonly"), mod_tpr_factors_fullyadjusted), ylim = c(-0.2, 1))
plot(Effect(c("synergistic"), mod_tpr_factors_fullyadjusted), ylim = c(-0.2, 1))

mod_tpr_interaction_fullyadjusted <- lmer(tpr_react_mc_z ~ stressconds * mindsetconds +
                                           s1.sex + s1.age + racefactor + pss.baseline + 
                                           stressmindset.baseline + fixedmindset.baseline +
                                           bothmindsets.baseline + selfesteem.baseline +
                                           testanxiety.baseline + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_interaction_fullyadjusted)
plot(mod_tpr_interaction_fullyadjusted)
confint(mod_tpr_interaction_fullyadjusted)
plot(Effect(c("stressconds", "mindsetconds"), mod_tpr_interaction_fullyadjusted))

eff_mod_tpr_interaction_fullyadjusted <- Effect(c("stressconds", "mindsetconds"), mod_tpr_interaction_fullyadjusted)
plot(eff_mod_tpr_interaction_fullyadjusted$fit[1:2], type = "n", lwd = 2, col = "orange",
     xlim = c(0.9, 2.1), ylim = c(0, 1), ylab = "TPR estimate (SD units)", frame.plot = F,
     xlab = "stress", xaxt = "n", yaxt = "n")
axis(1, at = c(1, 2), labels = c("-", "+"))
axis(2, at = c(0, 0.5, 1))
lines(x = c(0.95, 1.95), y =eff_mod_tpr_interaction_fullyadjusted$fit[1:2], lwd = 2, col = "orange")
lines(x = c(0.95, 0.95), y = c(eff_mod_tpr_interaction_fullyadjusted$lower[1], eff_mod_tpr_interaction_fullyadjusted$upper[1]), lwd = 2, col = "orange")
lines(x = c(1.95, 1.95), y = c(eff_mod_tpr_interaction_fullyadjusted$lower[2], eff_mod_tpr_interaction_fullyadjusted$upper[2]), lwd = 2, col = "orange")
lines(x = c(1.05, 2.05), y =eff_mod_tpr_interaction_fullyadjusted$fit[3:4], lwd = 2, col = "darkgreen")
lines(x = c(1.05, 1.05), y = c(eff_mod_tpr_interaction_fullyadjusted$lower[3], eff_mod_tpr_interaction_fullyadjusted$upper[3]), lwd = 2, col = "darkgreen")
lines(x = c(2.05, 2.05), y = c(eff_mod_tpr_interaction_fullyadjusted$lower[4], eff_mod_tpr_interaction_fullyadjusted$upper[4]), lwd = 2, col = "darkgreen")
legend("bottomleft", col = c("orange", "darkgreen"), lwd = 2, legend = c("mindset -", "mindset +"))

# Test heterogeneity of treatment effects: prior mindsets found to be predictive in study 2
mod_tpr_interaction_hte1 <- lmer(tpr_react_mc_z ~ stressconds * mindsetconds * bothmindsets.baseline +
                                            s1.sex + s1.age + racefactor + pss.baseline + 
                                            selfesteem.baseline +
                                            testanxiety.baseline + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_interaction_hte1)
plot(mod_tpr_interaction_hte1)
confint(mod_tpr_interaction_hte1)
plot(Effect(c("stressconds", "mindsetconds", "bothmindsets.baseline"), mod_tpr_interaction_hte1))



### Testing!

require(brms)

brm_tpr_factors_unadjusted <- brm(tpr_react_mc_z ~ stressonly + msonly + synergistic + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
plot(brm_tpr_factors_unadjusted)
summary(brm_tpr_factors_unadjusted)

draws_fit <- as_draws_array(brm_tpr_factors_unadjusted)

brm_tpr_interaction_unadjusted <- brm(tpr_react_mc_z ~ stressconds * mindsetconds + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
plot(brm_tpr_interaction_unadjusted)
summary(brm_tpr_interaction_unadjusted)
