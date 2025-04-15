# -------------------------------------------------------------------------
# Script adapted from Yeager et al (OSF: synergistic mindsets - Studies 1 to 6)
# Changes by Gustav Nilsonne (GN), Laura Caquelin (LC)
# All lines not marked with [Add by ...], [Modified by ...], or [Commented out by ...] are from the original script
# -------------------------------------------------------------------------

# [Add by GN] BEGIN:  Set working directory to source file directory ------------------------
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)
# [Add by GN] END ------------------------

# [Description added by GN] Require libraries
library(tidyr)
library(dplyr)
library(readr)
library(lme4)
# [Commented out by GN] Could not install from provided source "multibart_0.3.tar.gz" provided in the software folder.
# library(multibart) 
library(tidyverse)
library(effects)
require(brms)
library(lmerTest)
library(broom.mixed)

# [Description added by GN] Read data
dat_long <- read_csv("study4long.csv")

# [Description added by GN] Clean data and create variables for analyses

# [Modification by GN] BEGIN: Removed ------------------------
# dat_long_original = dat_long 
# dat_long_original$pnum = factor(dat_long_original$pnum)
# [Modification by GN] END ------------------------

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
# [Modification by GN] BEGIN: Removed ------------------------ 
#dat_long$sv_react <- dat_long$sv - dat_long$mean_base_sv
#dat_long$pep_react <- dat_long$pep - dat_long$mean_base_pep

### Subset data to the appropriate time in speech epoch - for TPR and SV
# dat_long = dat_long %>% filter(pid %in% goodid$pid, time==9|time==10|time==11|time==12|time==13) %>% 
#  group_by(pid) %>% summarize_all(mean, na.rm=T)
# [Modification by GN] END ------------------------

### Subset data to the appropriate time in recovery epoch - for PEP
#dat_long = dat_long %>% filter(pid %in% goodid$pid, time>13) %>%  # Just for PEP - subset to recovery
#  group_by(pid) %>% summarize_all(mean, na.rm=T)

# [Modification by GN] BEGIN: Removed ------------------------ 
#n = nrow(dat_long)
#t_train <- dat_long$time
# [Modification by GN] END ------------------------

# Correlations for Extended Data table
#cor(dat_long[which(dat_long$ctl==1), c("bothmindsets.baseline", "wellbeing", "task.appraisals")], use="pairwise.complete.obs")

# [Add by GN] BEGIN:  Make factors of continuous predictors that should be factors ------------------------
dat_long$stressconds <- as_factor(dat_long$stressconds)
dat_long$mindsetconds <- as_factor(dat_long$mindsetconds)
dat_long$stressonly <- as_factor(dat_long$stressonly)
dat_long$msonly <- as_factor(dat_long$msonly)
dat_long$synergistic <- as_factor(dat_long$synergistic)
dat_long$racefactor <- as_factor(dat_long$s1.race)
# [Add by GN] END ------------------------

# Set response variable
# [Modification by GN] BEGIN: Reversed the remove ------------------------
dat_long <- dat_long %>% filter(is.na(manipcheck)==F) # Toggle this depending on whether using the appraisal/well-being DVs
#dat_long <- dat_long %>% filter(is.na(wellbeing)==F) # Toggle this depending on whether using the appraisal/well-being DVs
# [Modification by GN] END ------------------------


#y = scale(dat_long$sv_react) 
#y = scale(dat_long$pep_react)
#y = scale(dat_long$tpr_react)
# [Modification by GN] BEGIN: Removed ------------------------
#y = scale(dat_long$task.appraisals)
#y = scale(dat_long$wellbeing)
# [Modification by GN] END ------------------------
#y = scale(dat_long$manipcheck)

# =======================================================================================
# [Replaced by GN and LC] The following section replaces the rest of the original script.
# A new analysis pipeline was implemented below, adapted to our replication
# ==============================

# Mean center outcome variable
dat_long$tpr_react_mc_z <- as.vector(scale(dat_long$tpr_react))

# Model 1: Unadjusted factors for TPR as the main outcome
mod_tpr_factors_unadjusted <- lmer(tpr_react_mc_z ~ stressonly + msonly + synergistic + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_factors_unadjusted)
confint(mod_tpr_factors_unadjusted)
confint(mod_tpr_factors_unadjusted, level = 0.98)
plot(mod_tpr_factors_unadjusted)
plot(Effect(c("stressonly"), mod_tpr_factors_unadjusted), ylim = c(-0.2, 1))
plot(Effect(c("msonly"), mod_tpr_factors_unadjusted), ylim = c(-0.2, 1))
plot(Effect(c("synergistic"), mod_tpr_factors_unadjusted), ylim = c(-0.2, 1))

# Model 2: Unadjusted interaction between stress and mindset
mod_tpr_interaction_unadjusted <- lmer(tpr_react_mc_z ~ stressconds * mindsetconds + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_interaction_unadjusted)
plot(mod_tpr_interaction_unadjusted)
plot(Effect(c("stressconds", "mindsetconds"), mod_tpr_interaction_unadjusted))

# Model 3: Fully adjusted model with additional covariates
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

# Model 4: Fully adjusted interaction model
mod_tpr_interaction_fullyadjusted <- lmer(tpr_react_mc_z ~ stressconds * mindsetconds +
                                           s1.sex + s1.age + racefactor + pss.baseline + 
                                           stressmindset.baseline + fixedmindset.baseline +
                                           bothmindsets.baseline + selfesteem.baseline +
                                           testanxiety.baseline + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_interaction_fullyadjusted)
plot(mod_tpr_interaction_fullyadjusted)
confint(mod_tpr_interaction_fullyadjusted)
plot(Effect(c("stressconds", "mindsetconds"), mod_tpr_interaction_fullyadjusted))

# Plotting effect size for interaction
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

### Plot mean TPR at baseline and after
## Mean tpr by condition and time
mean_tpr <- dat_long %>%
  filter(epoch %in% c("baseline", "speech")) %>% 
  group_by(cond.factor, epoch) %>%
  summarise(mean_TPR = mean(tpr, na.rm = TRUE),
            se_TPR = sd(tpr, na.rm = TRUE) / sqrt(n()))

## Plotting TPR by condition and time
plot1 <- ggplot(mean_tpr, aes(x = epoch, y = mean_TPR, color = cond.factor, group = cond.factor)) +
  geom_line(size = 2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_TPR - se_TPR, ymax = mean_TPR + se_TPR), width = 0.05) + 
  scale_x_discrete(labels = c("baseline" = "Baseline", "speech" = "Post")) + 
  scale_color_manual(values = c("#afafaf", "#f4a300", "#00407e",  "#a14040")) +  
  scale_fill_manual(values = c("#afafaf", "#f4a300", "#00407e", "#a14040")) +
  labs(x = "\nTime",
       y = "Mean TPR\n",
       color = "Condition", 
       fill = "Condition") +
  theme_minimal() + 
  theme(legend.position = "top",
        legend.box.just = "left",
        plot.margin = margin(20, 40, 10, 10),
        axis.title = element_text(size = 21),
        axis.text = element_text(size = 21),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        panel.grid = element_blank(), 
        axis.line = element_line(linewidth = 0.6),
        axis.ticks = element_line(size = 0.6))
  
ggsave(filename = "Plot_meanTPR_baselinevspost.png",plot = plot1, width = 8, height = 6, dpi = 300)

### Pairwise analysis
## Bayesian comparison
# Synergistic/Stress vs growth (Mindset)

dat_long <- dat_long %>%
  mutate(cond.factor = factor(cond.factor)) %>%
  mutate(cond.factor = relevel(cond.factor, ref = "Mindset"))

brm_tpr_factors_growth <- brm(tpr_react_mc_z ~ cond.factor + (1 | pid),
                              data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ],
                              iter = 4000, warmup = 2000)

#Addition of more iteration due to adjustment
brm_tpr_factors_growth_adjusted <- brm(tpr_react_mc_z ~ cond.factor + s1.sex + s1.age + racefactor + pss.baseline + 
                                         stressmindset.baseline + fixedmindset.baseline +
                                         bothmindsets.baseline + selfesteem.baseline +
                                         testanxiety.baseline + (1 | pid),
                                       data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ],
                                       iter = 6000, warmup = 3000)

summary(brm_tpr_factors_growth_adjusted)
resultsbrm_SynergisticStressvsgrowth <- posterior_summary(brm_tpr_factors_growth_adjusted)

# Synergistic vs Stress

dat_long <- dat_long %>% mutate(cond.factor = relevel(cond.factor, ref = "Stress"))

brm_tpr_factors_stress_adjusted <- brm(tpr_react_mc_z ~ cond.factor + s1.sex + s1.age + racefactor + pss.baseline + 
                                         stressmindset.baseline + fixedmindset.baseline +
                                         bothmindsets.baseline + selfesteem.baseline +
                                         testanxiety.baseline + (1 | pid),
                                       data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ],
                                       iter = 6000, warmup = 3000)

summary(brm_tpr_factors_stress_adjusted)
resultsbrm_Synergisticvsstress <- posterior_summary(brm_tpr_factors_stress_adjusted)

# Synergistic/Stress/Growth vs Control

dat_long <- dat_long %>% mutate(cond.factor = relevel(cond.factor, ref = "Control"))

brm_tpr_factors_control_adjusted <- brm(tpr_react_mc_z ~ cond.factor + s1.sex + s1.age + racefactor + pss.baseline + 
                                          stressmindset.baseline + fixedmindset.baseline +
                                          bothmindsets.baseline + selfesteem.baseline +
                                          testanxiety.baseline + (1 | pid),
                                        data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ],
                                        iter = 6000, warmup = 3000)

summary(brm_tpr_factors_control_adjusted)
resultsbrm_SynergisticStressGrowthvscontrol <- posterior_summary(brm_tpr_factors_control_adjusted)


# Table comparisons
final_table_bayesian   <- data.frame(
  Comparison = c("Synergistic vs growth", 
                 "Synergistic vs stress",
                 "Synergistic vs control",
                 "Stress vs control",
                 "Growth vs control", 
                 "Stress vs growth"),
  Estimate = c(round(resultsbrm_SynergisticStressvsgrowth["b_cond.factorSynergistic", "Estimate"], 2),
               round(resultsbrm_Synergisticvsstress["b_cond.factorSynergistic", "Estimate"], 2),
               round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorSynergistic", "Estimate"], 2),
               round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorStress", "Estimate"], 2),
               round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorMindset", "Estimate"], 2),
               round(resultsbrm_SynergisticStressvsgrowth["b_cond.factorStress", "Estimate"], 2)),
  CI_95 = c(paste0("[",round(resultsbrm_SynergisticStressvsgrowth["b_cond.factorSynergistic", "Q2.5"], 2), ";", round(resultsbrm_SynergisticStressvsgrowth["b_cond.factorSynergistic", "Q97.5"],2),"]"),
            paste0("[",round(resultsbrm_Synergisticvsstress["b_cond.factorSynergistic", "Q2.5"], 2), ";", round(resultsbrm_Synergisticvsstress["b_cond.factorSynergistic", "Q97.5"],2),"]"),
            paste0("[",round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorSynergistic", "Q2.5"], 2), ";", round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorSynergistic", "Q97.5"],2),"]"),
            paste0("[",round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorStress", "Q2.5"], 2), ";", round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorStress", "Q97.5"],2),"]"),
            paste0("[",round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorMindset", "Q2.5"], 2), ";", round(resultsbrm_SynergisticStressGrowthvscontrol["b_cond.factorMindset", "Q97.5"],2),"]"),
            paste0("[",round(resultsbrm_SynergisticStressvsgrowth["b_cond.factorStress", "Q2.5"], 2), ";", round(resultsbrm_SynergisticStressvsgrowth["b_cond.factorStress", "Q97.5"],2),"]")))

print(final_table_bayesian )

## Frequentist analysis
# Synergistic/Stress vs growth (Mindset)

dat_long <- dat_long %>% mutate(cond.factor = relevel(cond.factor, ref = "Mindset"))

lmer_tpr_factors_growth <- lmerTest::lmer(tpr_react_mc_z ~ cond.factor + s1.sex + s1.age + racefactor + pss.baseline + 
                                            stressmindset.baseline + fixedmindset.baseline +
                                            bothmindsets.baseline + selfesteem.baseline +
                                            testanxiety.baseline + (1 | pid),
                                          data = dat_long %>% filter(time %in% c(9, 10, 11, 12, 13)))

summary(lmer_tpr_factors_growth)

# Synergistic vs Stress

dat_long <- dat_long %>% mutate(cond.factor = relevel(cond.factor, ref = "Stress"))

lmer_tpr_factors_stress <- lmer(tpr_react_mc_z ~ cond.factor + s1.sex + s1.age + racefactor + pss.baseline + 
                                  stressmindset.baseline + fixedmindset.baseline +
                                  bothmindsets.baseline + selfesteem.baseline +
                                  testanxiety.baseline + (1 | pid),
                                data = dat_long %>% filter(time %in% c(9, 10, 11, 12, 13)))

summary(lmer_tpr_factors_stress)

# Synergistic/Stress/Growth vs Control

dat_long <- dat_long %>% mutate(cond.factor = relevel(cond.factor, ref = "Control"))

lmer_tpr_factors_control <- lmer(tpr_react_mc_z ~ cond.factor + s1.sex + s1.age + racefactor + pss.baseline + 
                                   stressmindset.baseline + fixedmindset.baseline +
                                   bothmindsets.baseline + selfesteem.baseline +
                                   testanxiety.baseline + (1 | pid),
                                 data = dat_long %>% filter(time %in% c(9, 10, 11, 12, 13)))

summary(lmer_tpr_factors_control)

# Table comparisons
# Function to extract results 

table_results <- function(model, comparison_label, target_comparison) {
  results <- tidy(model, conf.int = TRUE, conf.level = 0.95) %>%
    filter(term == target_comparison) %>% # Exclusion intercept
    mutate(
      Comparison = comparison_label,
      `Confidence interval` = paste0("[", round(conf.low, 2), "; ", round(conf.high, 2), "]"),
      Estimate = round(estimate, 2),
      pvalue = round(p.value, 3)) %>%
    select(Comparison, Estimate, `Confidence interval`, pvalue)
  
  return(results)
}

# Apply function
final_table_lmer <- bind_rows(
  table_results(lmer_tpr_factors_growth, "Synergistic vs Growth", "cond.factorSynergistic"),
  table_results(lmer_tpr_factors_stress, "Synergistic vs Stress", "cond.factorSynergistic"),
  table_results(lmer_tpr_factors_control, "Synergistic vs Control", "cond.factorSynergistic"),
  table_results(lmer_tpr_factors_control, "Stress vs Control", "cond.factorStress"),
  table_results(lmer_tpr_factors_control, "Growth vs Control", "cond.factorMindset"),
  table_results(lmer_tpr_factors_growth, "Stress vs Growth", "cond.factorStress"))

# Apply FDR correction to pvalue 
final_table_lmer <- final_table_lmer %>%
  mutate(p_adjusted = p.adjust(pvalue, method = "fdr"))

print(final_table_lmer)

### Interactions analysis
## Bayesian analysis
modbrm_tpr_interaction_fullyadjusted <- brm(tpr_react_mc_z ~ stressconds * mindsetconds +
                                            s1.sex + s1.age + racefactor + pss.baseline + 
                                            stressmindset.baseline + fixedmindset.baseline +
                                            bothmindsets.baseline + selfesteem.baseline +
                                            testanxiety.baseline + (1 | pid), 
                                         data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ],
                                         iter = 6000, warmup = 3000)
summary(modbrm_tpr_interaction_fullyadjusted)

resultsbrm_interaction <- posterior_summary(modbrm_tpr_interaction_fullyadjusted)

final_table_bayesian   <- data.frame(
  Comparison = c("Stress, main effect", 
                 "Growth, main effect",
                 "Interaction"),
  Estimate = c(round(resultsbrm_interaction["b_stressconds1", "Estimate"], 2),
               round(resultsbrm_interaction["b_mindsetconds1", "Estimate"], 2),
               round(resultsbrm_interaction["b_stressconds1:mindsetconds1", "Estimate"], 2)),
  CI_95 = c(paste0("[",round(resultsbrm_interaction["b_stressconds1", "Q2.5"], 2), ";", round(resultsbrm_interaction["b_stressconds1", "Q97.5"],2),"]"),
            paste0("[",round(resultsbrm_interaction["b_mindsetconds1", "Q2.5"], 2), ";", round(resultsbrm_interaction["b_mindsetconds1", "Q97.5"],2),"]"),
            paste0("[",round(resultsbrm_interaction["b_stressconds1:mindsetconds1", "Q2.5"], 2), ";", round(resultsbrm_interaction["b_stressconds1:mindsetconds1", "Q97.5"],2),"]")))

print(final_table_bayesian )

## Frequentist analysis

mod_tpr_interaction_fullyadjusted <- lmer(tpr_react_mc_z ~ stressconds * mindsetconds +
                                            s1.sex + s1.age + racefactor + pss.baseline + 
                                            stressmindset.baseline + fixedmindset.baseline +
                                            bothmindsets.baseline + selfesteem.baseline +
                                            testanxiety.baseline + (1 | pid), data = dat_long[dat_long$time %in% c(9, 10, 11, 12, 13), ])
summary(mod_tpr_interaction_fullyadjusted)

final_table_lmer_interaction <- bind_rows(
  table_results(mod_tpr_interaction_fullyadjusted, "Stress, main effect", "stressconds1"),
  table_results(mod_tpr_interaction_fullyadjusted, "Growth, main effect", "mindsetconds1"),
  table_results(mod_tpr_interaction_fullyadjusted, "Interaction", "stressconds1:mindsetconds1"))

print(final_table_lmer_interaction)

# =======================================================================================