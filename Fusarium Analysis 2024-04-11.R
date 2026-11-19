##### Housekeeping ####
# Load packages
library(readxl) 
library(tidyverse)
library(ggh4x)
library(tables)
library(nlme)
library(glmmTMB)
library(ordinal)
library(DHARMa)
library(ggResidpanel)
library(emmeans)
library(marginaleffects)
library(multcomp)
library(writexl)
library(rstatix)

# Set global options
#setwd("C:/Users/Simon/OneDrive - University of Florida/IFAS SCU/Silva de Moura, Vitor - IFAS SCU/Fusarium Wilt in Apiaceae/Alternative Analysis")
my_pal = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
options(ggplot2.discrete.colour = my_pal,
        ggplot2.discrete.fill = my_pal, 
        dplyr.width = Inf, pillar.print_max = 100, pillar.print_min = 100)
my_theme = theme_light() +
          theme(text = element_text(size = 14, color = 'black'), 
                strip.background = element_rect(color = 'black'), 
                strip.text = element_text(color = 'black'))
#theme_set(my_theme)
theme_set(theme_bw(base_size=14) + theme(strip.background = element_blank()))


# Define custom functions
Min = function(x){min(x, na.rm = T)}
Mean = function(x){mean(x, na.rm = T)}
Max = function(x){max(x, na.rm = T)}
N = function(x){sum(!is.na(x))}

# Read and prep data
fusarium = read.csv('Paper data_ Sep23_Fusarium experiment evaluation.csv') |> 
  mutate(across(c(1:6, 13:15),factor), 
         Discoloration_ord = ordered(ifelse(Rate_discoloration == 0, 1, Rate_discoloration)))
str(fusarium)

fusarium$Treatment <- factor(fusarium$Treatment, levels = c("KX31-inoculated", "Control"))

#### Mortality Analysis ####
# A zero-inflated gamma model separately estimates the probability of death
# as well as the effect on e.g. total weight among the living. Here
# Total_weight is used, but focus is on the zero-inflated part for now.
# Note that due to bad luck/insufficient data, it is not possible to
# estimate interaction effect of species or varieties and fusarium
# innoculation on mortality rate. This somewhat limits the estimates and
# comparisons we can make (for example, cannot compare innocated vs
# non-inoculated on mortality rate separately among varieties, because most
# varieties didn't exhibit any mortality among the non-inoculated). Also
# note the warnings can be ignored, they arise because Variety is nested
# within species
fusarium_mort_mod = glmmTMB(Total_weight ~ Block+Treatment*Plant*Variety, data = fusarium, 
                            family = ziGamma(link = 'inverse'),
                            ziformula = ~ 0+Treatment:Variety)

# Check model fit - when checking the residuals from a generalized linear model
# we need to use "simulated quantile residuals". By default these are made
# to follow a uniform distribution but I prefer to put them on the normal, so that
# they are more familiar/easier to diagnose. The code below does that, then
# plots the diagnostic plots
# Check model fit: couple of outliers, but not too bad.
fusarium_mort_res = simulateResiduals(fusarium_mort_mod) |> 
  residuals(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_mort_fit = fitted(fusarium_mort_mod)
resid_auxpanel(fusarium_mort_res, fusarium_mort_fit)

# to learn more about residuals for GLMs, see:
vignette('DHARMa', package = 'DHARMa')

# Calculate marginal means (aka least square means, aka adjusted means) for
# the zero inflation component. Using the "at" argument, we are specifying
# that we only want estimates for the inoculated treatments, since due to
# complete separation we don't have reliable estimates for some fo the
# uninoculated (no mortality, so incomparible)
(fusarium_mort_all_emm = emmeans(fusarium_mort_mod, ~ Treatment:Variety, 
                                     component = 'zi', type = 'response', 
                                     at = list(Treatment = 'KX31-inoculated')))

# To learn more about them, see the introduction vignette to the emmeans
# package, especially the section "Illustration: pigs experiment"
vignette('basics', package = 'emmeans')

# Just for illustration, here is the output for both treated and untreated,
# and we can see the two zero estimates for Mock C26 Celery and Mock
# Italian Parsly
emmeans(fusarium_mort_mod, ~ Treatment:Variety, component = 'zi', type = 'response')

# First, we want to test if there is an difference, on average, between
# parsley and celery. To do this, we need to define a contrast, which is a
# set of weights which define a comparison made of weighted estimates (the
# ordering of the weights should match the output from emmeans above)
fusarium_mort_plant_contr_list =  list('Parsley - Celery' = c(-0.5, -0.5, 0.5, 0.5))

# Mortality is higher in parsley than in celary 
(fusarium_mort_plant_contr = contrast(fusarium_mort_all_emm, fusarium_mort_plant_contr_list))

# Here we estimate the probability of death among inoculated parsley
# varieties. Risk ratio of 0.917/0.304 = 3.02 for the Curly variety...
(fusarium_mort_parsley_emm = emmeans(fusarium_mort_mod, ~ Treatment:Variety, 
                                     component = 'zi', type = 'response', 
                                     at = list(Variety = c('Curly', 'Italian'), 
                                               Treatment = 'KX31-inoculated')))
# ... and the difference is statistically significant (here it is reported
# as odds.ratio, which is closely related to risk ratio but the
# interpretation is a bit less intuitive)
(fusarium_mort_parsley_cld = cld(fusarium_mort_parsley_emm, reversed = T, Letters = letters))

# But no difference in the mortality rates among the celery varieties
(fusarium_mort_celery_emm = emmeans(fusarium_mort_mod, ~ Treatment:Variety, 
                                     component = 'zi', type = 'response', 
                                     at = list(Variety = c('C26', 'C42'), 
                                               Treatment = 'KX31-inoculated')))
# ... and the difference is statistically significant (here it is reported
# as odds.ratio, which is closely related to risk ratio but the
# interpretation is a bit less intuitive)
(fusarium_mort_celery_cld = cld(fusarium_mort_celery_emm, reversed = T, Letters = letters))

write_xlsx(x = list('Spp Avg Mort when Innoc' = as.data.frame(fusarium_mort_plant_contr), 
                    'Parsley Mortality when Innoc' = as.data.frame(fusarium_mort_parsley_cld), 
                    'Celery Mortality when Innoc' = as.data.frame(fusarium_mort_celery_cld)), 
           path = 'results/Mortality Analysis.xlsx')

#### Total Weight Analysis ####
# This uses the same model as above, the only difference is that we are now
# interested in the conditional mean weight (i.e. the fixed effects part,
# not the zero-inflation part)
#fusarium_ttlwt_mod = glmmTMB(Total_weight ~ Block + Treatment*Plant*Variety, data = fusarium, 
                            #family = ziGamma(link = 'inverse'),
                            #ziformula = ~ 0+Treatment:Variety)

##ziformula: modeling 0 and non 0 values

fusarium_ttlwt_mod = glmmTMB(Total_weight ~ Treatment*Plant*Variety + (1|Block), data = fusarium, 
                             family = tweedie(link = 'log'), 
                             dispformula = ~ Treatment:Variety)

##dispformula: model for dispersion (similar to varience)

# Model Checking
fusarium_ttlwt_res = simulateResiduals(fusarium_ttlwt_mod) |> 
  residuals(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_ttlwt_fit = fitted(fusarium_ttlwt_mod)
resid_auxpanel(fusarium_ttlwt_res, fusarium_ttlwt_fit)


# F-test indicates the story is in the treatment x species interaction
joint_tests(fusarium_ttlwt_mod, component = 'cond')

# Calculate marginal means. Note that component = 'response' means we are
# using estimates which include 0 weights (i.e. dead plants). If we want
# estimates only from among the surviving plants, use component = 'cmean'
# instead. Also note that in the formula, we are using Varieties (nested in
# plant) as a grouping variable via the "|". This means that when we
# subsequently perform any tests or comparisons, they are only made within
# each level of the grouping factor not between them.
(fusarium_ttlwt_emm = emmeans(fusarium_ttlwt_mod, ~ Treatment|Plant:Variety, type = 'response'))
(fusarium_ttlwt_cld = cld(fusarium_ttlwt_emm, reversed = T, Letters = letters))

(fusarium_ttlwt_contr = contrast(fusarium_ttlwt_emm, 'pairwise'))
fusarium_ttlwt_contr <- as.data.frame(fusarium_ttlwt_contr)

fusarium_ttlwt_contr$Treatment <- as.factor(paste0("Control"))
fusarium_ttlwt_contr$group2 <- as.factor(paste0("KX31-inoculated"))

fusarium_ttlwt_contr[fusarium_ttlwt_contr$p.value<=0.001,"p.adj.signif"]<-"***"
fusarium_ttlwt_contr[fusarium_ttlwt_contr$p.value>0.001 & fusarium_ttlwt_contr$p.value<=0.01,"p.adj.signif"]<-"**"
fusarium_ttlwt_contr[fusarium_ttlwt_contr$p.value>0.01 & fusarium_ttlwt_contr$p.value<=0.05,"p.adj.signif"]<-"*"
fusarium_ttlwt_contr[ fusarium_ttlwt_contr$p.value>0.05,"p.adj.signif"]<-"ns"

fusarium_ttlwt_contr <- fusarium_ttlwt_contr %>%
  mutate(y.position=c(39, 64, 28, 30))

fusarium_ttlwt_cld$Treatment <- factor(fusarium_ttlwt_cld$Treatment, levels = c("KX31-inoculated", "Control"))

# Plot the results
(total_plot <- fusarium_ttlwt_cld |> 
  as.data.frame() |> 
  #mutate(.group = gsub(' ', '', .group)) |> 
  ggplot(aes(x = Treatment, fill = Treatment)) +
  facet_grid(cols = vars(Variety),switch="both")+
  geom_col(aes(y = response)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = .2) +
  #geom_text(aes(y = asymp.UCL + 3, label = .group), size = 5) +
  theme_Publication() +
  labs(y = "Total plant weight (g)",x = "")  +
  scale_y_continuous(breaks = round(seq(0, 80, by = 10),1)) +
  guides(fill=guide_legend(title="Treatments")) +
  scale_fill_manual(values = c("#D3D3D3", "#000000"))+
  theme(axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        strip.background = element_blank(),strip.text.x = element_blank())+
        #strip.placement = "outside")+
  scale_x_discrete(label=c("", "","","", 
                           "", "","", ""))+
  stat_pvalue_manual(fusarium_ttlwt_contr))
total_plot

#### Root Weight Analysis####
fusarium_rootwt_mod = glmmTMB(Roots_weight ~ Treatment*(Plant/Variety)+ (1|Block), data = fusarium, 
                              family = tweedie(link = 'log'), 
                              dispformula = ~ Treatment)

print(fusarium_rootwt_mod)

# Check model fit: couple of outliers, but not too bad.
fusarium_root_res = simulateResiduals(fusarium_rootwt_mod) |> 
  residuals(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_root_fit = fitted(fusarium_rootwt_mod)
resid_auxpanel(fusarium_root_res, fusarium_root_fit)

# F-test indicates the story is in the treatment x species interaction
joint_tests(fusarium_rootwt_mod, nesting = NULL)

# Calculate marginal means. Note that component = 'response' means we are
# using estimates which include 0 weights (i.e. dead plants). If we want
# estimates only from among the surviving plants, use component = 'cmean'
# instead. Also note that in the formula, we are using Varieties (nested in
# plant) as a grouping variable via the "|". This means that when we
# subsequently perform any tests or comparisons, they are only made within
# each level of the grouping factor not between them.
(fusarium_rootwt_emm = emmeans(fusarium_rootwt_mod, ~ Treatment|Plant:Variety, 
                               type = 'response'))
(fusarium_rootwt_cld = cld(fusarium_rootwt_emm, reversed = T, Letters = letters))


fusarium_rootwt_cld[fusarium_rootwt_cld$Variety == 'Curly' &
                      fusarium_rootwt_cld$Treatment == 'KX31-inoculated', 
                    c('SE', 'asymp.LCL', 'asymp.UCL')] = 0 

fusarium_rootwt_cld[fusarium_rootwt_cld$Variety == 'Curly', 
                    c('.group')] = '' 


(fusarium_rootwt_contr = contrast(fusarium_rootwt_emm, 'pairwise'))
fusarium_rootwt_contr <- as.data.frame(fusarium_rootwt_contr)

fusarium_rootwt_contr$Treatment <- as.factor(paste0("Control"))
fusarium_rootwt_contr$group1 <- as.factor(paste0("Control"))
fusarium_rootwt_contr$group2 <- as.factor(paste0("KX31-inoculated"))

fusarium_rootwt_contr[fusarium_rootwt_contr$p.value<=0.001,"p.adj.signif"]<-"***"
fusarium_rootwt_contr[fusarium_rootwt_contr$p.value>0.001 & fusarium_rootwt_contr$p.value<=0.01,"p.adj.signif"]<-"**"
fusarium_rootwt_contr[fusarium_rootwt_contr$p.value>0.01 & fusarium_rootwt_contr$p.value<=0.05,"p.adj.signif"]<-"*"
fusarium_rootwt_contr[ fusarium_rootwt_contr$p.value>0.05,"p.adj.signif"]<-"*qc"

fusarium_rootwt_contr <- fusarium_rootwt_contr %>%
  mutate(y.position=c(17, 28, 9, 13))
root_plot<-root_plot 

# Plot the results
(root_plot  <- fusarium_rootwt_cld |> 
  as.data.frame() |> 
  mutate(.group = gsub(' ', '', .group)) |> 
  ggplot(aes(x = Treatment, fill = Treatment)) +
  facet_grid(cols = vars(Variety),switch="both")+
  geom_col(aes(y = response)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = .2) +
  #geom_text(aes(y = asymp.UCL + 3, label = .group), size = 5) +
  theme_Publication() +
  labs(y = "Root weight (g)",x = "")  +
  scale_y_continuous(breaks = round(seq(0, 80, by = 10),1)) +
  guides(fill=guide_legend(title="Treatments")) +
  scale_fill_manual(values = c("#D3D3D3", "#000000"))+
  theme(axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.2,0.1,1), 'lines'),
        strip.background = element_blank(),strip.text.x = element_blank())+
  #strip.placement = "outside")+
  scale_x_discrete(label=c("", "","","", 
                           "", "","", "")) + 
  stat_pvalue_manual(fusarium_rootwt_contr))
root_plot

#### Height Analysis####

fusarium_height_mod = glmmTMB(height_final ~ Treatment*Plant*Variety+ (1|Block), data = fusarium, 
                              family = tweedie(link = 'log'), 
                              dispformula = ~ Treatment:Variety)

# Check model fit: couple of outliers, but not too bad.
fusarium_height_res = simulateResiduals(fusarium_height_mod) |> 
  residuals(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_height_fit = fitted(fusarium_height_mod)
resid_auxpanel(fusarium_height_res, fusarium_height_fit)

# F-test indicates the story is in the treatment x species interaction
joint_tests(fusarium_height_mod, component = 'cond')

# Calculate marginal means. Note that type = 'response' means that the
# values are "inverse linked" (i.e. back-transformed) from the link (= log)
# scale to the original scale. Also note that in the formula, we are using
# Varieties (nested in plant) as a grouping variable via the "|". This
# means that when we subsequently perform any tests or comparisons, they
# are only made within each level of the grouping factor not between them.
(fusarium_height_emm = emmeans(fusarium_height_mod, ~ Treatment|Plant:Variety, 
                               type = 'response'))
(fusarium_height_cld = cld(fusarium_height_emm, reversed = T, Letters = letters))

(fusarium_height_contr = contrast(fusarium_height_emm, 'pairwise'))
fusarium_height_contr <- as.data.frame(fusarium_height_contr)

fusarium_height_contr$Treatment <- as.factor(paste0("Control"))
fusarium_height_contr$group1 <- as.factor(paste0("Control"))
fusarium_height_contr$group2 <- as.factor(paste0("KX31-inoculated"))

fusarium_height_contr[fusarium_height_contr$p.value<=0.001,"p.adj.signif"]<-"***"
fusarium_height_contr[fusarium_height_contr$p.value>0.001 & fusarium_height_contr$p.value<=0.01,"p.adj.signif"]<-"**"
fusarium_height_contr[fusarium_height_contr$p.value>0.01 & fusarium_height_contr$p.value<=0.05,"p.adj.signif"]<-"*"
fusarium_height_contr[ fusarium_height_contr$p.value>0.05,"p.adj.signif"]<-"qc"


fusarium_height_contr <- fusarium_height_contr %>%
  mutate(y.position=c(23, 26, 20, 22))
height_plot <- height_plot 


# Plot the results
(height_plot <- fusarium_height_cld |> 
  as.data.frame() |> 
  mutate(.group = gsub(' ', '', .group)) |> 
  ggplot(aes(x = Treatment, fill = Treatment)) +
  facet_nested(~ Plant + Variety, switch="both") +  
  geom_col(aes(y = response)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = .2) +
  #geom_text(aes(y = asymp.UCL + 3, label = .group), size = 5) +
  theme_Publication() +
  labs(y = "Plant height (cm)",x = "")  +
  scale_y_continuous(breaks = round(seq(0, 30, by = 10),1)) +
  guides(fill=guide_legend(title="Treatments")) +
  scale_fill_manual(values = c("#D3D3D3", "#000000"))+
  theme(axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        strip.background = element_blank(), strip.text.x = element_blank())+
  scale_x_discrete(label=c("", "","","", 
                           "", "","", "")) +
  stat_pvalue_manual(fusarium_height_contr))

ggsave("height_41.svg")

#### Height Analysis - 41 dai####

fusarium_height_mod = glmmTMB(height_middle ~ Treatment*Plant*Variety+ (1|Block), data = fusarium, 
                              family = tweedie(link = 'log'), 
                              dispformula = ~ Treatment:Variety)

# Check model fit: couple of outliers, but not too bad.
fusarium_height_res = simulateResiduals(fusarium_height_mod) |> 
  residuals(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_height_fit = fitted(fusarium_height_mod)
resid_auxpanel(fusarium_height_res, fusarium_height_fit)

# F-test indicates the story is in the treatment x species interaction
joint_tests(fusarium_height_mod, component = 'cond')

# Calculate marginal means. Note that type = 'response' means that the
# values are "inverse linked" (i.e. back-transformed) from the link (= log)
# scale to the original scale. Also note that in the formula, we are using
# Varieties (nested in plant) as a grouping variable via the "|". This
# means that when we subsequently perform any tests or comparisons, they
# are only made within each level of the grouping factor not between them.
(fusarium_height_emm = emmeans(fusarium_height_mod, ~ Treatment|Plant:Variety, 
                               type = 'response'))
(fusarium_height_cld = cld(fusarium_height_emm, reversed = T, Letters = letters))

(fusarium_height_contr = contrast(fusarium_height_emm, 'pairwise'))
fusarium_height_contr <- as.data.frame(fusarium_height_contr)

fusarium_height_contr$Treatment <- as.factor(paste0("Control"))
fusarium_height_contr$group1 <- as.factor(paste0("Control"))
fusarium_height_contr$group2 <- as.factor(paste0("KX31-inoculated"))

fusarium_height_contr[fusarium_height_contr$p.value<=0.001,"p.adj.signif"]<-"***"
fusarium_height_contr[fusarium_height_contr$p.value>0.001 & fusarium_height_contr$p.value<=0.01,"p.adj.signif"]<-"**"
fusarium_height_contr[fusarium_height_contr$p.value>0.01 & fusarium_height_contr$p.value<=0.05,"p.adj.signif"]<-"*"
fusarium_height_contr[ fusarium_height_contr$p.value>0.05,"p.adj.signif"]<-"qc"


fusarium_height_contr <- fusarium_height_contr %>%
  mutate(y.position=c(25, 29, 18, 23))

# Plot the results
(heightM_plot <- fusarium_height_cld |> 
    as.data.frame() |> 
    mutate(.group = gsub(' ', '', .group)) |> 
    ggplot(aes(x = Treatment, fill = Treatment)) +
    facet_grid(cols = vars(Variety),switch="both")+
    geom_col(aes(y = response)) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                  width = .2) +
    #geom_text(aes(y = asymp.UCL + 3, label = .group), size = 5) +
    theme_Publication() +
    labs(y = "Plant height (cm)",x = "")  +
    scale_y_continuous(breaks = round(seq(0, 30, by = 10),1)) +
    guides(fill=guide_legend(title="Treatments")) +
    scale_fill_manual(values = c("#D3D3D3", "#000000"))+
    theme(axis.ticks.x = element_blank(),
          plot.margin = unit(c(0,0.2,0,1), 'lines'),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          strip.placement = "outside")+
    scale_x_discrete(label=c("", "","","", 
                             "", "","", "")) +
    stat_pvalue_manual(fusarium_height_contr))


#### Discoloration Analysis ####
# This analysis will use the original percentage data (scaled to be a
# proportion by dividing by 100), and model the data as coming from the
# "ordinal beta distribution". The regular beta distribution is for
# modeling continuous proportions/indices bounded between 0 and 1, but
# doesn't allow observations of exactly 0 or 1. The ordinal beta
# distribution does allow observations of exactly 0 and 1, which is why
# we're using that one. Again, ignore the warnings because they are just
# because the variety is nested in species
fusarium_discol_mod = glmmTMB(Perc_discoloration/100 ~ Treatment*(Plant/Variety) + (1|Block), data = fusarium, 
                              family = ordbeta(link = 'logit'))

# Check the residuals - not perfect, but shold be ok
fusarium_discol_res = simulateResiduals(fusarium_discol_mod) |> 
  residuals(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_discol_fit = fitted(fusarium_discol_mod)
resid_auxpanel(fusarium_discol_res, fusarium_discol_fit)

# F-test: Treatment:Plant:Variety interaction is significant
joint_tests(fusarium_discol_mod)

# Calculate marginal means and perform comparisons, again using variety as
# a grouping variable since what we want is to test treated vs control
# within each variety
(fusarium_discol_emm = emmeans(fusarium_discol_mod, ~ Treatment|Plant:Variety, type = 'response'))
(fusarium_discol_cld = cld(fusarium_discol_emm, reversed = T, Letters = letters))

# Uninoculated C42 had no discoloration, causing a phenomenon called
# "quasi-complete separation" (its most famous in the context of logistic
# regression, but can happen in many types of GLMs). Basically it renders
# that estimate incomparible with other estimates, so I'm removing the SE, CI 
# letter code
fusarium_discol_cld[fusarium_discol_cld$Variety == 'C42' &
                      fusarium_discol_cld$Treatment == 'Control', 
                    c('SE', 'asymp.LCL', 'asymp.UCL')] = 0 

fusarium_discol_cld[fusarium_discol_cld$Variety == 'C42', 
                    c('.group')] = '' 

(fusarium_discol_contr = contrast(fusarium_discol_emm, 'pairwise'))
fusarium_discol_contr <- as.data.frame(fusarium_discol_contr)

fusarium_discol_contr$Treatment <- as.factor(paste0("Control"))
fusarium_discol_contr$group1 <- as.factor(paste0("Control"))
fusarium_discol_contr$group2 <- as.factor(paste0("KX31-inoculated"))

fusarium_discol_contr[fusarium_discol_contr$p.value<=0.001,"p.adj.signif"]<-"***"
fusarium_discol_contr[fusarium_discol_contr$p.value>0.001 & fusarium_discol_contr$p.value<=0.01,"p.adj.signif"]<-"**"
fusarium_discol_contr[fusarium_discol_contr$p.value>0.01 & fusarium_discol_contr$p.value<=0.05,"p.adj.signif"]<-"*"
fusarium_discol_contr[ fusarium_discol_contr$p.value>0.05,"p.adj.signif"]<-"ns"

fusarium_discol_contr[fusarium_discol_contr$Variety == 'C42', 
                    c('p.adj.signif')] = '*qc' 

fusarium_discol_contr <- fusarium_discol_contr %>%
  mutate(y.position=c(0.77, 0.79, 0.39, 0.30))

# plot results
(discol_plot <- fusarium_discol_cld |> 
  as.data.frame() |> 
  ggplot(aes(x = Treatment, fill = Treatment))+
  facet_nested(~ Plant + Variety, switch="both")+
  geom_col(aes(y = response))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = .2)+
  #geom_text(aes(y = asymp.UCL+.1, label = .group), size = 5)+
  theme_Publication() +
  labs(y = "Vascular discoloration (%)",x = "")  +
  scale_y_continuous(breaks = round(seq(0, 1, by = 0.1),1), labels = scales::percent) +
  guides(fill=guide_legend(title="Treatments")) +
  scale_fill_manual(values = c("#D3D3D3", "#000000"))+
  theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(),
        plot.margin = unit(c(0,0.2,0,1), 'lines'),
        strip.background = element_blank(),
        #strip.text.x = element_blank(), 
        strip.placement = "outside")+
  stat_pvalue_manual(fusarium_discol_contr))

ggsave("discoloration.svg")
ggsave("discoloration.tiff")

#### Discoloration vs. Total Weight ####
# Tweedie distribution solves a lot of the problems which were being encountered
# previously

# New and improved Tweedie version:
fusarium_colwt_mod = glmmTMB(Total_weight ~ (Plant/Variety)*Perc_discoloration, 
                             data = fusarium, family = tweedie(link = 'log'), 
                             dispformula = ~ Treatment)

# Not great, but this model is a challenging case
fusarium_colwt_res = simulateResiduals(fusarium_colwt_mod) |> 
  resid(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_colwt_fit = fitted(fusarium_colwt_mod)
resid_auxpanel(fusarium_colwt_res, fusarium_colwt_fit)

# McFadden's Psuedo-R2 and Adjusted Psuedo-R2
nullmod = update(fusarium_colwt_mod, formula. = .~1) # first get null model
# R2:
1 - (logLik(fusarium_colwt_mod)/logLik(nullmod)) 
# Adjusted R2:
1 - (logLik(fusarium_colwt_mod) - length(unlist(fixef(fusarium_colwt_mod))))/logLik(nullmod)

# If Perc_discoloration is significant, then there is a significant
# association between discoloration % and ttl. weight, If interaction terms
# Plant:Perc_Discoloration or Plant:Variety:Perc_discoloration are
# significant, then the associations between discoloration % and total
# weight are species- or variety-specific, where some may be significant
# and others not, or all significant but different from one another. Need
# emtrends() to find out what the story is
joint_tests(fusarium_colwt_mod)

# estimated slopes (on log scale because of the link function used)
(fusarium_colwt_emt = emtrends(fusarium_colwt_mod, ~ Plant, 
                               var = 'Perc_discoloration'))

# We can use the delta method to convert it to a more interpretable measure
# on the original scale: since we're using the log link, effects are
# multiplicative and we can express the "slope" as percent decline in total
# weight per 1% increase in discoloration. First we extract the model
# estimates and and variance-covariance matrix from our fitted model, then
# feed those and the appropriate formula into the detaMethod function:
# 
fe = fixef(fusarium_colwt_mod)$cond
keep = !is.na(fe) # some missing coefficients have NA values, we want to remove those
fe = fe[keep]
vc = vcov(fusarium_colwt_mod)$cond
vc = vc[keep, keep]

# Here are the estimates (along with SE and confidence limits) for percent
# decline (because negative) associated with 1% increase in discoloration.
# Note that you have to be a bit careful with the specification of the
# formula, because it depends on what contrast specifications R is using
# (this works for the default)
car::deltaMethod(object = fe,
                 g. = '100*(exp(Perc_discoloration)-1)',  # for celery, because that is the "reference level" for the plant factor
                 vcov = vc)
car::deltaMethod(object = fe,
                 g. = '100*(exp(Perc_discoloration + `PlantParsley:Perc_discoloration`)-1)', # for parsley, because we've added to the "reference level" the parley-effect
                 vcov = vc)

# Back to emtrends object to compare those "slopes"
contrast(fusarium_colwt_emt, 'pairwise') 

# generate predictions for each species over the 0-100% range of discoloration
fusarium_colwt_preds = emmeans(fusarium_colwt_mod, ~ Perc_discoloration|Plant, 
                               type = 'response', 
                               at = list(Perc_discoloration = 0:100))
# plot the results
(colwt_plot<-fusarium_colwt_preds |> 
  as.data.frame() |> 
  #mutate(Plant = ifelse(Variety %in% c('C26', 'C42'), 'Celery', 'Parsley')) |> 
  ggplot(aes(x = Perc_discoloration))+
  facet_grid(cols = vars(Plant))+
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), color = NA, alpha = .2)+
  geom_line(aes(y = response))+
  geom_point(aes(y = Total_weight), data = fusarium, alpha = .5) +
  #scale_y_continuous('Total Weight')+
  theme_Publication() +
  labs(x = "",y = "Total plant weight (g)") +
  theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(),
        strip.background = element_blank(),strip.text.x = element_blank(),
        plot.margin = unit(c(0,0.2,0.8,1), 'lines'))+
  xlim(0, 100)) 


colwt_plot

#### Discoloration vs. Root Weight ####
# This looks mostly reasonable to me - try the same model for height, I look
# again if the results are weird.
fusarium_colrt_mod = glmmTMB(Roots_weight ~ Plant*Perc_discoloration, 
                             data = fusarium, family = ziGamma(link = 'log'), 
                             ziformula = ~ Plant:Treatment)

fusarium_colrt_mod = glmmTMB(Roots_weight ~ (Plant/Variety)*Perc_discoloration, 
                             data = fusarium, family = tweedie(link = 'log'), 
                             dispformula = ~ Treatment)
joint_tests(fusarium_colrt_mod)

# Not great, but this model is a challenging case
fusarium_colrt_res = simulateResiduals(fusarium_colrt_mod) |> 
  resid(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_colrt_fit = fitted(fusarium_colrt_mod)
resid_auxpanel(fusarium_colrt_res, fusarium_colrt_fit)

# estimated slopes (on log scale because of the link function used)
(fusarium_colrt_emt = emtrends(fusarium_colrt_mod, ~ Variety|Plant, 
                               var = 'Perc_discoloration'))
# the difference is non-significant
contrast(fusarium_colrt_emt, 'pairwise') 

# generate predictions for each species over the 0-100% range of discoloration
fusarium_colrt_preds = emmeans(fusarium_colrt_mod, ~ Perc_discoloration|Plant,
                               type = 'response', 
                               at = list(Perc_discoloration = 0:100))
# plot the results
(colrt_plot <-fusarium_colrt_preds |> 
  as.data.frame() |> 
  ggplot(aes(x = Perc_discoloration))+
  facet_grid(cols = vars(Plant),switch="both") +  
  geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), color = NA, alpha = .2)+
  geom_line(aes(y = response))+
  geom_point(aes(y = Roots_weight), data = fusarium, alpha = .5) +
  theme_Publication() +
  labs(x = "Vascular discoloration (%)",y = "Root weight (g)") +
  theme(strip.background = element_blank(),strip.placement = "outside",
        strip.text.x = element_blank(),
        plot.margin = unit(c(0,0.2,0,1), 'lines'))+
  xlim(0, 100))

theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(),
strip.background = element_blank(),strip.text.x = element_blank(),
plot.margin = unit(c(0,0.2,0.8,1), 'lines'))

# We can use the delta method to convert it to a more interpretable measure
# on the original scale: since we're using the log link, effects are
# multiplicative and we can express the "slope" as percent decline in total
# weight per 1% increase in discoloration. First we extract the model
# estimates and and variance-covariance matrix from our fitted model, then
# feed those and the appropriate formula into the detaMethod function:
# 
fe = fixef(fusarium_colrt_mod)$cond
keep = !is.na(fe) # some missing coefficients have NA values, we want to remove those
fe = fe[keep]
vc = vcov(fusarium_colrt_mod)$cond
vc = vc[keep, keep]

# Here are the estimates (along with SE and confidence limits) for percent
# decline (because negative) associated with 1% increase in discoloration.
# Note that you have to be a bit careful with the specification of the
# formula, because it depends on what contrast specifications R is using
# (this works for the default)
car::deltaMethod(object = fe,
                 g. = '100*(exp(Perc_discoloration)-1)',  # for celery, because that is the "reference level" for the plant factor
                 vcov = vc)
car::deltaMethod(object = fe,
                 g. = '100*(exp(Perc_discoloration + `PlantParsley:Perc_discoloration`)-1)', # for parsley, because we've added to the "reference level" the parley-effect
                 vcov = vc)

#### Discoloration vs. Height ####
# Although the tweedie does help here, still a challenging model
fusarium_colht_mod = glmmTMB(height_final ~ (Plant/Variety)*Perc_discoloration,  
                             data = fusarium, family = tweedie(link = 'log'), 
                             dispformula = ~ Treatment:Variety)
joint_tests(fusarium_colht_mod)

# Not great, but this model is a challenging case
fusarium_colht_res = simulateResiduals(fusarium_colht_mod) |> 
  resid(quantileFunction = qnorm, outlierValues = c(-5, 5))
fusarium_colht_fit = fitted(fusarium_colht_mod)
resid_auxpanel(fusarium_colht_res, fusarium_colht_fit)

# estimated slopes (on log scale because of the link function used)
(fusarium_colht_emt = emtrends(fusarium_colht_mod, ~ Variety|Plant, 
                               var = 'Perc_discoloration'))
# the difference is non-significant
contrast(fusarium_colht_emt, 'pairwise') 

# Estimates and plots at the variety level (note lots of uncertainty regarding
# italian parsley) 
# generate predictions for each species over the 0-100% range of discoloration
fusarium_colht_preds = emmeans(fusarium_colht_mod, ~ Perc_discoloration|Variety,
                               type = 'response', 
                               at = list(Perc_discoloration = 0:100))
# plot the results
fusarium_colht_preds |> 
    as.data.frame() |> 
    ggplot(aes(x = Perc_discoloration, fill = Variety, color = Variety))+
    facet_grid(cols = vars(Plant))+
    geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), color = NA, alpha = .2)+
    geom_line(aes(y = response))+
    geom_point(aes(y = height_final), data = fusarium, alpha = .5) +
    theme_Publication() +
    labs(x = "Vascular discoloration (%)",y = "Plant height (cm)") 

# Alternatively, average over varieties and report species level trends
  fusarium_colht_preds = emmeans(fusarium_colht_mod, ~ Perc_discoloration|Plant,
                                 type = 'response', 
                                 at = list(Perc_discoloration = 0:100))
  # plot the results 
 (colht_plot <- fusarium_colht_preds |> 
    as.data.frame() |> 
    ggplot(aes(x = Perc_discoloration))+
    facet_grid(cols = vars(Plant))+
    geom_ribbon(aes(ymin = asymp.LCL, ymax = asymp.UCL), color = NA, alpha = .2)+
    geom_line(aes(y = response))+
    geom_point(aes(y = height_final), data = fusarium, alpha = .5) +
    theme_Publication() +
    labs(x = "",y = "Plant height (cm)") +
   theme(axis.ticks.x = element_blank(),axis.text.x=element_blank(),
         strip.background = element_blank(),strip.text.x = element_blank(),
         plot.margin = unit(c(0,0.2,0.8,1), 'lines'))+
   xlim(0, 100))
 
 colht_plot
  
 # We can use the delta method to convert it to a more interpretable measure
 # on the original scale: since we're using the log link, effects are
 # multiplicative and we can express the "slope" as percent decline in total
 # weight per 1% increase in discoloration. First we extract the model
 # estimates and and variance-covariance matrix from our fitted model, then
 # feed those and the appropriate formula into the detaMethod function:
 # 
 fe = fixef(fusarium_colht_mod)$cond
 keep = !is.na(fe) # some missing coefficients have NA values, we want to remove those
 fe = fe[keep]
 vc = vcov(fusarium_colht_mod)$cond
 vc = vc[keep, keep]
 
 # Here are the estimates (along with SE and confidence limits) for percent
 # decline (because negative) associated with 1% increase in discoloration.
 # Note that you have to be a bit careful with the specification of the
 # formula, because it depends on what contrast specifications R is using
 # (this works for the default)
 car::deltaMethod(object = fe,
                  g. = '100*(exp(Perc_discoloration)-1)',  # for celery, because that is the "reference level" for the plant factor
                  vcov = vc)
 car::deltaMethod(object = fe,
                  g. = '100*(exp(Perc_discoloration + `PlantParsley:Perc_discoloration`)-1)', # for parsley, because we've added to the "reference level" the parley-effect
                  vcov = vc)

####Figure####
library(ggpubr)
legend1 <- get_legend(height_plot)
legend1 <- as_ggplot(legend1)
legend1<- legend1 + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))

height_plot2 <- height_plot + theme(legend.position = "none")
total_plot2 <- total_plot + theme(legend.position = "none")
root_plot2 <- root_plot + theme(legend.position = "none")

p1<- ggarrange(height_plot2, total_plot2, root_plot2, labels = c("A", "B", "C"),
              ncol = 1, nrow = 3, common.legend = F)

legend2 <- get_legend(colht_plot)
legend2 <- as_ggplot(legend2)
legend2<- legend2 + theme(plot.margin = unit(c(0,0.2,0,1), 'lines'))


colht_plot2 <- colht_plot + theme(legend.position = "none")
colwt_plot2 <- colwt_plot + theme(legend.position = "none")
colrt_plot2 <- colrt_plot + theme(legend.position = "none")
p2<-ggarrange(colht_plot2,colwt_plot2, colrt_plot2, labels = c("D", "E", "F"),
          ncol = 1, nrow = 3, common.legend = F)

ggarrange(p1, p2, align = "h") +
  theme(plot.margin = margin(0.1,0.1,0.5,0.1, "cm")) 
ggarrange(legend1, legend2, ncol=1, nrow=2)
ggsave("fig6.svg", units = "px", width = 2800, height = 2100,
       device = "svg", dpi = 300)

ggsave("fig6_legend.svg", units = "px", width = 2800, height = 2100,
       device = "svg", dpi = 300, bg="white")


