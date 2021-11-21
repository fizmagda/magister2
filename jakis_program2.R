library(tidybayes)
library(rstan)
library(forcats)
library(brms)
library(rstanarm)
library(magrittr)
library(dplyr)
library(modelr)
library(ggdist)
library(ggplot2)
library(cowplot)
library(emmeans)
library(broom)
library(bayesplot)
library(MCMCglmm)
library(RColorBrewer)
library(withr)
library(rstantools)
library(HaDeX)

dat <- read_hdx(system.file(package = "HaDeX", 
                            "HaDeX/data/KD_190304_Nucb2_EDTA_CaCl2_test02_clusterdata.csv"))
View(head(dat,5))

#----------wyznaczenie masy peptydu: pepMass=z×(Center−protonMass)
?prepare_dataset
d<-prepare_dataset(dat)
dat['pepMass']

datprz <- read_hdx(system.file(package = "HaDeX", "HaDeX/data/KD_180110_CD160_HVEM.csv"))

# prepare dataset for states `CD160` and `CD160_HVEM` in given time parameters
datpr<-prepare_dataset(datprz,
                in_state_first = "CD160_0.001",
                chosen_state_first = "CD160_1",
                out_state_first = "CD160_1440",
                in_state_second = "CD160_HVEM_0.001",
                chosen_state_second = "CD160_HVEM_1",
                out_state_second = "CD160_HVEM_1440")
View(datpr)
comparison_plot(calc_dat = datpr,
                theoretical = TRUE,
                relative = TRUE,
                state_first = "Nucb2 Factor 1",
                state_second = "Nucb2 Factor 2") +
  labs(title = "Theoretical fraction exchanged in state comparison in 25 min time")
comparison_plot(calc_dat = datpr,
                theoretical = FALSE,
                relative = TRUE, 
                state_first = "Nucb2 Factor 1",
                state_second = "Nucb2 Factor 2") +
  labs(title = "Fraction exchanged in state comparison in 25 min time")
