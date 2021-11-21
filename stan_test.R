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


strsplit(Sys.getenv("PATH"), ";")


set.seed(5)
n = 10
n_condition = 5
ABC =
  data.frame(
    condition = rep(c("A","B","C","D","E"), n),
    response = rnorm(n * 5, c(0,1,2,1,-1), 0.5)
  )
head(ABC, 10)


ggplot(ABC,aes(x = response, y = fct_rev(condition))) +
  geom_point(alpha = 0.5) +
  ylab("condition")

ABC_stan<-stan_model(file="C:/Users/Magdalena/Documents/hdx_fixEf.stan")

print(compose_data(ABC))

m <- sampling(ABC_stan, data = compose_data(ABC), control = list(adapt_delta=0.99))

print(m, pars = c("overall_mean", "condition_mean_sd", "condition_mean", "response_sd"))
