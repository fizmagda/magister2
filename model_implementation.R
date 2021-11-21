
library(powerHDX)
library(dplyr)
library(mgcv)
library(stringr)
library(tidybayes)
library(rstan)
library(forcats)
library(brms)
library(rstanarm)
library(magrittr)
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

renv::restore()

prepare_data = function(data){
  data %>% 
    mutate(id = as.numeric(stringr::str_remove(id, " ")))
}

calculate_smth<- function(spectra)
  {
    #simulate one state
    noisy_curves = get_noisy_deuteration_curves(spectra, reference = 1,
                                                n_runs = 1,
                                                mass_deviations = 5,
                                                per_run_deviations = 0.1,
                                                compare_pairs = FALSE)
    res <- cbind(do.call(rbind, noisy_curves[[1]]), rep(1:100, each = 21*length(unique(spectra$PF))))
    data_1 <- res %>% 
      mutate(id = paste(V2, State)) %>% 
      dplyr::select(-V2)
    #simulate the same state
    noisy_curves_2 = get_noisy_deuteration_curves(spectra, reference = 1,
                                                  n_runs = 1,
                                                  mass_deviations = 5,
                                                  per_run_deviations = 0.1,
                                                  compare_pairs = FALSE)
    res_2 <- cbind(do.call(rbind, noisy_curves_2[[1]]), rep(1:(100*length(unique(spectra$PF))), each = 21))
    data_2 <- res_2 %>% 
      mutate(id = paste(V2, State)) %>% 
      dplyr::select(-V2) %>% 
      mutate(State = as.factor(paste0(State, "_2")))
    
    data <- prepare_data(rbind(data_1, data_2))
    pairs <- cbind(rbind(levels(data_1$State), 
                         levels(data_2$State)), 
                   combn(unique(levels(data_1$State)), 2))
    
data
    
  }  


new_pf <- c(10, 15, 20, 30, 40, 45, 50, seq(100, 2000, 100))
times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)

set.seed(17)
params <- readRDS("./all_params.RDS")
sequences <- sample(unique(params$sequence), 5, replace = FALSE)

library(dplyr)
all_params <- params %>% 
  dplyr::select(sequence, pH, step) %>% 
  unique() %>% 
  filter(sequence %in% sequences, pH == "6") %>% 
  group_by(sequence) %>% 
  slice(rep(1:n(), each = length(new_pf))) %>% 
  mutate(protection_factor = new_pf) %>% 
  mutate(charge = sample(1:5,length(new_pf),  replace = TRUE)) %>% 
  ungroup() %>% 
  data.frame()

theoretical_spectra <- do.call(rbind, lapply(1:nrow(all_params), function(ith_row) {
  simulate_theoretical_spectra(sequence = all_params[ith_row, "sequence"],
                               charge = all_params[ith_row, "charge"],
                               protection_factor = all_params[ith_row, "protection_factor"],
                               times = times,
                               pH = all_params[ith_row, "pH"],
                               temperature = 15,
                               n_molecules = 500,
                               time_step_const = all_params[ith_row, "step"],
                               use_markov = TRUE)
}))


spectra_by_seq <- split(theoretical_spectra, f = theoretical_spectra[, c("Sequence")])
spectra_by_seq$VGGEALGRL
d<-calculate_smth(spectra_by_seq[[1]])
ps<-as.data.frame(d)
#tests_results1 <- lapply(1:length(spectra_by_seq), 
 #                        function(ith_spectrum) {res_spectrum <- calculate_smth(spectra_by_seq[[ith_spectrum]])})

#------MODEL-------
#devtools::install_github("hadexversum/powerHDX")
library(powerHDX)
#library("tidybayes", lib.loc="C:/Users/Magdalena/Documents/mgr/magister/renv/library/R-4.0/x86_64-w64-mingw32")
#remotes::install_version("rstan",lib="C:/Users/Magdalena/Documents/mgr/magister/renv/library/R-4.0/x86_64-w64-mingw32")
#options(mc.cores = parallel::detectCores())
#ABC_stan<-stan_model(file="C:/Users/Magdalena/Documents/mgr/magister/fixEf.stan")
A_stan<-rstan::stan_model(file="C:/Users/Magdalena/Documents/mgr/magister/hdx_fixEf.stan")


print(tidybayes::compose_data(ps))
typeof(ps)

m <- rstan::sampling(A_stan, data = tidybayes::compose_data(ps), control = list(adapt_delta=0.99))





