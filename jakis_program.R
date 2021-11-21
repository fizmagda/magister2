app get ubuntu instalowanie pakietów
renv restore
library(powerHDX)
library(dplyr)
library(mgcv)
library(stringr)

renv::restore()
?readRDS
prepare_data = function(data){
  data %>% 
    mutate(id = as.numeric(stringr::str_remove(id, " ")))
}

deuteros = function(data, significance_level = 0.05) {
  States = unique(data$State)
  
  aic = loglik = Test_statistic = p_value = NA
  model = lm(Mass ~ Exposure*State, data = data)
  model_reduced = lm(Mass ~ Exposure, data = data)
  result = anova(model, model_reduced)
  aic= AIC(model)
  loglik = as.numeric(logLik(model))
  F_statistic = result$`F`[2]
  p_value = result$`Pr(>F)`[2]
  
  data.table::data.table(Test = "Deuteros lm",
                         State_1 = States[1],
                         State_2 = States[2],
                         Test_statistic = F_statistic,
                         P_value = p_value,
                         Significant_difference = (p_value <= significance_level),
                         AIC = aic,
                         logLik = loglik)
}

semiparametric5 <- function(data, significance_level = 0.05) {
  States = unique(data$State)
  
  aic = loglik = Test_statistic = p_value = NA
  model_reduced <- gam(Mass ~ s(Exposure) + Exposure,
                       data=data,
                       method="REML")
  model <- gam(Mass~s(Exposure) + Exposure*State,
               data = data,
               method="REML")
  
  result = anova.gam(model, model_reduced, test = "Chisq")
  aic = AIC(model)
  Test_statistic = result$Deviance[2]
  p_value = result$`Pr(>Chi)`[2]
  loglik = as.numeric(logLik(model))
  
  data.frame(Test = "semiparametric interaction",
             State_1 = States[1],
             State_2 = States[2],
             Test_statistic = Test_statistic,
             P_value = p_value,
             Significant_difference = (p_value <= significance_level),
             AIC = aic,
             logLik = loglik)
}

calculate_power_error <- function(spectra, n_runs) {
  lapply (1:n_runs, function(i) {
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
    
    n_pairs <- ncol(pairs)
    
    do.call(rbind, lapply(1:n_pairs, function(pair) {
      paired_data <- data %>% 
        filter(State %in% pairs[, pair])
      tryCatch({
        # deut <- deuteros(paired_data)
        # semi1 <- semiparametric1(paired_data)
        # semi2 <- semiparametric2(paired_data)
        # semi3 <- semiparametric3(paired_data)
        # semi4 <- semiparametric4(paired_data)
        # rbind(deut, semi1, semi2,
        #       semi3, semi4)
        deut <- deuteros(paired_data)
        semi5 <- semiparametric5(paired_data)
        rbind(deut, semi5)
      },
      error = function(e) {
        print(e)
        data.frame()
      })
    }))
  })
}

new_pf <- c(10, 15, 20, 30, 40, 45, 50, seq(100, 2000, 100))
times = c(5, 10, 20, 30, 40, 50, 60, 100, 300, 500, 900, 1200, 1500, 1800,
          2100, 2400, 3600, 7200, 21600, 43200)

set.seed(17)
params <- readRDS("./all_params.RDS")
sequences <- sample(unique(params$sequence), 50, replace = FALSE)


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

spectra_by_seq <- split(theoretical_spectra, 
                        f = theoretical_spectra[, c("Sequence")])

tests_results1 <- lapply(1:length(spectra_by_seq), function(ith_spectrum) {
  res_spectrum <- calculate_power_error(spectra_by_seq[[ith_spectrum]], n_runs = 100)
  saveRDS(res_spectrum, paste0("sequence_", ith_spectrum, "_rejection_rate.RDS"))
})
