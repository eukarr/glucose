library(tidyverse)
library(broom)
library(ggtext)
library(car)

#========this part of the script requires local files====================================
#========and should not be run if 'glucose_for_analysis.csv' is available================

setwd("C:/Users/Evgeny/Dropbox/_projects/CQD/glucose")

# this is a dataset of input and response variables collected manually
# id - sample label
# volume - reactor volume in mL
# initial mass - mass of the reaction mixture, g
# concentration - concentration of glucose in the reaction mixture, wt%
# duration - synthesis duration, h
# pH - initial pH of the reaction mixture
# temperature - syntesis temperature
# pH_final - pH of the filtrate
# yield_total - total yield of the product, %, including sediment
# conc_filtrate - concentration of the filtrate, mg/g
# yield_filtrate - yield of the soluble product, %
glucose <- read.csv("glucose_manual_data.csv", fileEncoding="UTF-8-BOM", na.strings = "n/a")

# reading absorbance spectra
samples_list_filepath <- paste0(getwd(), "/spectra/absorbance_citrate/samples_list.txt")
samples <- read.delim(file = samples_list_filepath)
colnames(samples) <- c('sample', 'files', 'comment')
samples$files <- strsplit(samples$files,", ")

spectrum_path <- paste0(getwd(), "/spectra/absorbance_citrate/")
file_list <- list.files(path=spectrum_path, pattern = "Absorbance_\\d+\\.txt")

citrate_spectra <- data.frame()

# each spectrum is corrected for background and cut
for (i in 1:length(file_list)){
  temp_data <- read.delim(file = paste0(spectrum_path, '\\', file_list[i]), skip = 14, dec = ',', header=FALSE)
  temp_data$spectrum <- sapply(strsplit(gsub(".txt", "", file_list[i]), "_"), function(x){x[2]})
  # baseline correction; lines 1336-1569 correspond to 800-900 nm range
  background <- mean(temp_data$V2[1336:1569])
  temp_data$V2 <- temp_data$V2 - background
  # lines 121-879 correspond to 255.2-599.7 nm range, 758 points per spectrum, at 0.454 nm spacing
  temp_data <- temp_data[121:879,]
  for (j in 1:length(samples$sample)){
    if (temp_data$spectrum[1] %in% samples[j, 2][[1]]) sample <- samples[j, 1]
  }
  temp_data$sample <- str_trim(sample)
  citrate_spectra <- rbind(citrate_spectra, temp_data)
}

rm(temp_data)
rm(samples)
rm(list = c("background", "file_list", "i", "j", "sample", "samples_list_filepath", "spectrum_path"))

colnames(citrate_spectra) = c("wavelength", "Abs", "spectrum", "sample_id")

# for each sample, a median spectrum of 2-3 is retained
citrate_median_spectra_by_sample <- citrate_spectra %>%
  group_by(sample_id, wavelength) %>%
  summarize(Abs = median(Abs))

# for each median spectrum, maximum Abs is saved, and the midpoint of the range corresponding to this Abs value
citrate_by_sample_abs <- citrate_median_spectra_by_sample %>%
  summarize(Abs_max = max(Abs),
            Abs_365 = Abs[round(wavelength, 1) == 364.9],
            Abs_365_norm = Abs_365 / Abs_max,
            Abs_315 = Abs[round(wavelength, 1) == 314.8],
            Abs_315_norm = Abs_315 / Abs_max,
            Abs_302 = Abs[round(wavelength, 1) == 302.4],
            Abs_302_norm = Abs_315 / Abs_max,
            wavelength_max = mean(wavelength[Abs == Abs_max]))

# adding spectral data to the main dataset
# Abs_max - the maximum Abs in baseline-corrected spectra
# wavelength_max - wavelength corresponding to the maximum Abs
glucose_full <- merge(glucose, citrate_by_sample_abs, by.x="id", by.y="sample_id", all.x = TRUE)      
write.csv(x = glucose_full, file = "glucose_full.csv")

glucose <- glucose_full[complete.cases(glucose_full), ]
write.csv(x = glucose, file = "glucose_for_analysis.csv")
#=======end of the script demanding local files=======================================



#=======START OF MAIN ANALYSIS========================================================

#loading data 
glucose <- read.csv("glucose_for_analysis.csv", fileEncoding="UTF-8-BOM", na.strings = "n/a")

# a function to encode factorial variables into the [-1, 1] range
encode_var = function(x) {
  coded <- 2 * (x - mean(x)) / (max(x) - min(x))
}

# preparing dataset with coded variables for ANOVA
glucose_coded <- glucose %>%
  mutate(concentration = encode_var(concentration),
         duration = encode_var(duration),
         pH = encode_var(trunc(pH)),
         temperature = encode_var(temperature)) 
write.csv(x = glucose_coded, file = "glucose_for_analysis_coded.csv")



#==============EFFECT OF VOLUME========================================
# v_eff is a dataset for testing the effect of volume
v_eff <- glucose %>%
  mutate(volume_test = NA)

# if for a given sample there is a 'paired' with a different volume
# and other conditions equal, mark both with TRUE
for (i in 1:nrow(v_eff)) {
  if (isTRUE(v_eff$volume_test[i])) next
  for (j in (i + 1):nrow(v_eff)) {
    if ((v_eff$concentration[i] == v_eff$concentration[j]) &
        ((v_eff$pH[i] == v_eff$pH[j])) &
        (v_eff$duration[i] == v_eff$duration[j]) &
        (v_eff$temperature[i] == v_eff$temperature[j]) &
        (v_eff$volume[i] != v_eff$volume[j])) {
      v_eff$volume_test[i] <- TRUE
      v_eff$volume_test[j] <- TRUE
    }
  }
}

# select the samples with at least one 'paired' counterpart
# and subset only necessary columns
v_eff <- v_eff[!(is.na(v_eff$volume_test)), ] %>%
  subset(select = c("id", "initial_mass", "temperature", "duration", "pH", "concentration", "volume", "yield_total", "yield_filtrate", "pH_final", "Abs_max")) %>%
  mutate(Abs_max_rel = Abs_max / (initial_mass * concentration),
         pH_final_rel = pH_final / (initial_mass * concentration))

# cross-join of the 'paired' samples differing only in volume
# '.x' corresponds to 5-mL reactors
# '.y' corresponds to 10-mL reactors
# and calculate the difference in the response variables in each pair
v_eff <- merge(v_eff[v_eff$volume == 5, ],
               v_eff[v_eff$volume == 10, ], 
                   by=c("temperature", "duration", "pH", "concentration")) %>%
  mutate(yield_total_diff = yield_total.x - yield_total.y,
         yield_filtrate_diff = yield_filtrate.x - yield_filtrate.y,
         pH_final_diff = pH_final.x - pH_final.y,
         pH_final_rel_diff = pH_final_rel.x - pH_final_rel.y,
         Abs_max_diff = Abs_max.x - Abs_max.y,
         Abs_max_rel_diff = Abs_max_rel.x - Abs_max_rel.y)

# long-format dataset of the variables scaled by SD for creating boxplot
v_eff_plot <- v_eff[, 23:28] %>%
  apply(MARGIN = 2, function(x) {x / sd(x)}) %>%
  as.data.frame() %>%
  gather(key = variable, value = value, factor_key=TRUE)

# boxplot for the distribution of the variables to assess the effect of reactor volume
# Fig. 2 in the paper
ggplot(data = v_eff_plot) +
  geom_hline(yintercept = 0, colour = "red") +
  geom_boxplot(aes(x=variable, y=value)) +
  ylab(label = "Z-score") +
  xlab(label = NULL) +
  scale_x_discrete(labels = c('Total yield',
                              'Filtrate yield',
                              'Final pH',
                              'Relative\nfinal pH',
                              'Maximum\nabsorbance',
                              'Relative\nmaximum\nabsorbance'))

v_eff_t <- glucose %>%
  subset(select = c("id", "concentration", "initial_mass", "volume", "yield_total", "yield_filtrate", "pH_final", "Abs_max")) %>%
  mutate(Abs_max_rel = Abs_max / (initial_mass * concentration),
         pH_final_rel = pH_final / (initial_mass * concentration)) %>%
  subset(select = c("id", "volume", "yield_total", "yield_filtrate", "pH_final", "pH_final_rel", "Abs_max", "Abs_max_rel")) %>%
  mutate(yield_total = (yield_total - mean(yield_total)) / sd(yield_total),
         yield_filtrate = (yield_filtrate - mean(yield_filtrate)) / sd(yield_filtrate),
         pH_final = (pH_final - mean(pH_final)) / sd(pH_final),
         pH_final_rel = (pH_final_rel - mean(pH_final_rel)) / sd(pH_final_rel),
         Abs_max = (Abs_max - mean(Abs_max)) / sd(Abs_max),
         Abs_max_rel = (Abs_max_rel - mean(Abs_max_rel)) / sd(Abs_max_rel))



# t-tests for significance of the effect of the reactor volume
t1 <- t.test(v_eff_t$yield_total[v_eff_t$volume == 5], v_eff_t$yield_total[v_eff_t$volume == 10])
t2 <- t.test(v_eff_t$yield_filtrate[v_eff_t$volume == 5], v_eff_t$yield_filtrate[v_eff_t$volume == 10])
t3 <- t.test(v_eff_t$pH_final[v_eff_t$volume == 5], v_eff_t$pH_final[v_eff_t$volume == 10])
t4 <- t.test(v_eff_t$pH_final_rel[v_eff_t$volume == 5], v_eff_t$pH_final_rel[v_eff_t$volume == 10])
t5 <- t.test(v_eff_t$Abs_max[v_eff_t$volume == 5], v_eff_t$Abs_max[v_eff_t$volume == 10])
t6 <- t.test(v_eff_t$Abs_max_rel[v_eff_t$volume == 5], v_eff_t$Abs_max_rel[v_eff_t$volume == 10])

# table with the results of the t-test
v_eff_result <- map_df(list(t1, t2, t3, t4, t5, t6), tidy) %>%
  subset(select = c("estimate", "statistic", "p.value", "conf.low", "conf.high")) %>%
  as.data.frame() %>%
  signif(3)

rownames(v_eff_result) <- c("Total yield",
                            "Filtrate yield",
                            "Final pH",
                            "Relative final pH",
                            "Absorbance in maximum",
                            "Relative absorbance in maximum")

v_eff_result
#==============END OF EFFECT OF VOLUME=================================



#==============ANALYSIS OF SOLIDS YIELD================================
# mean yields at 2 h treatment and 160C
glucose_coded %>%
  filter(temperature == -1 & duration == -1) %>%
  subset(select = c(yield_total, yield_filtrate)) %>%
  apply(MARGIN = 2, mean)

# mean yields at 8 h treatment and 180C
glucose_coded %>%
  filter(temperature == 1 & duration == 1) %>%
  subset(select = c(yield_total, yield_filtrate)) %>%
  apply(MARGIN = 2, mean)

# dataset with difference in the total yield and filtrate yield
yield_diff <- glucose_coded %>%
  subset(select = -c(pH_final, conc_filtrate, Abs_max, wavelength_max, initial_mass)) %>%
  mutate(yield_diff = (yield_total - yield_filtrate) / yield_filtrate)

ggplot(data = yield_diff) +
  geom_histogram(aes(x = yield_diff), bins = 20)

# selection of the outlier samples
yield_diff_limits <- quantile(yield_diff$yield_diff, c(0.15, 0.85))

yield_outliers <- yield_diff %>%
  filter(yield_diff < yield_diff_limits[1] | yield_diff > yield_diff_limits[2]) %>%
  arrange(abs(yield_diff))

# plot for the total and yields marking the outliers
# Fig. 3 in the paper
ggplot(data = glucose_coded, aes(x = yield_total, y = yield_filtrate)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1) +
  geom_hline(yintercept = 100, linetype = 'dashed') + 
  geom_vline(xintercept = 100, linetype = 'dashed') +
  geom_point(data = yield_outliers, aes(x = yield_total, y = yield_filtrate), color = "red") +
  xlab(label = "Total yield, %") +
  ylab(label = "Filtrate yield, %") +
  ylim(c(0, 100)) + xlim(c(0, 130))

# linear full model and its summary
model_yield_filtrate <- lm(data = glucose_coded, 
                           formula = yield_filtrate ~ concentration * duration * pH * temperature)
Anova(model_yield_filtrate, type = "II")
summary(model_yield_filtrate)


# model with only significant factors and its analysis         
model_yield_filtrate_reduced <- lm(data = glucose_coded, 
                                   formula = yield_filtrate ~ concentration + duration + pH + temperature +
                                     duration:temperature)
Anova(model_yield_filtrate_reduced, type = "II")
summary(model_yield_filtrate_reduced)



# checking the residuals of the reduced model
model_yield_filtrate_reduced_res <- data.frame(residual = model_yield_filtrate_reduced$residuals,
                                       measured = model_yield_filtrate_reduced[["model"]][["yield_filtrate"]],
                                       fitted = model_yield_filtrate_reduced$fitted.values)

ggplot(data = model_yield_filtrate_reduced_res) +
  geom_point(aes(x = measured, y = fitted)) +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1)

ggplot(data = model_yield_filtrate_reduced_res, aes(x = residual)) +
  geom_histogram(aes(y = stat(density)), bins = 15) +
  geom_density(adjust = 2, color = "red", size = 1)

ggplot(data = model_yield_filtrate_reduced_res) +
  geom_point(aes(x = measured, y = residual))


# reduced model excluding central point
model_yield_filtrate_reduced_no_central <- lm(data = glucose_coded %>%
                                             filter(concentration != 0), 
                                   formula = yield_filtrate ~ concentration + duration + pH + temperature +
                                     duration:temperature)

# mean filtrate yield in the central point
experiment_central_yield_filtrate <- glucose_coded %>%
  filter(concentration == 0) %>%
  summarize(mean(yield_filtrate, na.rm = TRUE)) %>%
  as.numeric(.)

# prediction of the filtrate yield in the central point using reduced linear model
predicted_central_yield_filtrate <- predict(model_yield_filtrate_reduced_no_central, 
                                            new = data.frame(concentration = 0,
                                                           duration = 0,
                                                           pH = 0,
                                                           temperature = 0), interval = "prediction")


# plot of the interaction of 'temperature' and 'duration' factors
ggplot(data = glucose %>%
         filter(concentration != 7.5) %>%
         group_by(temperature, duration) %>%
         summarize(yf_mean = mean(yield_filtrate),
                   upper = yf_mean + qt(0.975, df = length(yield_filtrate) - 1) * sd(yield_filtrate)/sqrt(length(yield_filtrate)),
                   lower = yf_mean - qt(0.975, df = length(yield_filtrate) - 1) * sd(yield_filtrate)/sqrt(length(yield_filtrate))),
       aes(x = temperature, y = yf_mean)) +
  geom_line(aes(color = factor(duration)), size = 1) +
  geom_point(aes(color = factor(duration)), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(duration)), width = 0.2) +
  xlab(label = "Temperature (°C)") +
  ylab(label = "Filtrate yield, %") +
  labs(color = "Duration (h)") +
  theme(legend.position = c(0.85, 0.85))
#==============END OF ANALYSIS OF YIELD================================



#==============ANALYSIS OF pH==========================================
glucose_10mL <- glucose_coded %>%
  filter(volume == 10)

glucose_5mL <- glucose_coded %>%
  filter(volume == 5)

cor.test(glucose_10mL$pH_final, glucose_10mL$yield_filtrate, method = 'spearman')
cor.test(glucose_5mL$pH_final, glucose_5mL$yield_filtrate, method = 'spearman')

# plot of general relationship between the filtrate pH and the gravimetric yield
# Fig. S2 in the paper
ggplot(data = glucose_coded, aes(x = yield_filtrate, y = pH_final, color = factor(volume))) +
  geom_point(size = 3) +
  geom_smooth(method = "loess", se = FALSE) +
  xlab(label = "Yield in filtrate (%)") +
  ylab(label = "Filtrate pH") +
  labs(color = "Reactor volume (mL)") +
  theme(legend.position = c(0.2, 0.85))
  

# linear full model and its summary
model_pH_final <- lm(data = glucose_10mL, 
                     formula = pH_final ~ concentration * duration * pH * temperature)
Anova(model_pH_final, type = "II")
summary(model_pH_final)

#linear model including only significant terms and its summary
model_pH_final_reduced <- lm(data = glucose_10mL, 
                     formula = pH_final ~ concentration + duration + pH + temperature +
                       concentration:pH + duration:temperature + concentration:pH:temperature)
Anova(model_pH_final_reduced, type = "II")
summary(model_pH_final_reduced)


# checking the residuals of the reduced model
model_pH_final_reduced_res <- data.frame(residual = model_pH_final_reduced$residuals,
                                               measured = model_pH_final_reduced[["model"]][["pH_final"]],
                                               fitted = model_pH_final_reduced$fitted.values)

ggplot(data = model_pH_final_reduced_res) +
  geom_point(aes(x = measured, y = fitted)) +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1)

ggplot(data = model_pH_final_reduced_res, aes(x = residual)) +
  geom_histogram(aes(y = stat(density)), bins = 15) +
  geom_density(adjust = 2, color = "red", size = 1)

ggplot(data = model_pH_final_reduced_res) +
  geom_point(aes(x = measured, y = residual))


# reduced model excluding central point
model_pH_final_reduced_no_central <- lm(data = glucose_10mL %>%
                                                filter(concentration != 0), 
                                        formula = pH_final ~ concentration + duration + pH + temperature +
                                          concentration:pH + duration:temperature + concentration:pH:temperature)

# mean filtrate pH in the central point
experiment_central_pH_final <- glucose_10mL %>%
  filter(concentration == 0) %>%
  summarize(mean(pH_final, na.rm = TRUE)) %>%
  as.numeric(.)

# prediction of the filtrate pH in the central point using reduced linear model
predicted_central_pH_final <- predict(model_pH_final_reduced_no_central, 
                                            new = data.frame(concentration = 0,
                                                             duration = 0,
                                                             pH = 0,
                                                             temperature = 0), interval = "prediction")

# plot of the interaction of 'temperature' and 'duration' factors
ggplot(data = glucose %>%
         filter(concentration != 7.5 & volume ==10) %>%
         group_by(temperature, duration) %>%
         summarize(pH_mean = mean(pH_final),
                   upper = pH_mean + qt(0.975, df = length(pH_final) - 1) * sd(pH_final)/sqrt(length(pH_final)),
                   lower = pH_mean - qt(0.975, df = length(pH_final) - 1) * sd(pH_final)/sqrt(length(pH_final))),
       aes(x = temperature, y = pH_mean)) +
  geom_line(aes(color = factor(duration)), size = 1) +
  geom_point(aes(color = factor(duration)), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(duration)), width = 0.2) +
  xlab(label = "Temperature (°C)") +
  ylab(label = "Filtrate pH") +
  labs(color = "Duration (h)") +
  theme(legend.position = c(0.85, 0.85))
#==============END OF ANALYSIS OF pH===================================



#==============ANALYSIS OF ABSORBANCE==================================
# choosing a sample spectrum at 2 h and 160°C 
glucose_coded %>%
  filter(duration == -1 & temperature == -1 & pH == 0 & volume == 10 & concentration == 1)
# sample 1410-3 chosen

# choosing a sample spectrum at 8 h and 180°C
glucose_coded %>%
  filter(duration == 1 & temperature == 1 & pH == 0 & volume == 10 & concentration == 1)
# sample 1610-8 chosen

# comparison of spectra of two products differing in conversion
ggplot(data = citrate_median_spectra_by_sample %>%
         filter(sample_id == "1410-3" | sample_id == "1610-8"), 
       aes(x = wavelength, y = Abs)) +
  geom_line(aes(color = sample_id)) +
  xlab(label = "Wavelength (nm)") +
  ylab(label = "Absorbance") +
  labs(color = "Treatment conditions") +
  scale_color_discrete(labels = c('2 h at 160 °C',
                              '8 h at 180 °C')) +
  theme(legend.position = c(0.85, 0.85))
  

# linear full model for absolute Abs_max in 10-mL reactors
#model_Abs_max <- lm(data = glucose_coded, 
#                    formula = Abs_max ~ concentration * duration * pH * temperature)
#Anova(model_Abs_max, type = "II")
#summary(model_Abs_max)


# linear full model for Abs_max normalized to glucose loading in all reactors
model_Abs_max <- lm(data = glucose_coded %>%
                      mutate(c_pct = ifelse(concentration == 1, 0.10, ifelse(concentration == -1, 0.05, 0.075))), 
                    formula = Abs_max / (initial_mass * c_pct) ~ concentration * duration * pH * temperature)
Anova(model_Abs_max, type = "II")
summary(model_Abs_max)


summary(glucose_coded$Abs_max)

# plot of the interaction of 'temperature' and 'duration' factors
ggplot(data = glucose %>%
         filter(concentration != 7.5) %>%
         mutate(Abs_max = Abs_max / (initial_mass * concentration)) %>%
         group_by(temperature, duration) %>%
         summarize(Abs_max_mean = mean(Abs_max),
                   upper = Abs_max_mean + qt(0.975, df = length(Abs_max) - 1) * sd(Abs_max)/sqrt(length(Abs_max)),
                   lower = Abs_max_mean - qt(0.975, df = length(Abs_max) - 1) * sd(Abs_max)/sqrt(length(Abs_max))),
       aes(x = temperature, y = Abs_max_mean)) +
  geom_line(aes(color = factor(duration)), size = 1) +
  geom_point(aes(color = factor(duration)), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(duration)), width = 0.2) +
  xlab(label = "Temperature (°C)") +
  ylab(label = "Absorbance at maximum, normalized") +
  labs(color = "Duration (h)") +
  theme(legend.position = c(0.45, 0.9))


# plot of the interaction of 'concentration' and 'duration' factors
ggplot(data = glucose %>%
         filter(concentration != 7.5) %>%
         mutate(Abs_max = Abs_max / (initial_mass * concentration)) %>%
         group_by(concentration, duration) %>%
         summarize(Abs_max_mean = mean(Abs_max),
                   upper = Abs_max_mean + qt(0.975, df = length(Abs_max) - 1) * sd(Abs_max)/sqrt(length(Abs_max)),
                   lower = Abs_max_mean - qt(0.975, df = length(Abs_max) - 1) * sd(Abs_max)/sqrt(length(Abs_max))),
       aes(x = concentration, y = Abs_max_mean)) +
  geom_line(aes(color = factor(duration)), size = 1) +
  geom_point(aes(color = factor(duration)), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(duration)), width = 0.2) +
  xlab(label = "Concentration (wt%)") +
  ylab(label = "Absorbance at maximum, normalized") +
  labs(color = "Duration (h)") +
  theme(legend.position = c(0.45, 0.9))


# a plot to reveal the difference in absorption between samples 1410-3 and 1610-8
ggplot(data = citrate_median_spectra_by_sample %>%
         filter(sample_id == "1410-3" | sample_id == "1610-8") %>%
         group_by(wavelength) %>%
         mutate(Abs_diff = (max(Abs) - min(Abs)) / mean(Abs)), 
       aes(x = wavelength, y = Abs_diff)) +
  geom_point() +
  xlab(label = "Wavelength (nm)") +
  ylab(label = "Absorbance difference")


# revealing the wavelength corresponding to highest between-groups difference in absorption
Abs_sd_within <- merge(citrate_median_spectra_by_sample, 
                     citrate_median_spectra_by_sample %>%
                       group_by(sample_id) %>%
                       summarize(peak = max(Abs)), 
                     by = 'sample_id') %>%
  mutate(Abs_norm = Abs / peak) %>%
  filter(grepl("0311", sample_id, fixed = TRUE)) %>%
  group_by(wavelength) %>%
  summarise(Abs_sd = sd(Abs_norm))

Abs_sd <- merge(citrate_median_spectra_by_sample, 
              citrate_median_spectra_by_sample %>%
                group_by(sample_id) %>%
                summarize(peak = max(Abs)), 
              by = 'sample_id') %>%
  mutate(Abs_norm = Abs / peak) %>%
  group_by(wavelength) %>%
  summarise(Abs_sd = sd(Abs_norm)) %>%
  mutate(sd_ratio = Abs_sd / Abs_sd_within$Abs_sd) %>%
  filter(wavelength > 295)

ggplot(data = Abs_sd, aes(x = wavelength, y = sd_ratio)) +
  geom_line()

print(c(Abs_sd$wavelength[Abs_sd$sd_ratio == max(Abs_sd$sd_ratio)], max(Abs_sd$sd_ratio)))




# full linear models for absorption at different wavelengths
model_Abs_302 <- lm(data = glucose_coded %>%
                      mutate(c_pct = ifelse(concentration == 1, 0.10, ifelse(concentration == -1, 0.05, 0.075))), 
                    formula = Abs_302 / (initial_mass * c_pct) ~ concentration * duration * pH * temperature)
Anova(model_Abs_302, type = "II")
summary(model_Abs_302)



model_Abs_315 <- lm(data = glucose_coded %>%
                      mutate(c_pct = ifelse(concentration == 1, 0.10, ifelse(concentration == -1, 0.05, 0.075))), 
                    formula = Abs_315 / (initial_mass * c_pct) ~ concentration * duration * pH * temperature)
Anova(model_Abs_315, type = "II")
summary(model_Abs_315)


model_Abs_365 <- lm(data = glucose_coded %>%
                      mutate(c_pct = ifelse(concentration == 1, 0.10, ifelse(concentration == -1, 0.05, 0.075))), 
                    formula = Abs_365 / (initial_mass * c_pct) ~ concentration * duration * pH * temperature)
Anova(model_Abs_365, type = "II")
summary(model_Abs_365)


# linear reduced model for Abs_max normalized to glucose loading in all reactors
model_Abs_max_red <- lm(data = glucose_coded %>%
                          mutate(c_pct = ifelse(concentration == 1, 0.10, ifelse(concentration == -1, 0.05, 0.075))), 
                        formula = Abs_max / (initial_mass * c_pct) ~ pH + concentration + duration * temperature + 
                          pH:temperature + concentration:temperature + concentration:duration)
Anova(model_Abs_max_red, type = "II")
summary(model_Abs_max_red)

# checking the residuals of the reduced model
model_Abs_max_red_res <- data.frame(residual = model_Abs_max_red$residuals,
                                         measured = model_Abs_max_red[["model"]][["Abs_max/(initial_mass * c_pct)"]],
                                         fitted = model_Abs_max_red$fitted.values)

ggplot(data = model_Abs_max_red_res) +
  geom_point(aes(x = measured, y = fitted)) +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1)

ggplot(data = model_Abs_max_red_res, aes(x = residual)) +
  geom_histogram(aes(y = stat(density)), bins = 15) +
  geom_density(adjust = 2, color = "red", size = 1)

ggplot(data = model_Abs_max_red_res) +
  geom_point(aes(x = measured, y = residual))


# reduced model excluding central point
model_Abs_max_red_no_central <- lm(data = glucose_coded %>%
                                          filter(concentration != 0) %>%
                                          mutate(c_pct = ifelse(concentration == 1, 0.10, 
                                                                ifelse(concentration == -1, 0.05, 0.075))), 
                                        formula = Abs_max / (initial_mass * c_pct) ~ pH + concentration + 
                                          duration * temperature + 
                                          pH:temperature + concentration:temperature + concentration:duration)

# mean Abs_max in the central point
experiment_central_Abs_max_red <- glucose_coded %>%
  filter(concentration == 0) %>%
  summarize(mean(Abs_max / (initial_mass * 0.075), na.rm = TRUE)) %>%
  as.numeric(.)

# prediction of the Abs_max in the central point using reduced linear model
predicted_central_Abs_max_red <- predict(model_Abs_max_red_no_central, 
                                      new = data.frame(concentration = 0,
                                                       duration = 0,
                                                       pH = 0,
                                                       temperature = 0), interval = "prediction")







# reduced linear model for absorption at 365 nm
model_Abs_365_red <- lm(data = glucose_coded %>%
                      mutate(c_pct = ifelse(concentration == 1, 0.10, ifelse(concentration == -1, 0.05, 0.075))), 
                    formula = Abs_365 / (initial_mass * c_pct) ~ concentration * duration * temperature +
                      pH + pH*temperature)
Anova(model_Abs_365_red, type = "II")
summary(model_Abs_365_red)

# checking the residuals of the reduced model
model_Abs_365_red_res <- data.frame(residual = model_Abs_365_red$residuals,
                                    measured = model_Abs_365_red[["model"]][["Abs_365/(initial_mass * c_pct)"]],
                                    fitted = model_Abs_365_red$fitted.values)

ggplot(data = model_Abs_365_red_res) +
  geom_point(aes(x = measured, y = fitted)) +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1)

ggplot(data = model_Abs_365_red_res, aes(x = residual)) +
  geom_histogram(aes(y = stat(density)), bins = 15) +
  geom_density(adjust = 2, color = "red", size = 1)

ggplot(data = model_Abs_365_red_res) +
  geom_point(aes(x = measured, y = residual))


# reduced model excluding central point
model_Abs_365_red_no_central <- lm(data = glucose_coded %>%
                                     filter(concentration != 0) %>%
                                     mutate(c_pct = ifelse(concentration == 1, 0.10, 
                                                           ifelse(concentration == -1, 0.05, 0.075))), 
                                   formula = Abs_365 / (initial_mass * c_pct) ~ pH + pH*temperature + 
                                     concentration * duration * temperature)
# mean Abs_max in the central point
experiment_central_Abs_365_red <- glucose_coded %>%
  filter(concentration == 0) %>%
  summarize(mean(Abs_365 / (initial_mass * 0.075), na.rm = TRUE)) %>%
  as.numeric(.)

# prediction of the Abs_max in the central point using reduced linear model
predicted_central_Abs_365_red <- predict(model_Abs_365_red_no_central, 
                                         new = data.frame(concentration = 0,
                                                          duration = 0,
                                                          pH = 0,
                                                          temperature = 0), interval = "prediction")









# plot of the interaction of 'concentration' and 'duration' factors for absorption at 365 nm
ggplot(data = glucose %>%
         filter(concentration != 7.5) %>%
         mutate(Abs_365 = Abs_365 / (initial_mass * concentration)) %>%
         group_by(temperature, duration) %>%
         summarize(Abs_365_mean = mean(Abs_365),
                   upper = Abs_365_mean + qt(0.975, df = length(Abs_365) - 1) * sd(Abs_365)/sqrt(length(Abs_365)),
                   lower = Abs_365_mean - qt(0.975, df = length(Abs_365) - 1) * sd(Abs_365)/sqrt(length(Abs_365))),
       aes(x = temperature, y = Abs_365_mean)) +
  geom_line(aes(color = factor(duration)), size = 1) +
  geom_point(aes(color = factor(duration)), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(duration)), width = 0.2) +
  xlab(label = "Temperature (C)") +
  ylab(label = "Absorbance at 365 nm, normalized") +
  labs(color = "Duration (h)") +
  theme(legend.position = c(0.45, 0.9))


# plot of the interaction of 'concentration' and 'duration' factors at low concentration for absorption at 365 nm
ggplot(data = glucose %>%
         filter(concentration != 7.5 & concentration == 5) %>%
         mutate(Abs_365 = Abs_365 / (initial_mass * concentration)) %>%
         group_by(temperature, duration) %>%
         summarize(Abs_365_mean = mean(Abs_365),
                   upper = Abs_365_mean + qt(0.975, df = length(Abs_365) - 1) * sd(Abs_365)/sqrt(length(Abs_365)),
                   lower = Abs_365_mean - qt(0.975, df = length(Abs_365) - 1) * sd(Abs_365)/sqrt(length(Abs_365))),
       aes(x = temperature, y = Abs_365_mean)) +
  geom_line(aes(color = factor(duration)), size = 1) +
  geom_point(aes(color = factor(duration)), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(duration)), width = 0.2) +
  xlab(label = "Temperature (C)") +
  ylab(label = "Absorbance at 365 nm, normalized") +
  labs(color = "Duration (h)") +
  theme(legend.position = c(0.45, 0.9))

# plot of the interaction of 'concentration' and 'duration' factors at high concentration for absorption at 365 nm
ggplot(data = glucose %>%
         filter(concentration != 7.5 & concentration == 10) %>%
         mutate(Abs_365 = Abs_365 / (initial_mass * concentration)) %>%
         group_by(temperature, duration) %>%
         summarize(Abs_365_mean = mean(Abs_365),
                   upper = Abs_365_mean + qt(0.975, df = length(Abs_365) - 1) * sd(Abs_365)/sqrt(length(Abs_365)),
                   lower = Abs_365_mean - qt(0.975, df = length(Abs_365) - 1) * sd(Abs_365)/sqrt(length(Abs_365))),
       aes(x = temperature, y = Abs_365_mean)) +
  geom_line(aes(color = factor(duration)), size = 1) +
  geom_point(aes(color = factor(duration)), size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(duration)), width = 0.2) +
  xlab(label = "Temperature (C)") +
  ylab(label = "Absorbance at 365 nm, normalized") +
  labs(color = "Duration (h)") +
  theme(legend.position = c(0.45, 0.9))


# contour plots of absorbance in response to synthesis conditions
time_axis = seq(-1.5, 1.5, 0.1)
temperature_axis = seq(-1.5 ,1.5, 0.1)

response_Abs_max = function(time, temperature, pH, conc){coef(model_Abs_max_red)["(Intercept)"] + 
    coef(model_Abs_max_red)["pH"]*pH + 
    coef(model_Abs_max_red)["concentration"]*conc + 
    coef(model_Abs_max_red)["duration"]*time + 
    coef(model_Abs_max_red)["temperature"]*temperature + 
    coef(model_Abs_max_red)["duration:temperature"]*time*temperature +
    coef(model_Abs_max_red)["pH:temperature"]*pH*temperature +
    coef(model_Abs_max_red)["concentration:temperature"]*conc*temperature +
    coef(model_Abs_max_red)["concentration:duration"]*conc*time}

Abs_max_predicted <- merge(x = time_axis, y = temperature_axis, all.x=TRUE, all.y=TRUE) %>%
  rename(duration = x, temperature = y) %>%
  merge(y = c(-1, 0, 1)) %>%
  rename(concentration = y) %>%
  merge(y = c(-1, 0, 1)) %>%
  rename(pH = y) %>%
  mutate(Abs_max = response_Abs_max(time = duration, temperature = temperature, pH = pH, conc = concentration))

ggplot(data = Abs_max_predicted, aes(x = duration, y = temperature, z = Abs_max)) +
  geom_density_2d() +
  geom_contour_filled() +
  facet_grid(pH ~ concentration, labeller = label_both)


response_Abs_365 = function(time, temperature, pH, conc){coef(model_Abs_365_red)["(Intercept)"] + 
    coef(model_Abs_365_red)["pH"]*pH + 
    coef(model_Abs_365_red)["concentration"]*conc + 
    coef(model_Abs_365_red)["duration"]*time + 
    coef(model_Abs_365_red)["temperature"]*temperature + 
    coef(model_Abs_365_red)["duration:temperature"]*time*temperature +
    coef(model_Abs_365_red)["temperature:pH"]*pH*temperature +
    coef(model_Abs_365_red)["concentration:temperature"]*conc*temperature +
    coef(model_Abs_365_red)["concentration:duration"]*conc*time +
    coef(model_Abs_365_red)["concentration:duration:temperature"]*conc*time*temperature}

Abs_365_predicted <- merge(x = time_axis, y = temperature_axis, all.x=TRUE, all.y=TRUE) %>%
  rename(duration = x, temperature = y) %>%
  merge(y = c(-1, 0, 1)) %>%
  rename(concentration = y) %>%
  merge(y = c(-1, 0, 1)) %>%
  rename(pH = y) %>%
  mutate(Abs_max = response_Abs_365(time = duration, temperature = temperature, pH = 0, conc = concentration))

ggplot(data = Abs_365_predicted, aes(x = duration, y = temperature, z = Abs_max)) +
  geom_density_2d() +
  geom_contour_filled() +
  facet_grid(pH ~ concentration, labeller = label_both)


# linear model with absorption band maximum position as response variable
model_wl_max <- lm(data = glucose_coded, 
                   formula = wavelength_max ~ concentration * duration * pH * temperature)
Anova(model_wl_max, type = "II")
summary(model_wl_max)
#==============END OF ANALYSIS OF ABSORBANCE===========================



#==============ANALYSIS OF DATASET SIZE================================
Abs_365_ratio_ref <- glucose_coded %>%
  select(id, concentration, duration, pH, temperature, Abs_365_norm) %>%
  rename(response = Abs_365_norm)

model_Abs_365_ratio_ref <- lm(data = Abs_365_ratio_ref, 
                        formula = response ~ concentration * duration * pH * temperature)
Anova(model_Abs_365_ratio_ref, type = "II")
summary(model_Abs_365_ratio_ref)

# this function creates a dataframe with 'repl' replicates at each experimental point 
# and 'repl_c' replicates of the central point
subset_empty <- function(repl, repl_c){
  single_replicate <- merge(x = c(-1, 1), y = c(-1, 1), all.x=TRUE, all.y=TRUE) %>%
    rename(temperature = x, duration = y) %>%
    merge(y = c(-1, 1), all.x=TRUE, all.y=TRUE) %>%
    rename(concentration = y) %>%
    merge(y = c(-1, 0, 1), all.x=TRUE, all.y=TRUE) %>%
    rename(pH = y) %>%
    mutate(response = NA)
  
  res <- NULL
  
  for (i in 1:repl){
    res <- rbind(res, single_replicate)
  }
  
  for (i in 1:repl_c){
    res <- rbind(res, data.frame(temperature = 0, duration = 0, concentration = 0, pH = 0, response = NA))
  }
  
  return(res)
}



# sampling the reference dataframe to fill the 'df' subset with response values
subset_sample <- function(df, ref, mult_sample = FALSE) {
  
  for (i in 1:nrow(df)){
    sample_row <- ref %>%
      filter(duration == df$duration[i] & 
               temperature == df$temperature[i] &
               pH == df$pH[i] &
               concentration == df$concentration[i]) %>%
      sample_n(., 1)
    df$response[i] <- sample_row$response[1]
    if (!mult_sample){
      ref <- ref %>%
        filter(id != sample_row$id[1])
    }
    
  }
  
  return(df)
}



# df_ref is a dataframe with 'reference' results. It should include the following columns:
# 'id', 'duration', 'temperature', 'concentration', 'pH', 'response'
# repl is a number of replicates of each design point (1:3)
# repl_c is a number of replicates of the central design point (1:8)
# runs is a number of differents subsets to test the partial models
test_subset <- function(df_ref, repl, repl_c, mult_sample = FALSE, runs, alpha = 0.05){
  set.seed(1)
  model <- lm(data = df_ref, 
              formula = response ~ concentration * duration * pH * temperature)
  
  res_signif <- NULL
  res_diff <- NULL
  res_diff_signif <- NULL
  signif_ref <- (Anova(model)["Pr(>F)"] < alpha)
  coef_ref <- coef(model)
  coef_ref_signif <- summary(model)[["coefficients"]][, "Pr(>|t|)"] < alpha
  
  for (i in 1:runs){
    subset_empty <- subset_empty(repl = repl, repl_c = repl_c)
    subset_filled <- subset_sample(df = subset_empty, ref = df_ref, mult_sample)
    
    model_sample <- lm(data = subset_filled, formula = response ~ concentration * duration * pH * temperature)
    
    signif_sample <- (Anova(model_sample)["Pr(>F)"] < alpha)
    
    res_signif <- c(res_signif, sum(signif_ref[-length(signif_ref)] == signif_sample[-length(signif_sample)]))
    
    coef_sample <- coef(model_sample)
    
    coef_diff <- mean(abs((coef_sample - coef_ref) / coef_ref))
    coef_diff_signif <- mean(abs((coef_sample[coef_ref_signif] - coef_ref[coef_ref_signif]) / coef_ref[coef_ref_signif]))
    
    res_diff <- c(res_diff, coef_diff)
    res_diff_signif <- c(res_diff_signif, coef_diff_signif)
    
  }
  # res_signif is a vector with number of factors correctly recognized as significant
  # res_diff is a vector of mean absolute difference between the model parameters
  # res_diff_signif is the same but includes only the coefficients considered significant in the reference model
  print(paste0("i = ", i))
  print(table(res_signif))
  print(median(res_diff))
  print(median(res_diff_signif))
  
  return(list(res_signif, res_diff, res_diff_signif))
}



Abs_365_ratio_ref <- glucose_coded %>%
  select(id, concentration, duration, pH, temperature, Abs_365_norm) %>%
  rename(response = Abs_365_norm)

model_Abs_365_ratio_ref <- lm(data = Abs_365_ratio_ref, 
                              formula = response ~ concentration * duration * pH * temperature)
Anova(model_Abs_365_ratio_ref, type = "II")
summary(model_Abs_365_ratio_ref)

set.seed(1)


test_ratio_3_8_01 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 3, repl_c = 8, mult_sample = TRUE, 
                                 runs = 1000, alpha = 0.01)

test_ratio_3_8_05 <- test_subset(df_ref = Abs_365_ratio_ref, 
                              repl = 3, repl_c = 8, mult_sample = TRUE, 
                              runs = 1000, alpha = 0.05)

test_ratio_3_8_10 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 3, repl_c = 8, mult_sample = TRUE, 
                                 runs = 1000, alpha = 0.10)



test_ratio_2_8_01 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 2, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_ratio_2_8_05 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 2, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_ratio_2_8_10 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 2, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)



test_ratio_1_8_01 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_ratio_1_8_05 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_ratio_1_8_10 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)




test_ratio_1_5_01 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 5, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_ratio_1_5_05 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 5, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_ratio_1_5_10 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 5, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)


test_ratio_1_2_01 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 2, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_ratio_1_2_05 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 2, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_ratio_1_2_10 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 2, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)


test_ratio_1_0_01 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 0, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_ratio_1_0_05 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 0, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_ratio_1_0_10 <- test_subset(df_ref = Abs_365_ratio_ref, 
                                 repl = 1, repl_c = 0, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)


Abs_365_ref <- glucose_coded %>%
  mutate(c_pct = ifelse(concentration == 1, 0.10, ifelse(concentration == -1, 0.05, 0.075)),
         Abs_365 = Abs_365 / (initial_mass * c_pct / 100)) %>%
  select(id, concentration, duration, pH, temperature, Abs_365) %>%
  rename(response = Abs_365)


test_abs_3_8_01 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 3, repl_c = 8, mult_sample = TRUE, 
                                 runs = 1000, alpha = 0.01)

test_abs_3_8_05 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 3, repl_c = 8, mult_sample = TRUE, 
                                 runs = 1000, alpha = 0.05)

test_abs_3_8_10 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 3, repl_c = 8, mult_sample = TRUE, 
                                 runs = 1000, alpha = 0.10)



test_abs_2_8_01 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 2, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_abs_2_8_05 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 2, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_abs_2_8_10 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 2, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)



test_abs_1_8_01 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_abs_1_8_05 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_abs_1_8_10 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 8, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)




test_abs_1_5_01 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 5, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_abs_1_5_05 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 5, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_abs_1_5_10 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 5, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)


test_abs_1_2_01 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 2, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_abs_1_2_05 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 2, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_abs_1_2_10 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 2, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)


test_abs_1_0_01 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 0, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.01)

test_abs_1_0_05 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 0, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.05)

test_abs_1_0_10 <- test_subset(df_ref = Abs_365_ref, 
                                 repl = 1, repl_c = 0, mult_sample = FALSE, 
                                 runs = 1000, alpha = 0.10)


test_abs_2_5_01 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 5, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.01)

test_abs_2_5_05 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 5, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.05)

test_abs_2_5_10 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 5, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.10)


test_abs_2_2_01 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 2, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.01)

test_abs_2_2_05 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 2, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.05)

test_abs_2_2_10 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 2, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.10)


test_abs_2_0_01 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 0, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.01)

test_abs_2_0_05 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 0, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.05)

test_abs_2_0_10 <- test_subset(df_ref = Abs_365_ref, 
                               repl = 2, repl_c = 0, mult_sample = FALSE, 
                               runs = 1000, alpha = 0.10)





temp <- test_subset(df_ref = Abs_365_ref, 
                               repl = 3, repl_c = 8, mult_sample = FALSE, 
                               runs = 10, alpha = 0.05)


correct_factors_abs <- data.frame(cbind(test_abs_3_8_01[[1]], test_abs_3_8_05[[1]], test_abs_3_8_10[[1]],
              test_abs_2_8_01[[1]], test_abs_2_8_05[[1]], test_abs_2_8_10[[1]],
              test_abs_1_8_01[[1]], test_abs_1_8_05[[1]], test_abs_1_8_10[[1]],
              test_abs_1_5_01[[1]], test_abs_1_5_05[[1]], test_abs_1_5_10[[1]],
              test_abs_1_2_01[[1]], test_abs_1_2_05[[1]], test_abs_1_2_10[[1]],
              test_abs_1_0_01[[1]], test_abs_1_0_05[[1]], test_abs_1_0_10[[1]],
              test_abs_2_5_01[[1]], test_abs_2_5_05[[1]], test_abs_2_5_10[[1]],
              test_abs_2_2_01[[1]], test_abs_2_2_05[[1]], test_abs_2_2_10[[1]],
              test_abs_2_0_01[[1]], test_abs_2_0_05[[1]], test_abs_2_0_10[[1]])) %>%
  `colnames<-`(c("abs_3_8_01", "abs_3_8_05", "abs_3_8_10",
                 "abs_2_8_01", "abs_2_8_05", "abs_2_8_10",
                 "abs_1_8_01", "abs_1_8_05", "abs_1_8_10",
                 "abs_1_5_01", "abs_1_5_05", "abs_1_5_10",
                 "abs_1_2_01", "abs_1_2_05", "abs_1_2_10",
                 "abs_1_0_01", "abs_1_0_05", "abs_1_0_10",
                 "abs_2_5_01", "abs_2_5_05", "abs_2_5_10",
                 "abs_2_2_01", "abs_2_2_05", "abs_2_2_10",
                 "abs_2_0_01", "abs_2_0_05", "abs_2_0_10")) %>%
  pivot_longer(cols = everything(), names_to = "test", values_to = "factors") %>%
  separate(test, c(NA, "repl", "repl_c", "alpha")) %>%
  mutate(repl = as.numeric(repl),
         repl_c = as.numeric(repl_c),
         alpha = as.numeric(alpha) / 100)


ggplot(data = correct_factors_abs %>% filter(repl_c == 8) %>%
         mutate(replicates = factor(repl, levels = c(3, 2, 1), labels = c("3", "2", "1")))) +
  facet_grid(rows = vars(replicates),
             labeller = label_both,
             scales = "free_y") +
  geom_bar(aes(x = factors, fill = factor(alpha)), 
           position = position_dodge(preserve = "single"), 
           width = 0.6 ) +
  xlim(c(6, 16)) +
  xlab("Correctly resolved factors") +
  theme(legend.position = c(0.15, 0.85), legend.title = element_markdown()) +
  labs(fill = "Threshold *p*-value")


ggplot(data = correct_factors_abs %>% filter(repl == 1) %>%
         mutate(replicates = factor(repl_c, levels = c(8, 5, 2, 0), labels = c("8", "5", "2", "0")))) +
  facet_grid(rows = vars(replicates),
             labeller = label_both,
             scales = "free_y") +
  geom_bar(aes(x = factors, fill = factor(alpha)), 
           position = position_dodge(preserve = "single"), 
           width = 0.6 ) +
  xlim(c(6, 16)) +
  xlab("Correctly resolved factors") +
  theme(legend.position = c(0.85, 0.9), legend.title = element_markdown()) +
  labs(fill = "Threshold *p*-value")


ggplot(data = correct_factors_abs %>% filter(repl == 2) %>%
         mutate(replicates = factor(repl_c, levels = c(8, 5, 2, 0), labels = c("8", "5", "2", "0")))) +
  facet_grid(rows = vars(replicates),
             labeller = label_both,
             scales = "free_y") +
  geom_bar(aes(x = factors, fill = factor(alpha)), 
           position = position_dodge(preserve = "single"), 
           width = 0.6 ) +
  xlim(c(8, 16)) +
  xlab("Correctly resolved factors") +
  theme(legend.position = c(0.15, 0.9), legend.title = element_markdown()) +
  labs(fill = "Threshold *p*-value")




mean_error_abs <- data.frame(cbind(test_abs_3_8_01[[2]], test_abs_3_8_05[[2]], test_abs_3_8_10[[2]],
                                   test_abs_2_8_01[[2]], test_abs_2_8_05[[2]], test_abs_2_8_10[[2]],
                                   test_abs_1_8_01[[2]], test_abs_1_8_05[[2]], test_abs_1_8_10[[2]],
                                   test_abs_1_5_01[[2]], test_abs_1_5_05[[2]], test_abs_1_5_10[[2]],
                                   test_abs_1_2_01[[2]], test_abs_1_2_05[[2]], test_abs_1_2_10[[2]],
                                   test_abs_1_0_01[[2]], test_abs_1_0_05[[2]], test_abs_1_0_10[[2]],
                                   test_abs_2_5_01[[2]], test_abs_2_5_05[[2]], test_abs_2_5_10[[2]],
                                   test_abs_2_2_01[[2]], test_abs_2_2_05[[2]], test_abs_2_2_10[[2]],
                                   test_abs_2_0_01[[2]], test_abs_2_0_05[[2]], test_abs_2_0_10[[2]])) %>%
  `colnames<-`(c("abs_3_8_01", "abs_3_8_05", "abs_3_8_10",
                 "abs_2_8_01", "abs_2_8_05", "abs_2_8_10",
                 "abs_1_8_01", "abs_1_8_05", "abs_1_8_10",
                 "abs_1_5_01", "abs_1_5_05", "abs_1_5_10",
                 "abs_1_2_01", "abs_1_2_05", "abs_1_2_10",
                 "abs_1_0_01", "abs_1_0_05", "abs_1_0_10",
                 "abs_2_5_01", "abs_2_5_05", "abs_2_5_10",
                 "abs_2_2_01", "abs_2_2_05", "abs_2_2_10",
                 "abs_2_0_01", "abs_2_0_05", "abs_2_0_10")) %>%
  pivot_longer(cols = everything(), names_to = "test", values_to = "error") %>%
  separate(test, c(NA, "repl", "repl_c", "alpha"), remove = FALSE) %>%
  mutate(repl = as.numeric(repl),
         repl_c = as.numeric(repl_c),
         alpha = as.numeric(alpha) / 100)


ggplot(data = mean_error_abs %>% 
         filter(alpha == 0.05 & (repl_c == 0 | repl_c == 8))) +
  geom_violin(aes(x = test, y = error)) +
  xlab("Subset: DR; CR") +
  ylab("Mean relative error") +
  scale_x_discrete(labels = c('1; 0',
                              '1; 8',
                              '2; 0', 
                              '2; 8',
                              '3; 8'))



mean_error_sign_abs <- data.frame(cbind(test_abs_3_8_01[[3]], test_abs_3_8_05[[3]], test_abs_3_8_10[[3]],
                                   test_abs_2_8_01[[3]], test_abs_2_8_05[[3]], test_abs_2_8_10[[3]],
                                   test_abs_1_8_01[[3]], test_abs_1_8_05[[3]], test_abs_1_8_10[[3]],
                                   test_abs_1_5_01[[3]], test_abs_1_5_05[[3]], test_abs_1_5_10[[3]],
                                   test_abs_1_2_01[[3]], test_abs_1_2_05[[3]], test_abs_1_2_10[[3]],
                                   test_abs_1_0_01[[3]], test_abs_1_0_05[[3]], test_abs_1_0_10[[3]],
                                   test_abs_2_5_01[[3]], test_abs_2_5_05[[3]], test_abs_2_5_10[[3]],
                                   test_abs_2_2_01[[3]], test_abs_2_2_05[[3]], test_abs_2_2_10[[3]],
                                   test_abs_2_0_01[[3]], test_abs_2_0_05[[3]], test_abs_2_0_10[[3]])) %>%
  `colnames<-`(c("abs_3_8_01", "abs_3_8_05", "abs_3_8_10",
                 "abs_2_8_01", "abs_2_8_05", "abs_2_8_10",
                 "abs_1_8_01", "abs_1_8_05", "abs_1_8_10",
                 "abs_1_5_01", "abs_1_5_05", "abs_1_5_10",
                 "abs_1_2_01", "abs_1_2_05", "abs_1_2_10",
                 "abs_1_0_01", "abs_1_0_05", "abs_1_0_10",
                 "abs_2_5_01", "abs_2_5_05", "abs_2_5_10",
                 "abs_2_2_01", "abs_2_2_05", "abs_2_2_10",
                 "abs_2_0_01", "abs_2_0_05", "abs_2_0_10")) %>%
  pivot_longer(cols = everything(), names_to = "test", values_to = "error") %>%
  separate(test, c(NA, "repl", "repl_c", "alpha"), remove = FALSE) %>%
  mutate(repl = as.numeric(repl),
         repl_c = as.numeric(repl_c),
         alpha = as.numeric(alpha) / 100)


ggplot(data = mean_error_sign_abs %>% 
         filter(alpha == 0.05 & (repl_c == 0 | repl_c == 8))) +
  geom_violin(aes(x = test, y = error)) +
  xlab("Subset: DR; CR") +
  ylab("Mean relative error") +
  scale_x_discrete(labels = c('1; 0',
                              '1; 8',
                              '2; 0', 
                              '2; 8',
                              '3; 8'))












correct_factors_ratio <- data.frame(cbind(test_ratio_3_8_01[[1]], test_ratio_3_8_05[[1]], test_ratio_3_8_10[[1]],
                                        test_ratio_2_8_01[[1]], test_ratio_2_8_05[[1]], test_ratio_2_8_10[[1]],
                                        test_ratio_1_8_01[[1]], test_ratio_1_8_05[[1]], test_ratio_1_8_10[[1]],
                                        test_ratio_1_5_01[[1]], test_ratio_1_5_05[[1]], test_ratio_1_5_10[[1]],
                                        test_ratio_1_2_01[[1]], test_ratio_1_2_05[[1]], test_ratio_1_2_10[[1]],
                                        test_ratio_1_0_01[[1]], test_ratio_1_0_05[[1]], test_ratio_1_0_10[[1]])) %>%
  `colnames<-`(c("ratio_3_8_01", "ratio_3_8_05", "ratio_3_8_10",
                 "ratio_2_8_01", "ratio_2_8_05", "ratio_2_8_10",
                 "ratio_1_8_01", "ratio_1_8_05", "ratio_1_8_10",
                 "ratio_1_5_01", "ratio_1_5_05", "ratio_1_5_10",
                 "ratio_1_2_01", "ratio_1_2_05", "ratio_1_2_10",
                 "ratio_1_0_01", "ratio_1_0_05", "ratio_1_0_10")) %>%
  pivot_longer(cols = everything(), names_to = "test", values_to = "factors") %>%
  separate(test, c(NA, "repl", "repl_c", "alpha")) %>%
  mutate(repl = as.numeric(repl),
         repl_c = as.numeric(repl_c),
         alpha = as.numeric(alpha) / 100)




ggplot(data = correct_factors_ratio %>% filter(repl_c == 8) %>%
         mutate(replicates = factor(repl, levels = c(3, 2, 1), labels = c("3", "2", "1")))) +
  facet_grid(rows = vars(replicates),
             labeller = label_both,
             scales = "free_y") +
  geom_bar(aes(x = factors, fill = factor(alpha)), 
           position = position_dodge(preserve = "single"), 
           width = 0.6 ) +
  xlim(c(8, 16)) +
  xlab("Correctly resolved factors") +
  theme(legend.position = c(0.15, 0.85), legend.title = element_markdown()) +
  labs(fill = "Threshold *p*-value")


ggplot(data = correct_factors_ratio %>% filter(repl == 1) %>%
         mutate(replicates = factor(repl_c, levels = c(8, 5, 2, 0), labels = c("8", "5", "2", "0")))) +
  facet_grid(rows = vars(replicates),
             labeller = label_both,
             scales = "free_y") +
  geom_bar(aes(x = factors, fill = factor(alpha)), 
           position = position_dodge(preserve = "single"), 
           width = 0.6 ) +
  xlim(c(10, 16)) +
  xlab("Correctly resolved factors") +
  theme(legend.position = c(0.15, 0.85), legend.title = element_markdown()) +
  labs(fill = "Threshold *p*-value")



mean_error_ratio <- data.frame(cbind(test_ratio_3_8_01[[2]], test_ratio_3_8_05[[2]], test_ratio_3_8_10[[2]],
                                   test_ratio_2_8_01[[2]], test_ratio_2_8_05[[2]], test_ratio_2_8_10[[2]],
                                   test_ratio_1_8_01[[2]], test_ratio_1_8_05[[2]], test_ratio_1_8_10[[2]],
                                   test_ratio_1_5_01[[2]], test_ratio_1_5_05[[2]], test_ratio_1_5_10[[2]],
                                   test_ratio_1_2_01[[2]], test_ratio_1_2_05[[2]], test_ratio_1_2_10[[2]],
                                   test_ratio_1_0_01[[2]], test_ratio_1_0_05[[2]], test_ratio_1_0_10[[2]])) %>%
  `colnames<-`(c("ratio_3_8_01", "ratio_3_8_05", "ratio_3_8_10",
                 "ratio_2_8_01", "ratio_2_8_05", "ratio_2_8_10",
                 "ratio_1_8_01", "ratio_1_8_05", "ratio_1_8_10",
                 "ratio_1_5_01", "ratio_1_5_05", "ratio_1_5_10",
                 "ratio_1_2_01", "ratio_1_2_05", "ratio_1_2_10",
                 "ratio_1_0_01", "ratio_1_0_05", "ratio_1_0_10")) %>%
  pivot_longer(cols = everything(), names_to = "test", values_to = "error") %>%
  separate(test, c(NA, "repl", "repl_c", "alpha"), remove = FALSE) %>%
  mutate(repl = as.numeric(repl),
         repl_c = as.numeric(repl_c),
         alpha = as.numeric(alpha) / 100)


ggplot(data = mean_error_ratio %>% 
         filter(alpha == 0.05 & (repl_c == 0 | repl_c == 8))) +
  geom_violin(aes(x = test, y = error)) +
  xlab("Subset: DR; CR") +
  ylab("Mean relative error") +
  scale_x_discrete(labels = c('1; 0',
                              '1; 8',
                              '2; 8',
                              '3; 8'))




mean_error_sign_ratio <- data.frame(cbind(test_ratio_3_8_01[[3]], test_ratio_3_8_05[[3]], test_ratio_3_8_10[[3]],
                                        test_ratio_2_8_01[[3]], test_ratio_2_8_05[[3]], test_ratio_2_8_10[[3]],
                                        test_ratio_1_8_01[[3]], test_ratio_1_8_05[[3]], test_ratio_1_8_10[[3]],
                                        test_ratio_1_5_01[[3]], test_ratio_1_5_05[[3]], test_ratio_1_5_10[[3]],
                                        test_ratio_1_2_01[[3]], test_ratio_1_2_05[[3]], test_ratio_1_2_10[[3]],
                                        test_ratio_1_0_01[[3]], test_ratio_1_0_05[[3]], test_ratio_1_0_10[[3]])) %>%
  `colnames<-`(c("ratio_3_8_01", "ratio_3_8_05", "ratio_3_8_10",
                 "ratio_2_8_01", "ratio_2_8_05", "ratio_2_8_10",
                 "ratio_1_8_01", "ratio_1_8_05", "ratio_1_8_10",
                 "ratio_1_5_01", "ratio_1_5_05", "ratio_1_5_10",
                 "ratio_1_2_01", "ratio_1_2_05", "ratio_1_2_10",
                 "ratio_1_0_01", "ratio_1_0_05", "ratio_1_0_10")) %>%
  pivot_longer(cols = everything(), names_to = "test", values_to = "error") %>%
  separate(test, c(NA, "repl", "repl_c", "alpha"), remove = FALSE) %>%
  mutate(repl = as.numeric(repl),
         repl_c = as.numeric(repl_c),
         alpha = as.numeric(alpha) / 100)

ggplot(data = mean_error_sign_ratio %>% 
         filter(alpha == 0.05 & (repl_c == 0 | repl_c == 8))) +
  geom_violin(aes(x = test, y = error)) +
  xlab("Subset: DR; CR") +
  ylab("Mean relative error") +
  scale_x_discrete(labels = c('1; 0',
                              '1; 8',
                              '2; 8',
                              '3; 8'))




















table(test_ratio_3_8[1])
test_ratio_3_8[2]
test_ratio_3_8[3]

test_ratio_2_8 <- test_subset(df_ref = Abs_365_ratio_ref, 
                              repl = 2, repl_c = 8, mult_sample = FALSE, 
                              runs = 100, alpha = 0.01)



temp_2 <- test_subset(df_ref = Abs_365_ratio_ref, repl = 1, repl_c = 0, runs = 100, alpha = 0.01)























#temp <- subset_empty(repl = 1, repl_c = 8)
#temp_2 <- subset_sample(df = temp, ref = Abs_365_ratio_ref, mult_sample = FALSE)

#model_Abs_365_ratio_sample <- lm(data = temp_2, 
#                           formula = response ~ concentration * duration * pH * temperature)
#Anova(model_Abs_365_ratio_sample, type = "II")
#summary(model_Abs_365_ratio_sample)


# signif_ref <- (Anova(model_Abs_365_ratio_ref)["Pr(>F)"] < 0.05)
# signif_sample <- (Anova(model_Abs_365_ratio_sample)["Pr(>F)"] < 0.05)

# sum(signif_ref[-length(signif_ref)] == signif_sample[-length(signif_sample)])

# sum((Anova(model_Abs_365_ratio_ref)["Pr(>F)"] < 0.05) == (Anova(model_Abs_365_ratio_sample)["Pr(>F)"] < 0.05))


res_signif <- NULL
res_diff <- NULL
res_diff_signif <- NULL
signif_ref <- (Anova(model_Abs_365_ratio_ref)["Pr(>F)"] < 0.05)
coef_ref <- coef(model_Abs_365_ratio_ref)
coef_ref_signif <- summary(model_Abs_365_ratio_ref)[["coefficients"]][, "Pr(>|t|)"] < 0.05

for (i in 1:10000){
  temp <- subset_empty(repl = 1, repl_c = 3)
  temp_2 <- subset_sample(df = temp, ref = Abs_365_ratio_ref, mult_sample = FALSE)
  
  model_Abs_365_ratio_sample <- lm(data = temp_2, 
                                   formula = response ~ concentration * duration * pH * temperature)
  
  signif_sample <- (Anova(model_Abs_365_ratio_sample)["Pr(>F)"] < 0.05)
  
  res_signif <- c(res_signif, sum(signif_ref[-length(signif_ref)] == signif_sample[-length(signif_sample)]))
  
  coef_sample <- coef(model_Abs_365_ratio_sample)
  
  coef_diff <- mean(abs((coef_sample - coef_ref) / coef_ref))
  coef_diff_signif <- mean(abs((coef_sample[coef_ref_signif] - coef_ref[coef_ref_signif]) / coef_ref[coef_ref_signif]))
  
  res_diff <- c(res_diff, coef_diff)
  res_diff_signif <- c(res_diff_signif, coef_diff_signif)
  
  if (!(i %% 50)){
    print(paste0("i = ", i))
    print(table(res_signif))
    print(median(res_diff))
    print(median(res_diff_signif))
  }
}

table(res_signif)


# we will consider the absorbance at 365 nm
Abs_365_ref <- glucose_coded %>%
  mutate(c_pct = ifelse(concentration == 1, 0.10, ifelse(concentration == -1, 0.05, 0.075)),
         Abs_365 = Abs_365 / (initial_mass * c_pct / 100)) %>%
  select(id, concentration, duration, pH, temperature, Abs_365) %>%
  rename(response = Abs_365)

model_Abs_365_ref <- lm(data = Abs_365_ref, 
                    formula = response ~ concentration * duration * pH * temperature)
Anova(model_Abs_365_ref, type = "II")
summary(model_Abs_365_ref)


res_signif <- NULL
res_diff <- NULL
res_diff_signif <- NULL
signif_ref <- (Anova(model_Abs_365_ref)["Pr(>F)"] < 0.05)
coef_ref <- coef(model_Abs_365_ref)
coef_ref_signif <- summary(model_Abs_365_ref)[["coefficients"]][, "Pr(>|t|)"] < 0.05

for (i in 1:10000){
  temp <- subset_empty(repl = 1, repl_c = 3)
  temp_2 <- subset_sample(df = temp, ref = Abs_365_ref, mult_sample = FALSE)
  
  model_Abs_365_sample <- lm(data = temp_2, 
                                   formula = response ~ concentration * duration * pH * temperature)
  
  signif_sample <- (Anova(model_Abs_365_sample)["Pr(>F)"] < 0.05)
  
  res_signif <- c(res_signif, sum(signif_ref[-length(signif_ref)] == signif_sample[-length(signif_sample)]))
  
  coef_sample <- coef(model_Abs_365_sample)
  
  coef_diff <- mean(abs((coef_sample - coef_ref) / coef_ref))
  coef_diff_signif <- mean(abs((coef_sample[coef_ref_signif] - coef_ref[coef_ref_signif]) / coef_ref[coef_ref_signif]))
  
  res_diff <- c(res_diff, coef_diff)
  res_diff_signif <- c(res_diff_signif, coef_diff_signif)
  
  if (!(i %% 50)){
    print(paste0("i = ", i))
    print(table(res_signif))
    print(median(res_diff))
    print(median(res_diff_signif))
  }
}

table(res_signif)













# this function creates a dataframe with 'repl' replicates at each experimental point 
# and 'repl_c' replicates of the central point
Abs_365_empty <- function(repl, repl_c){
  single_replicate <- merge(x = c(-1, 1), y = c(-1, 1), all.x=TRUE, all.y=TRUE) %>%
    rename(temperature = x, duration = y) %>%
    merge(y = c(-1, 1), all.x=TRUE, all.y=TRUE) %>%
    rename(concentration = y) %>%
    merge(y = c(-1, 0, 1), all.x=TRUE, all.y=TRUE) %>%
    rename(pH = y) %>%
    mutate(Abs_365 = NA)
  
  res <- NULL
  
  for (i in 1:repl){
    res <- rbind(res, single_replicate)
  }
  
  for (i in 1:repl_c){
    res <- rbind(res, data.frame(temperature = 0, duration = 0, concentration = 0, pH = 0, Abs_365 = NA))
  }

  return(res)
}

temp <- Abs_365_empty(repl = 1, repl_c = 3)

# sampling the reference dataframe to fill the 'df' subset with Abs_365 values
Abs_365_sample <- function(df, ref = Abs_365_ref, mult_sample = TRUE) {
  Abs_set <- ref
  
  for (i in 1:nrow(df)){
    sample_row <- Abs_set %>%
      filter(duration == df$duration[i] & 
               temperature == df$temperature[i] &
               pH == df$pH[i] &
               concentration == df$concentration[i]) %>%
      sample_n(., 1)
    df$Abs_365[i] <- sample_row$Abs_365[1]
    if (!mult_sample){
      Abs_set <- Abs_set %>%
        filter(id != sample_row$id[1])
    }
    
  }
  
  return(df)
}

temp_2 <- Abs_365_sample(df = temp, mult_sample = TRUE)

model_Abs_365_sample <- lm(data = temp_2, 
                        formula = Abs_365 ~ concentration * duration * pH * temperature)
Anova(model_Abs_365_sample, type = "II")
summary(model_Abs_365_sample)


(coef(model_Abs_365_ref) - coef(model_Abs_365_sample)) < 1e-12



sum(Abs_365_ref$Abs_365)
sum(temp_2$Abs_365)

#==============END OF ANALYSIS OF DATASET SIZE=========================




#=================== OLD STUFF








library(daewr)

effects <-coef(model_pH)
effects <- effects[c(2:length(effects))]
halfnorm(effects,names(effects),alpha=.15)

model_pH <- lm(data = glucose_pH, formula = pH_final ~ concentration + duration + temperature + pH + concentration:pH + duration:temperature + pH:temperature + concentration:pH:temperature)
Anova(model_pH, type = "II")
summary(model_pH)








model_pH <- lm(data = glucose_pH, formula = pH_final ~ concentration * duration * pH * temperature)

model_null <- lm(data = glucose_pH, formula = pH_final ~ 1)

model_pH <- step(object = model_null, scope = model_pH$call, direction = "both", trace = 100, k = 2)
Anova(model_pH, type = "II")
summary(model_pH)



residuals_pH <- data.frame(pH_final = model_pH$model$pH_final,
                           residuals = model_pH$residuals)

ggplot(residuals_pH, aes(x = pH_final, y = residuals)) + 
  geom_point() +
  geom_smooth()

qqPlot(residuals_pH$residuals)



predicted_central_pH <- predict(model_pH, new = data.frame(concentration = 0,
                                                           duration = 0,
                                                           pH = 0,
                                                           temperature = 0), interval = "prediction")

experiment_central_pH <- glucose_coded %>%
  filter(concentration == 0 & duration == 0 & pH == 0 & temperature == 0) %>%
  summarize(mean(pH_final, na.rm = TRUE)) %>%
  as.numeric(.)


model_pH <- lm(data = glucose_coded, formula = model_pH$call)
Anova(model_pH, type = "II")
summary(model_pH)

library(DAAG)

glucose_pH_cv <- cv.lm(data = model_pH$model, form.lm = formula(model_pH), m = 5, dots = 0, seed = 123, plotit = FALSE, printit = TRUE)

ggplot(data = glucose_pH_cv, aes(x = pH_final, y = cvpred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)




library(GGally)
ggpairs(glucose_full, columns = c(8:11, 14) , ggplot2::aes(colour=factor(volume)))



library(car)



library(daewr)

effects <-coef(model_pH)
effects <- effects[c(2:length(effects))]
halfnorm(effects,names(effects),alpha=.15)

model_pH <- lm(data = glucose_pH, formula = pH_final ~ concentration + duration + temperature + pH + concentration:pH + duration:temperature + pH:temperature + concentration:pH:temperature)
Anova(model_pH, type = "II")
summary(model_pH)


model_pH <- lm(data = glucose_pH, formula = pH_final ~ concentration * duration * pH * temperature)

model_null <- lm(data = glucose_pH, formula = pH_final ~ 1)

model_pH <- step(object = model_null, scope = model_pH$call, direction = "both", trace = 100, k = 2)
Anova(model_pH, type = "II")
summary(model_pH)



residuals_pH <- data.frame(pH_final = model_pH$model$pH_final,
                           residuals = model_pH$residuals)

ggplot(residuals_pH, aes(x = pH_final, y = residuals)) + 
  geom_point() +
  geom_smooth()

qqPlot(residuals_pH$residuals)



predicted_central_pH <- predict(model_pH, new = data.frame(concentration = 0,
                                                           duration = 0,
                                                           pH = 0,
                                                           temperature = 0), interval = "prediction")

experiment_central_pH <- glucose_coded %>%
  filter(concentration == 0 & duration == 0 & pH == 0 & temperature == 0) %>%
  summarize(mean(pH_final, na.rm = TRUE)) %>%
  as.numeric(.)


model_pH <- lm(data = glucose_coded, formula = model_pH$call)
Anova(model_pH, type = "II")
summary(model_pH)

library(DAAG)

glucose_pH_cv <- cv.lm(data = model_pH$model, form.lm = formula(model_pH), m = 5, dots = 0, seed = 123, plotit = FALSE, printit = TRUE)

ggplot(data = glucose_pH_cv, aes(x = pH_final, y = cvpred)) +
  geom_point() + 
  geom_abline(slope = 1, intercept = 0)








