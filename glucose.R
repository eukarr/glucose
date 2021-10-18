library(tidyverse)
library(broom)

setwd("C:/Users/Evgeny/Dropbox/_projects/CQD/glucose")

glucose <- read.csv("glucose_manual_data.csv", fileEncoding="UTF-8-BOM", na.strings = "n/a")
# this is a dataset of input and response variables collected manually
# id - sample label
# volume - reactor volume in mL
# initial mass - mass of the reaction mixture, g
# concentration - concentration of glucose in the reaction mixture, wt%
# duration - synthesis duration, h
# pH - initial pH of the reaction mixture
# temperature - syntehsis temperature
# pH_final - pH of the filtrate
# yield_total - total yield of the product, %, including sediment
# conc_filtrate - concentration of the filtrate, mg/g
# yield_filtrate - yield of the soluble product, %


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

colnames(citrate_spectra) = c("wavelength", "Abs", "spectrum", "sample_id")

# for each sample, a median spectrum of 2-3 is retained
citrate_median_spectra_by_sample <- citrate_spectra %>%
  group_by(sample_id, wavelength) %>%
  summarize(Abs = median(Abs))

# for each median spectrum, maximum Abs is saved, and the midpoint of the range corresponding to this Abs value
citrate_by_sample_abs <- citrate_median_spectra_by_sample %>%
  summarize(Abs_max = max(Abs), 
            wavelength_max = mean(wavelength[Abs == Abs_max]))

# adding spectral data to the main dataset
# Abs_max - the maximum Abs in baseline-corrected spectra
# wavelength_max - wavelength corresponding to the maximum Abs
glucose_full <- merge(glucose, citrate_by_sample_abs, by.x="id", by.y="sample_id", all.x = TRUE)      
write.csv(x = glucose_full, file = "glucose_full.csv")


glucose <- glucose_full[complete.cases(glucose_full), ]
write.csv(x = glucose, file = "glucose_for_analysis.csv")


# a function to encode factorial variables into the [-1, 1] range
encode_var = function(x) {
  coded <- 2 * (x - mean(x)) / (max(x) - min(x))
}

# a dataset of coded variables for ANOVA
glucose_coded <- glucose %>%
  mutate(concentration = encode_var(concentration),
         duration = encode_var(duration),
         pH = encode_var(pH),
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
  subset(select = c("id", "temperature", "duration", "pH", "concentration", "volume", "yield_total", "yield_filtrate", "pH_final", "conc_filtrate", "Abs_max"))


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
         conc_filtrate_diff = conc_filtrate.x - conc_filtrate.y,
         Abs_max_diff = Abs_max.x - Abs_max.y)
  

# long-format dataset of the variables scaled by SD for creating boxplot
v_eff_plot <- v_eff[, 19:23] %>%
  apply(MARGIN = 2, function(x) {x / sd(x)}) %>%
  as.data.frame() %>%
  gather(key = variable, value = value, factor_key=TRUE)

ggplot(data = v_eff_plot) +
  geom_hline(yintercept = 0, colour = "red") +
  geom_boxplot(aes(x=variable, y=value)) +
  ylab(label = "Z-score") +
  xlab(label = NULL) +
  scale_x_discrete(labels = c('Total yield',
                              'Filtrate yield',
                              'Final pH', 
                              'Filtrate\nconcentration', 
                              'Maximum\nabsorbance'))
  

summary(v_eff$yield_filtrate_diff)

t1 <- t.test(v_eff$yield_total.x, v_eff$yield_total.y)
t2 <- t.test(v_eff$yield_filtrate.x, v_eff$yield_filtrate.y)
t3 <- t.test(v_eff$pH_final.x, v_eff$pH_final.y)
t4 <- t.test(v_eff$conc_filtrate.x, v_eff$conc_filtrate.y)
t5 <- t.test(v_eff$Abs_max.x, v_eff$Abs_max.y)


v_eff_result <- map_df(list(t1, t2, t3, t4, t5), tidy) %>%
  subset(select = c("estimate", "statistic", "p.value", "conf.low", "conf.high")) %>%
  as.data.frame() %>%
  signif(3)

rownames(v_eff_result) <- c("Total yield",
                            "Filtrate yield",
                            "Final pH",
                            "Concentration in filtrate",
                            "Absorbance in maximum")

v_eff_result

#==============END OF EFFECT OF VOLUME=================================


#==============ANALYSIS OF YIELD=======================================
glucose_coded %>%
  filter(temperature == -1 & duration == -1) %>%
  subset(select = c(yield_total, yield_filtrate)) %>%
  apply(MARGIN = 2, mean)

glucose_coded %>%
  filter(temperature == 1 & duration == 1) %>%
  subset(select = c(yield_total, yield_filtrate)) %>%
  apply(MARGIN = 2, mean)


yield_diff <- glucose_coded %>%
  subset(select = -c(pH_final, conc_filtrate, Abs_max, wavelength_max, initial_mass)) %>%
  mutate(yield_diff = (yield_total - yield_filtrate) / yield_filtrate)

ggplot(data = yield_diff) +
  geom_histogram(aes(x = yield_diff), bins = 20)

yield_diff_limits <- quantile(yield_diff$yield_diff, c(0.15, 0.85))

yield_outliers <- yield_diff %>%
  filter(yield_diff < yield_diff_limits[1] | yield_diff > yield_diff_limits[2]) %>%
  arrange(abs(yield_diff))

ggplot(data = glucose_coded, aes(x = yield_total, y = yield_filtrate)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "blue", size = 1) +
  geom_hline(yintercept = 100, linetype = 'dashed') + 
  geom_vline(xintercept = 100, linetype = 'dashed') +
  geom_point(data = yield_outliers, aes(x = yield_total, y = yield_filtrate), color = "red") +
  xlab(label = "Total yield, %") +
  ylab(label = "Filtrate yield, %") +
  ylim(c(0, 100)) + xlim(c(0, 130))

model_yield_filtrate <- lm(data = glucose_coded, formula = yield_filtrate ~ concentration * duration * pH * temperature)
Anova(model_yield_filtrate, type = "II")
summary(model_yield_filtrate)


ggplot(data = glucose %>%
         filter(concentration != 7.5) %>%
         group_by(temperature, duration) %>%
         summarize(yield_filtrate = mean(yield_filtrate))) +
  geom_line(aes(x = temperature, y = yield_filtrate, color = factor(duration)), size = 1) +
  geom_point(aes(x = temperature, y = yield_filtrate, color = factor(duration)), size = 4) +
  xlab(label = "Temperature (°C)") +
  ylab(label = "Filtrate yield, %") +
  labs(color = "Duration (h)") +
  theme(legend.position = c(0.85, 0.85))
  

         
model_yield_filtrate_reduced <- lm(data = glucose_coded, 
                                   formula = yield_filtrate ~ concentration + duration + pH + temperature +
                                     duration:temperature)
Anova(model_yield_filtrate_reduced, type = "II")
summary(model_yield_filtrate_reduced)




model_yield_filtrate_reduced_no_central <- lm(data = glucose_coded %>%
                                             filter(concentration != 0), 
                                   formula = yield_filtrate ~ concentration + duration + pH + temperature +
                                     duration:temperature)

experiment_central_yield_filtrate <- glucose_coded %>%
  filter(concentration == 0) %>%
  summarize(mean(yield_filtrate, na.rm = TRUE)) %>%
  as.numeric(.)

predicted_central_yield_filtrate <- predict(model_yield_filtrate_reduced_no_central, 
                                            new = data.frame(concentration = 0,
                                                           duration = 0,
                                                           pH = 0,
                                                           temperature = 0), interval = "prediction")


#==============END OF ANALYSIS OF YIELD================================



#==============ANALYSIS OF pH==========================================
model_pH_final <- lm(data = glucose_coded%>%
                       filter(volume == 10), 
                     formula = pH_final ~ concentration * duration * pH * temperature)
Anova(model_pH_final, type = "II")
summary(model_pH_final)


model_pH_final <- lm(data = glucose_coded%>%
                       filter(volume == 10 & pH == 0), 
                     formula = pH_final ~ concentration * duration * temperature)
Anova(model_pH_final, type = "II")
summary(model_pH_final)


model_pH_final <- lm(data = glucose_coded%>%
                       filter(volume == 10 & pH == -1), 
                     formula = pH_final ~ concentration * duration * temperature)
Anova(model_pH_final, type = "II")
summary(model_pH_final)


model_pH_final <- lm(data = glucose_coded%>%
                       filter(volume == 10 & pH == 1), 
                     formula = pH_final ~ concentration * duration * temperature)
Anova(model_pH_final, type = "II")
summary(model_pH_final)
#==============END OF ANALYSIS OF pH===================================




#=================== OLD STUFF




model_yield_filtrate <- lm(data = glucose_coded %>%
                             filter(concentration != 0 & duration != 0 & pH != 0 & temperature != 0),
                           formula = yield_filtrate ~ concentration * duration * pH * temperature)

Anova(model_yield_filtrate, type = "II")
summary(model_yield_filtrate)




model_yield_total <- lm(data = glucose_coded, formula = yield_total ~ concentration * duration * pH * temperature)
Anova(model_yield_total, type = "II")
summary(model_yield_total)

model_yield_total <- lm(data = glucose_coded %>%
                          filter(concentration != 0 & duration != 0 & pH != 0 & temperature != 0),
                        formula = yield_total ~ concentration * duration * pH * temperature)

Anova(model_yield_total, type = "II")
summary(model_yield_total)




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



















ggplot() +
  geom_point(data=glucose_full, aes(x=conc_filtrate, y=Abs, color=factor(volume)))







temp <- merge(citrate_median_spectra_by_sample, glucose, by.x="sample_id", by.y="id", all.x = TRUE)
temp$Abs <- temp$Abs / (temp$initial_mass * temp$concentration)

citrate_median_spectra_by_composition <- temp %>%
  group_by(concentration, duration, pH, temperature, wavelength) %>%
  summarize(Abs = median(Abs))

citrate_median_spectra_by_composition$conditions <- paste0("pH_", citrate_median_spectra_by_composition$pH, " ",
                                                           "c_", citrate_median_spectra_by_composition$concentration, " ",
                                                           "t_", citrate_median_spectra_by_composition$duration, " ",
                                                           "T_", citrate_median_spectra_by_composition$temperature)

ggplot() +
  geom_line(data=citrate_median_spectra_by_composition, aes(x=wavelength, y=Abs, color=conditions))


ggplot() +
  geom_line(data=citrate_median_spectra_by_composition[citrate_median_spectra_by_composition$pH == 3, ], 
            aes(x=wavelength, y=Abs, color=conditions))


ggplot() +
  geom_line(data=citrate_median_spectra_by_composition[citrate_median_spectra_by_composition$pH == 6, ], 
            aes(x=wavelength, y=Abs, color=conditions))


citrate_samples_reduced_by_composition <- citrate_median_spectra_by_composition %>%
  group_by(conditions) %>%
  mutate(Abs = Abs / max(Abs))

ggplot() +
  geom_line(data=citrate_samples_reduced_by_composition, aes(x=wavelength, y=Abs, color=conditions))

ggplot() +
  geom_line(data=citrate_samples_reduced_by_composition[citrate_samples_reduced_by_composition$wavelength < 300, ],
            aes(x=wavelength, y=Abs, color=conditions))



ggplot() +
  geom_point(data=citrate_spectra, aes(x=wavelength, y=Abs, color=sample))

ggplot() +
  geom_line(data=citrate_median_spectra, aes(x=wavelength, y=Abs, color=sample))

citrate_by_sample_abs <- citrate_spectra_by_sample %>%
  summarize(Abs = max(Abs))

citrate_by_sample_abs








# glucose_full <- merge(glucose, citrate_by_sample_abs, by.x="id", by.y="sample", all.x = TRUE)      

# write.csv(x = glucose_full, file = "glucose_analysis_full.csv")

ggplot() +
  geom_point(data=glucose_full, aes(x=conc_filtrate, y=Abs, color=factor(volume)))


library(GGally)
ggpairs(glucose_full, columns = c(8:11, 14) , ggplot2::aes(colour=factor(volume)))

str(glucose_full)


# library(tidyverse)
library(car)
setwd("C:/Users/Evgeny/Dropbox/R projects/CQD")

glucose <- read.csv("glucose_analysis.csv", fileEncoding="UTF-8-BOM", na.strings = "n/a")

str(glucose[glucose$volume == 5, ])

small_reactor <- glucose[glucose$volume == 5, ]

str(glucose)

#encode_var = function(x) {
#  coded <- 2 * (x - mean(x)) / (max(x) - min(x))
#}

#glucose_coded <- glucose %>%
#  mutate(concentration = encode_var(concentration),
#         duration = encode_var(duration),
#         pH = encode_var(pH),
#         temperature = encode_var(temperature)) 

#write.csv(x = glucose_coded, file = "glucose_analysis_coded.csv")


glucose_pH <- glucose_coded %>%
  select(concentration, duration, pH, temperature, pH_final) %>%
  filter(concentration != 0 & duration != 0 & pH != 0 & temperature != 0)


model_pH <- lm(data = glucose_pH, formula = pH_final ~ concentration * duration * pH * temperature)
Anova(model_pH, type = "II")
summary(model_pH)


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








model_total_yield <- lm(data = glucose, formula = yield_total ~ concentration * duration * pH * temperature)
Anova(model_total_yield, type = "II")
summary(model_total_yield)

model_filtrate_yield <- lm(data = glucose, formula = yield_filtrate ~ concentration * duration * pH * temperature)
Anova(model_filtrate_yield, type = "II")
summary(model_filtrate_yield)

ggplot(model_pH) + 
  geom_point(aes(x = model_pH$model$pH_final, y = model_pH$residuals))

glucose %>% group_by(pH) %>% summarize_all(mean, na.rm = TRUE)


model_pH_main <- lm(data = glucose, formula = pH_final ~ concentration + duration + pH + temperature)
Anova(model_pH_main, type = "II")
summary(model_pH_main)

model_total_yield_main <- lm(data = glucose, formula = yield_total ~ concentration + duration + pH + temperature)
Anova(model_total_yield_main, type = "II")
summary(model_total_yield_main)

model_filtrate_yield_main <- lm(data = glucose, formula = yield_filtrate ~ concentration + duration + pH + temperature)
Anova(model_filtrate_yield_main, type = "II")
summary(model_filtrate_yield_main)

glucose_average <- glucose %>% group_by(concentration, duration, pH, temperature) %>%
  summarise(mean(pH_final), mean(yield_total), mean(yield_filtrate))

model_pH_average <- lm(data = glucose_average, formula = `mean(pH_final)` ~ concentration * duration * pH * temperature)
Anova(model_pH_average, type = "II")
summary(model_pH_average)

model_total_yield_average <- lm(data = glucose_average, formula = `mean(yield_total)` ~ concentration * duration * pH * temperature)
Anova(model_total_yield_average, type = "II")
summary(model_total_yield_average)

model_filtrate_yield_average <- lm(data = glucose_average, formula = `mean(yield_filtrate)` ~ concentration * duration * pH * temperature)
Anova(model_filtrate_yield_average, type = "II")
summary(model_filtrate_yield_average)


model_pH_average_main <- lm(data = glucose_average, formula = `mean(pH_final)` ~ concentration + duration + pH + temperature)
Anova(model_pH_average_main, type = "II")
summary(model_pH_average_main)

model_total_yield_average_main <- lm(data = glucose_average, formula = `mean(yield_total)` ~ concentration * duration * pH * temperature)
Anova(model_total_yield_average_main, type = "II")
summary(model_total_yield_average_main)

model_filtrate_yield_average_main <- lm(data = glucose_average, formula = `mean(yield_filtrate)` ~ concentration * duration * pH * temperature)
Anova(model_filtrate_yield_average_main, type = "II")
summary(model_filtrate_yield_average_main)


model_pH_2nd <- lm(data = glucose, formula = pH_final ~ (concentration + duration + pH + temperature)^2)
summary(model_pH_2nd)

model_filtrate_yield_2nd <- lm(data = glucose, formula = yield_filtrate ~ (concentration + duration + pH + temperature)^2)
summary(model_filtrate_yield_2nd)

model_filtrate_yield_final <- lm(data = glucose, formula = yield_filtrate ~ concentration + duration + temperature+ duration:temperature + concentration:duration + concentration:temperature)
summary(model_filtrate_yield_final)

glucose_coded <- glucose %>% mutate(pH = (pH - min(pH))/(max(pH) - min(pH)),
                                    temperature = (temperature - min(temperature))/(max(temperature) - min(temperature)),
                                    concentration = (concentration - min(concentration))/(max(concentration) - min(concentration)),
                                    duration = (duration - min(duration))/(max(duration) - min(duration)))

model_filtrate_yield_coded_final <- lm(data = glucose_coded, formula = yield_filtrate ~ concentration + duration + temperature+ duration:temperature + concentration:duration + concentration:temperature)
summary(model_filtrate_yield_coded_final)


model_pH_coded_final <- lm(data = glucose_coded, formula = pH_final ~ concentration + duration + temperature+ duration:temperature + concentration:duration + concentration:temperature)
summary(model_pH_coded_final)

model_total_yield_coded_final <- lm(data = glucose_coded, formula = yield_total ~ concentration + duration + temperature+ duration:temperature + concentration:duration + concentration:temperature)
summary(model_total_yield_coded_final)


model_filtrate_yield_test <- lm(data = glucose_coded, formula = yield_total ~ (concentration + duration + temperature)^2)
summary(model_filtrate_yield_test)









print(object.size(dataset), standard = "legacy", units = "Mb")


