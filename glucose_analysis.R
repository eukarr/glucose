library(tidyverse)
library(car)
setwd("C:/Users/Evgeny/Dropbox/R projects/CQD")

glucose <- read.csv("glucose_analysis.csv", fileEncoding="UTF-8-BOM", na.strings = "n/a")

str(glucose)

encode_var = function(x) {
  coded <- 2 * (x - mean(x)) / (max(x) - min(x))
}

glucose_coded <- glucose %>%
  mutate(concentration = encode_var(concentration),
         duration = encode_var(duration),
         pH = encode_var(pH),
         temperature = encode_var(temperature)) 

write.csv(x = glucose_coded, file = "glucose_analysis_coded.csv")


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
