library(dplyr)
library(brms)
library(ggplot2)

#linear model
lm <- lm(data=data, controlrate/Start_invert_mass ~ Common_name)
summary(lm)
plot(lm)

library(emmeans)
emm <- emmeans(lm(controlrate/Start_invert_mass ~ Common_name, data = data), ~ Common_name)
pairs(emm, adjust = "tukey")
species_list <- split(data, data$Common_name)

#bayesian model

#formula
formularough <- bf(logcontrolrate ~ 1 + me(logmass, log_mass_se) + (1 + me(logmass, log_mass_se) | Batch_id))

#per species bayesian models
fits <- lapply(names(species_list), function(sp) {
  brm(
    formula,
    data = species_list[[sp]],
    family = gaussian(),
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.999, max_treedepth = 15)
  )
})
names(fits) <- names(species_list)

#extract slopes with uncertainty
slopes <- lapply(fits, function(fit) {
  post <- as_draws_df(fit)
  tibble(
    mean = mean(post$b_me_logmasslog_mass_se),
    lower = quantile(post$b_me_logmasslog_mass_se, 0.025),
    upper = quantile(post$b_me_logmasslog_mass_se, 0.975)
  )
})

slopes_df <- bind_rows(slopes, .id = "Species")

#summaries
summary(fits[["Painted Woodlouse"]])
summary(fits[["Pill Woodlouse"]])
summary(fits[["Rough Woodlouse"]])
summary(fits[["Striped Woodlouse"]])


#bayesian model for just rough woodlice
fit_rough <- brm(formularough,
           data = subset(data, Common_name == "Rough Woodlouse"),
           family = gaussian(),
           prior = c(
             prior(normal(0, 5), class = "b"),
             prior(student_t(3, 0, 2.5), class = "sd"),
             prior(student_t(3, 0, 2.5), class = "sigma")
           ),
           chains = 4, iter = 4000, control = list(adapt_delta = 0.9999, max_treedepth = 15))

summary(fit_rough)


#hypothesis testing against 0.75

variables(fit_rough)
# Extract posterior draws
posterior <- as_draws_df(fit_rough)

# Probability slope > 0.75
mean(posterior$bsp_melogmasslog_mass_se > 0.75)

# Probability slope < 0.75
mean(posterior$bsp_melogmasslog_mass_se < 0.75)

# Posterior mean difference from 0.75
mean(posterior$bsp_melogmasslog_mass_se - 0.75)

# 95% credible interval for difference from 0.75
quantile(posterior$bsp_melogmasslog_mass_se - 0.75, probs = c(0.025, 0.975))


#hypothesis testing against 0.47


# Probability slope > 0.47
mean(posterior$bsp_melogmasslog_mass_se > 0.47)

# Probability slope < 0.47
mean(posterior$bsp_melogmasslog_mass_se < 0.47)

# Posterior mean difference from 0.47
mean(posterior$bsp_melogmasslog_mass_se - 0.47)

# 95% credible interval for difference from 0.47
quantile(posterior$bsp_melogmasslog_mass_se - 0.47, probs = c(0.025, 0.975))


#plotting


# Posterior predictions for plotting
newdata <- data.frame(
  logmass = seq(min(data$logmass[data$Common_name == "Rough Woodlouse"]),
                max(data$logmass[data$Common_name == "Rough Woodlouse"]),
                length.out = 100),
  log_mass_se = 0,               # set measurement error to zero for prediction
  Batch_id = NA                  # no grouping for marginal predictions
)

se_val <- min(data$log_mass_se[data$Common_name == "Rough Woodlouse" & data$log_mass_se > 0])

newdata$log_mass_se <- se_val

pp <- posterior_epred(fit_rough, newdata = newdata, re_formula = NA)
newdata$Estimate <- apply(pp, 2, mean)
newdata$Q2.5     <- apply(pp, 2, quantile, 0.025)
newdata$Q97.5    <- apply(pp, 2, quantile, 0.975)

# Observed data subset
roughdat <- subset(data, Common_name == "Rough Woodlouse")

# Plot

ggplot(roughdat, aes(x = logmass, y = logcontrolrate, color = Batch_id)) +
  geom_point(alpha = 0.6) +
  geom_line(data = newdata,
            aes(x = logmass, y = Estimate),
            color = "blue",
            inherit.aes = FALSE) +
  geom_ribbon(data = newdata,
              aes(x = logmass, ymin = Q2.5, ymax = Q97.5),
              alpha = 0.2, fill = "blue",
              inherit.aes = FALSE) +
  theme_minimal() +
  labs(
       x = log[10]~'Body Mass',
       y = log[10]~'Consumption Rate',
       color = 'Batch number')




