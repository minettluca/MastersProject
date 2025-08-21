#bayesian model

library(openxlsx)
library(dplyr)
library(tidyr)
library(segmented)
library(quantreg)
library(smatr)
library(brms)
library(ggplot2)
library(posterior)
library(patchwork)

setwd('C:/Users/chees/OneDrive/UNI one drive/Masters/final project/Data')

data <- read.xlsx('Woodlice_Consumption.xlsx', sheet = 'Data (2)')

#first, turn row 9 into character names
new_names <- as.character(unlist(data[9, ]))

#apply the new names to the data
names(data) <- new_names


#delete the first 9 rows (including the row that became headers)
data <- data[-(1:9), ]

#reset row names so they don’t look dumb
rownames(data) <- NULL

#filter only woodlice who survived
data <- filter(data, Alive == TRUE)

#convert mass to numeric to take log
data$Start_invert_mass<- as.numeric(data$Start_invert_mass)
data$Mass_difference <- as.numeric(data$Mass_difference)

#calculate daily rate
data$dailyrate <- data$Mass_difference/7


#add log mass and rate column
data <- mutate(data, logmass=log10(Start_invert_mass), lograte = log10(dailyrate))

#descriptive stats
mean(data$Start_invert_mass)
sd(data$Start_invert_mass)
min(data$Start_invert_mass)*1000
max(data$Start_invert_mass)*1000

mean(data$dailyrate)
sd(data$dailyrate)*1000
min(data$dailyrate)*1000
max(data$dailyrate)*1000

data$massmg <- data$Start_invert_mass*1000
data$ratemg <- data$controlrate*1000

#plots for descriptive stats

#hist for body mass
q1 <- ggplot(data, aes(x = massmg)) +
  geom_histogram(
    bins = 15,
    color = "black",
    fill = "steelblue",
    alpha = 0.7
  ) +
  labs(
    x = "Starting Invertebrate Mass (mg)",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)

#hist for rate
q2 <- ggplot(data, aes(x = ratemg)) +
  geom_histogram(
    bins = 15,
    color = "black",
    fill = "steelblue",
    alpha = 0.7
  ) +
  labs(
    x = expression(Rate~of~consumption~(mg~day^-1)),
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)

combined <- q1 / q2 + 
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.tag.position = c(0.9, 0.9),  # x, y from bottom-left (0,0) to top-right (1,1)
      plot.tag = element_text(size = 16, face = "bold")
    )
  )
combined

#linear model
model <- lm(lograte~logmass, data=data)

summary(model)

plot(x=data$Start_invert_mass, y= data$Leaf_mass_difference)
plot( x=data$logmass, y=data$controlrate)
abline(model)

lm_model <- lm(lograte ~ logmass, data = data)
summary(lm_model)

abline(lm_model)
plot(data$logmass, data$lograte, pch = 1, xlab = "logmass", ylab = "lograte")
plot(data$logmass, data$logcontrolrate)

#find variance in mass
datamass <- filter(data, !is.na(End_invert_mass))
datamass$End_invert_mass <- as.numeric(datamass$End_invert_mass)
datamass$logendmass <- log10(datamass$End_invert_mass)
datamass$difference <- datamass$logmass - datamass$logendmass

error <- var(datamass$difference)/2
error <- as.numeric(error)
error

totalvar <- var(data$logmass)
totalvar

errorestimate <- error/totalvar
errorestimate

#control
background <- 0.0005035
data$controlrate <- data$dailyrate - background
data$controlrate <- as.numeric(data$controlrate)
data$logcontrolrate <- log10(data$controlrate)
data <- data[data$controlrate >= 0, ]

#Bayesian model

data$log_mass_se <- sqrt(0.02386093)  # constant measurement error SD for all

# Model formula:
# me() for measurement error in predictor
# (1 + me(log_mass_obs, log_mass_se) | group) for random intercept and slope by group

formula <- bf(
  logcontrolrate ~ 1 + me(logmass, log_mass_se) + (1 | Common_name) + (1 | Batch_id)
)

formula <- bf(logcontrolrate ~ 1 + me(logmass, log_mass_se) + (1 | Common_name) + (1 + me(logmass, log_mass_se) | Batch_id))

formula2 <- bf(
  logcontrolrate ~ 1 + me(logmass, log_mass_se) + 
    (1 + me(logmass, log_mass_se) | Common_name) + 
    (1 | Batch_id)
)

# Fit Bayesian errors-in-variables mixed model
fit <- brm(formula,
           data = data,
           family = gaussian(),
           prior = c(
             prior(normal(0, 5), class = "b"),
             prior(student_t(3, 0, 2.5), class = "sd"),
             prior(student_t(3, 0, 2.5), class = "sigma")
           ),
           chains = 4, iter = 4000, control = list(adapt_delta = 0.9999, max_treedepth = 15))

fit2 <- brm(
  formula,
  data = data,
  family = gaussian(),
  prior = c(
    prior(normal(2/3, 0.2), class = "b", coef = "melogmasslog_mass_se"),
    prior(normal(-2.2, 0.5), class = "Intercept"),
    prior(normal(0, 0.25), class = "sd", group = "Batch_id", lb = 0),
    prior(normal(0, 0.25), class = "sd", group = "Common_name", lb = 0),
    prior(normal(0, 0.25), class = "sigma", lb = 0)
  ),
  chains = 4,
  iter = 4000,
  control = list(adapt_delta = 0.9999, max_treedepth = 15)
)



#compare to 0.75 scaling law

post <- as_draws_df(fit)
slope <- post$bsp_melogmasslog_mass_se 

mean(slope > 0.75)
quantile(slope - 0.75, c(0.025, 0.975))
mean(abs(slope - 0.75) < 0.05)

mean(slope > 0.55)
quantile(slope - 0.55, c(0.025, 0.975))
mean(abs(slope - 0.55) < 0.05)

#plot


# Generate predictions
newdata <- data.frame(
  logmass = seq(min(data$logmass), max(data$logmass), length.out = 100),
  log_mass_se = mean(data$log_mass_se),  # Use average measurement error
  Batch_id = data$Batch_id[1]  # Use first batch as reference
)

# Get predicted values with uncertainty
preds <- fitted(fit, newdata = newdata, allow_new_levels = TRUE)

# Combine for plotting
plot_data <- cbind(newdata, preds)

# Create the plot
ggplot() +
  geom_ribbon(data = plot_data, 
              aes(x = logmass, ymin = Q2.5, ymax = Q97.5), 
              alpha = 0.3, fill = "blue") +
  geom_line(data = plot_data, 
            aes(x = logmass, y = Estimate), 
            color = "blue", size = 1) +
  geom_point(data = data, 
             aes(x = logmass, y = logcontrolrate, color = Batch_id), 
             alpha = 0.6) +
  scale_color_discrete(name = "Batch number") +
  labs(x = expression(log[10]~"Body mass"), 
       y = expression(log[10]~"Consumption rate")) +
  theme_minimal()

#check outliers
pp_check(fit, type = "scatter_avg")
pp_check(fit, type = "intervals")

install.packages('broom.mixed')
library(broom.mixed)
res <- augment(fit, type.residuals = "pearson")
plot(res$.resid, res$log_mass)  # residuals vs predictor
hist(res$.resid)                # residual distribution
