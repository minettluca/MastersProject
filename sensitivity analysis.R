#prior sensitivity
library(brms); library(posterior); library(loo)

bf0 <- bf(logcontrolrate ~ 1 + me(logmass, log_mass_se) + (1|Common_name) + (1|Batch_id))

priors_list <- list(
  c(prior(normal(2/3,0.50), class="b", coef="melogmasslog_mass_se")),
  c(prior(normal(2/3,0.20), class="b", coef="melogmasslog_mass_se")),
  c(prior(normal(0.8,0.20),  class="b", coef="melogmasslog_mass_se")),
  c(prior(normal(0.5,0.20),  class="b", coef="melogmasslog_mass_se")),
  c(prior(student_t(3,2/3,0.25), class="b", coef="melogmasslog_mass_se"))
)

fits_prior <- lapply(priors_list, function(p) brm(bf0, data=data, family=gaussian(), prior=p, refresh=0))
slopes <- sapply(fits_prior, function(f) fixef(f)["melogmasslog_mass_se",])
slopes

#likelihood sensitivity
fit_gauss   <- brm(bf0, data=data, family=gaussian(), refresh=0)
fit_student <- brm(bf0, data=data, family=student(),  refresh=0)

bf_het <- bf(logcontrolrate ~ 1 + me(logmass, log_mass_se) + (1|Common_name) + (1|Batch_id),
             sigma ~ 1 + logmass)
fit_het <- brm(bf_het, data=data, family=gaussian(), refresh=0)

rbind(gauss = fixef(fit_gauss)["melogmasslog_mass_se",],
      student = fixef(fit_student)["melogmasslog_mass_se",],
      het = fixef(fit_het)["melogmasslog_mass_se",])

#random effect stucture sensitivity
fit_no_sp  <- brm(bf(logcontrolrate ~ 1 + me(logmass, log_mass_se) + (1|Batch_id)),
                  data=data, family=gaussian(), refresh=0)

fit_sp_slope <- brm(bf(logcontrolrate ~ 1 + me(logmass, log_mass_se) + 
                         (1 + me(logmass, log_mass_se)|Common_name) + (1|Batch_id)),
                    data=data, family=gaussian(), refresh=0)

rbind(no_species = fixef(fit_no_sp)["melogmasslog_mass_se",],
      baseline   = fixef(fit_gauss)["melogmasslog_mass_se",],
      sp_slope   = fixef(fit_sp_slope)["melogmasslog_mass_se",])

#measurment error senstivity 
scale_f <- c(0.5, 1, 1.5, 2)
fits_me <- lapply(scale_f, function(k){
  data_k <- transform(data, log_mass_se = k*log_mass_se)
  brm(bf0, data=data_k, family=gaussian(), refresh=0)
})
cbind(scale=scale_f,
      t(sapply(fits_me, function(f) fixef(f)["melogmasslog_mass_se", c("Estimate","l-95% CI","u-95% CI")])))



