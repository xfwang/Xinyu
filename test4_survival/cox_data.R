rm(list=ls())
library(dplyr)
library(survival)
library(gbm)
library(Rcpp)
library(rstan)

##########################################################################
# ----------------------   DATA LUNG (survial)   ------------------------
# inst: Institution code
# time: Survival time in days
# status: censoring status 1=censored, 2=dead
# age: Age in years
# sex: Male=1 Female=2
# ph.ecog: ECOG performance score (0=good 5=dead)
# ph.karno: Karnofsky performance score (bad=0-good=100) rated by physician
# pat.karno: Karnofsky performance score as rated by patient
# meal.cal: Calories consumed at meals
# wt.loss: Weight loss in last six months
##########################################################################
data(lung) # 228 10


## PACKAGE gbm
lung.mod <- gbm(Surv(time, status) ~ ., distribution="coxph", data=lung,
                n.trees=2000, shrinkage=0.01, cv.folds=5,verbose =FALSE)
relative.influence(lung.mod)

## PACKAGE survival
cox.fit <- coxph(Surv(time, status) ~ ., data = lung)
summary(cox.fit)$coefficients

## stan
lung_dat <- lung %>%
  subset(select = c("age", "ph.karno", "pat.karno", "meal.cal", "wt.loss")) %>%
  lapply(scale) %>%
  as.data.frame() %>%
  dplyr::mutate(time = lung[["time"]]) %>%
  dplyr::mutate(status = lung$status) %>%
  na.omit()
stan_file <- "test4_survival/coxph.stan"
idx.lung <- lung_dat$status == 2
stan_lung <- list(Nobs = sum(idx.lung),
                  Ncen = sum(!idx.lung),
                  yobs = lung_dat[idx.lung, "time"],
                  ycen = lung_dat[!idx.lung, "time"],
                  K = 5,
                  xobs = lung_dat[idx.lung, c("age", "ph.karno", "pat.karno", "meal.cal", "wt.loss")],
                  xcen = lung_dat[!idx.lung, c("age", "ph.karno", "pat.karno", "meal.cal", "wt.loss")]
                  )
stan.fit <- stan(stan_file,
                 data = stan_lung,
                 chains = 1,
                 pars = c("beta","alpha","lambda"), 
                 iter = 1000,
                 seed = 1328
                 )
