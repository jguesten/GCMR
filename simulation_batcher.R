rm(list=ls())

source("get_env.R")
source("libraries.R")
source("run_simulation.R")
library("tidyverse")
require("rstan")
source("modelling/stan_models.R")

options(mc.cores = 2)

# pipeline of parameter values

pipe_rasch = list(
  nreps = 2, # number of simulated data sets
  nsubj = 100, # number of subjects
  nitem = c(50), # number of items
  mu_gamma= c(0.5,0.75,0.9),
  mu_beta = c(-2,0,2),
  delta_beta = c(0),#c(0,0.5,1),
  mu_theta = c(0),
  sigma_beta = c(2),#c(0.1,0.5,2.5),
  sigma_theta = c(1),#c(0.1,0.5,2.5),
  missing = c(0),
  sets = F) %>%
  cross_df()

pipe_missing = list(
  nreps = 10, # number of simulated data sets
  nsubj = 100, # number of subjects
  nitem = c(50), # number of items
  mu_gamma= c(0.5,0.75,0.9),
  mu_beta = c(0),
  delta_beta = c(0),#c(0,0.5,1),
  mu_theta = c(0),
  sigma_beta = c(2),#c(0.1,0.5,2.5),
  sigma_theta = c(1),#c(0.1,0.5,2.5),
  missing = c(0,0.25,0.5),
  sets = F) %>%
  cross_df()

pipe_sets = list(
  nreps = 2, # number of simulated data sets
  nsubj = 300, # number of subjects
  nitem = c(150), # number of items
  mu_gamma= c(0.5,0.75,0.9),
  mu_beta = c(0),
  delta_beta = c(0),#c(0,0.5,1),
  mu_theta = c(0),
  sigma_beta = c(0.1,0.5,2.5),#c(0.1,0.5,2.5),
  sigma_theta = c(1),#c(0.1,0.5,2.5),
  missing = c(0),
  sets = T) %>%
  cross_df()

pipe_mturk = list(
  nreps = 10, # number of simulated data sets
  nsubj = 1550, # number of subjects
  nitem = c(3000), # number of items
  mu_gamma= c(0.75),
  mu_beta = c(0),
  delta_beta = c(0),#c(0,0.5,1),
  mu_theta = c(0),
  sigma_beta = c(1),#c(0.1,0.5,2.5),
  sigma_theta = c(1),#c(0.1,0.5,2.5),
  missing = c(0.9),
  sets = F) %>%
  cross_df()


# concatenate pipeline


pipe = list(pipe_rasch,
            pipe_missing
            #pipe_sets
            )

# specify stan models to estimate

stan_mods=list('one_pl_bc2_stan','rasch','Pr')
names(stan_mods)=stan_mods
names(stan_mods)[names(stan_mods)=="one_pl_bc2_stan"]="bc"

mod=lapply(stan_mods,function(x){stan_model(model_code = eval(parse(text=x)))}) #compile models and save to vector


# run pipeline
lapply(pipe,function(x)run_simulation(cp=x,stan_mods,mod))
