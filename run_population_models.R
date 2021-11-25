# import packages

require(rstan)
require(tidyverse)
require(qs)
source("stan_models.R")

# set directories here
data_dir = 
model_dir = 

# import mturk data
turk_trialdata <-
  readRDS(paste0(data_dir,"trialdata_prepared.rds"))

# code Domain
trialdata = rename(turk_trialdata, "Domain" = "Stim_Type(object=1, scene=2)")
trialdata$Domain = factor(trialdata$Domain, labels = (c("Object", "Scene")))

# add PR measure
trialdata = turk_trialdata   %>%
  as.data.frame()   %>%
  group_by(ID, Trial_Type)   %>%
  dplyr::summarize(pc = mean(result, na.rm = T))   %>%
  spread(Trial_Type, pc)   %>%
  dplyr::summarize(pr = `0` - (1 - `1`))   %>%
  merge(trialdata)

saveRDS(trialdata, paste0(data_dir, "/mturk_trialdata.rds"))

# convert data into stanformat and adjust variables

get_stan_dat = function(data){
  
# prepare format for stan
temp = cbind(ID = unique(data$ID), jj = 1:length(unique(data$ID)))
data = merge(data, temp)
temp = cbind(stim_id = unique(data$stim_id),
             ii = 1:length(unique(data$stim_id)))
data = merge(data, temp) %>% rename(c("y" = "result", "z" = "Trial_Type","domain" = "Domain","age" = "Age"))

temp=tibble(jj=data$jj,age=data$age)%>%unique
unique_jj=temp$jj
unique_age=temp$age

data$domain=data$domain%>%as.numeric()-1.5 #Object:-0.5,Scene:0.5
temp=tibble(ii=data$ii,domain=data$domain)%>%unique
unique_ii=temp$ii
unique_domain=temp$domain


d = list(
  N = nrow(data),
  J = length(unique(data$jj)),
  I = length(unique(data$ii)),
  ID = data$ID,
  domain=data$domain,
  age = scale(data$age,scale=FALSE)%>%.[,1],
  stim_id = data$stim_id,
  jj = data$jj,
  ii = data$ii,
  unique_jj=unique_jj,
  unique_age=unique_age,
  unique_ii=unique_ii,
  unique_domain=unique_domain,
  z = data$z,
  y = data$y
)
}

## run model
options(mc.cores = parallel::detectCores())

run_mod = function(data, model) {
  
  # generate and run model
  mod <- stan_model(model_code = eval(parse(text = model)))
  
  fit = rstan::sampling(mod,
                        data = data,
                        chains = 2,
                        iter = 2000)
  fit
}

## domain - wise model

model="GCMR" ## Set model here

domain_mod=list()
for (i in c("Object","Scene")) {
  tempdat = trialdata %>%
    subset(Domain == i) %>%
    dplyr::select(ID, Trial_Type, stim_id, result) %>%
    get_stan_dat()
  domain_mod[[i]] = run_mod(tempdat, model)
}

## full model
tempdat = trialdata %>% 
  dplyr::select(ID, Trial_Type, stim_id, result, Age, Domain) %>%
  get_stan_dat()
full_mod = run_mod(tempdat, model)

## Save models
for(i in domain_mod){
qsave(domain_mod[[i]],file=paste0(model_dir,"fit_mturk_",model,"_model_", i, "_.qs"))
}
qsave(full_mod,file=paste0(model_dir,"fit_mturk_full_",model,"_model.qs"))
