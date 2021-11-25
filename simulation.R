
run_simulation=function(cp,stan_mods,mod){

library("tidyverse")
library("sn")
library("rstan")
library("gridExtra")
library("doParallel")
library("qs")
source("stan_models.R")

# ensure this script returns the same results on each run
set.seed(8675309)

## specify alpha and beta pars for gamma so that variance remains constant
gamma_pars=list(
  '0.5' = c(10,10),
  '0.75' = c(11.0625,3.6875),
  '0.9' = c(5.904,0.656)
)

# FUNCTIONS

inv_logit=function(x)(exp(x)/(1+exp(x)))
logit=function(x)x/(1-x)

get_rmse=function(obs,fit){sqrt(mean((obs-fit)^2))} # get root of mean squared error

make_stanlist = function(stanlist){
  source('data.R')

  stanlist$ii = get_ii(stanlist$item_id)
  stanlist$jj = get_jj(stanlist$sub_id)
  stanlist$N = length(stanlist$y)
  stanlist$I = length(unique(stanlist$ii))
  stanlist$J = length(unique(stanlist$jj))
  stanlist
}


# set up the custom data simulation function
my_sim_data <- function(nsubj,nitem,mu_theta,mu_beta,sigma_beta,delta_beta,sigma_theta,alpha_gamma,beta_gamma,missing){ # residual (standard deviation)

    # simulate items
    # simulate a sample of items
    items <- data.frame(
      item_id = 1:sum(nitem),
      trial_type = rep(c("l", "r"), nitem),
      beta_m = c(rnorm(nitem[["l"]], mu_beta+delta_beta, sigma_beta),rnorm(nitem[["r"]], mu_beta-delta_beta,sigma_beta))
    )
    # effect code category
    items$z <- dplyr::recode(items$trial_type, "l" = 0, "r" = 1)

    # simulate subjects
    subjects <- data.frame(
      sub_id = 1:nsubj,
      theta = rnorm(nsubj, mu_theta, sigma_theta),
      gamma = rbeta(n=nsubj,alpha_gamma, beta_gamma, ncp = 0)
  )
  # simulate trials
  dat_sim <- crossing(sub_id = subjects$sub_id,
  item_id = items$item_id) %>%
  inner_join(subjects, "sub_id") %>%
  inner_join(items, "item_id") %>%
  #mutate(err = rnorm(nrow(.), mean = 0, sd = err_sd)) %>%
  mutate(prob_y = inv_logit(theta-beta_m)+(1-inv_logit(theta-beta_m))*abs(z-gamma)) %>%
  mutate(y = rbinom(nrow(.),1,prob_y))
  #mutate(y = rbinom(nrow(.),1,inv_logit(theta-beta_m)))
  #%>%
  #select(sub_id, item_id, category, cat, RT)
  dat_sim
}

run_mod=function(simdat,mod){
  # make stanlist
  stanlist = make_stanlist(simdat)

  # run model
  fit = sampling(mod, data = stanlist, chains = 2, iter = 1000)

  # extract parameters

  pars_sim=rstan::extract(fit)

  pars_sim
}


## RUN SIMULATION AND MODEL
# define and compile models to run

simdat=list()
fit_pars=list()

# Try in parallel

  cl <- makeForkCluster(nnodes = getOption("mc.cores", 2L)) ### use half of the available cores
  registerDoParallel(cl)
    simdat=foreach(i=1:nrow(cp)) %do% {
      alpha_gamma = gamma_pars[[as.character(cp$mu_gamma)[i]]][1]
      beta_gamma = gamma_pars[[as.character(cp$mu_gamma)[i]]][2]
      sigma_beta = cp[i,]$sigma_beta
      sigma_theta = cp[i,]$sigma_theta
      mu_theta = cp[i,]$mu_theta
      mu_beta = cp[i,]$mu_beta
      delta_beta = cp[i,]$delta_beta
      nreps = cp[i,]$nreps
      nsubj = cp[i,]$nsubj
      nitem = c(l=cp[i,]$nitem,r=cp[i,]$nitem)
      missing = cp[i,]$missing
      sets = cp[i,]$sets


      temp=replicate(nreps,my_sim_data(nsubj = nsubj, nitem = nitem, mu_theta = mu_theta, mu_beta = mu_beta, sigma_beta = sigma_beta, delta_beta = delta_beta, sigma_theta = sigma_theta, alpha_gamma = alpha_gamma, beta_gamma = beta_gamma))

    #}
    if(!missing==0){
      tempdat=apply(temp,2,function(x)
      {tempdat=as.data.frame(x)%>%mutate(rowid=1:nrow(.))
      missdat=tempdat%>%group_by(sub_id,trial_type)%>%sample_n(missing*(sum(nitem)/2))%>%mutate(y=NA)
      tempdat[missdat$rowid,]=missdat
      tempdat=dplyr::select(tempdat, -rowid)
      as.list(tempdat)
      })
      temp=do.call(cbind,tempdat)
    }

    if(sets==T){
      tempdat=apply(temp,2,function(x)
      {tempdat=as.data.frame(x)%>%
        arrange(sub_id,beta_m)%>%
        group_by(sub_id)%>%
        mutate(set=c(rep('a',length(sub_id)/4),rep('b',length(sub_id)/4),rep('c',length(sub_id)/4),rep('d',length(sub_id)/4)))%>%ungroup()
      tempdat=tempdat%>%mutate(group=c(rep('easy',nrow(.)/3),rep('medium',nrow(.)/3),rep('hard',nrow(.)/3)))
      tempdat$y[tempdat$group=='easy'&(tempdat$set=="c"|tempdat$set=="d")]=NA
      tempdat$y[tempdat$group=='medium'&(tempdat$set=="a"|tempdat$set=="d")]=NA
      tempdat$y[tempdat$group=='hard'&(tempdat$set=="a"|tempdat$set=="b")]=NA
      as.list(tempdat)
      })
      temp=do.call(cbind,tempdat)
    }
    simdat=temp

    pars.fit=foreach(s=names(stan_mods), .final = function(s) setNames(s, names(stan_mods))) %do% {
      # set up stan model
          apply(simdat,2,function(x){
          x=as.data.frame(x)
          x=x[complete.cases(x),]
          x=as.list(x)
          run_mod(x,mod[[s]])
          }) # run model for each dataset column and extract pars
        ## add extracted pars to corresponding grid parameter list
    }
    list(pars.fit,simdat)
    }
    stopCluster(cl)

fit_pars=lapply(simdat,function(x)x[[1]])
simdat=lapply(simdat,function(x)x[[2]])

## save model and data
cp_names=paste0(names(cp)[(sapply(cp,function(x)length(unique(x))>1))],collapse="_")
model_names=paste0(names(stan_mods),collapse="_")

cp_file=paste0(model_dir,'/','cp_',model_names,'_',nreps,'_reps_',cp_names,'_cp.qs')
fit_pars_file = paste0(model_dir,'/','fit_pars_',model_names,'_',nreps,'_reps_',cp_names,'_cp.qs')
simdat_file = paste0(data_dir,'/','simdat_',model_names,'_',nreps,'_reps_',cp_names,'_cp.qs')

qsave(fit_pars,fit_pars_file)
qsave(simdat,simdat_file)
qsave(cp,cp_file)

cp = qread(cp_file)
fit_pars = qread(fit_pars_file)
simdat = qread(simdat_file)

#combine simulated and estimated pars


result_tab =list()
#Structure:
  # cp level
      # model
          # reps

for(i in 1:nrow(cp)){
    theta_df = sapply(names(stan_mods),function(x){sapply(fit_pars[[i]][[x]],function(y)apply(y$theta,2,mean))},simplify = F,USE.NAMES = T)
    gamma_df = list(bc = sapply(fit_pars[[i]]$bc,function(x)apply(x$gamma,2,mean)),
    gamma_sim = apply(simdat[[i]],2,function(x)unique(x$gamma))) #%>% as.data.frame()

    result_tab[[i]] = list(theta_df,#beta_df,
                       gamma_df)#,y_df)
}


# add parameter value column from cp
for(i in 1:length(result_tab)){
  for(m in 1:length(result_tab[[i]])){
    sim_ind=sapply(sapply(c(1:ncol(result_tab[[i]][[m]]$bc)),function(x)rep(x,nrow(result_tab[[i]][[m]]$bc))),print)
    result_tab[[i]][[m]]=as.data.frame(sapply(result_tab[[i]][[m]],function(x)(sapply(x,print,USE.NAMES = TRUE))))
    lapply(names(cp),function(x)result_tab[[i]][[m]][[x]]<<-cp[[x]][i])
    result_tab[[i]][[m]]=cbind(result_tab[[i]][[m]],sim_ind)
  }
 }
}


