### IMPORT LIBRARIES ###

require(tidyverse)
require(qs)
require(lme4)
require(reghelper)
require(visreg)
require(lmerTest)
source("model.R")

## set directories
data_dir = 
model_dir = 

### Load Subject Data ###
sub_df=readRDS(paste0(data_dir,"subject_results.rds"))%>%
  rename("Domain"="domain")

### set variables ###
domain_list=c("Object","Scene")
mod_list=c("2HT","GCMR","Rasch")

### set LME Predictor tibble ###
LME_df=sub_df%>%dplyr::select(
  ID,
  Age,
  Gender,
  Handedness,
  Domain,
  Browser,
  Platform,
  Condition,
  Random_Seq,
  version_order) %>%
    mutate(ID=as.factor(ID),
          Age=as.numeric(Age),
          Condition=as.factor(Condition),
          Handedness=as.factor(Handedness),
          Gender=as.factor(Gender),
          Domain=as.factor(Domain),
          Browser=as.factor(Browser),
          Platform=as.factor(Platform),
          Condition=as.factor(Condition),
          Random_Seq=as.factor(Random_Seq),
          version_order=as.factor(version_order)
          )

### DOMAIN-WISE DATA FOR ALL MODS ###

domain_dat = sapply(c("Object",
  "Scene"), function(x)
    qread(paste0(
      data_dir, "/data_mturk_stanformat_", x, ".qs"
    ))%>%as_tibble()%>%mutate(Domain=x), simplify = F, USE.NAMES = TRUE)

### DEFINE HELPER FUNCTIONS ###
domod_wrapper = function(domain_list, mod_list, func, invisible = T) {
  out = sapply(domain_list, function(dom){
    sapply(mod_list, function(mod) {
      func(dom, mod)
    },simplify = F,USE.NAMES = T)
  },simplify = F,USE.NAMES = T)
  out
}
extract_par=function(domain,mod){
  # extract pars
  fit=qread(paste0(model_dir, "/fit_mturk_",mod,"_model_", domain, "_.qs"))
  
  if(!(mod=="Rasch")){
    pars = list(
      tibble(score = rstan::extract(fit, pars = "theta")$theta %>% apply(2, mean))%>%mutate(jj=1:nrow(.),Model=mod,Domain=domain,Par="theta"),
      tibble(score = rstan::extract(fit, pars = "gamma")$gamma %>% apply(2, mean))%>%mutate(jj=1:nrow(.),Model=mod,Domain=domain,Par="gamma")
    )
  }else{
    pars = list(
      tibble(score = rstan::extract(fit, pars = "theta")$theta %>% apply(2, mean))%>%mutate(jj=1:nrow(.),Model=mod,Domain=domain,Par="theta")
    )
  }
}

### EXTRACT PARAMETERS ###
pars_list=domod_wrapper(domain_list,mod_list,extract_par)

par_df=pars_list%>%
  lapply(.,function(x)do.call(rbind,x))%>%lapply(.,function(x)do.call(rbind,x))%>%do.call(rbind,.)

### PREPARE DATA ###
full_dat=do.call(rbind,domain_dat)%>%select(ID,jj)%>%unique()%>%merge(par_df)

final_dat=merge(full_dat,sub_df)
final_dat=final_dat[!duplicated(final_dat),] ## remove duplicated rasch scores
final_dat=final_dat%>%group_by(Par,Model)%>%mutate(scaled_score=as.vector(scale(score)))

### DEFINE VAR GRID ###
crosslist=list(
  mod_list=mod_list,
  var_list=c("theta","gamma")
)%>%cross()

## remove nonexistent rasch gamma element ##
for(i in 1:length(crosslist)){
  if(crosslist[[i]]$mod_list=="Rasch"&crosslist[[i]]$var_list=="gamma")
    crosslist=crosslist[-i]
  }

### RUN LME MODELS ###
for(i in 1:length(crosslist)){
    assign(paste0(crosslist[[i]][1],"_",crosslist[[i]][2]),final_dat%>%subset(Model==crosslist[[i]][1]&Par==crosslist[[i]][2])%>%lmer(scaled_score~Age*Domain+(1|ID),.))
}


