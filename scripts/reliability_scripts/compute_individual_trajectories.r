args = commandArgs(TRUE)

outdir=as.character(args[1])
model=as.character(args[2])
feature=as.character(args[3])
int=as.numeric(args[4])
seD=as.numeric(args[5])
meanD=as.numeric(args[6])
err=as.numeric(args[7])


compute_individual_trajectories = function(outdir, model, feature, int, seD, meanD, err){
  library(tidyverse)
  if ( meanD > 0) { meanD = - meanD}
  
  df.params = list(
    outdir = outdir, 
    model = model, 
    feature = feature, 
    int = int, 
    seD = seD, 
    meanD = meanD, 
    err = err
  )
  beta = meanD + c(-1,0,1) * seD
  n = 5000
  
  
  sim_params = crossing(
    tps = seq(3,9,2), 
    total_time = seq(2,12,2)
  )
  
  sim_ind = lapply(beta, function(x) {
    simres = pmap_dfr(sim_params, function(tps, total_time,beta) {
      tibble(
        id = seq_len(n),
        beta = int*0.01*x,
        Int = int
      ) %>% 
        uncount(tps) %>% 
        nest_by(id) %>% 
        mutate(
          time = list(seq(0, total_time, length.out = tps))
        ) %>% 
        unnest(cols = c(data,time)) %>% 
        ungroup() %>% 
        mutate(
          y = Int + beta*time + rnorm(nrow(.), 0, int*0.01*err)
        ) %>% 
        nest_by(id, .keep = TRUE) %>% 
        pmap_dfr(function(id, data) {
          mod = lm(y ~time, data = data)
          tibble(id = id, 
                 beta_true = unique(data$beta),
                 beta_est = coef(mod)[["time"]],
                 Int = unique(data$Int)
          ) 
        })%>% 
        mutate(tps = factor(tps), total_time = total_time)
    })
  })
  
  sim_ind =  data.table::rbindlist(sim_ind) %>% 
    mutate(beta_true = 100*beta_true/Int, 
           beta_est = 100*beta_est/Int)
  
  
  subset = 
    sim_ind[abs(sim_ind$beta_true - meanD) < 0.0001,] %>% 
    group_by(tps, total_time) %>% 
    summarise(sd_obs_beta = sd(beta_est))
  
  i = 1
  nsubs = 1000000
  
  df.extreme = lapply(1:dim(subset)[1], likelihood_outliers, nsubs = nsubs, subset = subset, meanD = meanD, seD = seD)
  names(df.extreme) =  paste(subset$tps, subset$total_time, sep = ".")
  df.extreme = 
    data.table::rbindlist(df.extreme, idcol = "grot") %>% 
    separate(grot, c("tps", "total_time"),  
             remove = T, 
             convert = T) 
  
  P_extreme = df.extreme %>% 
    group_by(tps, total_time) %>% 
    summarise(pct_extrem = n()/nsubs,
              probM_gvEx = sum(true_beta>0)/n(),
              probD_gvEx = sum(true_beta<meanD)/n()) 
  
  
  
  save(sim_ind,
       df.params, 
       df.extreme, 
       P_extreme,
       file = file.path(outdir, 
                        paste0(model,
                               "_",
                              feature, 
                              ".Rda")))

}


likelihood_outliers = function(i, nsubs, subset, meanD, seD) {
  error = subset$sd_obs_beta[i]
  value_normal = meanD + seD
  true_beta = rnorm(nsubs, meanD, seD)
  P = ecdf(true_beta)
  normal = P(value_normal)
  extreme = 1-P(0)
  observed = true_beta + rnorm(nsubs, sd = error)
  df = data.frame(true_beta, observed)
  
  df.extreme = df %>% filter(observed > 0)
  df.extreme = 
    df.extreme %>% 
    mutate(class = if_else(true_beta < value_normal,1,
                           if_else(true_beta > 0,-1,0)))
  
  return(df.extreme)
}

compute_individual_trajectories(outdir, model, feature, int, seD, meanD, err)
