args = commandArgs(TRUE)

outdir=as.character(args[1])
analysis=as.character(args[2])
feature=as.character(args[3])
oname=as.character(args[4])
meanF=as.numeric(args[5])
seF=as.numeric(args[6])
meanD=as.numeric(args[7])
seD=as.numeric(args[8])
meanE=as.numeric(args[9])
n_subjects=as.numeric(args[10])
n_icc=as.numeric(args[11])

wrapper_computational_reliability = function(outdir, analysis, feature, oname, meanF, seF, meadD, seD,meanE, n_subjects, n_icc) {
  library(tidyverse)
  library(psych)
  
  set.seed(123)
  idpdir = file.path(outdir, analysis, feature)
  try(dir.create(idpdir))

  sim_params <- crossing(
    tps = seq(from = 3, to = 9, by = 2),
    total_time = seq(from = 2, to = 12, by = 2)
  )
  
  simres <- pmap_dfr(sim_params, function(tps, total_time){
    tibble(
      id = seq_len(n_subjects),
      beta = rnorm(n_subjects, meanF*meanD*0.01, sd = meanF*seD*0.01),
      Int = rnorm(n_subjects, meanF, sd = seF)
    ) %>% 
      uncount(tps) %>% 
      nest_by(id) %>% 
      mutate(
        time = list(seq(from = 0, to = total_time, length.out = tps))
      ) %>% 
      unnest(cols = c(data, time)) %>% 
      ungroup() %>% 
      mutate(
        y = Int + beta * time + rnorm(nrow(.), 0, sd = meanF*meanE*0.01)
      ) %>% 
      nest_by(id, .keep = TRUE) %>% 
      pmap_dfr(function(id, data){
        mod <- lm(y ~ time, data = data)
        tibble(
          id = id,
          beta_true = unique(data$beta),
          beta_est = coef(mod)[["time"]],
          beta_se = sqrt(vcov(mod)[["time", "time"]]),
          Int_true = unique(data$Int),
        )
      }) %>% 
      mutate(tps = factor(tps), total_time = total_time)
  })

  
  simres.summary = simres %>% 
    group_by(tps, total_time) %>% 
    summarise(
      beta_rmse = sqrt(mean((beta_true - beta_est)^2)),
      avg_beta_se = mean(beta_se),
      seD.obs = sd(beta_est),
      seD.real = sd(beta_true)
    )
  simres.summary$seD = seD
  simres.summary$meanE = meanE
  simres.summary$meanF = meanF
  # simres.summary = 
  #   simres.summary %>% 
  #   mutate(valueseD = value /(meanF*seD*0.01))

  # mod = glm(beta_se ~ total_time*as.numeric(tps), data = simres, family = gaussian(link = "log"))
  # df.glm = broom.mixed::tidy(mod)
  
  # gs1 = ggplot(simres.summary, aes(x = total_time, y = valueseD, 
  #                                  group = tps, color = tps)) + 
  #   geom_point() + 
  #   geom_line() + 
  #   facet_wrap(vars(name), scales = "free")  
  
  df.out = lapply(1:n_icc, repeat_simres, simres = simres, meanE = meanE, meanF = meanF)
  names(df.out) = 1:n_icc
  df.out = data.table::rbindlist(df.out, idcol = "it") %>% 
    separate(idall, c("id", "total_time", "tps")) %>% 
    select(-c(beta_true,beta_se)) %>% 
    pivot_wider(names_from = "it", 
                values_from = "beta_est", 
                names_prefix = "it")
  
  df.icc = sim_params
  for (i in 1:dim(df.icc)[1]) {
    grot = 
      df.out %>% filter(tps == df.icc$tps[[i]] & total_time == df.icc$total_time[[i]]) %>% 
      select(-c("id", "total_time", "tps"))
    df.icc$icc[[i]] <- ICC(grot)     
  }
  
  df.icc$icc21 = sapply(df.icc$icc, function(x) {x$results$ICC[2]})
  df.icc$icc21.p = sapply(df.icc$icc, function(x) {x$results$p[2]})
  df.icc$icc21.ci = sapply(df.icc$icc, function(x) {x$results$ICC[2] - x$results$`lower bound`[2]})
  df.icc$icc2k = sapply(df.icc$icc, function(x) {x$results$ICC[5]})
  df.icc$icc2k.p = sapply(df.icc$icc, function(x) {x$results$p[5]})
  df.icc$icc2k.ci = sapply(df.icc$icc, function(x) {x$results$ICC[5] - x$results$`lower bound`[5]})
  df.icc = 
    df.icc %>% mutate(tps = as.factor(tps))
  
  
  # gs2 = ggplot(df.icc, aes(x = total_time, y = icc21, 
  #                                  group = tps, color = tps)) + 
  #   geom_point() + 
  #   geom_line()
  # 
  # gs3 = ggplot(df.icc, aes(x = total_time, y = icc2k, 
  #                          group = tps, color = tps)) + 
  #   geom_point() + 
  #   geom_line()
  
  
  df.out = list()
  df.out$icc = df.icc
  df.out$simres.summary = simres.summary
  # df.out$df.glm = df.glm
  # df.out$fig$gs1 = gs1
  # df.out$fig$gs2 = gs2
  # df.out$fig$gs3 = gs3
  
  save("df.out", 
       file = file.path(idpdir, 
                        paste0("comp.", oname, ".Rda")))
    
}


repeat_simres = function(seed = 1, simres, meanE, meanF) {
  set.seed(seed)
  xx = simres %>% select(id, beta_true, Int_true, tps, total_time) %>% 
    mutate(tps= as.numeric(as.character(tps)), 
           meanE = meanE) %>%  
    uncount(tps, .remove = F) %>% 
    group_by(id, total_time, tps) %>%
    mutate(numbering = row_number(), 
           time = total_time*(numbering-1)/(tps-1)) %>% 
    ungroup() %>% 
    mutate(
      y = Int_true + beta_true * time + rnorm(nrow(.), 0, meanF*0.01*meanE)
    ) %>% 
    mutate(idall = paste(id, total_time, tps)) %>% 
    nest_by(idall,.keep = TRUE) %>% 
    pmap_dfr(function(idall, data){
      mod <- lm(y ~ time, data = data)
      tibble(
        idall = idall,
        beta_true = unique(data$beta_true),
        beta_est = coef(mod)[["time"]],
        beta_se = sqrt(vcov(mod)[["time", "time"]]),
      )
    })
  return(xx)
}


wrapper_computational_reliability(outdir, analysis, feature, oname, meanF, seF, meadD, seD,meanE, n_subjects, n_icc)




