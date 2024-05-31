args = commandArgs(TRUE)

outdir=as.character(args[1])
deltafile=as.character(args[2])
errorfile=as.character(args[3])
intfile=as.character(args[4])
modelfile=as.character(args[5])
roifile=as.character(args[6])
n_subjects=as.numeric(args[7])
n_icc=as.numeric(args[8])

# outdir = "/ess/p274/cluster/projects/p039_image_brain_change/data_reliability_long/df_mri/all/computational/aggregated/modality.model1.thickness"
# deltafile = here('data-raw/tabulated/parameters_reliability/delta.change.rda')
# errorfile = here('data-raw/tabulated/parameters_reliability/s2c_error_matrix.Rda')
# intfile = here('data-raw/tabulated/parameters_reliability/s2c_Int_matrix.Rda')
# modelfile = here('data_reliability_long/df_mri/all/computational/model1.csv')
# roifile = file.path(outdir, "feature.roi.tsv")
# n_subjects = 250
# n_icc = 3

reliability_wrapper = function(outdir, deltafile, errorfile, intfile, modelfile, roifile, n_subjects, n_icc) {
  library(MASS)
  library(rio)
  library(tidyverse)
  library(psych)
  
  features.roi = read_tsv(roifile, col_names = F) %>% .$X1
  # deltfa matrix
  load(deltafile)
  # error matrix
  load(errorfile)
  # Int matrix
  load(intfile)
  # import df.model
  df.model = rio::import(modelfile)
  
  set.seed(1234)
  
  # select data of interest
  idx = match(features.roi, rownames(corr.n3))
  Mdelta = corr.n3[idx,idx]
  
  idx = match(features.roi, rownames(M.Int))
  Mmean = M.Int[idx,idx]
  
  idx = match(features.roi, rownames(M.error))
  Merror = M.error[idx,idx]
  
  idx = match(features.roi, df.model$feature)
  Mmodel = df.model[idx,]
  
    
  # set up simulation parameters
  sim_params <- crossing(
    tps = seq(from = 3, to = 9, by = 2),
    total_time = seq(from = 2, to = 12, by = 2)
  )
  
  sim_params.extended = 
    sim_params %>% 
    rowwise() %>% 
    mutate(
      time = list(seq(from = 0, 
                      to = total_time, 
                      length.out = tps))) %>% 
    unnest(time) %>% 
    unite("id", 
          c("tps", "total_time", "time"), 
          sep = "-", 
          remove = F)
  
   df = prepare_Int_Beta(n_subjects, Mdelta, Mmean, Mmodel)
  
   
   # prepare error
   error = prepare_error(sim_params.extended, n_subjects, Merror, Mmodel)
  
   # main routine
   x = main_routine(df,error)
   simres = x[[2]]
   
   # summary
   simres.summary = get_simres_summary(x[[2]], Mmodel)
   simres.summary.All = get_simres_summary(x[[1]], Mmodel, T)
   
    # lm
    mod = glm(std.error ~ total_time*as.numeric(tps), data = x[[2]], family = gaussian(link = "log"))
    df.glm = broom.mixed::tidy(mod)
  
    
    df.icc = wrapper_icc(n_icc, sim_params, n_subjects, Merror, Mmodel, df, sim_params.extended)
  
    
    # plotting
    # gs1 = ggplot(simres.summary, aes(x = total_time, y = valueseD, 
    #                                  group = tps, color = tps)) + 
    #   geom_point() + 
    #   geom_line() + 
    #   facet_wrap(vars(name), scales = "free") 
    #   
    # 
    #   
    # gs2 = ggplot(df.icc, aes(x = total_time, y = icc21, 
    #                          group = tps, color = tps)) + 
    #   geom_point() + 
    #   geom_line()
    # 
    # gs3 = ggplot(df.icc, aes(x = total_time, y = icc2k, 
    #                          group = tps, color = tps)) + 
    #   geom_point() + 
    #   geom_line()
    # 
    # prepare savings
    df.out = list()
    df.out$icc = df.icc[[1]]
    df.out$icc.all = df.icc[[2]]
    df.out$simres.summary = simres.summary
    df.out$simres.summary.All = simres.summary.All
    df.out$df.glm = df.glm
    # df.out$fig$gs1 = gs1
    # df.out$fig$gs2 = gs2
    # df.out$fig$gs3 = gs3
    # 
    save("df.out", 
         file = file.path(outdir, "reliability.Rda"))
}

make_error = function(x, n_subjects, Merror, Mmodel) {
  grot = mvrnorm(n_subjects, rep(0, dim(Merror)[1]), Merror)
  grot =scale(grot, scale = 1/(Mmodel$meanF*Mmodel$meanE/100))
  error = data.frame(
    id = seq_len(n_subjects),
    grot)
  return(error)
}

prepare_Int_Beta = function(n_subjects, Mdelta, Mmean, Mmodel) {
  # prepare Beta parameters
  beta = mvrnorm(n_subjects, rep(0, dim(Mdelta)[1]), Mdelta)
  beta =scale(beta, scale = 1/(Mmodel$meanF*Mmodel$seD/100))
  beta = sweep(beta, 2, (Mmodel$meanF*Mmodel$meanD/100), FUN = "+")
  beta = data.frame(
    id = seq_len(n_subjects),
    beta
  ) %>% 
    pivot_longer(-id, names_to = "features", values_to = "Beta")
  
  # prepare Int Parameters
  Int = mvrnorm(n_subjects, rep(0, dim(Mmean)[1]), Mmean)
  Int =scale(Int, scale = 1/Mmodel$seF)
  Int = sweep(Int, 2, Mmodel$meanF, FUN = "+")
  Int = data.frame(
    id = seq_len(n_subjects),
    Int
  ) %>% 
    pivot_longer(-id, names_to = "features", values_to = "Int")
  
  # merge int and beta
  df = left_join(Int, beta)
  return(df)
}

prepare_error = function(sim_params.extended, n_subjects, Merror, Mmodel) {
  error = lapply(1:dim(sim_params.extended)[1], 
                 make_error, 
                 n_subjects = n_subjects, 
                 Merror = Merror, 
                 Mmodel = Mmodel)
  names(error) = sim_params.extended$id
  error = data.table::rbindlist(error, idcol = "grot") %>% 
    separate(grot, c("tps", "total_time", "time"), remove = T, sep = "-", convert = T) %>% 
    pivot_longer(-c(id, tps, total_time, time), 
                 names_to = "features", 
                 values_to = "error")
  return(error)
}

main_routine = function(df, error) {
  x = full_join(df, error) %>% 
    mutate(
      y = Int + Beta*time + error) 
  
  xAll =
    x %>% group_by(id, tps, total_time, features, time)
  
  xMean =
    x %>% 
    group_by(id, tps, total_time, time) 
  
  xAll =
    xAll %>% 
    summarise(y = mean(y)) %>% 
    nest() %>%
    mutate(mod = map(data, ~ lm(y ~ time, data = .x)), 
           tidy = map(mod, ~broom::tidy(.x))) %>% 
    unnest(tidy) %>% 
    filter(term == "time") %>% 
    dplyr::select(estimate, std.error)
  
  xMean =
    xMean %>% 
    summarise(y = mean(y)) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ lm(y ~ time, data = .x)), 
           tidy = map(mod, ~broom::tidy(.x))) %>% 
    unnest(tidy) %>% 
    filter(term == "time") %>% 
    dplyr::select(estimate, std.error)
  
  x = list(xAll, 
       xMean)
  return(x)
}

get_simres_summary = function(simres, Mmodel,is.all = F) {
  
  if ( is.all == F) {
    simres.summary = simres %>% 
      group_by(tps, total_time)
  } else {
    simres.summary = simres %>% 
      group_by(tps, total_time, features)
  }
  
  simres.summary = 
    simres.summary %>% 
    summarise(
      avg_beta_se = mean(std.error)
    ) %>% 
    pivot_longer(cols = contains("beta"))
  simres.summary$seD = mean(Mmodel$seD)
  simres.summary = 
    simres.summary %>% 
    mutate(valueseD = value /seD)
  return(simres.summary)
}

wrapper_compute_y = function(seed = 1, sim_params.extended, n_subjects, Merror, Mmodel,df) {
  error = prepare_error(sim_params.extended, n_subjects, Merror, Mmodel)
  simres = main_routine(df,error)
  return(simres)
}


wrapper_icc = function(n_icc, sim_params, n_subjects, Merror, Mmodel, df, sim_params.extended) {
  df.out = lapply(1:n_icc, wrapper_compute_y, sim_params.extended = sim_params.extended, n_subjects = n_subjects, Merror = Merror,  Mmodel = Mmodel, df = df)
  names(df.out) = 1:n_icc
  df.out.mean = 
    lapply(df.out, function (x) x[[2]])  %>% 
    data.table::rbindlist(., idcol = "it") %>%  
    dplyr::select(-c(std.error)) %>% 
    pivot_wider(names_from = "it", 
                values_from = "estimate", 
                names_prefix = "it")
  
  df.out.all = 
    lapply(df.out, function (x) x[[1]])  %>% 
    data.table::rbindlist(., idcol = "it") %>%  
    dplyr::select(-c(std.error)) %>% 
    pivot_wider(names_from = "it", 
                values_from = "estimate", 
                names_prefix = "it")
  
  
  df.icc = sim_params
  for (i in 1:dim(df.icc)[1]) {
    grot = 
      df.out.mean %>% filter(tps == df.icc$tps[[i]] & total_time == df.icc$total_time[[i]]) %>% 
      dplyr::select(-c("id", "total_time", "tps"))
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
  
  df.icc.all = expand.grid(tps = unique(sim_params$tps), total_time = unique(sim_params$total_time), features = unique(df.out.all$features))
  
  for (i in 1:dim(df.icc.all)[1]) {
    grot = 
      df.out.all %>% filter(tps == df.icc.all$tps[[i]] & total_time == df.icc.all$total_time[[i]] & features == df.icc.all$features[[i]]) %>% 
      dplyr::select(-c("id", "total_time", "tps", "features"))
    df.icc.all$icc[[i]] <- ICC(grot)     
  }
  
  df.icc.all$icc21 = sapply(df.icc.all$icc, function(x) {x$results$ICC[2]})
  df.icc.all$icc21.p = sapply(df.icc.all$icc, function(x) {x$results$p[2]})
  df.icc.all$icc21.ci = sapply(df.icc.all$icc, function(x) {x$results$ICC[2] - x$results$`lower bound`[2]})
  df.icc.all$icc2k = sapply(df.icc.all$icc, function(x) {x$results$ICC[5]})
  df.icc.all$icc2k.p = sapply(df.icc.all$icc, function(x) {x$results$p[5]})
  df.icc.all$icc2k.ci = sapply(df.icc.all$icc, function(x) {x$results$ICC[5] - x$results$`lower bound`[5]})
  df.icc.all = 
    df.icc.all %>% mutate(tps = as.factor(tps))
  
  out = list(df.icc, 
             df.icc.all)
  return(out)
}




#
reliability_wrapper(outdir, deltafile, errorfile, intfile, modelfile, roifile, n_subjects, n_icc) 
