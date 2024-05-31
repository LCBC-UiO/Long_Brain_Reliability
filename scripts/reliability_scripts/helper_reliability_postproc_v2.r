squeue = function(user, job) {
  df.squeue = system(paste0("squeue --name=", job," -u ",user), intern = T) %>% 
    strsplit(., " +") %>% 
    simplify2array() %>% 
    t() %>% 
    as.data.frame()
  return(df.squeue)
  
}

retrieve_rm_vars_error = function() {
  rm_vars_error = c(
    "CSF",
    "SubCortGrayVol",
    "TotalGrayVol" ,
    "SupraTentorialVol",
    "SupraTentorialVolNotVent",
    "Right-Inf-Lat-Vent",
    "Right-Cerebellum-White-Matter",
    "Right-Cerebellum-Cortex",
    "Right-VentralDC",
    "Right-vessel",
    "Right-choroid-plexus",
    "Left-Inf-Lat-Vent",
    "Left-Cerebellum-White-Matter",
    "Left-Cerebellum-Cortex",
    "Left-VentralDC",
    "Left-vessel",
    "Left-choroid-plexus",
    "3rd-Ventricle",
    "4th-Ventricle",
    "Brain-Stem",
    "Left-Accumbens-area",
    "Right-Accumbens-area")
  return(rm_vars_error)
}

retrieve_rm_globalVars_error = function() {
  rm_globalvars_error = c(
    "MeanThickness",
    "CorticalArea",
    "CorticalVolume",
    "SubCorticalVolume",
    "SupraTentorialVolume")
  return(rm_globalvars_error)
}

select_parameter_reliability = function(dat) {
  rm_vars_error = retrieve_rm_vars_error()  
  rm_globalvars_error = retrieve_rm_globalVars_error()  
  
  dat = 
    dat %>% 
    filter(!feature %in% rm_vars_error) %>% 
    filter(!feature %in% rm_globalvars_error)
  return(dat)
}

wrapper_select_parameter_reliability = function(df.parameters, string.meas) {
  dat  = 
    df.parameters %>% 
    dplyr::select(feature, 
                  modality, 
                  contains(string.meas)) 
  names(dat) = gsub(string.meas, "", names(dat))
  if(string.meas == "pct.err.mean_") {
  dat = 
    dat %>% select(-all) }
  else if (string.meas == "sdDelta_") {
    dat = 
      dat %>% select(-c(Y, O, wayne))
  }
  dat = select_parameter_reliability(dat)
  return(dat)
}

pairwise.icc.parameters = function(dat) {
  x = dat %>% select(-c(feature, modality))
  c = 0
  icc.pw = ii = jj = c()
  for (i in 1:dim(x)[2]) {
    for (j in 1:dim(x)[2]) {
      if(i == j) next
      ss = ICC(x[,c(i,j)])
      c = c + 1
      icc.pw = c(icc.pw, ss$results$ICC[2])
      jj = c(jj, names(x)[j])
      ii = c(ii, names(x)[i])
    }
  }
  df = data.frame(dataset.1 = ii,
                  dataset.2 = jj,
                  icc = icc.pw)
  return(df)
}

icc_reliability_full = function(dat.me, dat.slopeV) {
  
  dat.me.long = 
    dat.me %>% 
    pivot_longer(-c(feature, modality), 
                 names_to = "dataset.me", 
                 values_to = "meanE")
  
  dat.slopeV.long = 
    dat.slopeV %>% 
    pivot_longer(-c(feature, modality), 
                 names_to = "dataset.slope", 
                 values_to = "seD")
  
  model = full_join(dat.me.long, dat.slopeV.long)
  model$model = "full"
  
  df.formula.icc = expand.grid(feature = unique(model$feature),
                               model = unique(model$model), 
                               t = seq(3,9,by = 2), 
                               total_time = seq(2,12, by =2)) %>% 
    left_join(model, .) %>% 
    mutate(ws = meanE^2*((total_time^2*t*(t+1))/ (12*(t-1)))^-1, 
           bs = seD^2, 
           icc = bs /(ws + bs))
  df.formula.icc = 
    df.formula.icc %>% 
    mutate(set = paste(dataset.me, dataset.slope, sep ="_")) 
  
  save(df.formula.icc, file = file.path(outreliability.comp, "icc.formula.between.rda"))
  return(df.formula.icc)
}

pairwise.reliability.icc = function(df.formula.icc) {
  me = unique(df.formula.icc$dataset.me)
  slope = unique(df.formula.icc$dataset.slope)
  me.combn = combn(me, 2)
  slope.combn = combn(slope, 2)
  out = list()
  c = 0
  l = dim(me.combn)[2] * dim(slope.combn)[2] * 2
  pb = txtProgressBar(0, l)
  for (i in 1:dim(me.combn)[2]) {
    for (j in 1:dim(slope.combn)[2]) {
      for (k in 1:2) {
        c = c + 1
        setTxtProgressBar(pb, c)
        if (k == 1) {
          dat = df.formula.icc %>% 
            filter(set == paste(me.combn[1,i], slope.combn[1,j], sep = "_") | 
                     set == paste(me.combn[2,i], slope.combn[2,j], sep = "_"))
        } else if (k == 2) {
          dat = df.formula.icc %>% 
            filter(set == paste(me.combn[1,i], slope.combn[2,j], sep = "_") | 
                     set == paste(me.combn[2,i], slope.combn[1,j], sep = "_"))
        }
        dat = 
          dat %>% 
          select(feature, modality, set, t, total_time, icc) %>% 
          pivot_wider(values_from = "icc", 
                      names_from = set)  
        
        icc.pw = ICC(dat[,5:6])$results$ICC[2]
        jj = names(dat)[5]
        ii = names(dat)[6]
        
        dat.sub =dat %>%  
          group_by(t, total_time) %>% 
          nest() %>% 
          mutate(mod = map(data, ~ ICC(.[, 3:4])), 
                 icc.sub = map_dbl(mod, ~ .$results$ICC[2])) %>% 
          select(-c(data, mod)) %>%
          mutate(jj = names(dat)[5], 
                 ii = names(dat)[6]) %>%  
          group_by(jj, ii) %>% 
          nest()
        
        df = data.frame(icc = icc.pw, ii = ii, jj = jj)
        df =left_join(df, dat.sub)
        out[[c]] = df
      }
    }
  }
  out = data.table::rbindlist(out)
  return(out)
}

icc2.reliability.icc = function(dat) {
  set.seed(123)
  n = 500
  me = unique(dat$dataset.me)
  slope = unique(dat$dataset.slope)
  slope.combn <- gtools::permutations(n = length(slope), r = 6, v = slope)
  
  nn = sample.int(dim(slope.combn)[1], n )
  out = list()
  pb = txtProgressBar(0, n)
  for (i in 1:n) {
    setTxtProgressBar(pb,i)
    idx = nn[i]
    grot = paste(me,slope.combn[idx,], sep = "_")
    x = dat %>% 
      filter(set %in% grot) %>% 
      select(feature, modality, set, t, total_time, icc) %>% 
      pivot_wider(values_from = "icc", 
                  names_from = set)  
    icc.pw = ICC(x[,5:10])$results$ICC[5]
    
    dat.sub =x %>%  
      mutate(i) %>% 
      group_by(t, total_time,i) %>% 
      nest() %>% 
      mutate(mod = map(data, ~ ICC(.[, 3:8])), 
             icc.sub = map_dbl(mod, ~ .$results$ICC[5])) %>% 
      select(-c(data, mod)) %>% 
      ungroup() %>% 
      group_by(i) %>% 
      nest()
    
    permn = list(grot)
    df = data.frame(icc = icc.pw, 
                    dat.sub, 
                    tibble(permn))
    
    out[[i]] = df
  }
  out = data.table::rbindlist(out)
  return(out)
}

one2all.icc.parameters = function(dat) {
  x = dat %>% select(-c(feature, modality))
  c = 0
  icc.pw = ii = c()
  for (i in 1:dim(x)[2]) {
    icc.row1 = x[,i]
    icc.row2 = rowMeans(x[,-1])
    ss = ICC(cbind(icc.row1, icc.row2))
    icc.pw = c(icc.pw, ss$results$ICC[2])
    ii = c(ii, names(x)[i])
  }
  df = data.frame(dataset.1 = ii,
                  icc = icc.pw)
  return(df)
}
one2all.icc.reliability  = function(dat) {
  out = list()
  grot = expand.grid(me = unique(dat$dataset.me), 
                     slope = unique(dat$dataset.slope))
  
  for (i in 1:dim(grot)[1]) {
    dat.1 = dat %>% filter(dataset.me == grot$me[i] & dataset.slope == grot$slope[i])
    dat.2 = dat %>% filter(!dataset.me == grot$me[i] & !dataset.slope == grot$slope[i]) %>% 
      group_by(feature, modality, t, total_time) %>% 
      summarise(icc.avg = mean(icc, na.rm = T))
    dat.icc = left_join(dat.1, dat.2)
    
    icc.pw = ICC(cbind(dat.icc$icc, dat.icc$icc.avg))$results$ICC[2]
    
    dat.sub =dat.icc %>%  
      group_by(t, total_time, set) %>% 
      nest() %>% 
      mutate(mod = map(data, ~ ICC(cbind(.$icc, .$icc.avg))), 
             icc.sub = map_dbl(mod, ~ .$results$ICC[2])) %>% 
      select(-c(data, mod)) %>% 
      ungroup() %>% 
      group_by(set) %>% 
      nest()
    
    df = data.frame(icc = icc.pw, 
                    dat.sub, 
                    dataset.me = grot$me[i], 
                    dataset.slope = grot$slope[i])
    
    out[[i]] = df
  }
  out = data.table::rbindlist(out)
  return(out)
}


scale_sd_change = function(se, db) {
  sd.delta = apply(delta[db$n >2 & db$time>2,], 2, sd, na.rm = T)
  df.se = sweep(se, 2, sd.delta, FUN = "/")
  df.se = 
    cbind(
      db %>% filter(n >2, time> 2),
      df.se)
  return(df.se)
}


prepare_se_data = function(se,db) {
  se = se[db$n > 2 & db$time>2, ]
  phenos = names(se)
  # remove extreme error measures
  mad.M = apply(se, 2,median, na.rm = T)
  mad.se = apply(se, 2,mad, na.rm = T)
  mad.se.idx = sweep(se, 2,mad.M, FUN = "-")
  mad.se.idx = sweep(mad.se.idx, 2,mad.se, FUN = "/")
  mad.se.idx = abs(mad.se.idx) > 5
  grot = table(mad.se.idx)
  100*grot[1]/sum(grot)  
  se[mad.se.idx] = NaN
  return(se)
}




glmer_empirical = function(df.se, phenos.se, outreliability, df.harmonize) {
  df.out = list()
  timengrid = 
    expand.grid(
      tps=c(3,5,7,9), 
      total_time = c(2,4,6,8,10,12))
  
  
  
  pb = txtProgressBar(0,length(phenos.se))
  for (i in 1:length(phenos.se)) {
    setTxtProgressBar(pb, i)
      try(df.out$glmer[[i]] <- glmer(get(phenos.se[i]) ~ total_time*tps + (1|dataset) , data = df.se, family = gaussian(link = "log")))    
  }
  names(df.out$glmer) = phenos.se
  
  # get estimates models glmer
  df.out.tidy.glmer = lapply(df.out$glmer, function(x) {broom.mixed::tidy(x)})
  names(df.out.tidy.glmer) = phenos.se
  df.out.tidy.glmer = data.table::rbindlist(df.out.tidy.glmer, idcol = "feature") %>% 
    filter(effect == "fixed") %>% 
    mutate(log10p = -log10(p.value), 
           log10psigned = if_else(statistic > 0, 
                                  log10p,
                                  -log10p))
  
  # extract predictions
  grot = lapply(df.out$glmer, 
                function(x) {cbind(timengrid,predict(x, timengrid, re.form=~0, type = "response"))})
  names(grot) = phenos.se
  df.out.grid.glmer = data.table::rbindlist(grot, idcol = "feature") 
  names(df.out.grid.glmer)[4] = "fit"
  df.out.grid.glmer = left_join(df.out.grid.glmer, df.harmonize)  
  
  # get summary for predictions
  df.emp.summary.glmer = 
    df.out.grid.glmer %>% 
    group_by(tps, total_time) %>% 
    summarize(mean.error = mean(fit), 
              sd.error = sd(fit)) %>% 
    mutate(t = as.factor(tps))
  

  
   
  df.emp = list()
  df.emp$phenos.se = phenos.se
  df.emp$df.se = df.se
  df.emp$df.out = df.out
  df.emp$df.out.tidy.glmer = df.out.tidy.glmer 
  df.emp$df.out.grid.glmer = df.out.grid.glmer
  df.emp$df.emp.summary.glmer = df.emp.summary.glmer
 
    save(df.emp,
         file = file.path(outreliability, "df.empirical.Rda"))
  
  return(df.emp)
}


gamlss_slope = function(df.delta, outreliability, phenos.delta) {
  library(gamlss)
  library(broom.mixed)
  
  
  tps_tp_df <- expand.grid(
    tps = seq(from = 3, to = 9, by = 2), 
    total_time = seq(2,12, by = 2))
  
  
  x = df.delta %>% 
    select(rid, 
           dataset, 
           total_time, 
           tps,
           all_of(phenos.delta)) %>% 
    pivot_longer(phenos.delta, 
                 names_to = "feature", 
                 values_to = "delta") %>% 
    drop_na() %>% 
    mutate(ds = as.factor(dataset)) %>% 
    group_by(feature) %>% 
    mutate(delta = delta - mean(delta, na.rm = T))  %>% 
    nest()
  
  l = length(x$feature)
  grot = list()
  for (i in 1:l) {
    grot$feature[i] = x$feature[i]  
    grot.dat = x$data[[i]]
    grot$mod1[[i]] =    gamlss(formula = delta~ 1 +random(ds),
                               sigma.formula = delta ~  cs(total_time, df = 1) + tps + cs(total_time, df = 1):tps,
                               data = grot.dat,
                               sigma.link = "log")
    grot$mod2[[i]] =   gamlss(formula = delta~ 1 ,
                              sigma.formula = delta ~  cs(total_time, df = 3) + tps + cs(total_time, df = 3):tps,
                              data = grot.dat)
    grot$mod3[[i]] =  gamlss(formula = delta~ 1 +random(ds),
                             sigma.formula = delta ~  cs(total_time, df = 1) + tps + cs(total_time, df = 1):tps,
                             data = grot.dat)
    # grot$mod4[[i]] =  gamlss(formula = delta~ 1,
    #                          sigma.formula = delta ~  cs(total_time, df = 2),
    #                          data = grot.dat)
    # grot$mod5[[i]] =  gamlss(formula = delta~ 1,
    #                          sigma.formula = delta ~  cs(total_time, df = 2) + cs(tps, df = 2),
    #                          data = grot.dat)
    # grot$mod6[[i]] =  gamlss(formula = delta~ 1,
    #                          sigma.formula = delta ~  cs(total_time, df = 3) + cs(tps, df =32) + cs(total_time, df = 3):cs(tps, df = 3),
    #                          data = grot.dat)
    grot$pred.mod1[[i]] =  exp(predict(grot$mod1[[i]], what = "sigma", newdata = tps_tp_df))
    grot$pred.mod2[[i]] =  exp(predict(grot$mod2[[i]], what = "sigma", newdata = tps_tp_df))
    grot$pred.mod3[[i]] =  exp(predict(grot$mod3[[i]], what = "sigma", newdata = tps_tp_df))
    # grot$pred.mod4[[i]] =  exp(predict(grot$mod4[[i]], what = "sigma", newdata = tps_tp_df))
    # grot$pred.mod5[[i]] =  exp(predict(grot$mod5[[i]], what = "sigma", newdata = tps_tp_df))
    #grot$pred.mod6[[i]] =  exp(predict(grot$mod6[[i]], what = "sigma", newdata = tps_tp_df))
    
  }
  
  
  df.bs = 
    tibble(feature = grot$feature,
           mod1 = grot$mod1,
           pred1 = grot$pred.mod1,
           mod2 = grot$mod2,
           pred2 = grot$pred.mod2,
           mod3 = grot$mod3, 
           pred3 = grot$pred.mod3) 
           #mod4 = grot$mod4, 
           # pred4 = grot$pred.mod4,
           # mod5 = grot$mod5, 
           # pred5 = grot$pred.mod5, 
           #mod6 = grot$mod6, 
           #pred6 = grot$pred.mod6)
  
  df.bs = 
    df.bs %>% 
    mutate(tidy.mod1 = map(mod1, ~broom.mixed::tidy(.)),
           tidy.mod2 = map(mod2, ~broom.mixed::tidy(.)), 
           tidy.mod3 = map(mod3, ~broom.mixed::tidy(.))) %>% 
           # tidy.mod4 = map(mod4, ~broom.mixed::tidy(.)),
           # tidy.mod5 = map(mod5, ~broom.mixed::tidy(.)),
           #tidy.mod6 = map(mod6, ~broom.mixed::tidy(.))) %>% 
    #select(-c(mod1, mod2, mod3, mod4, mod5, mod6))
    select(-c(mod1, mod2, mod3))
  
  
  df.bs = 
    df.bs %>% mutate(pred1 = map(pred1, ~cbind(.x, tps_tp_df)),
                     pred2 = map(pred2, ~cbind(.x, tps_tp_df)),
                     pred3 = map(pred3, ~cbind(.x, tps_tp_df)))
                     #pred4 = map(pred4, ~cbind(.x, tps_tp_df)), 
                     #pred5 = map(pred5, ~cbind(.x, tps_tp_df)), 
                     #pred6 = map(pred6, ~cbind(.x, tps_tp_df)))
  
  # df.bs = 
  #   df.bs %>% 
  #   mutate(pvalue_tps = map_dbl(tidy.mod1, ~.[4,6][[1]]))
  # df.bs$pfdr = p.adjust(df.bs$pvalue_tps, "fdr")         
  # 
  # df.bs = df.bs %>% 
  #   mutate(model.selected = if_else(pfdr < .05, "mod1", "mod2"), 
  #          pred = if_else(model.selected == "mod1", pred1, pred2))
  # 
  #df.slope  = df.bs  %>% select(feature, pred) %>% unnest(pred)
  df.slope  = df.bs  %>% select(feature, pred3) %>% unnest(pred3)
  names(df.slope)[2] = "pred_bs"
  
  save(df.bs,
       df.slope,
       file = file.path(outreliability, "df.empirical.slope.Rda"))
  
  return(df.slope)
}


rename_model = function(df.model, globalVars = F) {
  names(df.model) = c("feature", 
                      "modality", 
                      "meanF", 
                      "seF", 
                      "meanD", 
                      "seD", 
                      "meanE")
  rm_vars_error = retrieve_rm_vars_error()  
  rm_globalvars_error = retrieve_rm_globalVars_error()  
  
  if (globalVars == F) {
    df.model = 
      df.model %>% 
      filter(!feature %in% rm_vars_error) %>% 
      filter(!feature %in% rm_globalvars_error)
  } else if (globalVars == T) {
    df.model = 
      df.model %>% 
      filter(feature %in% rm_globalvars_error)
  }
  return(df.model)
}

wrapper_submit_model = function(df.model, modelname, nsubs, nicc, outdir, analysis) {
  scriptname =  here("scripts/reliability_scripts/submit_computational_reliability.sh")

    for (i in 1:dim(df.model)[1]) {
    system(paste(
      "module purge;",
      "sbatch",
      scriptname, 
      outdir, 
      analysis,
      paste0("'",df.model$feature[i],"'"),
      modelname, 
      df.model$meanF[i],
      df.model$seF[i],
      df.model$meanD[i],
      df.model$seD[i],
      df.model$meanE[i], 
      nsubs,
      nicc,
      sep = " "))
    
    df.squeue = squeue("p274-didacvp","computational_reliability")
    print(paste("script running on sbatch, N:", length(df.squeue$V1)-1, as.character(i), df.model$feature[i]), sep = " ")
  }
}

summary_computational_model = function(modelname, analysis, outreliability.comp) {
  df.comp = list()
  df.model = 
    rio::import(file.path(outreliability.comp, 
                          paste(modelname, 
                                "csv",
                                sep = ".")))  
  df.icc = df.simres.summary = df.glm.simres = list()
  pb = txtProgressBar(0,dim(df.model)[1])
  for (i in 1:dim(df.model)[1]) {
    setTxtProgressBar(pb,i)
    load(file.path(outreliability.comp, 
                   analysis, 
                   df.model$feature[i],
                   paste("comp", 
                         modelname, 
                         "Rda", 
                         sep = ".")))
    df.icc[[i]] = 
      df.out$icc %>% 
      select(tps, total_time, icc21, icc2k)
    
    df.simres.summary[[i]] = 
      df.out$simres.summary
  }
  names(df.icc) = names(df.simres.summary) = df.model$feature
  df.icc = data.table::rbindlist(df.icc, idcol = "feature")
  df.simres.summary = data.table::rbindlist(df.simres.summary, idcol = "feature")
  
  df.icc = left_join(df.icc, df.harmonize)
  df.simres.summary = left_join(df.simres.summary, df.harmonize)
  
  
  df.comp$icc = df.icc
  df.comp$simres.summary= df.simres.summary 
  
  save(df.comp,
       file = file.path(outreliability.comp, 
                        paste("df", 
                              modelname, 
                              "rda", 
                              sep = ".")))
  
  return(df.comp)
}

wrapper_individual_trajectories = function() {
  scriptname =  here("scripts/reliability_scripts/submit_compute_individual_trajectories.sh")
  indtrajdir = file.path(outreliability.comp, "individual_trajectories")
  try(dir.create(indtrajdir))
  
  for (i in 1:dim(model)[1]) {
    roi = model[i,]
    
    print(paste("script running on sbatch, N:", i , dim(model)[1], roi$feature, roi$model), sep = " ")
    
    system(paste(
      "module purge;",
      "sbatch",
      scriptname, 
      indtrajdir,
      roi$model,
      roi$feature,
      roi$meanF,
      roi$seD,
      roi$meanD,
      roi$meanE,
      sep = " "))
  }
}


figure_scheme_time = function() {
  library(tidyverse)
  library(patchwork)
  library(gghalves)
  
  set.seed(123)
  error.set = .2
  
  df = 
    rbind(
      expand.grid(
        tp = c(1,2), 
        it = c("s1","s2","s3"), 
        type = "short"),
      expand.grid(
        tp = c(1,2), 
        it = c("s1","s2","s3"), 
        type = "long")
    )
  
  df = 
    df %>% 
    mutate(decline = .2,
           time = if_else(tp == 1, 0,
                          if_else(tp == 2 & type == "short", 2,4)),
           error = if_else(it == "s1", 0, 
                           if_else(it == "s2" & tp == 1 | it == "s3" & tp == 2  , -error.set, error.set)),
           y = -time*decline + error)
  
  
  df.norm = 
    expand_grid(
      y = rnorm(1000, 0, .2),
      tp = c(1,2),
      type = c("short", "long")) %>% 
    mutate(decline = .2,
           time = if_else(tp == 1, 0,
                          if_else(tp == 2 & type == "short", 2,4)), 
           y = y - time*decline, 
           it = "s1")
  
  
  df.mod = 
    df %>% 
    group_by(type, it) %>% 
    nest() %>% 
    mutate(mod = map(data, ~lm(y ~time, data = .)), 
           tidy = map(mod, ~broom::tidy(.))) %>% 
    unnest(tidy) %>% 
    filter(term == "time") %>% 
    select(it, type, estimate)
  
  newlabs = c("short" = "2 yrs", "long" = "4 yrs")
  
  
  gs1 = ggplot(df %>% filter(type == "short"), 
               aes(x = time, y = y, group = it, color = it)) +
    geom_half_violin(data = df.norm %>% filter(type == "short"), 
                     mapping = aes(x = time, y = y, group = time), 
                     side = "r",
                     fill = "grey30",
                     alpha = .3,
                     size = 0) + 
    geom_line(aes(linetype = it),
              size = 2,
              position = position_dodge(width = .02)) +
    geom_point(size = 3, 
               position = position_dodge(width = .02),
               shape = 20) + 
    scale_color_manual(values = c("black", "azure2", "bisque")) + 
    #scale_color_manual(values = c("black", "darkolivegreen3", "coral2")) + 
    scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
    ylim(-1.3,.6) + 
    ylab("Observed Brain (a.u.)") +
    xlab("Years") + 
    ggtitle("2 years follow-up")
  
  
  
  
  
  gs2 = ggplot(df %>% filter(type == "long"), 
               aes(x = time, y = y, group = it, color = it)) +
    geom_half_violin(data = df.norm %>% filter(type == "long"), 
                     mapping = aes(x = time, y = y, group = time), 
                     side = "r",
                     fill = "grey30",
                     alpha = .3,
                     size = 0) + 
    geom_line(aes(linetype = it),
              size = 2,
              position = position_dodge(width = .02)) +
    geom_point(size = 3, 
               position = position_dodge(width = .02),
               shape = 20) + 
    scale_color_manual(values = c("black", "azure2", "bisque")) + 
    scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14, face = "bold"),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
    ylim(-1.3,.6) + 
    ylab("") +
    xlab("Years")+ 
    ggtitle("4 years follow-up")
  
  
  gs3 = 
    ggplot(df.mod, aes(x = it, y = estimate, color = it)) +
    geom_hline(yintercept = -.2, linetype = "dotdash", color = "grey50") + 
    geom_point(size = 8, shape = 20) + 
    facet_grid(cols = vars(type), labeller =labeller(.cols = newlabs)) + 
    scale_color_manual(values = c("black", "azure2", "bisque")) + 
    scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          strip.text = element_text(size=12),
          plot.background = element_rect(color = "grey50",
                                         size = 1)) + 
    ylab("Obs. yearly brain change")
  
  fname = here("data_reliability_long/results/figs/scheme/time.png")
  gs = gs1 + gs2 + inset_element(gs3, -0.15, 0, .33, .4)
  ggsave(file = fname, plot = gs)
  return(gs)
}

figure_scheme_obs = function() {
  library(tidyverse)
  library(patchwork)
  library(gghalves)
  
  set.seed(1)
  error.set = .2
  
  df = 
    rbind(
      expand.grid(
        tp = seq(1,2), 
        it = c("s1","s2","s3"), 
        type = "short"),
      expand.grid(
        tp = seq(1,7), 
        it = c("s1","s2","s3"), 
        type = "long")
    )
  
  df = 
    df %>% 
    rowwise() %>% mutate(error = rnorm(1,0,.2)) %>% 
    mutate(decline = .2,
           time = if_else(type == "short", (tp-1)*2,2*(tp-1)/6),
           error =if_else(it == "s1", 0, 
                          if_else(it == "s2" & time == 0 | it == "s3" & time == 2  , -error.set, 
                                  if_else( it == "s2" & time == 2 | it == "s3" & time == 0,error.set, 
                                           if_else(it == "s2",error - .07, error + .07)))),
           y = -time*decline + error)
  
  
  
  df.mod = 
    df %>% 
    group_by(type, it) %>% 
    nest() %>% 
    mutate(mod = map(data, ~lm(y ~time, data = .)), 
           tidy = map(mod, ~broom::tidy(.))) %>% 
    unnest(tidy) %>% 
    filter(term == "time") %>% 
    select(it, type, estimate)
  
  
  df.norm = 
    rbind(expand_grid(
      y = rnorm(1000, 0, .2),
      time = c(0,2),
      type = c("short")),
      expand_grid(
        y = rnorm(1000, 0, .2),
        time = seq(0,2, length.out = 7),
        type = c("long"))) %>% 
    mutate(decline = .2,
           y = y - time*decline, 
           it = "s1")
  
  
  
  
  gs1 = ggplot(df %>% filter(type == "short"), 
               aes(x = time, y = y, group = it, color = it)) +
    geom_half_violin(data = df.norm %>% filter(type == "short"), 
                     mapping = aes(x = time, y = y, group = time), 
                     side = "r",
                     fill = "grey30",
                     alpha = .3,
                     size = 0,
                     width = 2/7) + 
    geom_smooth(aes(linetype = it),
                method = "lm", 
                size = 2,
                position = position_dodge(width = .02)) +
    geom_point(size = 3, 
               position = position_dodge(width = .02),
               shape = 20) + 
    scale_color_manual(values = c("black", "azure2", "bisque")) + 
    #scale_color_manual(values = c("black", "darkolivegreen3", "coral2")) + 
    scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
    ylim(-1.3,.6) + 
    ylab("Observed Brain (a.u.)") +
    xlab("Years") + 
    ggtitle("2 observations")
  
  
  
  gs2 = ggplot(df %>% filter(type == "long"), 
               aes(x = time, y = y, group = it, color = it)) +
    geom_half_violin(data = df.norm %>% filter(type == "long"), 
                     mapping = aes(x = time, y = y, group = time), 
                     side = "r",
                     fill = "grey30",
                     alpha = .3,
                     size = 0, 
                     trim = F) + 
    geom_smooth(aes(linetype = it),
                method = "lm", 
                size = 2,
                position = position_dodge(width = .02),
                se = F) +
    geom_point(size = 3, 
               position = position_jitterdodge(jitter.width = .1, dodge.width = .0),
               shape = 20) + 
    scale_color_manual(values = c("black", "azure2", "bisque")) + 
    scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14, face = "bold"),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"))+
    ylim(-1.3,.6) + 
    ylab("") +
    xlab("Years") + 
    ggtitle("7 observations")
  
  
  newlabs = c("short" = "2 Obs", "long" = "7 Obs")
  gs3 = ggplot(df.mod, aes(x = it, y = estimate, color = it)) +
    geom_hline(yintercept = -.2, linetype = "dotdash", color = "grey50") + 
    geom_point(size = 8, shape = 20) + 
    facet_grid(cols = vars(type), labeller =labeller(.cols = newlabs)) + 
    scale_color_manual(values = c("black", "azure2", "bisque")) + 
    scale_linetype_manual(values=c("solid", "dashed", "dashed")) +
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 12),
          axis.title.x = element_blank(),
          plot.title = element_blank(),
          strip.text = element_text(size=12),
          plot.background = element_rect(color = "grey50",
                                         size = 1)) + 
    ylab("Obs. yearly brain change")
  
  fname = here("data_reliability_long/results/figs/scheme/observations.png")
  gs = gs1 + gs2 + inset_element(gs3, -0.15, 0, .33, .4)
  ggsave(file = fname, plot = gs)
  return(gs)
}

wrapper_figure_scheme = function() {
  library(here)
  library(ggplot2)
  options(bitmapType = "cairo")
  gs1 = figure_scheme_time()
  gs2 = figure_scheme_obs()
  gs = list(gs1, gs2)
  return(gs)
}


agegroupcomparison_wrapper = function(outdir, model, phenotypes, outreliability.comp) {
  library(tidyverse)
  library(broom)
  library(gamm4)
  load(file.path(outdir, "df.all.filt.Rda"))
  
  df.merge.pivot = 
    df.merge %>% 
    pivot_longer(-c(dataset, 
                    sub_id, 
                    rid, 
                    age, 
                    sex, 
                    site, 
                    sitenum,
                    sid), 
                 names_to = "features", 
                 values_to = "values")
  tmp = 
    df.merge.pivot %>% 
    filter(features == phenotypes & !is.na(values))
  
  
  mod <- gamm4(values ~ s(age, k = 20, bs = "cr") + sex, data = tmp, random = ~ (1|dataset) + (1|site) + (1|rid))
  rid.intercept = ranef(mod$mer)$rid
  grot.intercept = data.frame(rid = rownames(rid.intercept), 
                              meanZ = rid.intercept$`(Intercept)`)
  grot.intercept$meanZ = grot.intercept$meanZ
  
  
  #tmp$residuals =mod$gam$residuals %>% scale()
  tmp$residuals =mod$gam$residuals 
  df.merge.long = left_join(df.merge.long, tmp)
  
  mod.delta = 
    df.merge.long %>% 
    group_by(rid) %>% 
    mutate(xage = mean(age), 
           xtime = age - xage) %>% 
    dplyr::select(rid, xtime, residuals) %>% 
    filter(!is.na(residuals)) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ lm(residuals ~ xtime, data = .x)), 
           tidy = map(mod, ~broom::tidy(.x))) %>% 
    unnest(tidy) %>% 
    dplyr::select(rid, term, estimate, std.error) %>% 
    filter(term == "xtime") %>% 
    dplyr::select(-term) 
  names(mod.delta)[2] = "deltaZ"
  names(mod.delta)[3] = "seZ"
  
  # predict mean measure
  IndInt = 
    df.merge.long %>% 
    group_by(rid) %>% 
    mutate(xage = mean(age), 
           xtime = age - xage) %>% 
    summarise(age = first(xage), 
              sex = first(sex), 
              dataset  = first(dataset), 
              site = first(site))
  
  ## reference with respect to mean measure for individual!
  IndInt$Int = predict(mod$gam,IndInt)
  IndInt = 
    IndInt %>% 
    select(rid,Int)
  
  df.out = left_join(mod.delta, grot.intercept) %>% 
    left_join(.,IndInt)
  df.out = 
    df.out %>% 
    mutate(deltaZ = 100*deltaZ/Int, 
           seZ = 100*seZ/Int) %>% 
    select(-Int)
  names(df.out)[4] = "meanZ"
  
  
  
  subs = 
    df.merge.long%>% 
    group_by(rid)  %>% 
    mutate(time = max(age) - min(age)) %>% 
    summarise(n = n_distinct(age), 
              time = first(time), 
              dataset = first(dataset), 
              age = mean(age, na.rm = T))
  
  sss = subs %>% filter(time > 4 & n > 3, age < 60)
  var.Y = df.out %>% 
    filter(rid %in% sss$rid) %>% 
    .$deltaZ %>% 
    mad
  
  
  sss = subs %>% filter(time > 4 & n > 3, age >  60, age < 80)
  var.O = df.out %>% 
    filter(rid %in% sss$rid) %>% 
    .$deltaZ %>% 
    mad
  
  
  load(file.path(outreliability.comp, "icc.formula.rda"))
  df.agegroup = df.formula.icc %>% 
    filter(feature == phenotypes,
           model == "model1") 
  df.agegroup = expand_grid(df.agegroup, agegroup = c("young", "old")) %>% 
    mutate(seD = if_else(agegroup == "young", var.Y, var.O),
           bs = seD^2, 
           icc = bs /(ws + bs),
           sd_error = sqrt(ws), 
           mean1 = meanD + seD, 
           mean2 = meanD - seD)
  df.agegroup = compute_overlap(df.agegroup) %>% ungroup()
  save(df.agegroup, file = file.path(outreliability.comp, paste("icc.formula", phenotypes, "rda", sep = ".")))
}


compute_overlap = function(dat)  {
  custom_pdf1 <- function(x, mean1, std1) dnorm(x, mean = mean1, sd = std1)
  custom_pdf2 <- function(x, mean2, std2) dnorm(x, mean = mean2, sd = std2)
  
  result <- mapply(function(mean1, std1, mean2, std2) {
    integrate(function(x) sqrt(custom_pdf1(x, mean1, std1) * custom_pdf2(x, mean2, std2)), -Inf, Inf)$value
  }, dat$mean1, dat$sd_error, dat$mean2, dat$sd_error)
  # Bhattacharyya Coefficient
  dat$bc.main.v.decl <- result
  
  
  result <- mapply(function(mean1, std1, mean2, std2) {
    integrate(function(x) sqrt(custom_pdf1(x, mean1, std1) * custom_pdf2(x, mean2, std2)), -Inf, Inf)$value
  }, dat$mean1, dat$sd_error, dat$meanD, dat$sd_error)
  
  dat$bc.main.v.normal <- result
  
  dat = 
  dat %>% 
    mutate(p.main.v.decline = 1 - pnorm(mean1, mean = mean2, sd = sd_error),
           p.main.v.normal = 1 - pnorm(mean1, mean = meanD, sd = sd_error))
  return(dat)
}




pwr.r.test.wrapper = function(ES, alpha, sig.level, alternative) {
  n = 
    ceiling(pwr.r.test(r = ES, 
                       power = alpha, 
                       sig.level = sig.level, 
                       alternative = alternative)$n)
  return(n)
}

pwr.t.test.wrapper = function(ES, alpha, sig.level, alternative) {
  n = 
    ceiling(pwr.t.test(d = ES, 
                       power = alpha, 
                       sig.level = sig.level,
                       type = "two.sample",
                       alternative = alternative)$n)
  return(n)
}

pwr.anova.test.wrapper = function(ES, alpha, sig.level) {
  n = 
    ceiling(pwr.anova.test(k = 3, 
                           f = ES, 
                           power = alpha, 
                           sig.level = sig.level)$n)
  return(n)
}

sample_size_wrapper = function(df.formula.icc, df.statstest) {
  library(pwr)
  x =
    crossing(df.formula.icc, df.statstest) %>% 
    group_by(feature,
             model,
             t, 
             total_time) %>% 
    mutate(ObsES = sqrt(icc)* effectsize)
  
  x =
    x %>% 
    ungroup() %>% 
    rowwise() %>% 
    mutate(n = if_else(StatTest == "rpearson", 
                       pwr.r.test.wrapper(ObsES, power, 0.05, "two.sided"), 
                       if_else(StatTest == "two_sample",
                               pwr.t.test.wrapper(ObsES, power, 0.05, "two.sided"),
                               if_else(StatTest == "anova3l",
                                       pwr.anova.test.wrapper(ObsES, power, 0.05),NaN))))
  
  x = 
    x %>% group_by(feature,
                   model,
                   t, 
                   total_time) %>% 
    select(StatTest,
           effectsize,
           power, 
           n) %>% 
    nest()
  
  x = left_join(df.formula.icc, x)
  return(x)
}


fig_mean_icc = function(df.fig, oname) {

  pd=position_dodge(width = .4)
  
  dat = df.fig %>% 
    group_by(t,total_time) %>% 
    summarise(iccM = mean(icc, na.rm = T), iccSD = sd(icc, na.rm = T)) %>% 
    mutate(t = as.factor(t))
  
  gs = ggplot(dat, aes(x = total_time, y = iccM, group = t, color = t)) + 
    geom_errorbar(aes(ymin=iccM-iccSD, ymax=iccM+iccSD), colour="grey50", width=.4, position=pd, alpha = .4)  +
    geom_point(position = pd, size = 5, shape = 20) + 
    geom_smooth(position = pd, size = 2, se = F) + 
    scale_color_manual(name = "Obs.", values = c("coral2", "bisque", "azure2", "darkolivegreen3")) + 
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic()+ 
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=14, face = "bold"),
          legend.title = element_text(size=14, face = "bold")) +
    ylab("ICC") + 
    xlab("Follow-up time") + 
    ylim(0,1) + 
    guides(color = guide_legend(keywidth = 3, keyheight = 1))
  fname = here("data_reliability_long/results/figs/icc_main", paste(oname, "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 5, height = 7)
  return(gs)
}

fig_mean_icc_diff = function(df.fig, oname) {
  pd=position_dodge(width = .4)
  
  dat = df.fig %>% 
    group_by(t,total_time) %>% 
    summarise(iccM = mean(diff, na.rm = T), iccSD = sd(diff, na.rm = T)) %>% 
    mutate(t = as.factor(t))
  
  gs = ggplot(dat, aes(x = total_time, y = iccM, group = t, color = t)) + 
    geom_hline(yintercept = 0, linetype = 4, color = "grey50") + 
    geom_errorbar(aes(ymin=iccM-iccSD, ymax=iccM+iccSD), colour="grey50", width=.4, position=pd, alpha = .4)  +
    geom_point(position = pd, size = 5, shape = 20) + 
    geom_smooth(position = pd, size = 2, se = F) + 
    scale_color_manual(name = "Obs.", values = c("coral2", "bisque", "azure2", "darkolivegreen3")) + 
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic()+ 
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=14, face = "bold"),
          legend.title = element_text(size=14, face = "bold")) +
    ylab(expression(paste(Delta, " ICC"))) + 
    xlab("Follow-up time") + 
    ylim(-.5,.5) + 
    guides(color = guide_legend(keywidth = 3, keyheight = 1))
  fname = here("data_reliability_long/results/figs/icc_diffs", paste(oname, "png", sep = "."))
ggsave(file = fname, plot = gs, width = 5, height = 7)
}

fig_icc_regional = function(df,oname) {
  library(ggseg)
  library(viridis)
  meas = unique(df$modality)
  rICC = c(.1,.75)
  
  for (m in meas) {
    data = df %>% filter(modality == m)
    fname = here("data_reliability_long/results/figs/icc_main", paste("ggseg", oname, m, "png", sep = "."))
    
    if (!m == "svolume") {
      gs = ggplot(data) + 
        geom_brain(atlas = dk,
                   aes(fill = icc)) +
        theme_void() + 
        theme(legend.position = 'none') + 
        scale_fill_viridis(limits = rICC, oob = scales::squish)
      ggsave(fname, plot = gs)
    } else {
      gs = ggplot(data) + 
        geom_brain(atlas = aseg,
                   side = "coronal",
                   aes(fill = icc)) +
        theme_void() + 
        theme(legend.position = 'bottom',
              legend.text=element_text(size=22, face = "bold"),
              legend.title = element_text(size=26, face = "bold", vjust = .7),
              legend.key.size = unit(1.5,'cm')) + 
        scale_fill_viridis(name = "ICC", limits = rICC, breaks = c(.1, .3,.5, .7), oob = scales::squish) 
      ggsave(fname, plot = gs)
    }
  }
}

fig_icc_regional_difference = function(dat, altmodel) {
  df = dat
  library(ggseg)
  library(viridis)
  meas = unique(df$modality)
  rICC = c(-.2,.2)
  
  for (m in meas) {
    data = df %>% filter(modality == m)
    fname =  here("data_reliability_long/results/figs/icc_diffs", paste("ICC_rdiff", altmodel, m,  "png", sep = "."))
    
    if (!m == "svolume") {
      gs = ggplot(data) + 
        geom_brain(atlas = dk,
                   mapping = aes(fill = diff)) +
        theme_void() + 
        theme(legend.position = 'none') + 
        scale_fill_gradientn(limits = rICC, oob = scales::squish, colours = colorspace::diverge_hcl(15))
      ggsave(fname, plot = gs)
    } else {
      gs = ggplot(data) + 
        geom_brain(atlas = aseg,
                   side = "coronal",
                   aes(fill = diff)) +
        theme_void() + 
        theme(legend.position = 'bottom',
              legend.text=element_text(size=22, face = "bold"),
              legend.title = element_text(size=26, face = "bold", vjust = .7),
              legend.key.size = unit(1.5,'cm')) + 
        scale_fill_gradientn(limits = rICC, breaks = c(-.2,1, 0 ,1, .2), oob = scales::squish, colours = colorspace::diverge_hcl(15))
      ggsave(fname, plot = gs)
    }
  }
}

fig_icc_diffs = function(df.diff, oname) {
  df.diff$G = 1
  df.diff$modality <- factor(df.diff$modality, levels=c('area', 'thickness', 'volume', 'svolume'),
                             labels=c('area', 'thickness', 'volume', 'subcortical'))
  gs = 
    ggplot(df.diff, aes(G, diff, group = modality, fill = modality)) + 
    geom_violin(size = 1.5, alpha = .6, color = "grey50") + 
    geom_boxplot(color = "grey50", 
                 size = 1.5, 
                 fill = "white", 
                 alpha = 1,
                 width = .5, 
                 outlier.shape = NA) +
    scale_fill_viridis(option = "D", discrete = T) + 
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          legend.position = "none", 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 20, face = "bold"), 
          axis.line = element_line(size = 1.5, color = "grey50"), 
          axis.ticks = element_line(color = "grey50"), 
          axis.ticks.x = element_blank(), 
          axis.text = element_text(size = 15, "grey50"), 
          axis.text.x = element_blank(), 
          strip.text = element_text(size = 12.5, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50")) + 
    facet_wrap(modality~., nrow = 1, switch = "x") +
    labs(y = "\u0394ICC")
  fname = here("data_reliability_long/results/figs/icc_diffs", paste("ICC_diff", oname,  "png", sep = "."))
  ggsave(fname, plot = gs, width = 5, height = 7)  
}
 

fig_icc_diffs_sp = function(df.diff, oname) {
  library(ggpointdensity)
  library(viridis)
  if (oname == "model2") {
    df.diff$comparison = df.diff$model2
    ytitle = "ICC - Sim. Model2"
  } else if (oname == "empirical") {
    df.diff$comparison = df.diff$empirical
    ytitle = "ICC - Empirical Data"
  }
  
  gs = 
    ggplot(df.diff, aes(model1, comparison)) + 
    geom_pointdensity(adjust = 4) +
    geom_smooth(color = "grey50", size = 2) + 
    scale_color_viridis() + 
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          legend.position = "none", 
          axis.title = element_text(size = 20, face = "bold"),
          axis.line = element_line(size = 1.5, color = "grey50"), 
          axis.ticks = element_line(color = "grey50"), 
          axis.text = element_text(size = 15, "grey50"), 
          strip.text = element_text(size = 15, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50")) +
    labs(y = ytitle,
         x = "ICC - Main Sim. Model")
  
  fname = here("data_reliability_long/results/figs/icc_diffs", paste("ICC_scatterplot", oname,  "png", sep = "."))
  ggsave(fname, plot = gs, width = 5, height = 7)  
}

fig_icc_diffs_sp_set_design = function(df.diff, oname, time, tps) {
  if (oname == "model2") {
    df.diff$comparison = df.diff$model2
    ytitle = "ICC - Sim. Model2"
  } else if (oname == "empirical") {
    df.diff$comparison = df.diff$empirical
    ytitle = "ICC - Empirical Data"
  }
  df.diff$modality <- factor(df.diff$modality, levels=c('area', 'thickness', 'volume', 'svolume'),
                             labels=c('area', 'thickness', 'volume', 'subcortical'))
  
  dat = 
    df.diff %>% filter(total_time == time, t == tps)
  gs = 
    ggplot(dat, aes(model1, comparison, color = modality)) + 
    geom_smooth(method = "gam", color = "black", size = 2) +
    geom_point(size = 2) +
    scale_colour_viridis_d() +
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          legend.title = element_blank(),
          legend.position = "bottom", 
          legend.text = element_text(size=14),
          axis.title = element_text(size = 20, face = "bold"),
          axis.line = element_line(size = 1.5, color = "grey50"), 
          axis.ticks = element_line(color = "grey50"), 
          axis.text = element_text(size = 15, "grey50")) +
    labs(y = ytitle,
         x = "ICC - Main Sim. Model") +
    xlim(range(dat$model1, dat$comparison)) +
    ylim(range(dat$model1, dat$comparison))
  fname = here("data_reliability_long/results/figs/icc_diffs", paste("ICC_scatterplot_Set_design", oname, "fu", time, "tps",tps,  "png", sep = "."))
  ggsave(fname, plot = gs, width = 5, height = 7)  
  
  gs1 = 
    ggplot(dat, aes(model1, comparison, group = modality , color = modality)) + 
    geom_smooth(method = "gam", color = "black", size = 2) +
    geom_point(size = 2) +
    scale_colour_viridis_d() +
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          legend.title = element_blank(),
          legend.position = "bottom", 
          legend.text = element_text(size=14),
          axis.title = element_text(size = 20, face = "bold"),
          axis.line = element_line(size = 1.5, color = "grey50"), 
          axis.ticks = element_line(color = "grey50"), 
          axis.text = element_text(size = 15, "grey50"),
          strip.text = element_text(size = 12.5, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50")) +
    labs(y = ytitle,
         x = "ICC - Main Sim. Model") +
    facet_wrap(modality~., nrow = 1, switch = "x")+
    xlim(range(dat$model1, dat$comparison)) +
    ylim(range(dat$model1, dat$comparison))
  fname = here("data_reliability_long/results/figs/icc_diffs", paste("ICC_scatterplot_Set_design_facet",oname,  "fu", time, "tps",tps,  "png", sep = "."))
  ggsave(fname, plot = gs1, width = 10, height = 7)  
}


figure_icc_parameters = function(dat.me, dat.slopeV) {
  icc.me = 
    dat.me %>% 
    select(-c(feature, modality)) %>% 
    ICC(.)
  
  icc.slopeV = 
    dat.slopeV %>% 
    select(-c(feature, modality)) %>% 
    ICC(.)
  
  
  dat = 
    tribble(
      ~meas, ~icc.type, ~icc, ~icc_lb, ~icc_up,
      "slopeV", "ICC(2,1)", icc.slopeV$results$ICC[2], icc.slopeV$results$`lower bound`[2], icc.slopeV$results$`upper bound`[2],
      "slopeV", "ICC(2,k)", icc.slopeV$results$ICC[5], icc.slopeV$results$`lower bound`[5], icc.slopeV$results$`upper bound`[5],
      "me", "ICC(2,1)", icc.me$results$ICC[2], icc.me$results$`lower bound`[2], icc.me$results$`upper bound`[2],
      "me", "ICC(2,k)", icc.me$results$ICC[5], icc.me$results$`lower bound`[5], icc.me$results$`upper bound`[5])
  
  gs = 
    ggplot(dat, aes(x = meas, y = icc, group = icc.type, fill = icc.type)) +
    geom_errorbar(dat, mapping = aes(ymin = icc_lb, ymax = icc_up), width = 0.2, position=position_dodge(width = .5), linewidth = 2, color= "black") +
    geom_point(position=position_dodge(width = .5), shape = 21, size = 6, color = "black", alpha = .6) +
    scale_fill_manual(name = "", values = c("coral2", "darkolivegreen3")) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = .2)) + 
    scale_x_discrete(labels = c("error (%)", "slope (%)")) +
    theme_classic() + 
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=14, face = "bold"),
          legend.title = element_text(size=14, face = "bold")) +
    ylab("ICC") + 
    xlab("Parameter") +
    guides(color = guide_legend(keywidth = 3, keyheight = 1))
  fname = here("data_reliability_long/results/figs/model_parameter_comparison", paste("icc", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 4, height = 7)
  return(gs)
}


fig_pairwise.icc.parameters = function(dat, oname) {
  library(GGally)
  library(viridis)
  x = dat %>% 
    select(-c(feature, modality))
  
  
  my_fn <- function(data, mapping,  ...){
    
    # grab data
    x <- eval_data_col(data, mapping$x)
    y <- eval_data_col(data, mapping$y)
    
    # calculate correlation
    corr <- ICC(cbind(x, y))$results$ICC[2]
    
    fill <- viridis(100, 
                    direction = 1)[findInterval(corr, 
                                                seq(0, 
                                                    1, 
                                                    length=100))]
    
    ggally_text(
      label = as.character(round(corr, 2)),
      mapping = aes(),
      xP = 0.5, yP = 0.5,
      color = 'black',
      alpha = .5,
      size = 8, 
      fontface='bold',
      ...
    ) +
      theme_void() +
      theme(panel.background = element_rect(fill=fill, color = NA),
            panel.grid.major = element_blank())
    
  }
  
  my_dens <- function(data, mapping) {   
    ggplot(data = data, mapping = mapping) +
      geom_density(alpha = 0.5, color = "black", fill = "maroon", alpha = .4) +
      theme_classic()
  }
  
  my_fn2 <- function(data, mapping, method="lm", ...){
    p <- ggplot(data = data, mapping = mapping) + 
      geom_smooth(method=method, color ="black", size = 1, ...) +
      geom_point(color = "black",fill = "maroon", size = 2, shape = 21) + 
      theme_classic()
    p
  }
  
  
  gs =  ggpairs(
    x,
    lower = list(continuous = my_fn2),
    diag = list(continuous = my_dens),
    upper = list(continuous = my_fn)
  )
  
  
  fname = here("data_reliability_long/results/figs/model_parameter_comparison", 
               paste("icc.pairwise",oname, "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 8, height = 7)
  return(gs)
}

figure_icc_parameters_bymodality= function(dat.me, dat.slopeV) {
  
  icc.me = 
    dat.me %>% 
    select(-feature) %>% 
    group_by(modality) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ICC(.)),
           icc.21_M = map_dbl(mod, ~.$results$ICC[2]),
           icc.2_M = map_dbl(mod, ~.$results$ICC[5]),
           icc.21_lb = map_dbl(mod, ~.$results$`lower bound`[2]),
           icc.2_lb = map_dbl(mod, ~.$results$`lower bound`[5]),
           icc.21_ub = map_dbl(mod, ~.$results$`upper bound`[2]),
           icc.2_ub = map_dbl(mod, ~.$results$`upper bound`[5])) %>% 
    select(-c(mod, data)) %>% 
    pivot_longer(-modality, names_to = "names", values_to = "values") %>% 
    separate(names, c("icc.type", "grot.2"), sep = "_") %>% 
    pivot_wider(names_from = grot.2, 
                values_from = values)  %>% 
    mutate(meas = "me")
  
  
  icc.slopeV = 
    dat.slopeV %>% 
    select(-feature) %>% 
    group_by(modality) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ICC(.)),
           icc.21_M = map_dbl(mod, ~.$results$ICC[2]),
           icc.2_M = map_dbl(mod, ~.$results$ICC[5]),
           icc.21_lb = map_dbl(mod, ~.$results$`lower bound`[2]),
           icc.2_lb = map_dbl(mod, ~.$results$`lower bound`[5]),
           icc.21_ub = map_dbl(mod, ~.$results$`upper bound`[2]),
           icc.2_ub = map_dbl(mod, ~.$results$`upper bound`[5])) %>% 
    select(-c(mod, data)) %>% 
    pivot_longer(-modality, names_to = "names", values_to = "values") %>% 
    separate(names, c("icc.type", "grot.2"), sep = "_") %>% 
    pivot_wider(names_from = grot.2, 
                values_from = values) %>% 
    mutate(meas = "slope")
  
  
  dat = rbind(icc.me, icc.slopeV)
  dat$modality <- factor(dat$modality, 
                         levels=c('area', 'thickness', 'volume', 'svolume'),
                         labels=c('area', 'thickness', 'volume', 'subcortical'))
  
  gs = 
    ggplot(dat, aes(x = meas, y = M, group = icc.type, fill = icc.type)) +
    geom_errorbar(dat, mapping = aes(ymin = lb, ymax = ub), width = 0.2, position=position_dodge(width = .5), linewidth = 2, color= "black") +
    geom_point(position=position_dodge(width = .5), shape = 21, size = 6, color = "black", alpha = .6) +
    scale_fill_manual(name = "", labels = c("ICC(2,k)", "ICC(2,1)"), values = c("darkolivegreen3", "coral2")) + 
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1, by = .2)) + 
    scale_x_discrete(labels = c("error (%)", "slope (%)")) +
    theme_classic() + 
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=14, face = "bold"),
          legend.title = element_text(size=14, face = "bold")) +
    ylab("ICC") + 
    xlab("Parameter") +
    facet_wrap(modality~., nrow = 1) +
    guides(color = guide_legend(keywidth = 3, keyheight = 1))
  fname = here("data_reliability_long/results/figs/model_parameter_comparison", paste("icc_bymodality", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 10, height = 7)
  return(gs)
}

figure_icc_across_datasets = function(df.icc.reliability.icc2, df.icc.reliability.pairwise) {
  
  grot1 = df.icc.reliability.icc2 %>% 
    summarise(iccM = mean(icc), 
              sd = sd(icc),
              icc.type = "ICC(2,k)")
  
  grot2 = 
    df.icc.reliability.pairwise %>% 
    summarise(iccM = mean(icc), 
              sd = sd(icc),
              icc.type = "ICC(2,1)")
  
  dat = rbind(grot1, 
              grot2)
  
  
  gs = 
    ggplot(dat, aes(x = icc.type, y = iccM, group = icc.type, fill = icc.type)) +
    geom_errorbar(dat, mapping = aes(ymin = iccM - sd, ymax = iccM + sd), width = 0.2, linewidth = 2, color= "black") +
    geom_point(shape = 21, size = 6, color = "black", alpha = .6) +
    scale_fill_manual(name = "", values = c("coral2", "darkolivegreen3")) + 
    scale_y_continuous(name = "ICC", limits = c(0,1), breaks = seq(0,1, by = .2)) + 
    #scale_x_discrete(labels = c("error (%)", "slope (%)")) +
    theme_classic() + 
    theme(legend.position = 'none',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_blank(),
          legend.text=element_text(size=14, face = "bold"),
          legend.title = element_text(size=14, face = "bold"))
  fname = here("data_reliability_long/results/figs/icc_diffs", paste("icc_wholemodel_bw_datasets", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 3, height = 7)
  
  
  
  grot1 = 
    df.icc.reliability.icc2 %>% 
    unnest(data) %>% 
    select(t, total_time, icc.sub) %>% 
    mutate(icc.type = "ICC(2,k)")
  
  
  grot2 = 
    df.icc.reliability.pairwise %>% 
    unnest(data) %>% 
    select(t, total_time, icc.sub) %>% 
    mutate(icc.type = "ICC(2,1)") 
  
  dat = 
    rbind(grot1, 
          grot2)
  
  dat = 
    dat %>% 
    group_by(t, total_time, icc.type) %>% 
    summarise(icc = mean(icc.sub), sd.icc = sd(icc.sub))
  
  
  pd=position_dodge(width = .4)
  gs = 
    ggplot(dat %>% mutate(t = as.factor(t)), aes(x = total_time, y = icc, group = t, color = t)) + 
    geom_errorbar(aes(ymin=icc-sd.icc, ymax=icc+sd.icc), colour="grey50", width=.4, position=pd, alpha = .4)  +
    geom_point(position = pd, size = 5, shape = 20) + 
    geom_line(position = pd, size = 2) + 
    scale_color_manual(name = "Obs.", values = c("coral2", "bisque", "azure2", "darkolivegreen3")) + 
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic()+ 
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=14, face = "bold"),
          legend.title = element_text(size=14, face = "bold")) +
    ylab("ICC") + 
    xlab("Follow-up time") + 
    ylim(0,1) + 
    facet_wrap(icc.type~., nrow = 1) +
    guides(color = guide_legend(keywidth = 3, keyheight = 1))
  fname = here("data_reliability_long/results/figs/icc_diffs", paste("icc_byTFU_bw_datasets", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 5, height = 7)
}

# deprecated
fig_model_parameter_comparison = function(dat) {
  dat$group <- factor(dat$group, 
                      levels=c('area', 'thickness', 'volume', 'svolume'),
                      labels=c('area', 'thickness', 'volume', 'subcortical'))
  
  gs1 = ggplot(dat, aes(m1, m2, color = group)) +
    geom_point() + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_smooth(method = "lm", color = "blue") + 
    theme_classic() + 
    theme(legend.title = element_blank(),
          legend.position = "bottom", 
          axis.title = element_text(size = 16),
          legend.text = element_text(size=16)) +
    ylab("Replication Model") +
    xlab("Main Model")
  
  
  
  # gs2 = ggplot(dat, aes(m1, m2, color = group)) +
  #   geom_point() + 
  #   geom_abline(intercept = 0, slope = 1) + 
  #   geom_smooth(method = "lm", color = "blue") + 
  #   theme_classic() + 
  #   theme(legend.title = element_blank(),
  #         legend.position = "none", 
  #         axis.title = element_blank()) + 
  #   facet_grid(cols = vars(group))
  # 
  # gs2
  
  return = list(gs1)
}


get_simulated_results = function(outreliability.comp, modelname) {
  
  fname = file.path(outreliability.comp, 
                    paste("df", 
                          modelname, 
                          "rda", 
                          sep = "."))
  load(fname)
  dat = 
    df.comp$icc %>% select(-icc2k) %>% 
    mutate(model = modelname)
  
  # dat2 = 
  #   df.comp$simres.summary %>% 
  #   pivot_wider(names_from = "name", values_from = "value")
  
   dat2 = 
     df.comp$simres.summary
  
  dat = inner_join(dat,dat2)
  return(dat)
}

fig_corresponce_ws = function(dat) {
  gs1 = ggplot(dat, aes(beta_rmse^2, ws, color = tps)) + 
    geom_point() +
    geom_abline(intercept = 0, slope = 1) + 
    theme_classic() + 
    theme(legend.title = element_blank(),
          legend.position = "none", 
          axis.title = element_text(size = 16),
          legend.text = element_text(size=16)) +
    ylab("Var " ~ Delta) +
    xlab(bquote("RMSE " ~ Delta^"2"))
  
  gs2 =
    ggplot(dat, aes(avg_beta_se, beta_rmse)) +
    geom_point() + 
    geom_smooth(method = "lm", aes(color = tps), se = T) + 
    geom_abline(intercept = 0, slope = 1) + 
    facet_wrap(vars(t), nrow = 1) + 
    theme_classic() + 
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16),
          plot.title = element_blank(),
          strip.text = element_text(size=16),
          plot.background = element_blank()) + 
    xlab("SD " ~ Delta) +
    ylab(bquote("RMSE " ~ Delta))
  
  gs = gs1 + gs2 + plot_layout(widths = c(1,4))
  fname = here("data_reliability_long/results/figs/measure_equivalences/ws-v-rmseveta-v-betase.png")
  ggsave(file = fname, plot = gs,  width = 20, height = 20/5)
  return(gs)
}


fig_compare_icc_form_sim = function(dat.formula, dat.simulation) {
  
  common_vars = c("feature", "tps","total_time","modality", "model")
  
  dat1 = dat.formula %>% 
    rename("tps" = "t",
           "icc_formula" = "icc") %>% 
    mutate(tps = as.factor(tps)) %>% 
    select(c(all_of(common_vars), icc_formula))
  
  dat2 = 
    dat.simulation %>%  rename("icc_simulations" = "icc21") %>% 
    select(c(all_of(common_vars), icc_simulations))
  
  dat = inner_join(dat2,dat1)
  
  gs = ggplot(dat, aes(icc_formula, icc_simulations, color = model)) + 
    geom_point() + 
    geom_smooth(method = "lm", color = "black") + 
    #facet_wrap(vars(model)) +
    theme_classic() + 
    theme(legend.position = 'none',
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text= element_text(size = 16),
          axis.title = element_text(size = 16),
          plot.title = element_blank(),
          strip.text = element_text(size=16),
          plot.background = element_blank()) + 
    xlab("ICC Estimated") +
    ylab("ICC Simulated Dataset")
  
  fname = here("data_reliability_long/results/figs/measure_equivalences/icc.formula-v-simulation.png")
  ggsave(file = fname, plot = gs,  width = 6, height = 5)
  return(gs)
}


fig_individual_overlap = function(dat, pheno) {
  fname = here("data_reliability_long/results/figs/individual_overlap", 
               paste(pheno, "png", sep = "."))
  dat = 
    dat %>% 
    pivot_longer(c("mean1", "meanD", "mean2"),
                 names_to = "distr") %>% 
    mutate(mod = map2(value, sd_error, ~ rnorm(5000, .x,.y))) %>% 
    unnest(mod) %>% 
    mutate(distr = factor(distr, 
                          levels = c("mean1", "mean2", "meanD"), 
                          labels = c("Maintainer", "Decliner", "Normal")))
  
  
  gs= ggplot(dat, aes(x = mod, fill = distr)) +
    geom_density(alpha = 0.5) +
    labs(x = "% Annual Change", y = "") +
    scale_fill_manual(values = c("darkgreen", "darkred", "lightgray")) + 
    facet_grid(vars(total_time), vars(t)) + 
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          axis.title = element_text(size = 18, "grey50"), 
          axis.line = element_line(size = 1.5, color = "grey50"), 
          axis.ticks = element_line(color = "grey50"), 
          axis.text = element_text(size = 15, "grey50"), 
          strip.text = element_text(size = 15, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size=15))
  ggsave(file = fname, plot = gs,  width = 12, height = 10)
}


fig_bc_overlap = function() {
  
  
  load(file.path(outreliability.comp, "icc.formula.rda"))
  
  X = 
    df.formula.icc %>% 
    filter(model == "model1") %>% 
    group_by(total_time, t) %>% 
    summarise(Xbc.main.v.decl = mean(bc.main.v.decl, na.rm = T),
              SDbc.main.v.decl = sd(bc.main.v.decl, na.rm = T),
              Xbc.main.v.normal = mean(bc.main.v.normal, na.rm = T),
              SDbc.main.v.normal = sd(bc.main.v.normal, na.rm = T),
              Xp.main.v.decl = mean(p.main.v.decline, na.rm = T), 
              SDp.main.v.decl = sd(p.main.v.decline, na.rm = T),
              Xp.main.v.normal = mean(p.main.v.normal, na.rm = T), 
              SDp.main.v.normal = sd(p.main.v.normal, na.rm = T))
  X = 
    X %>% 
    pivot_longer(-c(total_time, t),
                 names_to = "measure", 
                 values_to = "value") %>% 
    separate(measure, 
             c("stat", "measure"), 
             sep="\\.",
             extra = "merge") %>% 
    pivot_wider(names_from = stat,
                values_from = value)
  
  
  pd=position_dodge(width = .4)
  X = X %>% mutate(t = as.factor(t))
  
  
  gs = ggplot(X, aes(x = total_time, y = Xbc, group = interaction(t, measure), color = t, linetype = measure)) +   
    geom_errorbar(aes(ymin=Xbc-SDbc, ymax=Xbc+SDbc), colour="grey50", width=.4, position=pd, alpha = .4)  +
    geom_point(position = pd, size = 5, shape = 20) + 
    geom_line(position = pd, size = 2) + 
    scale_color_manual(name = "Obs.", values = c("coral2", "bisque", "azure2", "darkolivegreen3")) + 
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic()+ 
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=12, face = "bold"),
          legend.title = element_text(size=12, face = "bold")) +
    ylab("Bhattacharyya Coefficient") + 
    xlab("Follow-up time") + 
    ylim(0,1) + 
    guides(color = guide_legend(keywidth = 3, keyheight = 1,nrow = 2, byrow =T),
           linetype = guide_legend(keywidth = 2.5, keyheight = 1,nrow = 2, byrow =T, title = "", labels =c("maintainer - decliner", "normal - decliner")))
  fname = here("data_reliability_long/results/figs/individual_overlap", 
               paste("Bhattacharyya", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 5, height = 7)
  
  
  
  gs = ggplot(X, aes(x = total_time, y = Xp, group = interaction(t, measure), color = t, linetype = measure)) +   
    geom_errorbar(aes(ymin=Xp-SDp, ymax=Xp+SDp), colour="grey50", width=.4, position=pd, alpha = .4)  +
    geom_point(position = pd, size = 5, shape = 20) + 
    geom_line(position = pd, size = 2) + 
    scale_color_manual(name = "Obs.", values = c("coral2", "bisque", "azure2", "darkolivegreen3")) + 
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic()+ 
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=12, face = "bold"),
          legend.title = element_text(size=12, face = "bold")) +
    ylab("probability missclassification") + 
    xlab("Follow-up time") + 
    ylim(0,.5) + 
    guides(color = guide_legend(keywidth = 3, keyheight = 1,nrow = 2, byrow =T),
           linetype = guide_legend(keywidth = 2.5, keyheight = 1,nrow = 2, byrow =T, title = "", labels =c("maintainer - decliner", "normal - decliner")))
  fname = here("data_reliability_long/results/figs/individual_overlap", 
               paste("p", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 5, height = 7)
  
}

fig_maintainer = function(model1, pheno) {
  
  indtrajdir = file.path(outreliability.comp, "individual_trajectories")
  load(file.path(indtrajdir,paste0("model1_", pheno, ".Rda")))
  
  fname = here("data_reliability_long/results/figs/maintenance", 
               paste(pheno, "png", sep = "."))
  
  densities <- df.extreme %>% 
    filter(tps %in% c(3,7), total_time %in% c(2,6,10)) %>%
    group_by(tps, total_time) %>%
    do(., dens = density(.$true_beta))
  
  
  grot = list()
  for (i in 1:length(densities$dens)) {
    grot[[i]] = data.frame(x = densities$dens[[i]]$x, 
                           y = densities$dens[[i]]$y, 
                           tps = densities$tps[i],
                           total_time = densities$total_time[i])
  }
  
  df.plot = data.table::rbindlist(grot) %>% 
    mutate(feature = pheno)
  df.plot = left_join(df.plot, model1)
  
  
  df.annot = P_extreme %>% 
    filter(tps %in% c(3,7), total_time %in% c(2,6,10)) %>% 
    mutate(an = as.character(paste0(round(100*pct_extrem,1), "%")))
  
  
  M = -abs(df.plot$meanD %>% unique())
  gs = ggplot(df.plot, aes(x, y, color = "black")) +
    geom_line(size = 1) + 
    theme_classic() +
    geom_area(data = df.plot %>% filter(x >0), fill = "darkgreen", alpha = .5) + 
    geom_area(data = df.plot %>% filter(x < M), fill = "darkred", alpha = .5) + 
    geom_hline(yintercept = 0) + 
    labs(x = "% True Annual Change", y = "") +
    facet_grid(vars(total_time),vars(tps)) +
    #scale_color_manual(values = c("darkred", "lightgray", "darkgreen")) + 
    facet_grid(vars(total_time), vars(tps)) + 
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          axis.title = element_text(size = 18, "grey50"), 
          axis.line = element_line(size = 1.5, color = "grey50"), 
          axis.ticks = element_line(color = "grey50"), 
          axis.text = element_text(size = 15, "grey50"), 
          strip.text = element_text(size = 15, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50"),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size=15)) + 
    geom_text(data = df.annot, aes(x = Inf , y = Inf, label = an), vjust = 1, hjust = 1, color = "grey50", size = 6)
  
  ggsave(file = fname, plot = gs, width = 7, height = 7)
  #1- pnorm(0, -3.753, 2.15)
}

fig_icc_globalvars = function() {
  outreliability.comp="/ess/p274/cluster/projects/p039_image_brain_change/data_reliability_long/df_mri/all/computational"
  tps = c(3,7)
  
  #  load and prepare global vars
  load(file.path(outreliability.comp, "icc.formula.global.rda"))
  df.formula.icc.global = df.formula.icc %>% 
    filter(t %in% tps,
           model == "model1") %>% 
    mutate(t = as.factor(t))
  df.formula.icc.global$modality <- factor(df.formula.icc.global$modality, levels=c('area', 'thickness', 'volume', 'svolume'),
                                           labels=c('area', 'thickness', 'volume', 'subcortical'))
  
  
  
  load(file.path(outreliability.comp, "icc.formula.rda"))
  dat = df.formula.icc %>% 
    filter(model == "model1",
           t %in% tps) %>% 
    group_by(t,total_time, modality) %>% 
    summarise(iccSD = sd(icc, na.rm = T),
              icc = mean(icc, na.rm = T)) %>% 
    mutate(t = as.factor(t), 
           iccmin = icc - iccSD,
           iccmax = icc + iccSD, 
           iccmax = if_else(iccmax >1,1,iccmax), 
           iccmin = if_else(iccmin < 0, 0, iccmin)) 
  
  dat$modality <- factor(dat$modality, levels=c('area', 'thickness', 'volume', 'svolume'),
                         labels=c('area', 'thickness', 'volume', 'subcortical'))
  
  df.formula.icc.global$feature <- factor(df.formula.icc.global$feature, 
                                          levels=c("CorticalArea","SupraTentorialVolume","MeanThickness", "CorticalVolume","SubCorticalVolume"),
                          labels = c("Cortical Area", "Supratentorial Volume", "Mean Thickness", "Cortical Volume", "Subcortical Volume"))
  
  
  
  smooth.ribbon = data.frame(total_time = seq(2,12, by = 0.1))
  x = 
    dat %>% 
    group_by(t, modality) %>% 
    nest() %>% 
    mutate(mod.min = map(data, ~loess(iccmin ~ total_time, data = .)), 
           mod.max = map(data, ~loess(iccmax ~ total_time, data = .)),
           mod.mean = map(data, ~loess(icc ~ total_time, data = .)),
           iccmin = map(mod.min, ~predict(., newdata = smooth.ribbon)), 
           iccmax = map(mod.max, ~predict(., newdata = smooth.ribbon)), 
           icc = map(mod.mean, ~predict(., newdata = smooth.ribbon)))
  
  
  xx = x %>% select(-c(data, mod.min, mod.max, mod.mean)) %>% unnest(iccmin, iccmax, icc)
  xx$total_time = rep(smooth.ribbon$total_time, x$t %>% length()) 
  gs = 
    ggplot(dat, aes(x = total_time, y = icc, group = interaction(t,modality))) + 
    geom_ribbon(data = xx, mapping = aes(ymin=iccmin, ymax=iccmax), colour="grey50", linewidth = 0, alpha = .3)  +
    geom_smooth(color = "black", size = .5, se = F) + 
    geom_smooth(df.formula.icc.global, mapping = aes(x = total_time, y = icc, group = interaction(t,feature), color = feature ), se = F, size = 1.2) + 
    facet_grid(modality ~t) +
    scale_color_manual(name = "Global Feature", values = c("coral2", "bisque", "azure2", "darkolivegreen3", "blueviolet")) + 
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic()+ 
    theme(legend.position = 'bottom',
          strip.text = element_text(size = 14, color = "grey50", face = "bold"),
          strip.background = element_rect(linewidth =  1.5, color = "grey50"),
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=12, face = "bold"),
          legend.title = element_text(size=12, face = "bold")) +
    ylab("ICC") + 
    xlab("Follow-up time") + 
    ylim(0,1) + 
    guides(color = guide_legend(keywidth = 2.5, keyheight = 1, nrow = 3, byrow =T))
  fname = here("data_reliability_long/results/figs/icc_globalVars", paste("icc_globalvars", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 7, height = 12)
  return(gs)
}


fig_parameters_globalvars = function() {
  outreliability.comp="/ess/p274/cluster/projects/p039_image_brain_change/data_reliability_long/df_mri/all/computational"
  
  #  load and prepare global vars
  load(file.path(outreliability.comp, "icc.formula.global.rda"))
  df.formula.icc.global = df.formula.icc %>%
    filter(model == "model1") %>% 
    group_by(feature, modality) %>% 
    summarise(seD = first(seD), 
              meanE = first(meanE)) %>% 
    mutate(G = as.factor("1"))
  
  df.formula.icc.global$modality <- factor(df.formula.icc.global$modality, levels=c('area', 'thickness', 'volume', 'svolume'),
                                           labels=c('area', 'thickness', 'volume', 'subcortical'))
  
  
  load(file.path(outreliability.comp, "icc.formula.rda"))
  dat = df.formula.icc %>% 
    filter(model == "model1") %>% 
    group_by(feature, modality)  %>% 
    summarise(seD = first(seD), 
              meanE = first(meanE)) %>% 
    mutate(G = as.factor("1"))
  dat$modality <- factor(dat$modality, levels=c('area', 'thickness', 'volume', 'svolume'),
                         labels=c('area', 'thickness', 'volume', 'subcortical'))
  
  gs1 = 
    ggplot(dat, aes(x = G, y = seD, group = modality, fill = modality )) +
    geom_violin(size = 1.5, alpha = .6, color = "grey50") + 
    geom_boxplot(color = "grey50", 
                 size = 1.5, 
                 fill = "white", 
                 alpha = 1,
                 width = .5, 
                 outlier.shape = NA) +
    geom_point(data = df.formula.icc.global, aes(x = G,y = seD, color = feature), size = 4 ) + 
    scale_fill_viridis(option = "D", discrete = T) + 
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          legend.position = "bottom", 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 16, face = "bold"), 
          axis.line = element_line(size = 1.5, color = "grey50"), 
          axis.ticks = element_line(color = "grey50"), 
          axis.ticks.x = element_blank(), 
          axis.text = element_text(size = 15, "grey50"), 
          axis.text.x = element_blank(), 
          strip.text = element_text(size = 12.5, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50"),
          legend.text=element_text(size=12, face = "bold"),
          legend.title = element_text(size=12, face = "bold")) + 
    facet_wrap(modality~., nrow = 1, switch = "x") + 
    scale_color_manual(name = "Global Feature", values = c("coral2", "bisque", "azure2", "darkolivegreen3", "blueviolet")) + 
    guides(color = guide_legend(keywidth = 2.5, keyheight = 1, nrow = 3, byrow =T)) + 
    guides(fill = 'none') + 
    ylab("standard deviation of the slopes (%)")
  
  
  gs2 = 
    ggplot(dat, aes(x = G, y = meanE, group = modality, fill = modality )) +
    geom_violin(size = 1.5, alpha = .6, color = "grey50") + 
    geom_boxplot(color = "grey50", 
                 size = 1.5, 
                 fill = "white", 
                 alpha = 1,
                 width = .5, 
                 outlier.shape = NA) +
    geom_point(data = df.formula.icc.global, aes(x = G,y = meanE, color = feature), size = 4 ) + 
    scale_fill_viridis(option = "D", discrete = T) + 
    theme_classic() + 
    theme(panel.background = element_rect(fill = "white"), 
          legend.position = "bottom", 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(size = 16, face = "bold"), 
          axis.line = element_line(size = 1.5, color = "grey50"), 
          axis.ticks = element_line(color = "grey50"), 
          axis.ticks.x = element_blank(), 
          axis.text = element_text(size = 15, "grey50"), 
          axis.text.x = element_blank(), 
          strip.text = element_text(size = 12.5, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50"),
          legend.text=element_text(size=12, face = "bold"),
          legend.title = element_text(size=12, face = "bold")) + 
    facet_wrap(modality~., nrow = 1, switch = "x") + 
    scale_color_manual(name = "Global Feature", values = c("coral2", "bisque", "azure2", "darkolivegreen3", "blueviolet")) + 
    guides(color = guide_legend(keywidth = 2.5, keyheight = 1, nrow = 3, byrow =T)) + 
    guides(fill = 'none') + 
    ylab(" Coefficient of variation (%)")
  
  gs = gs1 + gs2
  fname = here("data_reliability_long/results/figs/icc_globalVars", paste("parameters", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 8, height = 12)
  return(gs)
}

fig_agegroupcomparison = function(df.formula.icc, phenotypes, mdl) {
  pd=position_dodge(width = .4)
  
  dat = df.formula.icc %>% 
    filter(feature == phenotypes, model %in% mdl ) %>% 
    mutate(t = as.factor(t))
  
  gs = ggplot(dat %>% filter(t %in% c(3,7)), 
              aes(x = total_time, y = icc, group = model, color = model)) + 
    geom_point(size = 5, shape = 20) + 
    geom_smooth(size = 2,se = F) + 
    scale_color_manual(name = "Group", labels = c("Old", "Young"), values = c("coral2", "azure2")) + 
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic()+ 
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=14, face = "bold"),
          legend.title = element_text(size=14, face = "bold"),
          strip.text = element_text(size = 14, color = "grey50", face = "bold"),
          strip.background = element_rect(linewidth =  1.5, color = "grey50")) +
    ylab("ICC") + 
    xlab("Follow-up time") + 
    ylim(0,1) + 
    facet_wrap(t~., nrow = 1) + 
    guides(color = guide_legend(keywidth = 3, keyheight = 1))
  
  fname = here("data_reliability_long/results/figs/icc_agecomparison", paste(phenotypes, "model1", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 5, height = 5)
  return(gs)
}

fig_sample_size_mean = function(dat, oname) {
  library(pwr)
  r1 =  ceiling(pwr.r.test(r = .1, 
                           power = .8, 
                           sig.level = .05)$n)
  r2 =  ceiling(pwr.r.test(r = .3, 
                           power = .8, 
                           sig.level = .05)$n)
  r3 =  ceiling(pwr.r.test(r = .5, 
                           power = .8, 
                           sig.level = .05)$n)
  
  dat2 = tribble(
    ~effectsize, ~meaN,
    .1,   r1,
    .3,   r2,
    0.5,   r3
  )
  
  dat = dat %>% filter(model == oname) %>% 
    unnest(cols = c(data)) %>% 
    filter(StatTest == "rpearson", power == .8) %>% 
    group_by(total_time, 
             t, 
             model,
             effectsize) %>% 
    summarise(meanN = mean(n, na.rm = T), 
              sdN = sd(n, na.rm = T)) %>% 
    mutate(t = as.factor(t),
           effectsize = as.factor(effectsize))
  
  
  pd=position_dodge(width = .4)
  ybreaks = round(scales:::trans_breaks("sqrt", function(x) x^2)(c(60, 6000), n = 10),-1)
  gs = 
    ggplot(dat, aes(x = total_time, y = meanN, group = interaction(effectsize,t), color = effectsize, linetype = t)) + 
    geom_hline(data = dat2, mapping = aes(yintercept = meaN, group = effectsize), color = "grey50",linetype = 4, linewidth = 1) + 
    geom_errorbar(aes(ymin=meanN-sdN, ymax=meanN+sdN), colour="grey50", width=.4, position=pd, alpha = .4)  +
    geom_point(position = pd, size = 3, shape = 20) + 
    geom_smooth(se = F,position = pd, size = 1.2) + 
    scale_y_continuous(trans='sqrt', breaks = ybreaks)  +
    scale_color_manual(name = "r", values = c("coral2", "azure2", "darkolivegreen3")) + 
    scale_linetype(name = "Obs.") +
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic() +
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=12, face = "bold"),
          legend.title = element_text(size=12, face = "bold")) +
    ylab("Sample Size") + 
    xlab("Follow-up time") + 
    guides(color = guide_legend(keywidth = 2.5, keyheight = 1,nrow = 2, byrow =T),
           linetype = guide_legend(keywidth = 2.5, keyheight = 1,nrow = 2, byrow =T))
  fname = here("data_reliability_long/results/figs/samplesize", paste("mean", oname, "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 6, height = 8)
  return(gs)
}


fig_sample_size_region_tps = function(dat, ft, model, tps = c(3,7)) {
  
  library(metR)
  library(pwr)
  
  p =
    dat %>% 
    ungroup() %>% 
    dplyr::filter(feature == ft & model == model) %>% 
    summarize(feature = first(feature), 
              modality = first(modality), 
              seD = first(seD), 
              meanE = first(meanE), 
              model = first(model))
  
  p = expand.grid(feature = unique(p$feature), 
                  t = seq(3,9,by = 2), 
                  total_time = seq(2,12, by =.1),
                  r = seq(0.05, 0.5, by = 0.01), 
                  alpha = .8) %>% 
    left_join(p, .) %>% 
    mutate(ws = meanE^2*((total_time^2*t*(t+1))/ (12*(t-1)))^-1, 
           bs = seD^2, 
           icc = bs /(ws + bs),
           obsr = r*sqrt(icc)) %>% 
    rowwise() %>% 
    mutate(n = pwr.r.test.wrapper(obsr,alpha, .05, "two.sided"))
  
  
  fillbreaks =  c(50, 
                  100, 
                  200,
                  300,
                  400,
                  500, 
                  750, 
                  1000,
                  1500, 
                  2000, 
                  2500,
                  3000, 
                  4000, 
                  5000,
                  6000, 
                  8000, 
                  10000)
  
  
  dat = p %>% filter(t %in% tps)
  dat = 
    dat %>% 
    filter(total_time <= 8, 
           r <= .4, 
           total_time >= 2, 
           r >=.05)
  
  gs = 
    ggplot(dat, aes(x = total_time, y = r, z = n, fill = n)) +
    geom_raster(interpolate = T) + #interpolate for success 
    geom_contour(aes(z = n), color = "white", linetype = 4, size = 1, breaks = fillbreaks) + 
    geom_text_contour(aes(z = n), 
                      color = "white", 
                      breaks = fillbreaks,
                      size = 4, 
                      skip = 0,
                      label.placer = label_placer_n(1), 
                      nudge_x = .005, 
                      nudge_y = .005) +
    scale_fill_gradientn(trans = "log10", limits = c(50, 6000), oob = scales::squish,colours=rainbow(90)) + 
    #scale_fill_viridis(trans = "log10", limits = c(100, 5000), oob = scales::squish, direction = -1) + 
    theme_classic() + 
    theme(legend.position = "bottom", 
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=16, face = "bold"),
          legend.title = element_text(size=16, face = "bold"),
          strip.text = element_text(size = 16, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50"),
          legend.key.width = unit(1.5, "cm"),
          panel.spacing = unit(1, "lines")) + 
    scale_x_continuous(expand = c(0,0), name = "Follow-up Time") + 
    scale_y_continuous(expand = c(0,0), name = "Correlation (r)",breaks = seq(.05,.4, by = .05))  + 
    coord_cartesian(xlim = c(2, 8), ylim = c(0.05, .4), expand = FALSE) +
    facet_wrap(as.factor(t)~.)
  
  fname = here("data_reliability_long/results/figs/samplesize", paste(model, ft, "png", sep = "."))
  ggsave(plot = gs, file = fname, width = 6, height = 6)
  return(gs)
}


fig_mean_icc_modality = function(df.fig, oname) {
  pd=position_dodge(width = .4)
  
  dat = df.fig %>% 
    group_by(t,total_time, modality) %>% 
    summarise(iccM = mean(icc, na.rm = T), iccSD = sd(icc, na.rm = T)) %>% 
    mutate(t = as.factor(t))
  
  dat$modality <- factor(dat$modality, levels=c('area', 'thickness', 'volume', 'svolume'),
                             labels=c('area', 'thickness', 'volume', 'subcortical'))
  
  gs = ggplot(dat, aes(x = total_time, y = iccM, group = interaction(modality, t), color = modality)) + 
    geom_errorbar(aes(ymin=iccM-iccSD, ymax=iccM+iccSD), colour="grey50", width=.4, position=pd, alpha = .4)  +
    geom_point(position = pd, size = 5, shape = 20) + 
    geom_smooth(position = pd, size = 2, se = F) + 
    scale_color_manual(name = "Modality", values = c("coral2", "bisque", "azure2", "darkolivegreen3")) + 
    scale_x_continuous(breaks = seq(2,12, by = 2)) + 
    theme_classic()+ 
    facet_wrap(t~., nrow = 1) +
    theme(legend.position = 'bottom',
          axis.text = element_text(size = 14, face = "bold"), 
          axis.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=14, face = "bold"),
          legend.title = element_text(size=14, face = "bold"),
          strip.text = element_text(size = 12.5, color = "grey50", face = "bold"), 
          strip.background = element_rect(size = 1.5, color = "grey50")) +
    ylab("ICC") + 
    xlab("Follow-up time") + 
    ylim(0,1) + 
    guides(color = guide_legend(keywidth = 3, keyheight = 1))
  
  
  fname = here("data_reliability_long/results/figs/icc_main", paste(oname ,"modality", "png", sep = "."))
  ggsave(file = fname, plot = gs, width = 12, height = 7)
  return(gs)
}


fig_linechart = function(infile, oname) {
  x = infile %>% 
    dplyr::select(dataset, rid, age, sub_id ) %>% 
    group_by(rid) %>% 
    mutate(minAge = min(age),
           maxAge = max(age)) %>% 
    ungroup() %>% 
    group_by(dataset) %>% 
    mutate(nrid = n_distinct(rid)) %>% 
    ungroup() %>% 
    arrange(minAge, maxAge)
  
  xx = data.frame(rid = unique(x$rid) )
  xx = left_join(xx,
                 x %>% dplyr::select(rid, dataset) %>% group_by(rid,dataset) %>% tally())
  xx = xx %>% group_by(dataset) %>% mutate(i = row_number())
  x = left_join(x,xx)
  
  gs = 
    ggplot(x, aes(x = age,
                  y = i, 
                  group = as.character(i), 
                  color = dataset)) + 
    geom_line(size = .1) + 
    geom_point(size =.1) + 
    #scale_color_viridis_c() +
    theme_classic() + 
    theme(legend.position = 'none') +
    facet_wrap(.~dataset, scales = "free_y") +
    ylab("N. of Participants")
  ggsave(plot = gs,  width = 5, height = 7, file = here("data_reliability_long/results/figs/demographics", paste("linechart", oname, "png", sep = ".")))
  return(gs)
} 


fig_raincloud = function(infile, oname){
  # plot raincloud-like plot
  gs = 
    ggplot(infile, aes(x = age, y = dataset, fill = stat(x))) +
    geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "age", option = "C") +
    scale_y_discrete(limits = unique(rev(infile$dataset))) +
    theme_ridges(font_size = 13, grid = TRUE) + 
    theme(legend.position = 'none',
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0.5))
  ggsave(plot = gs, width = 5, height = 7,  file = here("data_reliability_long/results/figs/demographics", paste("raincloud", oname, "png", sep = ".")))
  return(gs)
}

fig_hist_tps = function(dat, oname){
  dat = 
    dat %>% 
    mutate(n = if_else(n > 10, 10,n))
  
  gs = 
    ggplot(dat, aes(n, fill = dataset)) + 
    geom_bar() + 
    scale_x_continuous(breaks = c(seq(2,10,by = 1 ))) + 
    theme_classic()
  ggsave(plot = gs ,
         width = 5, 
         height = 7,
         file = here("data_reliability_long/results/figs/demographics", paste("hist","tps", oname, "png", sep = ".")))
  return(gs)
}

wrapper_table_aov = function(df.formula.icc, oname) {
  
  mod = 
    df.formula.icc %>% 
    filter(model == oname) %>% 
    mutate(modality = as.factor(modality),
           tps = as.ordered(t), 
           total_time = as.ordered(total_time)) %>% 
    lm(icc ~ modality*tps*total_time, data = .)
  
  # anova (time, fu, modality)
  mod.aov = aov(mod)
  
  k1 = broom::tidy(mod.aov) %>% 
    kable(digits = 3) %>% 
    kable_styling(full_width = FALSE) 
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("aov",oname, "html", sep = "."))
  save_kable(k1, fname)
  
  # effect sizes
  k2 = effectsize::eta_squared(mod.aov) %>% 
    kable(digits = 3) %>% 
    kable_styling(full_width = FALSE) 
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("eta",oname, "html", sep = "."))
  save_kable(k2, fname)
  
  
  
  # main effect of follow-up point - descriptive
  dat = 
    df.formula.icc %>% 
    filter(model == oname) %>% 
    group_by(total_time) %>% 
    summarise(iccM = mean(icc, na.rm = T), 
              iccSD = sd(icc, na.rm = T))
  
  k3 = kable(dat, digits = 2) %>% 
    kable_styling(full_width = FALSE) 
  
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("me_fu",oname, "html", sep = "."))
  save_kable(k3, fname)
  
  
  # main effect of time - descrptive
  dat = 
    df.formula.icc %>% 
    filter(model == oname) %>% 
    group_by(t) %>% 
    summarise(iccM = mean(icc, na.rm = T), 
              iccSD = sd(icc, na.rm = T))
  
  k4 = kable(dat, digits = 2) %>% 
    kable_styling(full_width = FALSE) 
  
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("me_nobs",oname, "html", sep = "."))
  save_kable(k4, fname)
  
  
  k5 = paste(round(abs((dat$iccM[1] - dat$iccM[4])/6),3),
             "Increase in reliability for each number of time point", 
             sep = " ")
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("delta_nobs",oname, "html", sep = "."))
  writeLines(k5, fname) 
  
  # main effect modality
  Tukey_modality =  emmeans(mod.aov, pairwise ~  modality)
  k6 = broom::tidy(Tukey_modality$emmeans) %>% 
    kable(digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("me_modality",oname, "html", sep = "."))
  save_kable(k6, fname)
  
  k7 = broom::tidy(Tukey_modality$contrasts)  %>% 
    kable(digits = 3) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("me_modality_posthoc",oname, "html", sep = "."))
  save_kable(k7, fname)
  
  #single effects, time, t
  dat = 
    df.formula.icc %>% 
    filter(model == oname) %>% 
    group_by(t,total_time) %>% 
    summarise(iccM = mean(icc, na.rm = T), 
              iccSD = sd(icc, na.rm = T)) %>% 
    pivot_wider(names_from = t,
                values_from = c(iccM, iccSD))
  
  k8 = dat  %>% 
    kable(digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("me_fu_nbs",oname, "html", sep = "."))
  save_kable(k8, fname)
  
  
  
  dat = 
    df.formula.icc %>% 
    filter(model == oname) %>% 
    group_by(t,total_time, modality) %>% 
    summarise(iccM = mean(icc, na.rm = T), 
              iccSD = sd(icc, na.rm = T)) %>% 
    pivot_wider(names_from = t,
                values_from = c(iccM, iccSD)) %>% 
    arrange(modality) 
  
  k9 = dat  %>% 
    kable(digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/icc_aovs", 
               paste("me_fu_nbs_modality",oname, "html", sep = "."))
  save_kable(k9, fname)
  
  # # only visualization - ordered factors make interpretation more easy 
  # g = emmip(mod.aov, total_time ~modality, plotit = FALSE)
  # gs = emmip_ggplot(g)
  # fname = here("data_reliability_long/results/stats_and_tables/icc_main", paste("anova_modality_tps_fu","Tukey","modality","totaltime","png", sep = "."))
  # ggsave(plot = gs, file = fname)
  # 
  # g = emmip(mod.aov, total_time ~tps,plotit = FALSE)
  # gs = emmip_ggplot(g)
  # fname = here("data_reliability_long/results/stats_and_tables/icc_main", paste("anova_modality_tps_fu","Tukey","tps","totaltime","png", sep = "."))
  # ggsave(plot = gs, file = fname)
  
}

wrapper_table_parameter_comparison = function() {
  load(file.path(outreliability.comp, "icc_reliability_between_samples.Rda"))
  
  icc.me = 
    dat.me %>% 
    select(-c(feature, modality)) %>% 
    ICC(.) %>%
    .$results %>% 
    filter(type %in% c("ICC2", "ICC2k")) %>% 
    mutate(parameter = "error")
  
  icc.slopeV = 
    dat.slopeV %>% 
    select(-c(feature, modality)) %>% 
    ICC(.) %>% 
    .$results %>% 
    filter(type %in% c("ICC2", "ICC2k")) %>% 
    mutate(parameter = "slope")
  
  icc.global = rbind(icc.me,
                     icc.slopeV)
  
  k1 =icc.global %>% 
    kable(digits = 2) %>% 
    kable_styling(full_width = FALSE) 
  fname = here("data_reliability_long/results/stats_and_tables/model_comparison", 
               paste("global_agreement","html", sep = "."))
  save_kable(k1, fname)
  
  
  icc.me = 
    dat.me %>% 
    select(-feature) %>% 
    group_by(modality) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ICC(.)),
           icc.21_M = map_dbl(mod, ~.$results$ICC[2]),
           icc.2_M = map_dbl(mod, ~.$results$ICC[5]),
           icc.21_lb = map_dbl(mod, ~.$results$`lower bound`[2]),
           icc.2_lb = map_dbl(mod, ~.$results$`lower bound`[5]),
           icc.21_ub = map_dbl(mod, ~.$results$`upper bound`[2]),
           icc.2_ub = map_dbl(mod, ~.$results$`upper bound`[5])) %>% 
    select(-c(mod, data)) %>% 
    pivot_longer(-modality, names_to = "names", values_to = "values") %>% 
    separate(names, c("icc.type", "grot.2"), sep = "_") %>% 
    pivot_wider(names_from = grot.2, 
                values_from = values)  %>% 
    mutate(meas = "me")
  
  
  icc.slopeV = 
    dat.slopeV %>% 
    select(-feature) %>% 
    group_by(modality) %>% 
    nest() %>% 
    mutate(mod = map(data, ~ICC(.)),
           icc.21_M = map_dbl(mod, ~.$results$ICC[2]),
           icc.2_M = map_dbl(mod, ~.$results$ICC[5]),
           icc.21_lb = map_dbl(mod, ~.$results$`lower bound`[2]),
           icc.2_lb = map_dbl(mod, ~.$results$`lower bound`[5]),
           icc.21_ub = map_dbl(mod, ~.$results$`upper bound`[2]),
           icc.2_ub = map_dbl(mod, ~.$results$`upper bound`[5])) %>% 
    select(-c(mod, data)) %>% 
    pivot_longer(-modality, names_to = "names", values_to = "values") %>% 
    separate(names, c("icc.type", "grot.2"), sep = "_") %>% 
    pivot_wider(names_from = grot.2, 
                values_from = values) %>% 
    mutate(meas = "slope")
  
  
  icc.by.modality = 
    rbind(icc.me, 
          icc.slopeV)  
  
  k2 = icc.by.modality %>% 
    kable(digits = 2) %>% 
    kable_styling(full_width = FALSE) 
  fname = here("data_reliability_long/results/stats_and_tables/model_comparison", 
               paste("agreement_by_modality","html", sep = "."))
  save_kable(k2, fname)
} 

wrapper_table_effect_size = function(df.formula.icc, oname) {
  
  dat = 
    df.formula.icc  %>% 
    filter(model == oname) %>% 
    unnest(data) %>% 
    filter(StatTest == "rpearson")
  
  mod = 
    dat %>% 
    group_by(t, total_time, effectsize) %>% 
    summarise(xn = ceiling(mean(n)), 
              sdn = ceiling(sd(n)))
  mod = 
    mod %>% pivot_wider(names_from = "t", values_from = c("xn", "sdn")) %>%  arrange(effectsize)
  k1 = kable(mod, digits = 1) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/sample_size", paste(oname,"mean","html", sep = "."))
  save_kable(k1, fname)
  
  
  fs = c("lh_entorhinal_thickness",
         "Left-Hippocampus")
  
  mod.1 = 
    dat %>% filter(feature %in% fs ) %>% 
    select(t, total_time, effectsize, n, feature) %>% 
    pivot_wider(names_from = "t", values_from = n) %>%  arrange(feature, effectsize)
  k2= kable(mod.1, digits = 1) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/sample_size", paste(oname,"fs","html", sep = "."))
  save_kable(k2, fname)
}

wrapper_table_individual_overlap = function(df.formula.icc, oname) {
  X = 
    df.formula.icc %>% 
    filter(model == oname) %>% 
    group_by(total_time, t) %>% 
    summarise(Xbc.main.v.decl = mean(bc.main.v.decl, na.rm = T),
              SDbc.main.v.decl = sd(bc.main.v.decl, na.rm = T),
              Xbc.main.v.normal = mean(bc.main.v.normal, na.rm = T),
              SDbc.main.v.normal = sd(bc.main.v.normal, na.rm = T),
              Xp.main.v.decl = mean(p.main.v.decline, na.rm = T), 
              SDp.main.v.decl = sd(p.main.v.decline, na.rm = T),
              Xp.main.v.normal = mean(p.main.v.normal, na.rm = T), 
              SDp.main.v.normal = sd(p.main.v.normal, na.rm = T))
  X = 
    X %>% 
    pivot_longer(-c(total_time, t),
                 names_to = "measure", 
                 values_to = "value") %>% 
    pivot_wider(names_from = measure,
                values_from = value) 
  
  k1 = 
    X %>% 
    kable( digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/individual_overlap", paste(oname,"html", sep = "."))
  save_kable(k1, fname)
  
  
  
  fs = c("lh_entorhinal_thickness",
         "Left-Hippocampus")
  
  
  k2 = 
    df.formula.icc %>% 
    filter(model == "model1") %>%
    filter(feature %in% fs ) %>% 
    select(feature, t, total_time, starts_with("bc."), starts_with("p.")) %>% 
    pivot_wider(names_from = feature,
                values_from = -c(feature, t, total_time)) %>% 
    kable( digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/individual_overlap", paste(oname,"fs", "html", sep = "."))
  save_kable(k2, fname)
}

wrapper_table_agegroup = function(df.formula.icc, oname) {
  phenotypes.agegroup = 
    c("Left-Hippocampus", 
      "lh_entorhinal_thickness")
  
  dat = 
    df.formula.icc %>% 
    filter(feature %in% phenotypes.agegroup, 
           startsWith(model, oname)) %>% 
    select(feature, model, total_time, t, icc)
  
  dat = dat %>% pivot_wider(names_from = c(feature,t),  
                            values_from = icc)
  
  k1 = 
    dat %>% 
    kable( digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/agecomparisons", paste(oname,"html", sep = "."))
  save_kable(k1, fname)
}

wrapper_table_maintenance = function(oname) {
  phenotypes.m = 
    c("Left-Hippocampus", 
      "lh_entorhinal_thickness")
  
  
  grot = list()
  for (pheno in phenotypes.m) {
    indtrajdir = file.path(outreliability.comp, "individual_trajectories")
    load(file.path(indtrajdir,paste0(oname, "_", pheno, ".Rda")))
    grot[[pheno]] = 
      P_extreme %>% mutate(feature = pheno)
  }
  
  grot = data.table::rbindlist(grot) 
  grot = grot %>% 
    pivot_wider(names_from = tps, 
                values_from = c(pct_extrem, probM_gvEx, probD_gvEx))
  
  k1 = 
    grot %>% 
    kable( digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/maintenance", paste(oname,"html", sep = "."))
  save_kable(k1, fname)
}

wrapper_table_descriptive_test_retest = function(oname) {
  fname = here("data-raw/tabulated/parameters_reliability/demog/all.rda" )
  load(fname)
  
  k1 = 
    dat.demog %>% 
    group_by(dataset) %>% 
    summarise(nID = n_distinct(id), 
              m = sum(sex),
              nObs = first(n),
              ageM = mean(age, na.rm = T),
              ageSD = sd(age, na.rm = T), 
              agemin = min(age, na.rm = T),
              agemax = max(age, na.rm = T),
              delayM = mean(delay, na.rm = T),
              delaySD = sd(delay, na.rm = T), 
              delaymin = min(delay, na.rm = T),
              delaymax = max(delay, na.rm = T))
  k1 = 
    k1 %>% 
    kable( digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/descriptives", paste(oname,"html", sep = "."))
  save_kable(k1, fname)
}


wrapper_table_descriptive_long_dataset = function(dat, oname) {
  mod.group = 
    dat %>% 
    group_by(dataset) %>% 
    summarise(subs = n_distinct(rid), 
              males = sum(sex), 
              obs = sum(n), 
              ageX = mean(xage), 
              ageSD = sd(xage), 
              agemin = min(xage), 
              agemax = max(xage),
              timeX = mean(time), 
              timeSD = sd(time), 
              timemin = min(time), 
              timemax = max(time), 
              obsX = mean(tps), 
              obsSD = sd(tps),
              obsmin = min(tps),
              obsmax = max(tps),
              nsite = n_distinct(site))
  
  mod.all = 
    dat%>% 
    ungroup() %>% 
    summarise(subs = n_distinct(rid), 
              males = sum(sex), 
              obs = sum(n),
              ageX = mean(xage), 
              ageSD = sd(xage),
              agemin = min(xage), 
              agemax = max(xage),
              timeX = mean(time), 
              timeSD = sd(time), 
              timemin = min(time), 
              timemax = max(time), 
              obsX = mean(tps), 
              obsSD = sd(tps),
              obsmin = min(tps),
              obsmax = max(tps),
              nsite = n_distinct(site)) %>% 
    mutate(dataset = "all")
  
  k1 = bind_rows(mod.group, mod.all)
  
  k1 = 
    k1 %>% 
    kable( digits = 2) %>%
    kable_styling(full_width = F)
  fname = here("data_reliability_long/results/stats_and_tables/descriptives", paste(oname,"html", sep = "."))
  save_kable(k1, fname)
  
}