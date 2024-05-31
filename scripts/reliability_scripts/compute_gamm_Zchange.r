args = commandArgs(TRUE)

outdir=as.character(args[1])
model=as.character(args[2])
phenotypes=as.character(args[3])
global=as.logical(args[4])
cross=as.logical(args[5])


compute_gamm_Zchange = function(outdir, model, phenotypes, global = F, cross = F) {
  library(tidyverse)
  library(broom)
  library(gamm4)
  if (global == F) {
    load(file.path(outdir, "df.all.filt.Rda"))
  } else if (global == T) {
    load(file.path(outdir, "df.all.global.filt.Rda"))
  }  else if (cross == T) {
    load(file.path(outdir, "df.all.cross.filt.Rda"))
  }
  
  idpdir = file.path(outdir, model, phenotypes)
  try(dir.create(idpdir))
  
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
  if (cross == F) {
  save("df.out", 
       file = file.path(idpdir, "Z_predict_extended.Rda"))
  }
  
  out = get_data_empirical(df.merge.long, mod, rid.intercept, df.out)
  
  if (cross == T) {
    fname = file.path(idpdir, "model_parameters.cross.csv")
  } else {
    fname = file.path(idpdir, "model_parameters.csv")
  }
  write.table(out,
            quote = F,
            row.names = F, 
            col.names = T,
            file =fname)
}

get_data_empirical = function(df.merge.long, mod, rid.intercept, df.out) {
   library(gratia)
  # get data for empirical model 
   x = expand.grid(age = c(60,65,70,75,80), sex = c(0.5))
   x$fit = predict(mod$gam, x)
   sd1 = rid.intercept$`(Intercept)` %>% sd()
   ch1 = 100*(x$fit[5] - x$fit[1])/x$fit[1]/20
   ch2 = 100*(x$fit[5] - x$fit[2])/x$fit[2]/15
   ch3 = 100*(x$fit[5] - x$fit[3])/x$fit[3]/10
   xx = derivatives(mod$gam)
   s = xx %>% filter(between(data,60,80))
   ch4= 100*mean(s$derivative)/x$fit[3]
   s = xx %>% filter(between(data,70,80))
   ch5= 100*mean(s$derivative)/x$fit[4]
   
   fit.df = data.frame(t(x$fit))
   names(fit.df) = c("mean_60y", "mean_65y", "mean_70y", "mean_75y","mean_80y")
   ch.df = data.frame(sd1, ch1,ch2,ch3,ch4,ch5)
   
  subs = 
    df.merge.long%>% 
    group_by(rid)  %>% 
    mutate(time = max(age) - min(age)) %>% 
    summarise(n = n_distinct(age), 
              time = first(time), 
              dataset = first(dataset), 
              age = mean(age))
  
  dat = inner_join(df.out,subs)
  
  dat = 
    dat %>% 
    filter(time > 4 & n > 3)
  
  
  dat.all =
    dat %>% 
    ungroup() %>% 
    mutate(meanV = mean(deltaZ), 
           madV = mad(deltaZ)) %>% 
    mutate(new.deltaZ = if_else(between(deltaZ, meanV - madV*5, meanV + madV*5), deltaZ, NaN)) %>% 
    summarise(madDelta = mad(deltaZ), 
              sdDelta = sd(new.deltaZ, na.rm = T))
  
  
  dat.age =
    dat %>% 
    ungroup() %>% 
    mutate(ageG = if_else(age < 60, "Y", "O")) %>% 
    group_by(ageG) %>% 
    mutate(meanV = mean(deltaZ), 
           madV = mad(deltaZ)) %>% 
    ungroup() %>% 
    mutate(new.deltaZ = if_else(between(deltaZ, meanV - madV*5, meanV + madV*5), deltaZ, NaN)) %>% 
    group_by(ageG) %>% 
    summarise(madDelta = mad(deltaZ), 
              sdDelta = sd(new.deltaZ, na.rm = T)) %>% 
    pivot_wider(names_from = ageG, 
                values_from = -ageG)
  
  
  
  dat.ds =
    dat %>% 
    ungroup() %>% 
    group_by(dataset) %>% 
    mutate(meanV = mean(deltaZ), 
           madV = mad(deltaZ)) %>% 
    ungroup() %>% 
    mutate(new.deltaZ = if_else(between(deltaZ, meanV - madV*5, meanV + madV*5), deltaZ, NaN)) %>% 
    group_by(dataset) %>% 
    summarise(madDelta = mad(deltaZ), 
              sdDelta = sd(new.deltaZ, na.rm = T)) %>% 
    pivot_wider(names_from = dataset, 
                values_from = -dataset)
  
 
  
  out = data.frame(fit.df, ch.df, dat.all, dat.age, dat.ds)
  return(out)
}

compute_gamm_Zchange(outdir, model, phenotypes, global, cross)





