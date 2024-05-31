fetch_aparc_data_s2c = function(idir) {
  lfiles = list.files(idir)
  subs = 
    strsplit(lfiles,"sub-") %>% 
    reduce(rbind) %>% 
    .[,2] %>% 
    substr(., 1,7) %>% 
    unique()
  
  df.mri = lapply(lfiles, function(x) { import(file.path(idir, x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  names(df.mri) = lfiles
  
  df.mri.all = list()
  for (i in 1:length(subs)) {
    idx = lfiles[grep(subs[i], lfiles)]
    x = df.mri[idx]
    df.mri.all[[i]] = x %>% reduce(cbind)
  }
  names(df.mri.all) = subs
  df.mri = data.table::rbindlist(df.mri.all, idcol = "subs")
  return(df.mri)
}

add_global_variables = function(df, df.pheno) {
  for (i in 1:length(df.pheno$phenos)) {
    df = 
      df %>% 
      mutate(!!df.pheno$phenos[[i]] := rowMeans(select(., all_of(df.pheno$orig[[i]]))))
  }
  return(df)
}

save_measurment_error_data = function(odir, site, dat, phenotypes, df.pheno, is.cross = F) {
  
  if (is.cross == T) {
    odir = paste0(odir, "_cross")
    try(dir.create(odir))
  }
  save(
    dat,
    phenotypes,
    df.pheno,
    file = file.path(odir, 
                     paste(site, 
                           "data.rda", 
                           sep = "_")))
  
  
  idx = match(phenotypes,dat$avg.error$features)
  df.out = dat$avg.error[idx, ]
  
  write.table(
    df.out,
    quote = F,
    row.names = F,
    file = file.path(odir, 
                     paste(site, 
                           "data.csv", 
                           sep = "_")))
  
  idx = match(df.pheno$phenos,dat$avg.error$features)
  df.out = dat$avg.error[idx, ]
  write.table(
    df.out,
    quote = F,
    row.names = F,
    file = file.path(odir, 
                     paste(site, 
                           "data",
                           "globalVars.csv", 
                           sep = "_")))
}



fetch_mri_standard = function(sitedir, mrifiles, is.cross = F) {
  
  if (is.cross == T) {
    df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", "cross", x))})
  } else {
    df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", "desikan", x))})
  }
  
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri = 
    df.mri %>% 
    separate(input, c("subs",
                      "grot2"),
             sep = "ses", 
             remove = F) %>% 
    mutate(visit = substr(grot2, 1,2) %>% as.numeric(),
           subs = gsub("./", "", subs)) %>% 
    select(-grot2)
  return(df.mri)
}

merge_data = function(lfiles, odir, ofile, phenotypes, is.cross = F) {
  
  if (is.cross == T) {
    odir = paste0(odir, "_cross")
    }
  df = lapply(lfiles, function(x) { import(file.path(odir, x))})
  names(df) = sites
  df = data.table::rbindlist(df, idcol = "site")
  
  all = 
    df %>% 
    group_by(features) %>% 
    summarise_if(is.numeric, mean)
  all$site = "all"
  
  df = rbind(df, all)
  df= 
    df %>% 
    pivot_wider(values_from = -c(site, features), 
                names_from = site)
  
  df = df[match(phenotypes, df$features), ]
  
  save(df, 
       file = file.path(odir, 
                        paste(ofile, 
                              "rda", 
                              sep = ".")))
  
  write.table(
    df,
    quote = F,
    row.names = F,
    file = file.path(odir, 
                     paste(ofile, 
                           "csv", 
                           sep = ".")))
}

test_retest_error_2tp = function(df) {
  x = 
    df %>% 
    group_by(subs) %>% 
    pivot_longer(-c(subs, input, visit), 
                 names_to = "features", 
                 values_to = "values") %>% 
    arrange(subs,visit)
  
  df.error =
    x %>% 
    group_by(subs, features) %>% 
    mutate(v1 = first(values),
           v2 = last(values),
           M = mean(values),
           errM = 100*abs(v1-v2)/M) %>% 
    ungroup() %>% 
    group_by(features) %>% 
    mutate(pctM = median(errM, na.rm = T),
           pctmad = mad(errM, na.rm = T), 
           err.R = if_else(errM > pctM +  5*pctmad,NaN,errM))
  
  df.error.avg = 
    df.error %>% group_by(features) %>% 
    summarise(pct.err.mean.raw = mean(errM, na.rm =T),
              pct.err.mean = mean(err.R, na.rm =T),
              pct.err.median = median(errM, na.rm =T),
              meanV = mean(M), 
              sdV = sd(M))
  
  df.out = list()
  df.out$df = x
  df.out$subj.error = df.error
  df.out$avg.error = df.error.avg
  return(df.out)
}

test_retest_error_within = function(df) {
  x = 
    df %>% 
    group_by(subs) %>% 
    pivot_longer(-c(subs, input, visit), 
                 names_to = "features", 
                 values_to = "values") %>% 
    arrange(subs,visit)
  
  df.error =
    x %>% 
    group_by(subs, features) %>% 
    mutate(meanV = mean(values), 
           medianV = median(values), 
           madV = mad(values), 
           sdV = sd(values)) %>%
    ungroup() %>% 
    mutate(newvals = if_else(between(values, medianV - 5*madV, meanV + 5*madV), values, NaN)) %>% 
    group_by(subs, features) %>% 
    summarise(meanV = mean(values), 
              sdV = sd(values),
              meanR = mean(newvals, na.rm = T), 
              sdR = sd(newvals, na.rm = T),
              medianV = median(values), 
              madV = mad(values)) %>% 
    mutate(cv.robust = 100*madV/medianV, 
           cv.outliers = 100*sdR/meanR, 
           cv.raw = 100*sdV/meanV)  
  
  
  
  df.error.avg = 
    df.error %>% 
    ungroup() %>% 
    group_by(features) %>% 
    summarise(pct.err.mean.raw = mean(cv.raw, na.rm =T),
              pct.err.mean = mean(cv.outliers, na.rm =T),
              pct.err.median = mean(cv.robust, na.rm =T),
              meanV = mean(meanR), 
              sdV = mean(sdR))
  
  df.out = list()
  df.out$df = x
  df.out$subj.error = df.error
  df.out$avg.error = df.error.avg
  
  return(df.out)  
}



wrapper_demog_s2c = function() {
  # open subject. keep mri from those with demographic data
  dat.demog = import(file.path(sitedir,"S2C_data.csv")) %>% 
    rename("subs" = "subject_id", 
           "visit" = "wave_code") %>% 
    mutate(subs = as.character(subs))
  
  grot = 
    dat$df %>% 
    group_by(subs, visit) %>% 
    tally()
  
  dat.demog = inner_join(dat.demog , grot) %>% 
    group_by(subs) %>% 
    summarise(sex = first(subject_sex), 
            delay = 365.25*(max(visit_age) - min(visit_age)), 
            age = first(visit_age))
  
  dat.demog = 
    dat.demog %>% 
    mutate(n = 2,
           sex = if_else(sex =="Male",1 ,0)) %>% 
    select(subs, 
           sex, 
           age, 
           delay, 
           n) 
  names(dat.demog) = c("id", 
                        "sex", 
                        "age", 
                        "delay", 
                        "n")
  
  save(dat.demog, file = file.path(odir, "demog", paste(site, "rda", sep =".")))
}


wrapper_demog_maclaren = function() {
  dat.demog = 
    tribble(
      ~ID, ~Sex, ~Age, ~N, ~ISI,
      "S1", "M", 26,40, 31,
      "S2", "M", 31,40,31,
      "S3", "F", 30,40,  31
    )
  
  dat.demog = 
    dat.demog %>% 
    mutate(Sex = if_else(Sex =="M",1 ,0)) %>% 
    select(ID, 
           Sex, 
           Age, 
           ISI, 
           N) 
  names(dat.demog) = c("id", 
                       "sex", 
                       "age", 
                       "delay", 
                       "n")
  
  save(dat.demog, file = file.path(odir, "demog", paste(site, "rda", sep =".")))
}


wrapper_demog_oasis = function() {
  dat.demog  = import(file.path(sitedir, "oasis_cross-sectional.csv")) %>% 
    mutate(Delay = if_else(Delay == "N/A", NA, as.numeric(Delay))) %>% 
    filter(!is.na(Delay))
  
  dat.demog = 
    dat.demog %>% 
    mutate(n = 2,
           sex = if_else("M/F" =="M",1 ,0)) %>% 
    select(ID, 
           sex, 
           Age, 
           Delay, 
           n) 
  names(dat.demog) = c("id", 
                       "sex", 
                       "age", 
                       "delay", 
                       "n")
  
  save(dat.demog, file = file.path(odir, "demog", paste(site, "rda", sep =".")))
}

wrapper_demog_preventad= function() {
  dat.demog = read.csv(file.path(sitedir, "UM_phenotypic_data.csv"), na.strings = "#") %>% 
    group_by(SUBID) %>% 
    summarise(AGE_AT_SCAN_1 = max(AGE_AT_SCAN_1, na.rm = T),
              RETEST_DURATION = max(RETEST_DURATION, na.rm = T),
              SEX = first(SEX, na_rm = T))
  
  dat.demog = 
    dat.demog %>% 
    mutate(n = 2,
           sex = if_else(SEX==1,1 ,0)) %>% 
    select(SUBID, 
           sex, 
           AGE_AT_SCAN_1, 
           RETEST_DURATION, 
           n) 
  names(dat.demog) = c("id", 
                       "sex", 
                       "age", 
                       "delay", 
                       "n")
  
  save(dat.demog, file = file.path(odir, "demog", paste(site, "rda", sep =".")))
}

wrapper_demog_hnu= function() {
  dat.demog = read.csv(file.path(sitedir, "HNU_1_phenotypic_data.csv"), na.strings = "#") %>% 
    group_by(SUBID) %>% 
    summarise(AGE_AT_SCAN_1 = max(AGE_AT_SCAN_1, na.rm = T),
              RETEST_DURATION = max(RETEST_DURATION, na.rm = T),
              SEX = first(SEX, na_rm = T))
  
  dat.demog = 
    dat.demog %>% 
    mutate(n = 10,
           sex = if_else(SEX==1,1 ,0)) %>% 
    select(SUBID, 
           sex, 
           AGE_AT_SCAN_1, 
           RETEST_DURATION, 
           n) 
  names(dat.demog) = c("id", 
                       "sex", 
                       "age", 
                       "delay", 
                       "n")
  
  save(dat.demog, file = file.path(odir, "demog", paste(site, "rda", sep =".")))
}

wrapper_demog_gsp = function() {
  dat.demog  = import(file.path(sitedir, "GSP_retest_140630.csv")) %>% 
    separate(Subject_ID, c("ID", "Sess")) %>% 
    group_by(ID) %>% 
    summarise(Age = max(Age_Bin, na.rm = T),
              Delay= max(Delay, na.rm = T),
              SEX = first(Sex, na_rm = T))
  
  dat.demog = 
    dat.demog %>% 
    mutate(n = 2,
           sex = if_else(SEX=="M",1 ,0)) %>% 
    select(ID, 
           sex, 
           Age, 
           Delay, 
           n) 
  names(dat.demog) = c("id", 
                       "sex", 
                       "age", 
                       "delay", 
                       "n")
  
  save(dat.demog, file = file.path(odir, "demog", paste(site, "rda", sep =".")))
}

