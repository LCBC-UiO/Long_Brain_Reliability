set_links = function(folder, is.cross = F) {
  sitedir = here('data_reliability_long','df_mri', folder)
  if (is.cross == T) {
    load(file.path(sitedir, "df.cross.rda"))  
  }else {
    load(file.path(sitedir, "df.rda"))  
  }
  
  return(df.out)
}

rename_variables = function(df.out, site) { 
  grot = df.out %>% 
    ungroup() %>% 
    select(df.harmonize[[site]])
  names(grot) = df.harmonize$normative_modelling
  sex.f = sex.rename$female[grep(site,sex.rename$dataset)]
  sex.m = sex.rename$male[grep(site,sex.rename$dataset)]
  sex.fn = sex.rename$female[grep("normative",sex.rename$dataset)] %>% as.numeric()
  sex.mn = sex.rename$male[grep("normative",sex.rename$dataset)] %>% as.numeric()
  grot = 
    grot %>% 
    mutate(sex= 
             if_else(sex == sex.f, sex.fn,
                     if_else(sex == sex.m,sex.mn,NaN)))
  return(grot)
}


squeue = function(user, job) {
  df.squeue = system(paste0("squeue --name=", job," -u ",user), intern = T) %>% 
    strsplit(., " +") %>% 
    simplify2array() %>% 
    t() %>% 
    as.data.frame()
  return(df.squeue)
  
}

prepare.uio.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, scanner, sep = "_")) 
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("uio data renamed. ready to be merged")
}


prepare.ucam.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("ucam data renamed. ready to be merged")
}


prepare.mpib.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("mpib data renamed. ready to be merged")
}

prepare.umu.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("umu data renamed. ready to be merged")
}


prepare.ub.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("ub data renamed. ready to be merged")
}

prepare.vumc.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, as.character(scanloc), sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("vumc data renamed. ready to be merged")
}

prepare.habs.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "01", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("habs data renamed. ready to be merged")
}

prepare.ous.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, scanner, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("ous data renamed. ready to be merged")
}

prepare.aibl.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, model, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("aibl data renamed. ready to be merged")
}


prepare.adni.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    select(-site) %>% 
    mutate(site = paste(site, paste(Site, model,sep ="-"), sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           str_pad(df.out$site %>% 
               as.factor() %>% 
               as.numeric() %>% 
               as.character(),
             3, 
             pad = "0")) %>% 
    as.numeric()
  
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("adni data renamed. ready to be merged")
}

prepare.preventad.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, "1", sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("preventad data renamed. ready to be merged")
}



prepare.wayne.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, scanner, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("wayne data renamed. ready to be merged")
}

prepare.oasis3.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(site = paste(site, scanner, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum)
  return(df.out)
  
  cat("oasis3 data renamed. ready to be merged")
}

prepare.ukb.brainchart = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(input = paste(eid, tp, sep = "-"),
           site = paste(site, uk_biobank_assessment_centre_f54, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = rename_variables(df.out, site) %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum) %>% 
    filter(!is.na(EstimatedTotalIntraCranialVol))
  return(df.out)
  
  cat("ukb data renamed. ready to be merged")
}

prepare.ukb.brainchart_link = function(df.out, site, codesite) {
  df.out = 
    df.out %>% 
    mutate(input = paste(eid, tp, sep = "-"),
           site = paste(site, uk_biobank_assessment_centre_f54, sep = "_"))
  
  df.out$sitenum = 
    paste0(codesite, 
           "00", 
           df.out$site %>% 
             as.factor() %>% 
             as.numeric() %>% 
             as.character()) %>% 
    as.numeric()
  
  df.out = 
    df.out %>% 
    mutate(vol_tempole_lh = NaN, 
           vol_tempole_rh = NaN, 
           cth_tempole_lh = NaN, 
           cth_tempole_rh = NaN, 
           area_tempole_lh = NaN, 
           area_tempole_rh = NaN)
  df.out = rename_variables(df.out, "ukb_link") %>% 
    drop_na(sub_id, rid, age, sex, site, sitenum) %>% 
    filter(!is.na(EstimatedTotalIntraCranialVol))
  return(df.out)
  
  cat("ukb data renamed. ready to be merged")
}

stats.overview = function(df.merge.long, df.merge) {
  mod = list()
  df.merge.long = 
    df.merge.long %>% 
    group_by(rid) %>% 
    mutate(futime = max(age)- min(age),
           mage = mean(age))  
  
  df.merge = 
    df.merge %>% 
    group_by(rid) %>% 
    mutate(mage = mean(age))
  
  df.long.dmg = 
    df.merge.long %>% 
    group_by(rid) %>% 
    summarise(sex = first(sex), 
              mage = first(mage), 
              futime = first(futime),
              n = n_distinct(sub_id), 
              dataset = first(dataset))
  
  df.all.dmg = 
    df.merge %>% 
    group_by(rid) %>% 
    summarise(sex = first(sex), 
              mage = first(mage), 
              dataset = first(dataset))
  
  # n information
  mod$n$long.data = df.merge.long$sub_id %>% length()
  
  mod$n$all.data = df.merge$sub_id %>% length()
  
  mod$n$long.data.dataset = 
    df.merge.long %>% 
    group_by(dataset) %>% 
    tally()
  
  mod$n$all.data.dataset = 
    df.merge %>% 
    group_by(dataset) %>% 
    tally()
  
  mod$n$long.data.site = 
    df.merge.long %>% 
    group_by(site) %>% 
    tally()
  
  mod$n$all.data.site = 
    df.merge %>% 
    group_by(site) %>% 
    tally()
  
  mod$n$long.unique.subs = df.long.dmg$rid %>% length()
  mod$n$all.unique.subs = df.all.dmg$rid %>% length()
  
  mod$n$long.unique.subs.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    tally()
  
  mod$n$all.unique.subs.dataset = 
    df.all.dmg %>% 
    group_by(dataset) %>% 
    tally()
  
  # age information  
  mod$age$all.mean.subs.age.dataset = 
    df.all.dmg %>% 
    group_by(dataset) %>% 
    summarise(mean = mean(mage), sd = sd(mage))
  
  mod$age$long.mean.subs.age.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    summarise(mean = mean(mage), sd = sd(mage))
  
  mod$age$all.mean.subs.age = 
    df.all.dmg %>% 
    summarise(mean = mean(mage), sd = sd(mage))
  
  mod$age$long.mean.subs.age = 
    df.long.dmg %>%
    summarise(mean = mean(mage), sd = sd(mage))
  
  mod$n$sites = 
    df.merge %>% 
    summarise(sites = n_distinct(site))
  
  mod$n$sites.dataset = 
    df.merge %>% 
    group_by(dataset) %>% 
    summarise(sites = n_distinct(sitenum))
  
  # SEX
  mod$sex$all.sex.dataset = 
    df.all.dmg %>% 
    group_by(dataset) %>% 
    summarise(female = length(sex[sex == 0]),
              male = length(sex[sex == 1]))
  
  mod$sex$long.sex.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    summarise(female = length(sex[sex == 0]),
              male = length(sex[sex == 1]))
  
  mod$sex$all.sex = 
    df.all.dmg %>% 
    summarise(female = length(sex[sex == 0]),
              male = length(sex[sex == 1]))
  
  mod$sex$long.sex = 
    df.long.dmg %>%
    summarise(female = length(sex[sex == 0]),
              male = length(sex[sex == 1]))
  
  # follow up time and observations
  mod$follow.up$long.follow.up.span.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    summarise(mean = mean(futime), sd = sd(futime))
  
  mod$follow.up$long.follow.up.span = 
    df.long.dmg %>%
    summarise(mean = mean(futime), sd = sd(futime))
  
  mod$follow.up$long.observations.per.participant.dataset = 
    df.long.dmg %>% 
    group_by(dataset) %>% 
    summarise(mean = mean(n), sd = sd(n))
  
  mod$follow.up$long.observations.per.participant = 
    df.long.dmg %>%
    summarise(mean = mean(n), sd = sd(n))
  
  return(mod)
} 

linechart = function(infile, outdir, filename) {
  x = infile %>% 
    select(dataset, rid, age, sub_id ) %>% 
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
                 x %>% select(rid, dataset) %>% group_by(rid,dataset) %>% tally())
  xx = xx %>% group_by(dataset) %>% mutate(i = row_number())
  x = left_join(x,xx)
  
  ggplot(x, aes(x = age,
                y = i, 
                group = as.character(i), 
                color = dataset)) + 
    geom_line(size = .1) + 
    geom_point(size =.1) + 
    #scale_color_viridis_c() +
    theme_classic() + 
    theme(legend.position = 'none') +
    facet_wrap(.~dataset, scales = "free_y")
  ggsave( file = file.path(outdir, filename))
} 
  
  
filter_sites_with_low_n = function(df.merge.long,df.merge, thr.obs) {
  i = 0
  c = 0
  # step filters scans with low numbers of subs, and the removes subs with n == 1 obserbation
  while (i == 0) {
    c = c +1
    print(c)
    grot.long = 
      df.merge.long %>% 
      group_by(site) %>% 
      mutate(nscans = n_distinct(sub_id)) %>% 
      filter(nscans >= thr.obs) %>% 
      ungroup() %>% 
      group_by(rid) %>% 
      mutate(nsub = n_distinct(age)) %>% 
      filter(nsub > 1) %>% 
      ungroup()
    
    # removes scanners where n given tp == 1  < 25
    grot.all = 
      df.merge %>% 
      group_by(rid) %>% 
      mutate(minage = min(age)) %>% 
      filter(age == minage) %>% 
      ungroup() %>% 
      group_by(site) %>% 
      mutate(nscans = n_distinct(sub_id)) %>% 
      ungroup() %>% 
      filter(nscans >= thr.obs) 
    
    # filter df based on both condition
    sites = intersect(grot.all$site, grot.long$site)
    
    grot.long = 
      grot.long %>% 
      filter(site %in% sites)
    
    grot.all2 = 
      df.merge %>% 
      filter(site %in% sites)
    
    if  (is_empty(setdiff(unique(df.merge.long$site),sites))) {
      df.merge <<- grot.all2
      df.merge.long <<- grot.long
      i = 1
    } else {
      df.merge = grot.all2
      df.merge.long = grot.long
    }
  }
}


merge_normative_mri_data = function(outdir, analysis, phenotypes, thr.time, thr.mad, thr.miss=.2, globalVars = F) {
  sub_vars = c("sub_id", 
               "rid", 
               "dataset", 
               "age", 
               "sex", 
               "site", 
               "sitenum")
  
  load(file.path(outdir, "df.all.filt.Rda"))
  df = df.merge.long %>% 
    select(sub_vars) %>% 
    group_by(rid) %>% 
    mutate(xage = mean(age),
           agebsl = min(age),
           time = max(age)- min(age)) %>% 
    filter(time >thr.time)
  
  df.rid= 
    df %>% 
    group_by(rid) %>% 
    summarise(n = n_distinct(age),
              sex = first(sex),
              dataset = first(dataset),
              site = first(site),
              sites = n_distinct(site),
              time = first(time),
              xage = first(age))
  
  
  
  df.p = lapply(phenotypes.gamm, function(x) {import(file.path(outdir, analysis, x, "Z_predict_extended.Rda"))})
  names(df.p) = phenotypes.gamm
  df.p = df.p %>% data.table::rbindlist(idcol = "features")
  
  deltaZ = df.p %>% 
    dplyr::select(-c(meanZ, seZ)) %>% 
    pivot_wider(rid, 
                names_from = features, 
                values_from = deltaZ)
  
  meanZ = df.p %>% 
    dplyr::select(-c(deltaZ, seZ)) %>% 
    pivot_wider(rid, 
                names_from = features, 
                values_from = meanZ)
  
  
  seZ = df.p %>% 
    dplyr::select(-c(deltaZ, meanZ)) %>% 
    pivot_wider(rid, 
                names_from = features, 
                values_from = seZ)
  
  
  deltaZ = 
    left_join(
      df.rid %>% dplyr::select(rid), 
      deltaZ) %>% 
    dplyr::select(-rid)
  
  meanZ = 
    left_join(
      df.rid %>% dplyr::select(rid), 
      meanZ) %>% 
    dplyr::select(-rid)
  
  seZ = 
    left_join(
      df.rid %>% dplyr::select(rid), 
      seZ) %>% 
    dplyr::select(-rid)
  
  MD = colMeans(deltaZ, na.rm = T)
  madD = apply(deltaZ, 2, mad, na.rm = T)
  grot.delta = sweep(deltaZ, 2, MD, FUN = "-")
  grot.delta = sweep(grot.delta, 2, madD, FUN = "/")
  
  
  MM = colMeans(meanZ, na.rm = T)
  madM = apply(meanZ, 2, mad, na.rm = T)
  grot.mean = sweep(meanZ, 2, MM, FUN = "-")
  grot.mean = sweep(grot.mean, 2, madM, FUN = "/")
  
  
  delta.outlier = abs(grot.delta) < thr.mad
  mean.outlier = abs(grot.mean) < thr.mad
  
  
  delta.outlier[is.na(delta.outlier)] = F
  mean.outlier[is.na( mean.outlier)] = F
  # create output dataframe
  
  deltaZ[delta.outlier == F] = NaN
  meanZ[mean.outlier == F] = NaN
  
  df = list()
  df$outliers$delta$idx = delta.outlier
  df$outliers$mean$idx = mean.outlier
  df$outliers$delta$roi = colSums(delta.outlier == F)
  df$outliers$delta$rid = rowSums(delta.outlier == F)
  df$outliers$mean$roi = colSums(mean.outlier == F)
  df$outliers$mean$rid = rowSums(mean.outlier == F)
  df$outliers$delta$rmsubs.miss.data = df$outliers$delta$rid < length(phenotypes)*thr.miss
  df$outliers$mean$rmsubs.miss.data = df$outliers$mean$rid < length(phenotypes)*thr.miss
  df$outliers$thr.miss = thr.miss
  
  df$df$delta$df = deltaZ
  df$df$mean$df = meanZ
  df$df$base$df = df.rid
  df$df$delta$se = seZ
  
  df$misc$phenotypes = phenotypes
  if (globalVars == F) {
    df.ggeg = rename_variables_ggseg(phenotypes, df.harmonize)
    df$misc$phenotypes.ggseg = df.ggeg$phenotypes
    df$misc$phenotypes.atlas = df.ggeg$atlas
  }
  return(df)
}


rename_variables_ggseg = function(phenotypes, df.harmonize) { 
  library(ggseg)
  df.p = data.frame(phenotypes = phenotypes)
  gg.dest = unique(dk$data$label)
  gg.aseg = unique(aseg$data$label) 
  
  for (i in 1:length(phenotypes)) {
    library(ggseg)
    idx = grep(paste0("\\b",phenotypes[i],"\\b"),df.harmonize$normative_modelling)
    df.p$phenotypes[i] = df.harmonize$ggseg[idx]
    df.p$atlas[i] = ifelse(df.p$phenotypes[i] %in% gg.dest, "dk", 
                           ifelse(df.p$phenotypes[i] %in% gg.aseg, "aseg", "global"))
  }
  
  return(df.p)
}

