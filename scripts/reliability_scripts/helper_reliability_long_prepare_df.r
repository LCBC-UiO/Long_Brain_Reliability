mutate_ub = function(cohort, is.long = F) {
  if (is.long == T) {
    cohort =
      cohort %>% 
      mutate(Sess = if_else(`Round_ id` =="R01",1, 
                            if_else(`Round_ id` == "R02",2,3))) 
  } else {
    cohort = 
      cohort %>% 
      mutate(Sess = 1)   
  }
  return(cohort)
}


df.ub.change.names = function(df.mri) {
  df.mri$Subject_id[ df.mri$Subject_id == "hc001"] = "hc001_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc002"] = "hc002_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc003"] = "hc003_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc005"] = "hc005_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc010"] = "hc010_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc011"] = "hc011_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc014"] = "hc014_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc015"] = "hc015_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc016"] = "hc016_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc017"] = "hc017_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc018"] = "hc018_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc019"] = "hc019_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc020"] = "hc020_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc021"] = "hc021_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc023"] = "hc023_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc024"] = "hc024_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc027"] = "hc027_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc028"] = "hc028_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc029"] = "hc029_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc031"] = "hc031_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc032"] = "hc032_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc033"] = "hc033_nc059"
  df.mri$Subject_id[ df.mri$Subject_id == "hc035"] = "hc035_p"
  df.mri$Subject_id[ df.mri$Subject_id == "hc036"] = "hc036_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc039"] = "hc039_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc040"] = "hc040_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc041"] = "hc041_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc042"] = "hc042_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc043"] = "hc043_m"
  df.mri$Subject_id[ df.mri$Subject_id == "hc044"] = "hc044_nc060"
  df.mri$Subject_id[ df.mri$Subject_id == "nc052"] = "nc052_2"  
  df.mri$Subject_id[ df.mri$Subject_id == "001_MGS"] = "001_MSG"
  df.mri$Subject_id[ df.mri$Subject_id == "013_RFS"] = "013_RFJ"
  return(df.mri)
}

load_scan_info_adni = function(fname, field.strength) {
  
  mri.info = read.csv(file.path(sitedir, fname ))
  names(mri.info)[10] = "Image_ID"
  
  tmp =mri.info %>% 
    filter(Type =="Original") %>% 
    group_by(Subject.ID, 
             Visit,
             Imaging.Protocol) %>% 
    tally()
  
  mri.info =
    left_join(mri.info %>% 
                select(-Imaging.Protocol),
              tmp) %>% 
    mutate(field.strength = field.strength)
  return(mri.info)  
}

              "Mfg Model=Allegra"                        
#[6] "Mfg Model=Ingenuity"         "Mfg Model=Ingenia Elition X" "Mfg Model=MAGNETOM Vida"  
df.adni.set.scans = function() {
  triotim = c("Mfg Model=TrioTim",
              "Mfg Model=Trio")
  prisma = c("Mfg Model=Prisma_fit",
             "Mfg Model=Prisma")
  verio = c("Mfg Model=Verio")
  skyra = c("Mfg Model=Skyra",
            "Mfg Model=Skyra|DicomCleaner",
            "Mfg Model=Skyra_fit")
  achieva = c("Mfg Model=Achieva",
              "Mfg Model=Achieva dStream")
  ingenia = c("Mfg Model=Ingenia")
  signa = c("Mfg Model=Signa HDxt",
            "Mfg Model=SIGNA HDx",
            "Mfg Model=GENESIS_SIGNA",
            "Mfg Model=SIGNA EXCITE",
            "Mfg Model=SIGNA HDx",
            "Mfg Model=Signa HDxt",
            "Mfg Model=SIGNA Premier")
  discovery = c("Mfg Model=DISCOVERY MR750",
                "Mfg Model=DISCOVERY MR750w")
  intera = c("Mfg Model=Intera",
             "Mfg Model=Gyroscan Intera",
             "Mfg Model=Gyroscan NT",
             "Mfg Model=Intera",
             "Mfg Model=Intera Achieva")
  gemini = c("Mfg Model=GEMINI",
             " Mfg Model=Ingenuity")
  biographmMR = c("Mfg Model=Biograph_mMR")
  avanto=c("Mfg Model=Avanto")
  sonata=c("Mfg Model=Sonata",
           "Mfg Model=SonataVision",
           "Mfg Model=Espree")
  symphony=c("Mfg Model=Symphony",
             "Mfg Model=SymphonyTim",
             "Mfg Model=NUMARIS/4")
  allegra = c("Mfg Model=Allegra")
  models = list("triotim" = triotim,
                "prisma" = prisma,
                "verio" = verio,
                "skyra" = skyra, 
                "achieva" = achieva, 
                "ingenia" = ingenia, 
                "signa" = signa, 
                "discovery" = discovery, 
                "intera" = intera, 
                "gemini" = gemini, 
                "biographmMR" = biographmMR,
                "avanto" = avanto, 
                "sonata" = sonata, 
                "symphony" = symphony,
                "allegra" = allegra)
  
  return(models)
}
set_links = function(folder) {
  tabulateddata <<- here('data-raw/tabulated')
  sitedir <<- file.path(tabulateddata, folder)
  outdir <<- here('data_reliability_long','df_mri', folder)
  if (!dir.exists(outdir)) {dir.create(outdir)}
}

check.uio.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  df.uio = import(file.path(sitedir, "noas_query_2022-11-29_21-26-23_041b29e_b0f5e7b.csv"))
  df.uio = df.uio %>% filter(subject_shareable == 1)
  suppressWarnings({
    df.uio = 
      df.uio %>% 
      group_by(subject_id) %>% 
      mutate(mean_age = mean(visit_age),
             max_mms_score = max(mms_score, na.rm = T)) %>% 
      filter(mean_age > 20 & (max_mms_score > 25 | is.na(max_mms_score)))
  })
  
  grot.scanner = import(file.path(sitedir, "noas_query_2022-11-29_21-28-27_94ea888_b0f5e7b.csv"))
  grot.scanner.1 = grot.scanner %>% group_by(mri_info_folder) %>% summarise_if(is.numeric, mean)
  grot.scanner.2 = grot.scanner %>% group_by(mri_info_folder) %>% summarise_if(is.character, first)
  grot.scanner.3 = grot.scanner %>% group_by(mri_info_folder) %>% summarise(mri_info_date = first(mri_info_date))
  df.scanner = 
    left_join(grot.scanner.1, grot.scanner.2) %>% 
    left_join(., grot.scanner.3)
  
  df.scanner = 
    df.scanner %>% filter(subject_id %in% df.uio$subject_id)
  
  if (is.cross == T ) {
    idir = file.path(sitedir, "mri_fs7", "cross")
    lfiles = list.files(idir)
    subs =
      gsub("lh.", "", lfiles) %>% 
      gsub("rh.", "", .) %>% 
      strsplit(.,"\\.")  %>% 
      sapply(., `[`, 1) %>% unique()  
    
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
    df.mri.all = data.table::rbindlist(df.mri.all, idcol = "sub_id", fill = T) %>% 
      data.frame()
    
    df.mri = 
    df.mri.all %>% 
      mutate(subject_id = substr(input,1,7) %>% as.numeric(),
             visit_number = substr(input,9,10) %>% as.numeric(), 
             project = substr(input,12,13) %>% as.numeric(),
             wave = substr(input,15,16) %>% as.numeric(),
             scannerN = substr(input,18,19) %>% as.numeric(),
             scanner = if_else(scannerN == 11, "ousAvanto", 
                               if_else(scannerN == 12, "ousSkyra", 
                                       if_else(scannerN == 13, "ousPrisma", 
                                               if_else(scannerN == 20, "ntnuAvanto",NA ))))) %>% 
    filter(!is.na(scanner))
             
    } else {
    df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", "desikan", x))})
  
  
    for (df in 1:length(df.mri)) {
      names(df.mri[[df]])[1] <- "input"
    }
    df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
    
    df.mri =
      df.mri %>% 
      mutate(subject_id = substr(input, 5,11) %>% 
               as.numeric(),
             visit_number = substr(input, 16,16) %>% 
               as.numeric(),
             scanner = substr(input,17,30)) %>% 
      separate(scanner, into = "scanner", sep = ".lo")
  }
  
  df.out = inner_join(df.uio, df.mri)
  
  # create mri-to-database file
  if (is.cross == T ) {
    save(df.out, 
         df.scanner, 
         df.uio, 
         file = file.path(outdir, "df.cross.rda")) 
  } else {
    save(df.out, 
         df.scanner, 
         df.uio, 
         file = file.path(outdir, "df.rda")) 
  }
  
  cat("uio data checked. rda file save with valid data only")
}

check.ucam.fs = function(sitedir, outdir, mrifiles) {
  ## load non-IMAGING data
  tmp.tp1 = read_tsv(file.path(sitedir, "identifier_ses-01_p029.tsv"), na = c("", "NA", "NaN") )
  tmp.tp2 = read_tsv(file.path(sitedir, "identifier_ses-02_p029.tsv"), na = c("", "NA", "NaN") )
  df.camcan = rbind(tmp.tp1, tmp.tp2) %>% filter(!is.na(Age))
  df.camcan = 
    df.camcan %>% 
    group_by(Subject_id) %>% 
    mutate(min_Age = min(Age)) %>% 
    filter(min_Age > 18 & Flag_hasFreeSurfer == 1) 
  
  grot.df = import(file.path(sitedir, "subject.tsv")) %>% 
    rename("Subject_id" = "subject_id")
  
  df.camcan = inner_join(df.camcan, grot.df)
  
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", "desikan", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri =
    df.mri %>% 
    mutate(Subject_id = substr(input, 5,12),
           SubjectRound_id = substr(input, 17,17) %>% as.numeric())
  # missing subjects
  df.out = inner_join(df.camcan, df.mri)
 
  save(df.out, 
       df.camcan,
       file = file.path(outdir, "df.rda")) 
  cat("ucam data checked. rda file save with valid data only")
  
}


check.mpib.fs = function(sitedir, outdir, mrifiles) {
  ## load non-IMAGING data
  df.MPIB = read.csv(file.path(sitedir,"Data_BASEII_YlvaEniko211011.csv"))
  df.id.mpib = df.MPIB %>%
    mutate(ses01 = if_else(!is.na(MRT12_age), 1,NaN),
           ses02 = if_else(!is.na(MRT16_age), 2,NaN)) %>% 
    #select(ID, ses01, ses02) %>% 
    pivot_longer(cols = c(ses01,ses02), 
                 names_to = "grot1",
                 values_to = "SubjectRound_id",
                 values_drop_na = T) %>% 
    mutate(visit_age = if_else(SubjectRound_id ==1, MRT12_age, MRT16_age)) %>% 
    select(-c(Age, 
              MPIB_2012_Date1, 
              MRT12_age, 
              MRT16_age)) %>% 
    filter(MMSE > 25)
  
  
  # check if fs data
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", "desikan", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  
  df.mri =
    df.mri %>% 
    mutate(ID = substr(input, 5,11),
           SubjectRound_id = substr(input, 15,16) %>% as.numeric())
  
  # missing subjects
  df.out = left_join(df.id.mpib, df.mri)
  
  save(df.out, 
       df.id.mpib, 
       file = file.path(outdir, "df.rda"))
  
  
  cat("mpib data checked. rda file save with valid data only")
  
}

check.umu.fs = function(sitedir, outdir, mrifiles) {
  # load databases and get mri_filepaths
  
  grot.t5 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T5.xlsx"))
  grot.t6 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T6.xlsx"))
  grot.t7 = readxl::read_xlsx(file.path(sitedir,  "B_LB_p029_new_T7.xlsx"))
  
  db.exclusion = 
    import(file.path(sitedir, "Exclusions_Betula_Imaging_Sample_cleaned.xlsx"))
  
  db.exclusion = 
    c(paste0(db.exclusion  %>% 
               filter(`Exclude T5` == 1) %>% 
               .$`Subject code`, 
             "_T5"),
      paste0(db.exclusion  %>% 
               filter(`Exclude T6` == 1) %>% 
               .$`Subject code`, 
             "6"),
      paste0(db.exclusion  %>% 
               filter(`Exclude T7` == 1) %>% 
               .$`Subject code`, 
             "7"))
  
  
  df.Umea = data.table::rbindlist(list(grot.t5, 
                                  grot.t6, 
                                  grot.t7), 
                             use.names = F) %>% 
    rowwise() %>% 
    mutate(age = mean(c(rounded_age_HT, rounded_age_MT), na.rm = T),
           input = paste(Lifebrain_Subject_id, Lifebrain_StudyRound_id, sep = "_T")) %>% 
    filter(!`Demens::dementiaStatusAtEvaluationDate_binary` == 1 & !is.na(age) & !input %in% db.exclusion)
  
  # check if fs data
  df.mri = readxl::read_xlsx(file.path(sitedir, 
                            "mri_fs7_JamesProc", 
                            "FreeSurf711_LongPipe_VolumeAreaThickness_20210423.xlsx"))
  
  df.mri =
    df.mri %>% 
    rename("Lifebrain_Subject_id" = "subject") %>% 
    mutate(Lifebrain_StudyRound_id = if_else(TimePoint == "T5", 5, 
                                             if_else(TimePoint == "T6", 6,7)),
           input = paste(Lifebrain_Subject_id, Lifebrain_StudyRound_id, sep = "_T")) 
  df.mri = 
    df.mri %>% 
    filter(!is.na(`aseg:Left-Lateral-Ventricle.Volume_mm3`))
  # missing subjects
  df.out = inner_join(df.Umea, df.mri)
  

  
  save(df.out, 
       df.Umea,
       file = file.path(outdir, "df.rda")) 
  cat("umu data checked. rda file save with valid data only")
}

check.ub.fs = function(sitedir, outdir, mrifiles, is.cross = T) {
  # load databases and get mri_filepaths
  grot.PDcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                    sheet= "PDcohort")
  grot.WAHAcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                      sheet= "WAHAcohort")
  grot.MSAcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                     sheet= "MSAcohort")
  grot.CRcohort = readxl::read_xlsx(file.path(sitedir,  "UB_Depression_late_Eniko_May2021.xlsx"), 
                                    sheet= "CRcohort")
  grot.GABAcohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                      sheet= "GABAcohort")
  grot.iTBScohort = readxl::read_xlsx(file.path(sitedir, "UB_Depression_late_Eniko_May2021.xlsx"), 
                                      sheet= "iTBScohort")

  grot.WAHAcohort = grot.WAHAcohort %>% 
    rename(c("calculated_age_MRI" = "calculated_age",
             "Sex" = `Sex (1=woman)`,
             "Other_Praxis_Executive_ROCF_Copy" = "Praxis_Executive_ROCF_Copy")) %>% 
    mutate(Subject_id = as.character(Subject_id))
  grot.PDcohort = grot.PDcohort %>% 
    rename(c("calculated_age_MRI" = "calculated_age",
             "Sex" = `Sex (1=woman)`))

    
  grot.PDcohort  = mutate_ub(grot.PDcohort, T)
  grot.WAHAcohort  = mutate_ub(grot.WAHAcohort, T)
  grot.MSAcohort  = mutate_ub(grot.MSAcohort)
  grot.CRcohort = mutate_ub(grot.CRcohort)
  grot.GABAcohort = mutate_ub(grot.GABAcohort)
  grot.iTBScohort = mutate_ub(grot.iTBScohort)
  
  suppressMessages({
  df.UB = 
  list(grot.PDcohort, 
       grot.WAHAcohort, 
       grot.MSAcohort, 
       grot.CRcohort, 
       grot.GABAcohort, 
       grot.iTBScohort) %>% reduce(full_join)
  })
  df.UB = df.UB %>% 
    mutate(MMSE = if_else(MMSE < 0,NaN, MMSE)) %>% 
  drop_na(any_of(c("Subject_id",
               "Sex", 
               "Sess", 
               "calculated_age_MRI")))
  

  # # check if fs data
  if (is.cross == T ) {
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
    separate(input, 
             c("Subject_id", "sess"), 
             sep = "ses0", 
             remove = F) %>% 
    mutate(Subject_id = substr(Subject_id, 5, 30),
           Sess = substr(sess, 1,1) %>% as.numeric()) 
  
  # change ID strings for CR cohort
  stringi::stri_sub(df.mri$Subject_id[df.mri$Subject_id %>% startsWith("CTR")],7,6) <- "_"
  stringi::stri_sub(df.mri$Subject_id[grepl("^[[:digit:]]+", df.mri$Subject_id) & !grepl("[[:digit:]]+$", df.mri$Subject_id)],4,3)  <- "_"
  df.mri = df.ub.change.names(df.mri)
  
  
  df.out = inner_join(df.UB, df.mri)
 
  if (is.cross == T ) {
    save(df.out, 
         df.UB,
         file = file.path(outdir, "df.cross.rda")) 
  } else {
    save(df.out, 
         df.UB,
         file = file.path(outdir, "df.rda")) 
  }
  
  cat("ub data checked. linkage file created")
}

check.vumc.fs = function(sitedir, outdir, mrifiles) {
  df.nesda = import(file.path(sitedir, "220208_NESDA_Lifebrain_lateonsetdepr_long_feb22.xlsx"))
  df.nesda = 
    df.nesda %>%  
    group_by(Subject_ID) %>% 
    mutate(mean_age = mean(age)) %>% 
    filter(mean_age > 20 & !Subject_ID == 120441)
  #         scan_avail == 1) # check quality of scans
  
  
  # # check if fs data
  df.mri = lapply(mrifiles, function(x) { import(file.path(sitedir, "mri_fs7", "desikan", x))})
  for (df in 1:length(df.mri)) {
    names(df.mri[[df]])[1] <- "input"
  }
  df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
  
  df.mri =
    df.mri %>% 
    mutate(Subject_ID = substr(input, 5,10) %>% as.numeric(),
           timepoint = substr(input, 15,15) %>% 
             as.numeric()) 
  
  
  df.out = inner_join(df.nesda, df.mri) %>% 
    filter(scan_avail == 1)
  
  save(df.out, 
       df.nesda,
       file = file.path(outdir, "df.rda")) 
  
  cat("vumc data checked. linkage file created")
}


check.habs.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  ## load non-IMAGING data
  df.habs = import(file.path(sitedir, "Demographics_HABS_DataRelease_2.0.csv"))
  names(df.habs)[1] = "SubjIDshort"
  df.grot1 =import(file.path(sitedir, "ClinicalMeasures_HABS_DataRelease_2.0.csv"))
  df.grot2 =import(file.path(sitedir, "Cognition_HABS_DataRelease_2.0.csv"))
  df.grot3 =import(file.path(sitedir, "PACC_HABS_DataRelease_2.0.csv"))
  #df.grot4 =import(file.path(sitedir, "ADNI_MRI_FS6_XSec_HABS_DataRelease_2.0.csv"))
  x = left_join(df.habs, df.grot1)
  xx = left_join(df.grot2, df.grot3)
  
  df.habs = left_join(x, xx) %>% 
    rename("SubjID" = "SubjIDshort")
    
  df.habs = 
  df.habs %>% 
    group_by(SubjID) %>% 
    mutate(AgeBsl = min(NP_Age),
           SubjectRound_id = if_else(StudyArc == "HAB_1.0", 0,
                                     if_else(StudyArc == "HAB_2.0", 12,
                                             if_else(StudyArc == "HAB_3.0", 24,
                                                     if_else(StudyArc == "HAB_4.0", 36,
                                                             if_else(StudyArc == "HAB_5.0", 48,60)))))) 
  
  df.habs = 
    df.habs %>% 
    filter(HABS_DX == "CN")
  
  
  if (is.cross == T ) {
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
    mutate(input = gsub("./", "", input),
           SubjID = paste0("P_", substr(input, 9,14)),
           SubjectRound_id = substr(input, 19,20) %>% as.numeric(), 
           SubjectRound_id = if_else(SubjectRound_id == 18, 24,SubjectRound_id)) 
  # use for merging purposes use MR age for modeling. 
  
 # missing subjects
  df.out = inner_join(df.habs, df.mri)
  
  if (is.cross == T ) {
    save(df.out, 
         df.habs,
         file = file.path(outdir, "df.cross.rda")) 
  } else {
    save(df.out, 
         df.habs,
         file = file.path(outdir, "df.rda")) 
  }
  
  
  cat("habs data checked. rda file save with valid data only")
}


check.ous.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  
  load(file.path(sitedir, "cognNBM_updated.Rda"))
  df.scan.date = import(file.path(sitedir, "scan.date.csv")) %>% 
    separate(acq_time, c("scan_date","tod"), sep = " ") %>% 
    mutate(Subject_Timepoint = substr(ses, 6,6) %>% as.character(),
           scanner = substring(ses, 7), 
           ID = substring(sub,5)) 
  df.scan.date = 
    df.scan.date %>% group_by(sub, ses, Subject_Timepoint, scanner, ID) %>% 
    summarise(scan_date = first(scan_date), tod = first(tod)) %>% ungroup()
  
  df.out = 
    cognNBM_updated %>% 
    group_by(ID) %>% 
    summarise(Birth_Date = first(Birth_Date), 
              Sex = first(Sex), 
              Edu_Years = first(Edu_Years))
  df.out = left_join(df.out, df.scan.date)
  
  if (is.cross == T ) {
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
    mutate(input = gsub("./", "", input), 
           ID = substr(input, 5,11) %>% 
             as.character(),
           Subject_Timepoint = substr(input, 16,16) %>% 
             as.character(),
           scanner = substr(input,17,30)) %>% 
    separate(scanner, into = "scanner", sep = ".l")
  
  # missing subjects
  df.out = inner_join(df.out, df.mri)
  df.ous = cognNBM_updated
  
  df.out = 
    df.out %>% mutate(MRI_Age = (as.Date(scan_date)- as.Date(Birth_Date))/365.25 %>% as.numeric())
  class(df.out$MRI_Age) = "numeric"
  
  
  
  if (is.cross == T ) {
    save(df.out, 
         df.ous,
         file = file.path(outdir, "df.cross.rda")) 
  } else {
    save(df.out, 
         df.ous,
         file = file.path(outdir, "df.rda")) 
  }
  
  cat("ous data checked. rda file save with valid data only")
}

check.aibl.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  ## load non-IMAGING data
  
  df.aibl.demog =import(file.path(sitedir, "aibl_ptdemog_01-Jun-2018.csv")) %>% 
    select(-VISCODE)
  df.aibl.mmse =import(file.path(sitedir, "aibl_mmse_01-Jun-2018.csv"))  %>% 
    select(-EXAMDATE)
  df.aibl.cdr =import(file.path(sitedir, "aibl_cdr_01-Jun-2018.csv"))  %>% 
    select(-EXAMDATE)
  df.aibl.cogn =import(file.path(sitedir, "aibl_neurobat_01-Jun-2018.csv"))  %>% 
    select(-EXAMDATE)
  df.aibl.dx = import(file.path(sitedir, "aibl_pdxconv_01-Jun-2018.csv"))  
  
  #df.aibl.visits = import(file.path(sitedir, "aibl_visits_01-Jun-2018.csv"))
  df.grot.mri =import(file.path(sitedir, "aibl_mri3meta_01-Jun-2018.csv")) %>% 
    mutate(Tscan = 3)
  df.grot.mri3 =import(file.path(sitedir, "aibl_mrimeta_01-Jun-2018.csv")) %>% 
    mutate(Tscan = 1.5)
  df.aibl.mri =
    rbind(df.grot.mri, 
          df.grot.mri3) %>% 
    filter(MMSMPRAGE == 1) %>% 
    select(-EXAMDATE)
  
  df.aibl = purrr::reduce(
    list(df.aibl.demog, 
         df.aibl.mmse, 
         df.aibl.cdr, 
         df.aibl.cogn, 
         df.aibl.dx,
         df.aibl.mri), 
    left_join)
  
  df.aibl = df.aibl %>% 
    mutate(visit = if_else(VISCODE == "bl", 0, 
                           if_else(VISCODE == "m18", 18, 
                                   if_else(VISCODE == "m36", 36, 
                                           if_else(VISCODE == "m54", 54, 
                                                   if_else(VISCODE == "m72", 72, NaN))))))
  
  df.aibl = df.aibl %>%
    group_by(RID) %>% 
    mutate(time.norm =if_else(DXCURREN == 1, visit, NaN), 
           time.norm = max(time.norm, na.rm = T)) %>%
    ungroup() %>% 
    mutate(crit.max.normal = if_else(visit <= time.norm & !is.nan(time.norm) & !is.infinite(time.norm), 1,0))  %>% 
    filter(crit.max.normal == 1)
  
  # df.aibl = df.aibl %>% 
  #   filter(DXCURREN == 1 & MMSMPRAGE == 1)
  
  df.aibl.idaSearch = import(file.path(sitedir, "idaSearch_12_05_2022.csv")) %>% 
    separate(`Imaging Protocol`, 
             c("grot1","field_strength", "grot3", "model", "grot4", "slice_thicknes"), 
             sep="=|;") %>% 
    mutate(Session_ID = if_else(Visit =="Baseline", "bl",
                                if_else(Visit =="18 Month follow-up", "m18",
                                        if_else(Visit =="36 Month follow-up", "m36",
                                                if_else(Visit =="54 Month follow-up", "m54",
                                                        if_else(Visit =="72 Month follow-up", "m72","")))))) %>% 
    rename("Subjects_ID" = "Subject ID") %>% 
    select(-c(contains("grot"),
              Description, 
              Sex,
              Visit)) 
  
  df.aibl.idaSearch = 
    df.aibl.idaSearch %>% 
    group_by(Subjects_ID, Session_ID) %>% 
    summarise(field_strength = first(field_strength), 
              slice_thicknes = first(slice_thicknes), 
              Session_ID = first(Session_ID), 
              model = first(model), 
              Age = mean(Age, na.rm = T))
  
  df.aibl.paths = import(file.path(sitedir, "t1_paths_aibl.tsv")) %>% 
    mutate(Session_ID = if_else(Session_ID == "M00", "bl", Session_ID)) %>% 
    group_by(Subjects_ID, Session_ID) %>% 
    tally()
  
  
  df.aibl.paths = 
    inner_join(df.aibl.paths, df.aibl.idaSearch) %>% 
    rename(c("VISCODE" = "Session_ID", 
             "RID" = "Subjects_ID")) 
  
  df.aibl = left_join(df.aibl, df.aibl.paths)
  
  if (is.cross == T ) {
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
    mutate(grot1 = gsub("sub-AIBL", "", input), 
           grot1 = gsub("./", "", grot1)) %>% 
    separate(grot1, c("RID",
                      "grot2"),
             sep = "sesM") %>% 
    mutate(visit = substr(grot2, 1,2) %>% as.numeric(),
           RID = as.numeric(RID))
  
  
  # missing subjects
  df.out = inner_join(df.aibl, df.mri)
  
  if (is.cross == T) {
    save(df.out, 
         df.aibl,
         file = file.path(outdir, "df.cross.rda")) 
    
  } else {
    save(df.out, 
         df.aibl,
         file = file.path(outdir, "df.rda")) 
  }
  cat("aibl data checked. rda file save with valid data only")
  
}

check.adni.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  #ADNI.stripped = read.csv(file.path(sitedir, "ADNI_stripped.csv"))
  df.adni.participants =import(file.path(sitedir, "participants.tsv")) %>% 
    rename("RID" = "adni_rid") %>%
    mutate(RID = as.numeric(RID))
  df.adni.conversioninfo =import(file.path(sitedir, "t1_paths_new.tsv")) 
  
  load(file.path(sitedir, "adni.sessions.Rda"))
  df.adni.sessions[df.adni.sessions == "n/a"] = NA 
  idx = !apply (is.na(df.adni.sessions), 2, all)
  df.adni.sessions = df.adni.sessions[, ..idx ]
  df.adni = inner_join(df.adni.sessions, df.adni.participants)
  
  mri.info15 = load_scan_info_adni("idaSearch_1_25_2022_1.5T.csv", "1.5")
  mri.info3 = load_scan_info_adni("idaSearch_1_25_2022_3T.csv", "3.0")
  
  mri.info = rbind(mri.info15, 
                   mri.info3)
  
  df.adni.mri = 
    inner_join(df.adni.conversioninfo,mri.info) %>% 
    mutate(Subject.ID = as.character(Subject.ID)) %>% 
    separate(Subject.ID, c("Site", "ID"), sep = "_S_", remove = F)
  
  df.adni.mri$model = ""
  models = df.adni.set.scans()
  
  for (x in 1:length(models)) {
    xx = names(models)[[x]]
    df.adni.mri = 
      df.adni.mri %>% 
      mutate(model = if_else(Imaging.Protocol %in% models[[x]], xx, model))
  }
  
  df.adni.mri = 
    df.adni.mri %>% 
    mutate(RID = ID %>% as.numeric(),
           session_id = gsub("m", "ses-M", VISCODE),
           session_id = if_else(session_id =="bl", "ses-M00", session_id))
  
  
  df.adni = left_join(df.adni,df.adni.mri)
  
  
  
  if (is.cross == T ) {
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
    mutate(grot1 = gsub("sub-ADNI", "", input),
           grot1 = gsub("./", "", grot1)) %>% 
    separate(grot1, c("grot4",
                      "grot2"),
             sep = "ses") %>% 
    separate(grot2, 
             c("session_id", "grot3"), 
             sep = ".long") %>%
    separate(grot4, 
             c("grot5", "RID"), 
             sep = "S") %>% 
    mutate(RID = as.numeric(RID),
           session_id = paste("ses", session_id, sep = "-")) %>% 
    select(-starts_with("grot"))
  
  df.out = inner_join(df.adni, df.mri) %>% 
    filter(diagnosis == "CN")
  
  if (is.cross == T ) {
    save(df.out, 
         df.adni,
         file = file.path(outdir, "df.cross.rda")) 
    
  } else {
    save(df.out, 
         df.adni,
         file = file.path(outdir, "df.rda")) 
    
  }
  cat("adni data checked. rda file save with valid data only")
}

check.preventad.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  datadir = file.path(sitedir, "2022_11_23", "tabular")
  
  df.preventAD.main = import(file.path(sitedir, "data-2022-12-11T20_48_55.373Z.csv")) %>% 
    rename(c("CONP_ID" = "PSCID",
             "CONP_CandID" = "DCCID",
             "Study_visit_label" = "Visit Label"))
  
  df.preventAD.demographics = import(file.path(datadir, "Demographics_Registered_PREVENTAD.csv"))
  df.preventAD.genetics = import(file.path(datadir, "Genetics_Registered_PREVENTAD.csv"))
  df.preventAD.DKEFS = import(file.path(datadir, "Neuropsych_DKEFS-CWIT_Registered_PREVENTAD.csv"))
  df.preventAD.RAVLT = import(file.path(datadir, "Neuropsych_RAVLT_Registered_PREVENTAD.csv"))
  df.preventAD.TMT = import(file.path(datadir, "Neuropsych_TMT_Registered_PREVENTAD.csv" ))
  df.preventAD.RBANS = import(file.path(datadir, "RBANS_Registered_PREVENTAD.csv" ))
  df.preventAD.AD8 = import(file.path(datadir, "AD8_Registered_PREVENTAD.csv" ))
  df.preventAD.APS = import(file.path(datadir, "APS_Registered_PREVENTAD.csv" ))
  df.preventAD.CDR = import(file.path(datadir, "EL_CDR_MoCA_Registered_PREVENTAD.csv" ))
  df.preventAD.MedHist = import(file.path(datadir, "EL_Medical_history_Registered_PREVENTAD.csv"))                                                   
  df.preventAD.Reg = import(file.path(datadir, "Lab_Registered_PREVENTAD.csv"))
  #df.preventAD.siblings = import(file.path(datadir, "List_of_participants_with_only_1_sibling.txt"))
  #df.preventAD.Switch = import(file.path(datadir, "List_of_participants_switched_back_to_cohort.txt"))
  
  
  df.preventAD.main = 
    inner_join(df.preventAD.demographics, 
               df.preventAD.main)
  df.preventAD.main$age_at_mci = NaN
  df.preventAD.main$age_at_mci[df.preventAD.main$probable_MCI_visit == df.preventAD.main$Study_visit_label ] = df.preventAD.main$`Age At MRI In Months`
  
  
  df.preventAD.main = 
    df.preventAD.main %>% 
    group_by(CONP_ID) %>% 
    mutate(age_at_mci = min(age_at_mci, na.rm = T), 
           age = `Age At MRI In Months`/12,
           mci_visit = if_else(age_at_mci <= `Age At MRI In Months`, 1, 0))
  
  
  df.preventAD = full_join(df.preventAD.Reg, df.preventAD.main)
  
  df.preventAD = 
    reduce(
      list(df.preventAD, 
           df.preventAD.AD8, 
           df.preventAD.APS, 
           df.preventAD.CDR, 
           df.preventAD.DKEFS,
           #df.preventAD.genetics, 
           df.preventAD.MedHist, 
           df.preventAD.RAVLT, 
           df.preventAD.RBANS, 
           df.preventAD.TMT),
      left_join,
      by = c("CONP_ID", "CONP_CandID",  "Visit_label")
    )
  
  df.preventAD = 
    left_join(df.preventAD,
              df.preventAD.genetics)
  
  
  if (is.cross == T ) {
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
    mutate(grot1 = gsub("sub-", "", input),
           grot1 = gsub("./", "", grot1)) %>% 
    separate(grot1, c("CONP_CandID",
                      "grot2"),
             sep = "ses") %>% 
    separate(grot2, 
             c("Study_visit_label.x", "grot3"), 
             sep = ".long") %>% 
    mutate(CONP_CandID = as.numeric(CONP_CandID)) %>% 
    select(-starts_with("grot"))
  
  df.out = inner_join(df.preventAD, df.mri) 
  
  if (is.cross == T ) {
    save(df.out, 
         df.preventAD,
         file = file.path(outdir, "df.cross.rda")) 
    } else {
    save(df.out, 
         df.preventAD,
         file = file.path(outdir, "df.rda")) 
  }
  
  cat("preventad data checked. rda file save with valid data only")
}


check.wayne.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  
  df.wayne.list = list()
  df.wayne.list[[1]] = import(file.path(sitedir, "Wayne_longitudinal_4wave_1.5T_Siemens_Study10.xlsx"), sheet = "wave1")
  df.wayne.list[[2]] = import(file.path(sitedir, "Wayne_longitudinal_4wave_1.5T_Siemens_Study10.xlsx"), sheet = "wave2") %>% 
    mutate(AGE = as.numeric(AGE))
  df.wayne.list[[3]] = import(file.path(sitedir, "Wayne_longitudinal_4wave_1.5T_Siemens_Study10.xlsx"), sheet = "wave3")
  df.wayne.list[[4]] = import(file.path(sitedir, "Wayne_longitudinal_4wave_1.5T_Siemens_Study10.xlsx"), sheet = "wave4")
  df.wayne.1 = reduce(df.wayne.list, full_join) %>% 
    mutate(SEX = if_else(!is.na(`SEX (1=M, 0=F)`), `SEX (1=M, 0=F)`, `SEX(1=M, 0=F)`),
           scanner = "1.5T") %>% 
    select(-c(`SEX (1=M, 0=F)`, `SEX(1=M, 0=F)`))
  
  df.wayne.list = list()
  df.wayne.list[[1]] = import(file.path(sitedir, "Wayne_Longitudinal_2wave_4T_Siemens_Study11.xlsx"), sheet = "wave1")
  df.wayne.list[[2]] = import(file.path(sitedir, "Wayne_Longitudinal_2wave_4T_Siemens_Study11.xlsx"), sheet = "wave2") 
  df.wayne.2 = reduce(df.wayne.list, full_join) %>% 
    rename("SEX" = "SEX (1=M, 0=F)") %>% 
    mutate(scanner = "4.0T") 
  
  df.wayne = rbind(df.wayne.1, df.wayne.2)
  
  
  
  if (is.cross == T ) {
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
    separate(input, c("PATID",
                      "grot2"),
             sep = "sesWave",
             remove = F) %>% 
    separate(grot2, 
             c("wave", "grot3"), 
             sep = ".long") %>% 
    select(-starts_with("grot")) %>% 
    mutate(wave = as.numeric(wave), 
           PATID = gsub("./", "", PATID))
  
  
  # missing subjects
  df.out = inner_join(df.wayne, df.mri) 
  
  if (is.cross == T ) {
    save(df.out, 
         df.wayne,
         file = file.path(outdir, "df.cross.rda")) 
  } else {
    save(df.out, 
         df.wayne,
         file = file.path(outdir, "df.rda")) 
  }
  
  cat("wayne data checked. rda file save with valid data only")
}


check.ukb.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  # load databases and get mri_filepaths
  #list.files(file.path(sitedir, "45249", "normative_modelling"))
  dfdir=file.path(sitedir, "45249", "normative_modelling")
  
  load(file.path(dfdir, "basic_demographics.Rda"))
  #load(file.path(dfdir, "mri_data_tp23.Rda"))
  load(file.path(dfdir, "demographics_expanded.Rda"))
  
  df.ukb = df.dmg %>% filter(!is.na(age_3)) %>% 
    group_by(eid) %>% 
    select(-(c(age_0, age_1))) %>% 
    mutate(agebsl = min(c(age_2, age_3)))
  
  df.ukb = 
    left_join(df.ukb, 
              df.dmg.extra %>% 
                select(
                  eid,
                  uk_biobank_assessment_centre_f54_2_0,
                  uk_biobank_assessment_centre_f54_3_0
                ))
  grot.idx = grepl("age_",names(df.ukb))
  names(df.ukb)[grot.idx] = paste0(names(df.ukb)[grot.idx], "_0")
  
  df.ukb = 
    df.ukb %>% 
    mutate(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0),
           t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0),
           uk_biobank_assessment_centre_f54_2_0 = as.numeric(uk_biobank_assessment_centre_f54_2_0),
           uk_biobank_assessment_centre_f54_3_0 = as.numeric(uk_biobank_assessment_centre_f54_3_0)) %>% 
    pivot_longer(-c(eid, Sex, agebsl), 
                 names_to = "names", 
                 values_to = "values") %>% 
    mutate(grot1 = stringi::stri_reverse(names)) %>% 
    separate(grot1, 
             c("grot2", "tp", "grot3"),
             sep ="_",
             extra = "merge") %>% 
    mutate(names = stringi::stri_reverse(grot3)) %>% 
    select(-starts_with("grot")) %>% 
    pivot_wider(names_from = names, 
                values_from = values)
  
 
  if (is.cross == T ) {
    df.mri = lapply(mrifiles, function(x) { import(file.path(dfdir, "tsd_processed_fs.7.1.0", "cross", x))})
    for (df in 1:length(df.mri)) {
      names(df.mri[[df]])[1] <- "input"
    }
    df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
    
    df.mri =
      df.mri %>% 
      mutate(input = gsub("./", "", input),
             eid = substr(input, 5,11) %>% as.integer(),
             tp = substr(input, 15,15)) 
    
  } else {
    df.mri = lapply(mrifiles, function(x) { import(file.path(dfdir, "tsd_processed_fs.7.1.0", "desikan", x))})
    for (df in 1:length(df.mri)) {
      names(df.mri[[df]])[1] <- "input"
    }
    df.mri = purrr::reduce(df.mri, dplyr::left_join, by = "input")
    
    df.mri =
      df.mri %>% 
      mutate(eid = substr(input, 5,11) %>% as.integer(),
             tp = substr(input, 15,15)) 
    
  }
  
  #df.mri = lapply(mrifiles, function(x) { import(file.path(dfdir, "tsd_processed_fs.7.1.0", "desikan", x))})

  df.out = inner_join(df.ukb, df.mri)
  
  
  if (is.cross == T ) {
    save(df.out, 
         df.ukb,
         file = file.path(outdir, "df.cross.rda")) 
  } else {
    save(df.out, 
         df.ukb,
         file = file.path(outdir, "df.rda")) 
  }
  
  cat("ukb data checked. rda file save with valid data only")
}

check.ukb.fs_showcase = function(sitedir, outdir, mrifiles) {
#deprecated - due to use of cross-sectional only pipeline
# load databases and get mri_filepaths
#list.files(file.path(sitedir, "45249", "normative_modelling"))
dfdir=file.path(sitedir, "45249", "normative_modelling")
  
load(file.path(dfdir, "basic_demographics.Rda"))
load(file.path(dfdir, "mri_data_tp23.Rda"))
load(file.path(dfdir, "demographics_expanded.Rda"))

df.ukb = df.dmg %>% filter(!is.na(age_3)) %>% 
  group_by(eid) %>% 
  select(-(c(age_0, age_1))) %>% 
  mutate(agebsl = min(c(age_2, age_3)))

df.ukb = 
  left_join(df.ukb, 
            df.dmg.extra %>% 
              select(
                eid,
                uk_biobank_assessment_centre_f54_2_0,
                uk_biobank_assessment_centre_f54_3_0
              ))
grot.idx = grepl("age_",names(df.ukb))
names(df.ukb)[grot.idx] = paste0(names(df.ukb)[grot.idx], "_0")

df.ukb = 
  df.ukb %>% 
  mutate(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0),
         t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0),
         uk_biobank_assessment_centre_f54_2_0 = as.numeric(uk_biobank_assessment_centre_f54_2_0),
         uk_biobank_assessment_centre_f54_3_0 = as.numeric(uk_biobank_assessment_centre_f54_3_0)) %>% 
  pivot_longer(-c(eid, Sex, agebsl), 
               names_to = "names", 
               values_to = "values") %>% 
  mutate(grot1 = stringi::stri_reverse(names)) %>% 
  separate(grot1, 
           c("grot2", "tp", "grot3"),
           sep ="_",
           extra = "merge") %>% 
  mutate(names = stringi::stri_reverse(grot3)) %>% 
  select(-starts_with("grot")) %>% 
  pivot_wider(names_from = names, 
              values_from = values)

df.mri = 
  df.mri %>% 
  filter(eid %in% df.ukb$eid)
df.mri = 
  df.mri %>% 
  pivot_longer(-eid, 
               names_to = "names", 
               values_to = "values") %>% 
  mutate(grot1 = stringi::stri_reverse(names)) %>% 
  separate(grot1, 
           c("grot2", "tp", "grot3"),
           sep ="_",
           extra = "merge") %>% 
  mutate(names = stringi::stri_reverse(grot3)) %>% 
  select(-starts_with("grot")) %>% 
  pivot_wider(names_from = names, 
              values_from = values)


df.out = inner_join(df.ukb, df.mri)

save(df.out, 
     df.ukb,
     file = file.path(outdir, "df.rda")) 
cat("ukb data checked. rda file save with valid data only")
}


check.ukb.fs_showcase_link = function(sitedir, outdir, mrifiles) {
  #deprecated - due to use of cross-sectional only pipeline
  # load databases and get mri_filepaths
  #list.files(file.path(sitedir, "45249", "normative_modelling"))
  dfdir=file.path(sitedir, "45249", "normative_modelling")
  
  load(file.path(dfdir, "basic_demographics.Rda"))
  load(file.path(dfdir, "mri_data_tp23.Rda"))
  load(file.path(dfdir, "demographics_expanded.Rda"))
  
  df.ukb = df.dmg %>% filter(!is.na(age_3)) %>% 
    group_by(eid) %>% 
    select(-(c(age_0, age_1))) %>% 
    mutate(agebsl = min(c(age_2, age_3)))
  
  df.ukb = 
    left_join(df.ukb, 
              df.dmg.extra %>% 
                select(
                  eid,
                  uk_biobank_assessment_centre_f54_2_0,
                  uk_biobank_assessment_centre_f54_3_0
                ))
  grot.idx = grepl("age_",names(df.ukb))
  names(df.ukb)[grot.idx] = paste0(names(df.ukb)[grot.idx], "_0")
  
  df.ukb = 
    df.ukb %>% 
    mutate(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_2_0),
           t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0 = as.numeric(t2flair_used_in_addition_to_t1_to_run_freesurfer_f26500_3_0),
           uk_biobank_assessment_centre_f54_2_0 = as.numeric(uk_biobank_assessment_centre_f54_2_0),
           uk_biobank_assessment_centre_f54_3_0 = as.numeric(uk_biobank_assessment_centre_f54_3_0)) %>% 
    pivot_longer(-c(eid, Sex, agebsl), 
                 names_to = "names", 
                 values_to = "values") %>% 
    mutate(grot1 = stringi::stri_reverse(names)) %>% 
    separate(grot1, 
             c("grot2", "tp", "grot3"),
             sep ="_",
             extra = "merge") %>% 
    mutate(names = stringi::stri_reverse(grot3)) %>% 
    select(-starts_with("grot")) %>% 
    pivot_wider(names_from = names, 
                values_from = values)
  
  df.mri = 
    df.mri %>% 
    filter(eid %in% df.ukb$eid)
  df.mri = 
    df.mri %>% 
    pivot_longer(-eid, 
                 names_to = "names", 
                 values_to = "values") %>% 
    mutate(grot1 = stringi::stri_reverse(names)) %>% 
    separate(grot1, 
             c("grot2", "tp", "grot3"),
             sep ="_",
             extra = "merge") %>% 
    mutate(names = stringi::stri_reverse(grot3)) %>% 
    select(-starts_with("grot")) %>% 
    pivot_wider(names_from = names, 
                values_from = values)
  
  
  grot.df.out = inner_join(df.ukb, df.mri)
  
  load(file.path(outdir, "df.rda"))
  df.out = df.out %>% select(eid, tp)
  df.out = inner_join(df.out, grot.df.out)
  
  save(df.out, 
       df.ukb,
       file = file.path(outdir, "df.cross.rda")) 
  cat("ukb data checked. rda file save with valid data only")
}

check.oasis.fs = function(sitedir, outdir, mrifiles, is.cross = F) {
  datadir = file.path(sitedir, "data", "data_files")
  
  #load(file.path(sitedir,"sessions.info.extended.Rda"))
  load(file.path(sitedir,"sessions.info.Rda"))
  
  sessions = sessions %>% 
    mutate(Subject = substr(subs, 5,12), 
           days_to_visit = gsub("ses-d", "", ses) %>% as.numeric(), 
           scanner = paste(ManufacturersModelName, DeviceSerialNumber, sep = "_"))
  
  
  df.oasis.clinical = import(file.path(datadir, "UDSb4/csv/OASIS3_UDSb4_cdr.csv"))
  df.oasis.demo = import(file.path(datadir, "OASIS3_demographics.csv"))
  df.oasis.mri.info = import(file.path(datadir, "MRI_info.csv"))
  df.oasis.psych = import(file.path(datadir, "pychometrics/csv/OASIS3_UDSc1_cognitive_assessments.csv"))
  df.oasis.cognorm = import(file.path(sitedir,"data/cohort_files/CogNorm/csv/OASIS3_unchanged_CDR_cognitively_healthy.csv"))
  
  
  idx = !colSums(is.na(df.oasis.clinical)) == length(df.oasis.clinical$OASISID)
  df.oasis.clinical = df.oasis.clinical[,idx] %>% 
    mutate(session_clinical = OASIS_session_label) %>% 
    select(-OASIS_session_label) 
  
  idx = !colSums(is.na(df.oasis.demo)) == length(df.oasis.demo$OASISID)
  df.oasis.demo = df.oasis.demo[,idx]
  
  idx = !colSums(is.na(df.oasis.psych)) == length(df.oasis.psych$OASISID)
  df.oasis.psych = df.oasis.psych[,idx] %>% 
    mutate(session_psych = OASIS_session_label) %>% 
    select(-OASIS_session_label) 
  
  
  df.oasis.mri.info = 
    df.oasis.mri.info %>% 
    separate(`MR ID`, c("OASISID", "days_to_visit"), sep = "_MR_d", remove = F) %>% 
    filter(!is.na(days_to_visit)) %>% 
    mutate(days_to_visit = as.numeric(days_to_visit),
           `age at visit` = Age)
  
  df.oasis.cognorm= 
    df.oasis.cognorm %>%  
    mutate(OASISID = OASIS3_id, 
           cognorm = 1) %>% 
    select(OASISID, cognorm)
  
  
  
  df.oasis3 =
    reduce(list(
      df.oasis.psych, 
      df.oasis.mri.info,
      df.oasis.clinical,
      df.oasis.demo, 
      df.oasis.cognorm),
      full_join
    )
  
  df.oasis3 = 
    df.oasis3 %>% 
    mutate(time = days_to_visit/366.25, 
           AGE = AgeatEntry + time)
  
  df.oasis3 = 
    df.oasis3 %>% 
    mutate(time.dx1.normal =if_else(dx1 %in% c("Cognitively normal","No dementia"), time, NaN),
           time.dx1.exclude =if_else(!dx1 %in% c("Cognitively normal","No dementia") & !is.na(dx1), time, NaN)) %>% 
    group_by(OASISID) %>% 
    mutate(time.dx1.normal.max = max(time.dx1.normal, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(time.dx1.exclude.filt = if_else(time.dx1.exclude > time.dx1.normal.max, time.dx1.exclude, NaN)) %>% 
    group_by(OASISID) %>% 
    mutate(time.dx1.exclude.min = min(time.dx1.exclude.filt, na.rm = T)) %>% 
    ungroup() %>% 
    mutate(time.dx1.criteria = (time.dx1.exclude.min + time.dx1.normal.max)/2,
           crit.dx1.normal.max = if_else(time <= time.dx1.criteria, 1,0 ))
  
  
  if (is.cross == T ) {
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
    mutate(input = gsub("./", "", input), 
           Subject = substr(input, 5,12),
           days_to_visit = substr(input, 17,20) %>% 
             as.numeric())
  
  df.mri = inner_join(df.mri, sessions)
  
  df.out = inner_join(df.mri, df.oasis3)
  
  df.out = 
    df.out %>% filter(crit.dx1.normal.max == 1)
  
  
  
  if (is.cross == T ) {
    save(df.out, 
         df.oasis3,
         file = file.path(outdir, "df.cross.rda")) 
  } else {
    save(df.out, 
         df.oasis3,
         file = file.path(outdir, "df.rda")) 
  }
  
  cat("oasis3 data checked. rda file save with valid data only")
}
