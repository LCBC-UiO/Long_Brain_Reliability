args = commandArgs(TRUE)

outdir=as.character(args[1])
value=as.character(args[2])
m=as.numeric(args[3])
global= as.logical(args[4])

impute_outlier_data = function(outdir,value,m, global) {
  library(mice)
  set.seed(123)
  if (global == F) {
    filename = file.path(outdir, "df.merged.Rda")
  } else if (global == T) {
    filename = file.path(outdir, "df.merged.global.Rda")
  }
  
  load(filename)
  if(value == "delta") {
    db = df$df$delta$df
    idx = df$outliers$delta$idx
    rm.subs = df$outliers$delta$rmsubs.miss.data
  } else if (value == "mean") {
    db = df$df$mean$df
    idx = df$outliers$mean$idx
    rm.subs = df$outliers$mean$rmsubs.miss.data
  }
  db[!idx] = NaN
  db = db[rm.subs,]
  names(db) = paste0("v", 1:length(names(db)))
  imputed = mice(db, m=m)
  
  load(filename)
  if(value == "delta") {
    df$df$imputed.delta = imputed
  } else if (value == "mean") {
    df$df$imputed.mean = imputed
  }
  save(df, file = filename)
}


impute_outlier_data(outdir,value,m, global) 

