get_roi_data <- function(roi = 1, cope = 'fearGTcalm', ncpus = 1, tsv_dir = 'tsv'){
  f <- dir(tsv_dir, pattern = sprintf('\\d{4}.*%s.tsv', cope), full.names = TRUE)
  fregex <- '.*/(\\d{4})_(\\d{2})_(\\w+).tsv'
  
  sch400HO <- fread(file.path(tsv_dir, 'Sch400_HO_combined.tsv'), sep = ' ',  col.names = c('i', 'j', 'k', 'ROI'))[ROI == roi, c('i', 'j', 'k')]
  sch400HO[, v := 1:.N]
  
  cl <- parallel::makePSOCKcluster(ncpus)
  parallel::clusterExport(cl = cl, varlist = c('fregex', 'sch400HO'), envir = rlang::env(fregex = fregex, sch400HO = sch400HO))
  nada <- parallel::clusterEvalQ(cl, {library(data.table); setDTthreads(1)})
  f_split <- split(f, 1:ncpus)
  sea_data_list <- parallel::parLapply(cl = cl, X = f_split, fun = function(f_chunk){
    data_list <- lapply(1:length(f_chunk), function(i){
      d <- fread(f_chunk[[i]], sep = ' ',  col.names = c('i', 'j', 'k', 'y'))[sch400HO, on = c('i', 'j', 'k')][, c('v','y')]
      id <- gsub(fregex, '\\1', f_chunk[[i]])
      sess <- gsub(fregex, '\\2', f_chunk[[i]])
      cope <- gsub(fregex, '\\3', f_chunk[[i]])
      d[, c('id', 'sess', 'cope') := list(id, sess, cope)]
      return(d)
    })
    
    r <- rbindlist(data_list)
    return(r)
  })
  parallel::stopCluster(cl)
  
  sea_data_dt <- rbindlist(sea_data_list)
}
