library(data.table)
library(future)
library(future.apply)
data.table::setDTthreads(1)
plan(tweak(multisession, workers = future::availableCores()))
ma_files <- dir('fits/', pattern = '.*meta_sum.rds', full.names = TRUE)
names(ma_files) <- ma_files
ma_dt_list <- future_lapply(ma_files, readRDS)
data.table::setDTthreads(future::availableCores())
ma_dt <- rbindlist(ma_dt_list, idcol = 'file')
gsub_regex <- 'fits//alpha-fearGT(calm|scramble)-ROI([0-9]{1,3})(-more_\\d{2})*(-sub)*(-ns)*(-(?:brain)*pred(-alt)*_\\w+)*(-cntrlreg)*(-network)*-iccma-meta_sum.rds'
ma_dt[, c('contrast', 'roi', 'more', 'sub', 'ns', 'pred', 'alt', 'cntrlreg', 'network') := 
        lapply(list('\\1', '\\2', '\\3', '\\4', '\\5', '\\6', '\\7', '\\8', '\\9'), gsub, pattern = gsub_regex, x = file)]
ma_dt[, rn := 'icc_ma_id']
ma_dt[, file := NULL]
ma_dt[, sess := NULL]
ma_dt[, ROI := NULL]
