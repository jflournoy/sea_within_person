library(data.table)
d_ <- data.table(expand.grid(mode = c('pa', 'sleep'), summary = c('daily', 'weekly', 'id')))

d <- data.table::rbindlist(lapply(split(d_, d_[, c('mode', 'summary')]), function(x){
  if (x['summary'] == 'daily'){
    x <- cbind(x, data.table(time = c(7, 14, 21)))
  } else if (x['summary'] == 'weekly') {
    x <- cbind(x, data.table(time = c(1, 2, 3)))
  } else {
    x
  }
}), fill = TRUE)
d[is.na(time), time := 0]
d

readr::write_tsv(d, file = '~/code/abcd_ml/arguments_file.txt')
