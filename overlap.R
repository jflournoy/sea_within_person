d <- data.table::melt(data.table::fread('~/sea_icc_stuff/within_person/overlap.tsv'))[, 3:4]
s <- data.table::melt(data.table::fread('~/sea_icc_stuff/within_person/parcel_sizes.tsv'), value.name = 'total')[, 3:4]
d <- s[d, on = 'variable']
d[, p := value / total]
data.table::setorder(d, p)
# plot(1:dim(d)[[1]], d$p)
# abline(h = .2)
cat(paste(sort(c(410, 420, as.numeric(gsub('NZcount_', '', d[p > .2, variable])))), collapse = ','))
