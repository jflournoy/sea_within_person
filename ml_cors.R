library(psych)
psych_d <- readr::read_csv('./psychdata/stress_psych_sleep-midpoint5.csv')
psych_d$sess <- sprintf('%02d', psych_d$time)
psych_d$id <- as.character(psych_d$idnum)

ml_cors <- psych::statsBy(psych_d[, c('CHRONICSEV', 'EPISODICTOT', 'PANASNEG', 'sleep_dur.mean', 'idnum')], 
                          group = 'idnum', method = 'spearman')
print(ml_cors, short = FALSE)
knitr::kable(ml_cors$ci.bg$r.ci, format = 'simple', digits = 2)
knitr::kable(ml_cors$ci.wg$r.ci, format = 'simple', digits = 2)


