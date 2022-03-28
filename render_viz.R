library(argparse)

parser <- ArgumentParser(description='Run within-person alpha')
parser$add_argument('--subsample', action = 'store_true', help = 'Visualize results from subsample models.')
args <- parser$parse_args()
SUBSAMP <- args$subsample

setwd('~/sea_icc_stuff/within_person/')
ROI <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

message("ROI: ", ROI)
message("SUBSAMP: ", SUBSAMP)

report_fn <- sprintf("reports/report_ROI_%03d%s.html", ROI, ifelse(SUBSAMP, '-sub', ''))
if(!file.exists(report_fn)) rmarkdown::render('viz.Rmd', output_file = report_fn, params = list(ROI = ROI, SUBSAMP = SUBSAMP))
