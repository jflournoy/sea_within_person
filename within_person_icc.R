library(data.table)
library(brms)
library(ggplot2)
library(argparse)
CPUS_DEFAULT <- 1
parser <- ArgumentParser(description='Run within-person alpha')
parser$add_argument('--happy', action = 'store_true', help = 'happy versus calm instead of fear?')
parser$add_argument('--moresamples', action = 'store_true', help = 'add samples?')
parser$add_argument('--subsample', action = 'store_true', help = 'run the model with 15 randomly selected voxels (sampled within person-sessions)')
parser$add_argument('--chains', type = 'integer', help = 'Number of chains')
parser$add_argument('--threads', type = 'integer', help = 'Number of threads per chain')
parser$add_argument('--adaptdelta', type = 'double', help = 'Value for adapt delta setting', default = .999)
parser$add_argument('--nosmooth', action = 'store_true', help = 'use unsmoothed data?')
parser$add_argument('--predictors', action = 'store_true', help = 'compute variance explained by predictors')

# args <- parser$parse_args(c('--chains', '1', '--threads', '1', '--adaptdelta', '.999', '--nosmooth', '--subsample', '--predictors'))
args <- parser$parse_args()

print(args)

iterations <- 5000
warmup <- 1000
chains <- args$chains
threads <- args$threads
npersubsamp <- 15
adaptdelta <- args$adaptdelta

source('get_roi_data.R')

task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
cpus_total <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

if(is.na(task_id)){
  message('Not running as slurm array, setting ROI to some value')
  message('Assuming cpus total is at least ', chains*threads)
  this_ROI <- 421 #TO BE SET BY SLURM ENV VAR
  cpus_total <- CPUS_DEFAULT
} else {
  this_ROI <- task_id
  message('Running as slurm array with task id=ROI: ', task_id, '=', this_ROI)
  message('Total CPUs: ', cpus_total)
  message('Total CPUs needed: ', chains*threads)
  if(chains*threads > cpus_total){
    stop('Too few CPUs')
  }
}
setDTthreads(cpus_total)
flags <- args[c('happy', 'moresamples', 'subsample', 'nosmooth', 'predictors')]
HAPPY <- args$happy
MORESAMP <- args$moresamples
SUBSAMP <- args$subsample
NOSMOOTH <- args$nosmooth
PREDICTORS <- args$predictors

if(any(unlist(flags))){
  suffix <- paste0('-', paste(paste(c('happy', 'more_01', 'sub' ,'ns' , 'pred')[unlist(flags)], collapse = '-')))  
} else {
  suffix <- ''
}

cope <- ifelse(HAPPY, 'happyGTcalm', 'fearGTcalm')
outname <- sprintf('fits/alpha-%s-ROI%03d%s', cope, this_ROI, suffix)

if(MORESAMP){
  reg <- gsub('01', '(.*)', paste0(outname, '.rds'))
  morefiles <- dir(dirname(reg), pattern = basename(reg), full.names = TRUE)
  if(length(morefiles) > 0){
    max_morefiles <- max(as.numeric(gsub(reg, '\\1', morefiles)))
    outname <- gsub('01', sprintf('%02d', max_morefiles + 1), outname)
  }
}
#outname <- 'fits/alpha-fearGTcalm-ROI410_TEST_THREAD'
message('Outfile is: ', outname)

tsv_dir <- 'tsv'
if(NOSMOOTH){
  tsv_dir <- 'tsv_nosmooth'
}

roi_d <- get_roi_data(roi = this_ROI, cope = cope, ncpus = cpus_total, tsv_dir = tsv_dir)
roi_d[, y := scale(y)]
if(!dim(roi_d)[[1]] == max(roi_d$v) * 292){
  stop('Data length is not as expected')
}

psych_d <- readr::read_csv('/ncf/mclaughlin/stressdevlab/stress_pipeline/Group/NewNeuropoint/data/')

source('lme4_Addendum.R', print.eval = TRUE)

if(SUBSAMP){
  nvox <- max(roi_d$v)
  keep_vec <- c(rep(TRUE, npersubsamp), rep(FALSE, max(roi_d$v) - npersubsamp))
  roi_d <- roi_d[, keep := sample(x = keep_vec, size = .N), by = c('id', 'sess', 'cope')][keep == TRUE]
  if(dim(roi_d)[[1]] != npersubsamp*292){
    stop('Error with subsampling')
  }
  source('lme4_Addendum.R', print.eval = TRUE)
}
roi_d[, sess_fac := factor(sess)]

message(sprintf('Total voxels: %d', max(roi_d[, .N, by = c('id', 'sess')][, N])))

#  icc_form <- bf(y ~ 1 + sess + (1 | id))
alpha_model_form <- bf(y ~ 1 + sess_fac + (1 | id/sess))

dir.create('fits')
get_prior(alpha_model_form, data = roi_d)
priors <- prior(normal(0, 0.25), class = sd)

if(threads == 1){
  alpha_model_comp <- brm(alpha_model_form, data = roi_d, 
                          family = 'student',
                          chains = 1, cores = 1,
                          iter = 10, prior = priors,
                          file = 'fits/alpha_compiled')
  prior_summary(alpha_model_comp)
  alpha_model_fit <- update(alpha_model_comp, formula. = alpha_model_form, newdata = roi_d, 
                            family = 'student', prior = priors,
                            chains = chains, cores = chains, 
                            iter = iterations, warmup = warmup, control = list(adapt_delta = adaptdelta, max_treedepth = 15), 
                            file = outname)
} else if(threads > 1){
  message('Running with multiple threads, compiling every time.')
  threader <- threading(threads = threads, grainsize = NULL)
  alpha_model_fit <- brm(alpha_model_form, data = roi_d, 
                         family = 'student', prior = priors,
                         backend = 'cmdstan', threads = threader,
                         chains = chains, cores = chains, 
                         iter = iterations, warmup = warmup, control = list(adapt_delta = adaptdelta, max_treedepth = 15), 
                         file = outname)
} else {
  stop('Threads argument misspecified...')
}


