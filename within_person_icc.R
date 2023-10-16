.libPaths(c('/ncf/mclaughlin/users/jflournoy/R/x86_64-pc-linux-gnu-library/verse-4.1.1',
            '/ncf/mclaughlin/users/jflournoy/R/x86_64-pc-linux-gnu-library/4.1',
            "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
#library(ggplot2)
library(data.table)
library(brms)
library(argparse)
CPUS_DEFAULT <- 4
parser <- ArgumentParser(description='Run within-person alpha')
parser$add_argument('--happy', action = 'store_true', help = 'happy versus calm instead of fear?')
parser$add_argument('--brainaspred', action = 'store_true', help = 'Use Brain data to predict?')
parser$add_argument('--moresamples', action = 'store_true', help = 'add samples?')
parser$add_argument('--subsample', action = 'store_true', help = 'run the model with 15 randomly selected voxels (sampled within person-sessions)')
parser$add_argument('--chains', type = 'integer', help = 'Number of chains')
parser$add_argument('--threads', type = 'integer', default = 1, help = 'Number of threads per chain')
parser$add_argument('--adaptdelta', type = 'double', help = 'Value for adapt delta setting', default = .999)
parser$add_argument('--nosmooth', action = 'store_true', help = 'use unsmoothed data?')
parser$add_argument('--predictors', action = 'store_true', help = 'compute variance explained by predictors')
parser$add_argument('--altpred', action = 'store_true', help = 'use alternative specification for predictor models')
parser$add_argument('--scramble', action = 'store_true', help = 'fear (or happy) versus scramble instead of calm?')
parser$add_argument('--cntrl_regions', action = 'store_true', help = 'Use size-matched control regions?')
parser$add_argument('--networks', action = 'store_true', help = 'Network-level analysis')
parser$add_argument('--networkrand', action = 'store_true', help = 'Network-level analysis for random sets of parcels')
parser$add_argument('--iccma', action = 'store_true', help = 'Between-person ICC and meta-analysis')


# args <- parser$parse_args(c('--brainaspred', '--chains', '1', '--adaptdelta', '.999', '--nosmooth'))
# args <- parser$parse_args(c('--brainaspred', '--chains', '4', '--threads', '2', '--adaptdelta', '.999', '--nosmooth', '--subsample', '--predictors'))
# args <- parser$parse_args(c('--chains', '4', '--altpred', '--threads', '1', '--adaptdelta', '.999', '--nosmooth', '--subsample', '--predictors', '--scramble'))
# args <- parser$parse_args(c('--chains', '4', '--adaptdelta', '.999', '--nosmooth', '--iccma'))
args <- parser$parse_args()

print(dput(args))

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
  this_ROI <- 1 #TO BE SET BY SLURM ENV VAR
  cpus_total <- CPUS_DEFAULT
  message('Not running as slurm array, setting ROI to ', this_ROI)
  message('Assuming cpus total is at least ', chains*threads)
  iterations <- 100
  warmup <- 50
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

if(args$networkrand){
  args$networks <- args$networkrand
  this_ROI <- this_ROI + 100
}
flags <- args[c('happy', 'moresamples', 'subsample', 'nosmooth', 'predictors', 'altpred', 'brainaspred', 'cntrl_regions', 'networks', 'iccma')]
NETWORK <- args$networks  
HAPPY <- args$happy
MORESAMP <- args$moresamples
SUBSAMP <- args$subsample
NOSMOOTH <- args$nosmooth
PREDICTORS <- args$predictors
ALTPRED <- args$altpred
BRAINASPRED <- args$brainaspred
SCRAMBLE <- args$scramble
CNTRL_REGIONS <- args$cntrl_regions
ICCMA <- args$iccma

if(ALTPRED & !PREDICTORS){
  stop('Alternative predictor specifications implies predictors.')
}

if(BRAINASPRED & PREDICTORS){
  stop('--brainaspred and --predictors are exclusive')
}

if(any(unlist(flags))){
  suffix <- paste0('-', paste(paste(c('happy', 'more_01', 'sub' ,'ns' , 'pred', 'alt', 'brainpred', 'cntrlreg', 'network', 'iccma')[unlist(flags)], collapse = '-')))  
} else {
  suffix <- ''
}
if(NETWORK | ICCMA){
  NETWORK_NUMBER <- this_ROI
  if(SUBSAMP){
    message('Subsampling is ignored for network-level and between-person ICC analysis.')
  }
}
if(ICCMA && cpus_total != 9){
  warning(sprintf('Number of CPUs should be 9 but is %d', cpus_total))
}
cope <- ifelse(HAPPY, ifelse(SCRAMBLE, 'happyGTscramble', 'happyGTcalm'), ifelse(SCRAMBLE, 'fearGTscramble', 'fearGTcalm'))
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

if(NOSMOOTH){
  tsv_dir <- 'tsv_nosmooth'
} else {
  tsv_dir <- 'tsv'  
}

psych_d <- readr::read_csv('./psychdata/stress_psych_sleep-midpoint5.csv')
psych_d$sess <- sprintf('%02d', psych_d$time)
psych_d$id <- as.character(psych_d$idnum)

if(NETWORK){
  network_mapping <- data.table::fread('Schaefer2018_400Parcels_7Networks_order.txt',
                                       header = FALSE, col.names = c('ROI', 'ROI_Label'),
                                       select = 1:2)
  network_regex <- '^7Networks_(RH|LH)_(([A-Za-z]+)_*([A-Za-z]+)*)_(\\d{1,2})$'
  network_mapping[, c('hem', 'network', 'subnet', 'netroinum') :=
                    lapply(c('\\1', '\\3', '\\4', '\\5'), gsub, pattern = network_regex, x = ROI_Label)]
  network_mapping[, net_num := as.numeric(factor(network, levels = unique(network)))]
  if(!args$networkrand){
    this_ROI <- network_mapping[net_num == NETWORK_NUMBER, ROI]  
  } else {
    set.seed(547524) #from random.org
    network_sizes <- network_mapping[, .N, by = net_num]
    network_size <- network_sizes[net_num == NETWORK_NUMBER - 100, N]
    this_ROI <- sample(1:400, network_size)
  }
}

roi_d <- get_roi_data(roi = this_ROI, cope = cope, ncpus = cpus_total, tsv_dir = tsv_dir, cntrl_regions = CNTRL_REGIONS)
roi_d[, c('y_orig', 'y') := list(y, scale(y))]
if(ICCMA){
  roi_d <- roi_d[, .(y = mean(y, na.rm = TRUE)), by = c('ROI', 'id', 'sess', 'cope')]
}

roi_d[, sess_fac := factor(sess)]
if(NETWORK){
  roi_d <- roi_d[, .(y = mean(y, na.rm = TRUE),
                     y_se = sd(y, na.rm = TRUE) / sqrt(sum(!is.na(y)))),
                 by = c('id', 'sess', 'sess_fac', 'ROI')]
  roi_d[, roi_fac := factor(ROI)]
} 

if(!NETWORK && !ICCMA && !dim(roi_d)[[1]] == max(roi_d$v) * 292){
  stop('Data length is not as expected')
}

if(!(NETWORK | ICCMA) && SUBSAMP){
  nvox <- max(roi_d$v)
  keep_vec <- c(rep(TRUE, npersubsamp), rep(FALSE, max(roi_d$v) - npersubsamp))
  roi_d <- roi_d[, keep := sample(x = keep_vec, size = .N), by = c('id', 'sess', 'cope')][keep == TRUE]
  if(dim(roi_d)[[1]] != npersubsamp*292){
    stop('Error with subsampling')
  }
  #source('lme4_Addendum.R', print.eval = TRUE)
}

message(sprintf('Total observations per id/sess: %d', max(roi_d[, .N, by = c('id', 'sess')][, N])))

if(PREDICTORS){
  roi_d <- merge(roi_d, psych_d, by = c('id', 'sess'))
  if(ALTPRED){
    if(NETWORK){
      model_form <- list(eptot = bf(y | mi(sdy = y_se) ~ 1 + TIMECENTER + GCEN_EPISODICTOT + WCEN_EPISODICTOT + (1 | id/sess_fac), family = brms::student()),
                         chsev = bf(y | mi(sdy = y_se) ~ 1 + TIMECENTER + GCEN_CHRONICSEV + WCEN_CHRONICSEV + (1 | id/sess_fac), family = brms::student()),
                         panasneg = bf(y | mi(sdy = y_se) ~ 1 + TIMECENTER + GCEN_PANASNEG + WCEN_PANASNEG + (1 | id/sess_fac), family = brms::student()),
                         sleep = bf(y | mi(sdy = y_se) ~ 1 + TIMECENTER + GCEN_SLEEPDUR + WCEN_SLEEPDUR + (1 | id/sess_fac), family = brms::student()),
                         time = bf(y | mi(sdy = y_se) ~ 1 + TIMECENTER + (1 | id/sess_fac), family = brms::student()))
    } else {
      model_form <- list(eptot = bf(y ~ 1 + TIMECENTER + GCEN_EPISODICTOT + WCEN_EPISODICTOT + (1 | id/sess)),
                         chsev = bf(y ~ 1 + TIMECENTER + GCEN_CHRONICSEV + WCEN_CHRONICSEV + (1 | id/sess)),
                         panasneg = bf(y ~ 1 + TIMECENTER + GCEN_PANASNEG + WCEN_PANASNEG + (1 | id/sess)),
                         sleep = bf(y ~ 1 + TIMECENTER + GCEN_SLEEPDUR + WCEN_SLEEPDUR + (1 | id/sess)),
                         time = bf(y ~ 1 + TIMECENTER + (1 | id/sess)),
                         timefac = bf(y ~ 1 + (1 | id) + (1 | sess) + (1 | id:sess)))
    }
  } else {
    model_form <- list(eptot = bf(y ~ 1 + TIMECENTER + GCEN_EPISODICTOT + WCEN_EPISODICTOT + (1 | id)),
                       chsev = bf(y ~ 1 + TIMECENTER + GCEN_CHRONICSEV + WCEN_CHRONICSEV + (1 | id)),
                       panasneg = bf(y ~ 1 + TIMECENTER + GCEN_PANASNEG + WCEN_PANASNEG + (1 | id)),
                       sleep = bf(y ~ 1 + TIMECENTER + GCEN_SLEEPDUR + WCEN_SLEEPDUR + (1 | id)),
                       time = bf(y ~ 1 + TIMECENTER + (1 | id)),
                       all = bf(y ~ 1 + sess_fac + GCEN_EPISODICTOT + WCEN_EPISODICTOT + 
                                  GCEN_CHRONICSEV + WCEN_CHRONICSEV + 
                                  GCEN_PANASNEG + WCEN_PANASNEG + 
                                  GCEN_SLEEPDUR + WCEN_SLEEPDUR + (1 | id/sess)))
  }
  outname <- paste0(outname, '_', names(model_form))
} else if(BRAINASPRED){
  #Center data and get measurement error
  roi_d_reduced <- roi_d[, list(id_sess_mean = mean(y_orig, na.rm = TRUE), 
                                id_sess_se = sd(y_orig, na.rm = TRUE) / sqrt(sum(!is.na(y_orig)))), 
                         by = c('id', 'sess')]
  
  roi_d_reduced[, c('id_mean') := list(mean(id_sess_mean)), by = 'id']
  roi_d_reduced[, c('WCEN_BRAIN', 'GCEN_BRAIN', 'WB_SE') := 
                        list(id_sess_mean - id_mean, id_mean - mean(id_sess_mean), id_sess_se)]
  
  roi_d <- merge(roi_d_reduced[, c('id', 'sess', 'WCEN_BRAIN', 'GCEN_BRAIN', 'WB_SE')], psych_d, by = c('id', 'sess'))
  
  model_form <- list(PHQ9 = bf(PHQ9_TOT ~ 1 + TIMECENTER + GCEN_BRAIN + me(WCEN_BRAIN, WB_SE) + (1 | id)),
                     GAD7 = bf(GAD7_TOT ~ 1 + TIMECENTER + GCEN_BRAIN + me(WCEN_BRAIN, WB_SE) + (1 | id)))
  outname <- paste0(outname, '_', names(model_form))
    
} else if(NETWORK){
  model_form <- bf(y | mi(sdy = y_se) ~ 1 + sess_fac + (1 | id/sess_fac))
} else if(ICCMA){
  model_form <- bf(y ~ 1 + sess_fac + (1 | id))
} else {
  #ICC and alpha
  model_form <- bf(y ~ 1 + sess_fac + (1 | id/sess))  
}

message(paste(c('Outfile(s):', sprintf('%s', outname)), collapse = '\n'))

fit_model <- function(this_model_form, this_outname, data, chains, cores, future = FALSE, iter, warmup, threads, adaptdelta){
  if(!dir.exists('fits')){
    dir.create('fits')
  } 
  
  if(future){
    nbrOfWorkers()
  }

  if(!BRAINASPRED){ #bad practice
    priors <- prior(normal(0, 0.25), class = sd)
  } else {
    priors <- prior(student_t(3, 0, 4.4), class = sd) 
  }
  
  if(is.null(threads) || threads == 1){
    model_fit <- brm(this_model_form, data = data,
                     family = 'student', prior = priors,
                     backend = 'cmdstan',
                     chains = chains, cores = cores, future = FALSE,
                     iter = iterations, warmup = warmup, control = list(adapt_delta = adaptdelta, max_treedepth = 15),
                     file = this_outname)
  } else if(threads > 1){
    message('Running with multiple threads')

    if(future){
      warning('NOT SURE HOW THREADS AND FUTURE PLAY TOGETHER SO WATCH OUT!')
    }
    
    threader <- threading(threads = threads, grainsize = NULL)
    model_fit <- brm(this_model_form, data = data,
                     family = 'student', prior = priors,
                     backend = 'cmdstan', threads = threader,
                     chains = chains, cores = cores, future = future,
                     iter = iterations, warmup = warmup, control = list(adapt_delta = adaptdelta, max_treedepth = 15),
                     file = this_outname)
  } else {
    stop('Threads argument misspecified...')
  }
  return(model_fit)
}

message('Running!')

print(model_form)
print(outname)
message('Data size: ')
dput(dim(roi_d))

if(ICCMA){
  message('Estimating pairwise ICCs...')
  simultaneous_models <- cpus_total
  unique_sessions <- sort(unique(roi_d$sess))
  session_pairs <- mapply(c, unique_sessions[-length(unique_sessions)], unique_sessions[-1], SIMPLIFY = FALSE)
  if(length(session_pairs) != cpus_total){
    warning('Session pairs is not equal to number of CPUs!')
  }
  cl <- parallel::makeCluster(simultaneous_models)
  nada <- parallel::clusterExport(cl = cl, varlist = ls())  
  #(length(unique_sessions) - 1)
  session_pair_iccs_list <- parallel::parLapply(cl = cl, X = session_pairs, fun = function(session_pair){
    library(data.table)
    library(brms)
    twosess_d <- roi_d[sess %in% session_pair]
    sess_tag <- paste(session_pair, collapse = '_')
    this_outname <- paste(outname, sess_tag, sep = '_')
    a_fit <- fit_model(this_model_form = model_form, this_outname = this_outname, data = twosess_d,
                       chains = 1, cores = 1, future = FALSE, iter = iterations, warmup = warmup,
                       threads = threads, adaptdelta = adaptdelta)
    
    vc <- VarCorr(a_fit, summary = FALSE)
    
    sigma2 <- vc$residual__$sd^2
    tau_pi <- vc$id$sd^2
    total_var <- tau_pi + sigma2
    pvar_id <- tau_pi / total_var
    post_df <- as.data.table(posterior_summary(pvar_id, robust = TRUE))
    post_df[, sess := sess_tag]
    post_df[, est_logit := median(qlogis(pvar_id))]
    post_df[, sd_logit := sd(qlogis(pvar_id))]
    return(post_df)
  })
  parallel::stopCluster(cl)

  message('Estimating meta-analytic average...')
  meta_icc_dt <- rbindlist(session_pair_iccs_list)
  meta_icc_form <- bf(est_logit | se(sd_logit, sigma = TRUE) ~ 1)
  meta_fit <- brm(formula = meta_icc_form, data = meta_icc_dt,
                  family = 'gaussian',
                  chains = chains, cores = chains, iter = iterations, 
                  control = list(adapt_delta = 0.999, max_treedepth = 15),
                  file = outname, file_refit = 'on_change')
  
  meta_sum <- data.table(posterior_summary(plogis(as_draws_matrix(meta_fit, variable = 'b_Intercept'))), robust = TRUE)
  meta_sum[, sess := 'meta']
  
  meta_sum[, ROI := this_ROI]
  saveRDS(meta_sum, paste0(outname, '-meta_sum.rds'), compress = FALSE)
} else {
  fits <- mapply(fit_model, this_model_form = model_form, this_outname = outname, 
                 MoreArgs = list(data = roi_d, chains = chains, cores = chains, 
                                 iter = iterations, warmup = warmup, threads = threads, 
                                 adaptdelta = adaptdelta))  
}

sessionInfo()
