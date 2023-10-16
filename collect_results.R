library(data.table)
library(brms)
CPUS_DEFAULT <- 4
cpus_total <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

files_to_include_expression <- NULL
data_rows_to_keep_expression <- NULL

if(is.na(cpus_total)){
  message('Not running as slurm array, setting CPUs to some value: ', CPUS_DEFAULT)
  cpus_total <- CPUS_DEFAULT
  INTERACTIVE <- TRUE
} else {
  message('Running as slurm array with CPUs: ', cpus_total)
  INTERACTIVE <- FALSE
}
setDTthreads(cpus_total)

#collect
read_fit_files <- function(fit_files_dt){
  setDTthreads(1)
  fits <- lapply(file.path('fits', fit_files_dt$file), readRDS)
  cfit <- brms::combine_models(mlist = fits)
  return(cfit)
}
compute_reliability <- function(fit_files_dt){
  #fit_files_dt <- filenames[[1]][[1]]
  setDTthreads(1)
  print(fit_files_dt$file)
  cfit <- read_fit_files(fit_files_dt)

  model_info <- unique(fit_files_dt[, c('contrast', 'roi', 'sub', 'ns', 'pred', 'alt', 'cntrlreg', 'network')])

  #We don't compute our proportions of variance or alpha from predictor models,
  #except for the model in which we include all predictors.
  if(grepl('-(?:brain)*pred', model_info$pred) & model_info$pred != '-pred_all'){
    if(grepl('-pred(-alt)*_time', model_info$pred)){
      variable <- 'b_TIMECENTER'
    } else {
      variable <- grep('b(sp_me|_)WCEN', variables(cfit), value = TRUE)
    }
    post_samples <- brms::as_draws_matrix(cfit, variable = variable, inc_warmup = FALSE)
    post_sum <- data.table(
      brms::posterior_summary(post_samples, probs = c(.025, .975,
                                              round(abs(c(0, -1) + .05/421/2), 5)),
                              robust = TRUE, variable = variable),
      keep.rownames = TRUE)
    sign_est <- sign(post_sum$Estimate[[1]])
    p_sign_prop <- posterior::summarise_draws(post_samples, 
                                              list(p_SSEst = \(x) mean( (x / sign_est) > 0 ),
                                                   p_nSSEst = \(x) mean( (x / sign_est) < 0 ),
                                                   draws = length))
    setnames(p_sign_prop, 'variable', 'rn')
    post_sum <- merge(post_sum, p_sign_prop, by = c('rn'))
    post_sum[, c('contrast', 'roi', 'sub', 'ns', 'pred', 'alt', 'cntrlreg', 'network') := c(model_info)]
    
    d_df <- as.data.table(density(post_samples[, variable], n = 100)[c('x', 'y')])
    d_df$name <- variable
    
    r <- list(post_sum = post_sum, density_df = d_df)
  } else {
    vc <- VarCorr(cfit, summary = FALSE)
    tau_beta <- vc$id$sd^2
    sigma2 <- vc$residual__$sd^2
    
    if(model_info$network == '-network'){
      tau_pi <- vc$`id:sess_fac`$sd^2
      nvoxel <- length(cfit$data$y)/dim(unique(cfit$data[, c('sess_fac', 'id:sess_fac')]))[[1]]
    } else {
      tau_pi <- vc$`id:sess`$sd^2
      nvoxel <- length(cfit$data$y)/dim(unique(cfit$data[, c('sess', 'id:sess')]))[[1]]
    }
    alpha_post <- tau_pi / (tau_pi + sigma2/nvoxel)
    total_var <- tau_pi + tau_beta + sigma2
    pvar <- tau_pi / total_var
    pvar_id <- tau_beta / total_var
    post_df <- data.frame(alpha = alpha_post[,1],
                          prop_var = pvar[,1],
                          prop_var_id = pvar_id[,1])
    #Do additional work to look at different amounts of variance explained by
    #predictors.
    if(model_info$pred == '-pred_all'){
      null_model_files <- data.table(file = dir('fits',
                                                pattern = unique(gsub('-pred_all', '', fit_files_dt$file))))
      null_model <- read_fit_files(null_model_files)
      vc_null <- VarCorr(null_model, summary = FALSE)
      tau_pi_delta <- vc_null$`id:sess`$sd^2 - tau_pi
      tau_beta_delta <- vc_null$id$sd^2 - tau_beta
      sigma2_delta <- vc_null$residual__$sd^2 - sigma2
      tau_pi_delta_prop <- tau_pi_delta / vc_null$`id:sess`$sd^2
      tau_beta_delta_prop <- tau_beta_delta / vc_null$id$sd^2
      sigma2_delta_prop <- sigma2_delta / vc_null$residual__$sd^2
      post_df <- cbind(post_df,
                       data.frame(tau_pi_delta = tau_pi_delta[,1],
                                  tau_beta_delta = tau_beta_delta[,1],
                                  sigma2_delta = sigma2_delta[,1],
                                  tau_pi_delta_prop = tau_pi_delta_prop[,1],
                                  tau_beta_delta_prop = tau_beta_delta_prop[,1],
                                  sigma2_delta_prop = sigma2_delta_prop[,1]))
    }
    if(model_info$sub == '-sub' & model_info$ns == '-ns'){
      #save a sample of the full posterior to create an overall mean later
      fn <- file.path('post_samples', gsub('.rds$', '-psamp.rds', unique(fit_files_dt[, c('file')])))
      if(!file.exists(fn)) saveRDS(post_df[sample(1:dim(post_df)[[1]], 1000), ], file = fn)
    }
    post_df$idsess_right <- post_df$prop_var_id + post_df$prop_var
    d_df <- do.call(rbind, lapply(c('prop_var_id', 'idsess_right'), function(n){
      p <- post_df[, n]
      d <- as.data.frame(density(p, n = 100)[c('x', 'y')])
      d$name <- n
      return(d)
    }))
    post_sum <- data.table(posterior_summary(post_df, robust = TRUE), keep.rownames = TRUE)
    post_sum[, c('nvoxel', 'contrast', 'roi', 'sub', 'ns', 'pred', 'alt', 'cntrlreg', 'network') := c(nvoxel, model_info)]
    r <- list(post_sum = post_sum, density_df = d_df)
  }

  return(r)
}

collect_smoothness <- function(data_file_source = '/ncf/mclaughlin/stressdevlab/stress_pipeline/Group/FaceReactivity/NewNeuropoint/datafiles/setfilenames_fearGTcalm.txt'){
  data_files <- fread(data_file_source, header = FALSE, col.names = 'file')
  data_files[, c('ID', 'sess') := transpose(stringi::stri_match_all(file, regex = '.*/(\\d{4})/month(\\d{2})/.*'))[2:3]]
  data_files[, file := NULL]
  smoothness_data_files <- dir('smoothness', pattern = '\\d{1,3}.*txt', full.names = TRUE)
  smoothness_subcort_data_files <- dir('smoothness_subcort', pattern = '\\d{1,3}.*txt', full.names = TRUE)
  smoothness_cntrl_data_files <- dir('smoothness_cntrl', pattern = '\\d{1,3}.*txt', full.names = TRUE)
  names(smoothness_data_files) <- stringi::stri_match_first_regex(smoothness_data_files, '\\w+/(\\d{1,3})-.*')[, 2]
  names(smoothness_subcort_data_files) <- 400 + as.numeric(stringi::stri_match_first_regex(smoothness_subcort_data_files, '\\w+/(\\d{1,3})-.*')[, 2])
  names(smoothness_cntrl_data_files) <- stringi::stri_match_first_regex(smoothness_cntrl_data_files, '\\w+/(\\d{1,3})-.*')[, 2]

  smoothness_data_list <- lapply(list(cort = smoothness_data_files, subcort = smoothness_subcort_data_files, cntrl = smoothness_cntrl_data_files), 
                                 \(x){
                                   d_list <- lapply(x, \(y) tryCatch(cbind(data_files, 
                                                                           fread(y, header = TRUE, colClasses = 'numeric', col.names = 'smoothness')),
                                                                     error = \(y) return(NULL)))
                                   return(rbindlist(d_list, idcol = 'roi'))
                                 })
  smoothness_data <- rbindlist(smoothness_data_list, idcol = 'anatomy')
  return(smoothness_data)
}


rel_fn <- 'reliability_results.rds'
bd_fn <- 'boundary_density.rds'
smoothness_fn <- 'smoothness_results.rds'

#SET THESE TO NULL USUALLY
#files_to_include_expression <- "contrast == 'scramble' & pred == ''"
#files_to_include_expression <- "cntrlreg == '-cntrlreg' & contrast == 'scramble'"
#files_to_include_expression <- "network == '-network' & roi > 100"
#data_rows_to_keep_expression <- "cntrlreg != '-cntrlreg'"
if(!file.exists(rel_fn) | !is.null(files_to_include_expression)){
  gsub_regex <- 'alpha-fearGT(calm|scramble)-ROI([0-9]{1,3})(-more_\\d{2})*(-sub)*(-ns)*(-(?:brain)*pred(-alt)*_\\w+)*(-cntrlreg)*(-network)*.rds'
  
  alpha_fit_files <- data.table(file = dir('fits', pattern = 'alpha-fearGT.*rds'))
  alpha_fit_files[, c('contrast', 'roi', 'more', 'sub', 'ns', 'pred', 'alt', 'cntrlreg', 'network') :=
                    lapply(list('\\1', '\\2', '\\3', '\\4', '\\5', '\\6', '\\7', '\\8', '\\9'), gsub, pattern = gsub_regex, x = file)]
  
  if(!is.null(files_to_include_expression)){
    alpha_fit_files <- subset(alpha_fit_files, eval(parse(text = files_to_include_expression)))
    if(file.exists(rel_fn)){
      reliability_dt_existing <- readRDS(rel_fn)
      if(!is.null(data_rows_to_keep_expression)){
        reliability_dt_existing <- subset(reliability_dt_existing, eval(parse(text = data_rows_to_keep_expression)))
      }
    } else {
      reliability_dt_existing <- NULL
    }
  }
  
  filenames <- split(split(alpha_fit_files,
                           by = c('contrast', 'roi', 'sub', 'ns', 'pred', 'alt', 'cntrlreg', 'network'),
                           drop = TRUE), 1:cpus_total)
  
  cl <- parallel::makeCluster(cpus_total)

  parallel::clusterExport(cl, c('compute_reliability', 'read_fit_files'))
  nada <- parallel::clusterEvalQ(cl, {library(data.table); library(brms); setDTthreads(1)})
  system.time({reliability_list <- parallel::parLapply(cl = cl, filenames, function(some_fns){
    lapply(some_fns, \(x) tryCatch(compute_reliability(x), error = NULL))
  })})

  parallel::stopCluster(cl)

  #posterior summaries of the various reliability metrics
  reliability_list_ul <- unlist(reliability_list, recursive = FALSE)
  reliability_dt <- rbindlist(lapply(reliability_list_ul,
                                     `[[`, 'post_sum'),
                              fill = TRUE)
  #computed densities of the reliability metrices posteriors
  # boundary_density_dt <- rbindlist(lapply(reliability_list_ul,
  #                                         function(rel){
  #                                           rel$density_df$roi <- unique(rel$post_sum$roi)[[1]]
  #                                           rel$density_df$sub <- unique(rel$post_sum$sub)[[1]]
  #                                           rel$density_df$ns <- unique(rel$post_sum$ns)[[1]]
  #                                           return(rel$density_df)
  #                                         }), fill = TRUE)
  if(!is.null(reliability_dt_existing)){
    reliability_dt <- rbindlist(list(reliability_dt_existing, reliability_dt), fill = TRUE)
  }
  #reliability_dt[is.na(network), network := '']
  # ggplot(boundary_density_dt[sub == '' & ns == '-ns'], aes(x = x, y = y, group = roi)) +
  #   geom_line(alpha = .2) +
  #   facet_wrap(~ name, scales = 'free')
  saveRDS(reliability_dt, rel_fn)
  # saveRDS(boundary_density_dt, bd_fn)
} else {
  reliability_dt <- readRDS(rel_fn)
  boundary_density_dt <- readRDS(bd_fn)
}

if(!file.exists(smoothness_fn)){
  smoothness_data <- collect_smoothness()
  smoothness_data[smoothness < 2.71, smoothness := NA]
  saveRDS(smoothness_data, smoothness_fn)
} else {
  smoothness_data <- readRDS(smoothness_fn)
}
