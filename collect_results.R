library(data.table)
library(brms)
CPUS_DEFAULT <- 7
cpus_total <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

if(is.na(cpus_total)){
  message('Not running as slurm array, setting ROI to some value: ', CPUS_DEFAULT)
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
  setDTthreads(1)
  print(fit_files_dt$file)
  cfit <- read_fit_files(fit_files_dt)
  
  model_info <- unique(fit_files_dt[, c('roi', 'sub', 'ns', 'pred')])
  
  #We don't compute our proportions of variance or alpha from predictor models,
  #except for the model in which we include all predictors.
  if(grepl('-pred', model_info$pred) & model_info$pred != '-pred_all'){
    if(grepl('-pred_time', model_info$pred)){
      variable <- 'b_TIMECENTER'
    } else {
      variable <- grep('WCEN', variables(cfit), value = TRUE)  
    }
    post_sum <- data.table(
      brms::posterior_summary(cfit, probs = c(.025, .975, 
                                              round(abs(c(0, -1) + .05/421/2), 5)), 
                              robust = TRUE, variable = variable),
      keep.rownames = TRUE)
    post_sum[, c('nvoxel', 'roi', 'sub', 'ns', 'pred') := c(NA, model_info)]
    r <- list(post_sum = post_sum, density_df = NULL)
  } else {
    vc <- VarCorr(cfit, summary = FALSE)
    tau_pi <- vc$`id:sess`$sd^2
    tau_beta <- vc$id$sd^2
    sigma2 <- vc$residual__$sd^2
    
    nvoxel <- length(cfit$data$y)/dim(unique(cfit$data[, c('sess', 'id:sess')]))[[1]]
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
    post_sum[, c('nvoxel', 'roi', 'sub', 'ns', 'pred') := c(nvoxel, model_info)]
    r <- list(post_sum = post_sum, density_df = d_df)
  }
  
  return(r)
}

gsub_regex <- 'alpha-fearGTcalm-ROI([0-9]{1,3})(-more_\\d{2})*(-sub)*(-ns)*(-pred_\\w+)*.rds'

alpha_fit_files <- data.table(file = dir('fits', pattern = 'alpha-fearGTcalm.*rds'))
alpha_fit_files[, c('roi', 'more', 'sub', 'ns', 'pred') := 
                  lapply(list('\\1', '\\2', '\\3', '\\4', '\\5'), gsub, pattern = gsub_regex, x = file)]
filenames <- split(split(alpha_fit_files, 
                         by = c('roi', 'sub', 'ns', 'pred'),
                         drop = TRUE), 1:cpus_total)
rel_fn <- 'reliability_results.rds'
bd_fn <- 'boundary_density.rds'
if(!file.exists(rel_fn)){
  cl <- parallel::makeCluster(cpus_total)
  
  parallel::clusterExport(cl, c('compute_reliability', 'read_fit_files'))
  nada <- parallel::clusterEvalQ(cl, {library(data.table); library(brms); setDTthreads(1)})
  system.time({reliability_list <- parallel::parLapply(cl = cl, filenames, function(some_fns){
    lapply(some_fns, compute_reliability)
  })})
  
  parallel::stopCluster(cl)
  
  #posterior summaries of the various reliability metrics
  reliability_list_ul <- unlist(reliability_list, recursive = FALSE)
  reliability_dt <- rbindlist(lapply(reliability_list_ul, 
                                     `[[`, 'post_sum'), 
                              fill = TRUE)
  #computed densities of the reliability metrices posteriors
  boundary_density_dt <- rbindlist(lapply(reliability_list_ul, 
                                          function(rel){
                                            rel$density_df$roi <- unique(rel$post_sum$roi)[[1]]
                                            rel$density_df$sub <- unique(rel$post_sum$sub)[[1]]
                                            rel$density_df$ns <- unique(rel$post_sum$ns)[[1]]
                                            return(rel$density_df)
                                          }), fill = TRUE)
  # ggplot(boundary_density_dt[sub == '' & ns == '-ns'], aes(x = x, y = y, group = roi)) + 
  #   geom_line(alpha = .2) + 
  #   facet_wrap(~ name, scales = 'free')
  saveRDS(reliability_dt, rel_fn)
  saveRDS(boundary_density_dt, bd_fn)
} else {
  reliability_dt <- readRDS(rel_fn)
  boundary_density_dt <- readRDS(bd_fn)
}

if(INTERACTIVE){
  #average voxels per parcel
  reliability_dt[sub == '' & ns == '-ns' & rn == 'prop_var' & pred == ''][, .(nvoxel_median = quantile(nvoxel, .5),
                                                                              nvoxel_25 = quantile(nvoxel, .25),
                                                                              nvoxel_75 = quantile(nvoxel, .75))]
  
  gsub_regex <- 'alpha-fearGTcalm-ROI([0-9]{1,3})(-more_\\d{2})*(-sub)*(-ns)*(-pred_\\w+)*-psamp.rds'
  
  postsamp_files <- data.table(file = dir('post_samples', pattern = 'alpha-fearGTcalm.*rds'))
  postsamp_files[, c('roi', 'more', 'sub', 'ns', 'pred') := 
                   lapply(list('\\1', '\\2', '\\3', '\\4', '\\5'), gsub, pattern = gsub_regex, x = file)]
  
  postsamp_filenames <- split(postsamp_files, 1:cpus_total)
  
  prop_var_post_samps_fn <- 'prop_var_post_samps.rds'
  if(!file.exists(prop_var_post_samps_fn)){
    cl <- parallel::makeCluster(cpus_total)
    
    nada <- parallel::clusterEvalQ(cl, {library(data.table); setDTthreads(1)})
    system.time({post_samps_list <- parallel::parLapply(cl = cl, postsamp_filenames, function(some_fns){
      lapply(1:dim(some_fns)[[1]], function(i){
        d <- data.table(readRDS(file.path('post_samples', some_fns[i, file])))
        d[, roi := unique(some_fns[i, roi])]
        return(d)
      })
    })})
    
    parallel::stopCluster(cl)
    post_samps_dt <- rbindlist(unlist(post_samps_list, recursive = FALSE))
    saveRDS(post_samps_dt, prop_var_post_samps_fn)
  } else {
    post_samps_dt <- readRDS(prop_var_post_samps_fn)
  }
  
  setnames(post_samps_dt, 
           old = c('prop_var', 'prop_var_id'),
           new = c('sess', 'id'))
  post_samps_dt[, id_sess := id + sess]
  post_samps_dt[, group := fifelse(roi <= 400, 'cortical',
                                   fifelse(roi %in% c(401, 403, 412, 414), 'reference',
                                           fifelse(roi %in% c(404:407, 409:411, 415:421), 'subcortial',
                                                   NA_character_)))]
  
  post_samps_dt_split <- split(post_samps_dt, post_samps_dt$group)
  
  brain_group_prop_var <- rbindlist(lapply(post_samps_dt_split, function(adt){
    rbindlist(lapply(c('id', 'id_sess'), function(p){
      d <- data.table(posterior_summary(adt[, get(p)], robust = TRUE))
      d[, param:= p]
      return(d)
    }))
  }), id = 'group')
  brain_group_prop_var_dens <- rbindlist(lapply(post_samps_dt_split, function(adt){
    rbindlist(lapply(c('id', 'id_sess'), function(p){
      d <- as.data.table(density(adt[, get(p)], n = 512)[c('x','y')])
      d[, param:= p]
      return(d)
    }))
  }), id = 'group')
  brain_group_prop_var[, Est.Error := NULL]
  brain_group_prop_var_CI <- brain_group_prop_var[, c('group', 'Q2.5', 'Q97.5', 'param')]
  brain_group_prop_var <-  brain_group_prop_var[, c('group', 'Estimate', 'param')]
  brain_group_prop_var_w <- dcast(brain_group_prop_var, group ~ param, value.var = 'Estimate')
  setnames(brain_group_prop_var_w, c('group', 'id_xmax', 'sess_xmax'))
  brain_group_prop_var_w[, c('id_xmin', 'sess_xmin', 'err_xmin', 'err_xmax') :=
                           list(0, id_xmax, sess_xmax, 1)]
  brain_group_prop_var_w_l <- melt(brain_group_prop_var_w, id.vars = c('group'))
  brain_group_prop_var_w_l[, c('source', 'stat', 'variable') :=
                             list(gsub('(\\w+)_(\\w+)', '\\1', variable),
                                  gsub('(\\w+)_(\\w+)', '\\2', variable),
                                  NULL)]
  brain_group_prop_var_forplot <- dcast(brain_group_prop_var_w_l, ... ~ stat, value.var = 'value')
  region_groups <- c('cortical' = 3,
                     'reference' = 1,
                     'subcortial' = 2)
  brain_group_prop_var_forplot[, group_num := region_groups[group]]
  brain_group_prop_var_CI[, group_num := region_groups[group]]
  brain_group_prop_var_CI[, param := factor(param, levels = c('id', 'id_sess'), labels = c('id', 'sess'))]
  brain_group_prop_var_dens[, param := factor(param, levels = c('id', 'id_sess'), labels = c('id', 'sess'))]
  brain_group_prop_var_dens[, group_num := region_groups[group]]
  reliability_dt[, roi_num := as.numeric(roi)]
  boundary_density_dt[, roi_num := as.numeric(roi)]
  #Count
  
  #should be 0
  reliability_dt[
    rn == 'prop_var' & sub == '-sub'& ns == '-ns' & pred == '', ][
      data.table(roi = sprintf('%03d', 1:421)), on = 'roi'][
        is.na(Estimate), roi]
  
  #should be 421
  reliability_dt[
    rn == 'prop_var' & sub == '-sub' & ns == '-ns' & pred == '', ][
      data.table(roi = sprintf('%03d', 1:421)), on = 'roi'][
        !is.na(Estimate), .N]
  
  reliability_dt[
    rn == 'prop_var' & sub == '' & ns == '-ns' & pred == '', ][
      data.table(roi = sprintf('%03d', 1:421)), on = 'roi'][
        !is.na(Estimate), .N]
  
  
  reliability_dt[rn == 'prop_var' & sub == '-sub' & ns == '-ns'  & pred == '' & nvoxel < 15]
  
  #Plot 
  
  ho_sub_labels <- fread('HO_sub.csv', col.names = c('label', 'roi_num'))
  ho_sub_labels[, roi_num := roi_num + 400]
  reliability_dt_labels <- ho_sub_labels[reliability_dt, on = 'roi_num']
  reliability_dt_labels[, label_text := fifelse((Estimate < .22 & rn == 'prop_var') |
                                                  (Estimate < .75 & rn == 'alpha'), label, NA_character_)]
  
  reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & pred == '' & roi > 400, c('label', 'label_text', 'Estimate')]
  
  reliability_dt_labels[, roi_num_reverse := -1*(roi_num - 401)]
  
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  
  reliability_dt_labels_prop_w <- dcast(reliability_dt_labels[rn %in% c('prop_var_id', 'idsess_right') & 
                                                                ns == '-ns' & sub == '-sub' & pred == '' & 
                                                                !roi_num %in% c(402, 413), #exclude cerebral cortex from HO
                                                              c('rn', 'roi_num', 'Estimate', 'roi_num_reverse')],
                                        ... ~ rn, value.var = 'Estimate')
  
  setnames(reliability_dt_labels_prop_w, 
           old = c('idsess_right', 'prop_var_id'), 
           new = c('sess_xmax', 'id_xmax'))
  reliability_dt_labels_prop_w[, c('id_xmin', 'sess_xmin', 'err_xmin', 'err_xmax') :=
                                 list(0, id_xmax, sess_xmax, 1)]
  reliability_dt_labels_prop_l <- melt(reliability_dt_labels_prop_w, id.vars = c('roi_num', 'roi_num_reverse'))
  reliability_dt_labels_prop_l[, c('source', 'stat', 'variable') :=
                                 list(gsub('(\\w+)_(\\w+)', '\\1', variable),
                                      gsub('(\\w+)_(\\w+)', '\\2', variable),
                                      NULL)]
  reliability_dt_labels_prop_l_w <- dcast(reliability_dt_labels_prop_l, ... ~ stat, value.var = 'value')
  # Variance proportions
  
  source_order <- c(2, 3, 1)
  proportion_by_ROI_plot <- ggplot(reliability_dt_labels_prop_l_w, aes(x = NULL, y = roi_num_reverse)) + 
    geom_ribbon(aes(xmin = xmin, xmax = xmax, fill = source, group = source), 
                alpha = 1) + 
    annotate('rect', ymin = -19, ymax = 0, xmin = -.1, xmax = 1.1, fill = 'black', alpha = .1) + 
    coord_cartesian(xlim = c(0,1)) +
    labs(x = 'Cumulative Proportion', y = 'Parcel / Region') + 
    scale_fill_manual(breaks = c('err', 'id', 'sess')[source_order], 
                      values = c('#FFF8F3', '#EDC12C', '#71A1D4')[source_order],
                      labels = c('Error', 'Participant means', 'Within-session means')[source_order],
                      name = 'Source of variance') + 
    scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 1, seq(80, 320, 80), 'Schaefer-400')) +
    theme_minimal()
  print(proportion_by_ROI_plot)
  ggsave('prop_id_sess_err_variance_by_ROI.svg', plot = proportion_by_ROI_plot, width = 8, height = 5, dpi = 150)
  
  # Variance proportions by region-groups
  brain_group_prop_var_dens[, y_pct := y/max(y), by = c('group', 'param')]
  brain_group_prop_var_dens[, x := fifelse(x < 0, NA_real_, x)]
  grouped_prop_var <- ggplot(brain_group_prop_var_forplot[!is.na(group)], aes(x = NULL, y = group_num)) + 
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = group_num - .5, ymax = group_num + .5, group = source, fill = source),
              alpha = .85) + 
    # geom_rect(data = brain_group_prop_var_CI,
    #           aes(xmin = Q2.5, xmax = Q97.5, ymin = group_num - .45, ymax = group_num + .45, group = param, fill = param),
    #           alpha = .5) +
    geom_segment(data = brain_group_prop_var_CI,
                 aes(x = Q2.5, xend = Q97.5, y = group_num, yend = group_num, color = param),
                 size = 1) +
    # geom_segment(data = brain_group_prop_var_dens,
    #              aes(x = x, xend = x, y = group_num + .5, yend = group_num - .5, alpha = y_pct,
    #                  group = param, color = param)) +
    scale_fill_manual(breaks = c('err', 'id', 'sess')[source_order],
                      values = c('#FFF8F3', '#EDC12C', '#71A1D4')[source_order],
                      labels = c('Error', 'Participant means', 'Within-session means')[source_order],
                      name = 'Source of variance', aesthetics = c('fill', 'color')) +
    scale_y_continuous(breaks = c(1, 2, 3), labels = names(region_groups)[order(region_groups)]) + 
    scale_alpha_continuous(range = c(0, 1)) + 
    coord_cartesian(xlim = c(0, 1)) + 
    theme_minimal() +
    theme(axis.line = element_blank(), panel.grid = element_blank()) + 
    labs(x = 'Cumulative proportion variance', y = '')
  print(grouped_prop_var)
  ggsave('prop_id_sess_err_variance_grouped.svg', plot = grouped_prop_var, width = 8, height = 5, dpi = 150)
  
  ggplot(brain_group_prop_var_dens, aes(x = x, y = y_pct)) + 
    geom_line() + 
    facet_wrap(~ group + param, ncol = 2)
  
  ggplot(post_samps_dt, aes(x = id)) + 
    geom_density() + 
    facet_wrap(~ group, ncol = 2)
  
  ggplot(post_samps_dt, aes(x = id_sess)) + 
    geom_density() + 
    facet_wrap(~ group, ncol = 2)
  
  
  # Mean (with CI) proportions by regions
  
  tmp <- data.table(readRDS('post_samples/alpha-fearGTcalm-ROI001-sub-ns-psamp.rds'))
  
  # Between ICC 
  
  between_icc_plot <- ggplot(reliability_dt_labels[rn == 'prop_var_id' & sub == '-sub' & 
                                                     ns == '-ns' & pred == '' &
                                                     ! grepl('Cerebral Cortex', label)], 
                             aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
    geom_rect(ymin = -19, ymax = 0, xmin = 0, xmax = 1, fill = 'lightgrey') + 
    geom_errorbar(alpha = .5, width = 0, size = .25) +
    geom_point(size = .5) + 
    scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 1, seq(80, 320, 80), 'Schaefer-400')) +
    theme_minimal() + 
    coord_cartesian(x = c(0, 1)) + 
    labs(x = 'Proportion variance due to participant-session', y = 'Parcel')
  print(between_icc_plot)
  
  ns_icc_range <- range(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & ! grepl('Cerebral Cortex', label), Estimate])
  
  within_icc_plot <- ggplot(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' &
                                                  ns == '-ns' & pred == '' &
                                                    ! grepl('Cerebral Cortex', label)], 
                            aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
    geom_rect(ymin = -19, ymax = 0, xmin = 0, xmax = 1, fill = 'lightgrey') + 
    geom_errorbar(alpha = .5, width = 0, size = .10) +
    geom_point(size = 1, aes(color = Estimate)) + 
    scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 1, seq(80, 320, 80), 'Schaefer-400')) +
    theme_minimal() + 
    coord_cartesian(x = c(0, 1)) + 
    geom_label_repel(aes(label = label_text), box.padding = 1, max.overlaps = Inf, min.segment.length = 0, 
                     nudge_x = -1, nudge_y = 75, fill = alpha('white', .65)) + 
    scale_color_gradientn(colors = c('#ff0000', '#ff0000', '#ffff00', '#ffff00'), 
                          values = c(0, 0.10, 0.80, 1),
                          limits = c(0, 1),
                          rescaler = function(x, from = y) x) + 
    labs(x = 'Proportion variance due to participant-session', y = 'Parcel')
  print(within_icc_plot)
  # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ggsave('prop_var_id-session-ns.svg', width = 8, height = 5, dpi = 150)
  
  within_icc_plot <- ggplot(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' &
                                                    ns == '-ns' & pred == '' & 
                                                    ! grepl('Cerebral Cortex', label)], 
                            aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
    geom_rect(ymin = -19, ymax = 0, xmin = 0, xmax = 1, fill = 'lightgrey') + 
    geom_errorbar(alpha = .5, width = 0, size = .25) +
    geom_point(size = .5) + 
    scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 1, seq(80, 320, 80), 'Schaefer-400')) +
    theme_minimal() + 
    coord_cartesian(x = c(0, 1)) + 
    geom_label_repel(aes(label = label_text), box.padding = 1, max.overlaps = Inf, min.segment.length = 0, 
                     nudge_x = -1, nudge_y = 75, fill = alpha('white', .65)) + 
    labs(x = 'Proportion variance due to participant-session', y = 'Parcel')
  print(within_icc_plot)
  # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ggsave('prop_var_id-session.svg', width = 8, height = 5, dpi = 150)
  
  within_icc_plot <- ggplot(reliability_dt_labels[rn == 'alpha' & sub == '-sub' &
                                                    ns == '-ns' & pred == '' &
                                                    ! grepl('Cerebral Cortex', label)], 
                            aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
    geom_rect(ymin = -19, ymax = 0, xmin = 0, xmax = 1, fill = 'lightgrey') + 
    geom_errorbar(alpha = .5, width = 0, size = .25) +
    geom_point(size = .5) + 
    scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 1, seq(80, 320, 80), 'Schaefer-400')) +
    theme_minimal() + 
    coord_cartesian(x = c(0, 1)) + 
    geom_label_repel(aes(label = label_text), box.padding = 1, max.overlaps = Inf, min.segment.length = 0, 
                     nudge_x = -1, nudge_y = 75, fill = alpha('white', .65)) + 
    labs(x = 'Proportion variance due to participant-session', y = 'Parcel')
  print(within_icc_plot)
  # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ggsave('alpha_id-session-ns.svg', width = 8, height = 5, dpi = 150)
  
  within_icc_plot <- ggplot(reliability_dt_labels[rn == 'alpha' & sub == '' &
                                                    ns == '-ns' & pred == '' &
                                                    ! grepl('Cerebral Cortex', label)], 
                            aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
    geom_rect(ymin = -19, ymax = 0, xmin = 0, xmax = 1, fill = 'lightgrey') + 
    geom_errorbar(alpha = .5, width = 0, size = .25) +
    geom_point(size = .5) + 
    scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 1, seq(80, 320, 80), 'Schaefer-400')) +
    theme_minimal() + 
    coord_cartesian(x = c(0, 1)) + 
    geom_label_repel(aes(label = label_text), box.padding = 1, max.overlaps = Inf, min.segment.length = 0, 
                     nudge_x = -1, nudge_y = 75, fill = alpha('white', .65)) + 
    labs(x = 'Proportion variance due to participant-session', y = 'Parcel')
  print(within_icc_plot)
  # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ggsave('alpha-allvox_id-session-ns.svg', width = 8, height = 5, dpi = 150)

  reliability_dt_labels[, credible := (sign(Q0.006) * sign(Q99.994)) == 1]
  reliability_dt_labels[, credible_95 := (sign(Q2.5) * sign(Q97.5)) == 1]
  
  tau_pi_delta_prop_plot <- ggplot(reliability_dt_labels[rn == 'tau_pi_delta_prop' & sub == '-sub' &
                                                    ns == '-ns' & pred == '-pred_all' &
                                                    ! grepl('Cerebral Cortex', label)], 
                            aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
    geom_rect(ymin = -19, ymax = 0, xmin = -1, xmax = 1, fill = 'lightgrey') + 
    geom_errorbar(alpha = .5, width = 0, size = .25, aes(color = credible_95)) +
    geom_point(size = .5, aes(color = credible_95)) + 
    scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 1, seq(80, 320, 80), 'Schaefer-400')) +
    theme_minimal() + 
    coord_cartesian(x = c(-1, 1)) + 
    scale_color_manual(breaks = c(TRUE, FALSE), values = c('black', 'grey')) +
    geom_label_repel(aes(label = label_text), box.padding = 1, max.overlaps = Inf, min.segment.length = 0, 
                     nudge_x = -1, nudge_y = 75, fill = alpha('white', .65)) + 
    labs(x = 'Proportion variance due to participant-session', y = 'Parcel')
  print(tau_pi_delta_prop_plot)
  
  # which(! 1:400 %in% reliability_dt_labels[rn == 'tau_pi_delta_prop' & sub == '-sub' &
  #                                            ns == '-ns' & pred == '-pred_all' &
  #                                            roi_num <= 400, sort(roi_num)])
  # which(! 1:400 %in% reliability_dt_labels[rn == 'b_TIMECENTER' & sub == '-sub' &
  #                                            ns == '-ns' &
  #                                            roi_num <= 400, sort(roi_num)])
  
  plot_preds <- function(x, variable, outfile){
    plot <- ggplot(x[rn == variable & sub == '-sub' &
                       ns == '-ns' &
                       ! grepl('Cerebral Cortex', label)], 
                   aes(x = Estimate, xmin = Q0.006, xmax = Q99.994, y = roi_num_reverse)) +
      geom_rect(ymin = -19, ymax = 0, xmin = -1, xmax = 1, fill = 'lightgrey') + 
      geom_errorbar(alpha = .5, width = 0, size = .25, aes(color = credible)) +
      geom_point(size = .5, aes(color = credible)) + 
      scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 1, seq(80, 320, 80), 'Schaefer-400')) +
      theme_minimal() + 
      scale_color_manual(breaks = c(TRUE, FALSE), values = c('black', 'grey')) +
      labs(x = 'Estimate (99.988% CI)', y = 'Parcel') + 
      theme(legend.position = 'none')
    ggsave(outfile, width = 4, height = 5, dpi = 150)
    return(plot)    
  }
  #unique(reliability_dt_labels$rn)
  plot_preds(x = reliability_dt_labels, variable = 'b_TIMECENTER', outfile = 'time-ns-sub.svg')
  plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', outfile = 'episodic_stress-ns-sub.svg')
  plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', outfile = 'chronic_stress-ns-sub.svg')
  plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', outfile = 'negative_affect-ns-sub.svg')
  plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR', outfile = 'sleep-ns-sub.svg')
  # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  write_cifti_text <- function(x, variable, outfile){
    #x <- copy(reliability_dt_labels)
    x <- x[rn == variable & sub == '-sub' &
             ns == '-ns' &
             ! grepl('Cerebral Cortex', label), 
           c('roi_num', 'Estimate', 'credible')]
    setkey(x, roi_num)
    rois <- data.table(roi_num = 1:400)
    out_data <- x[rois]
    out_data[! credible %in% TRUE, Estimate := 0]
    write.table(out_data[, Estimate], file = outfile, row.names = FALSE, col.names = FALSE)
  }
  
  write_cifti_text(x = reliability_dt_labels, variable = 'b_TIMECENTER', outfile = 'time-ns-sub.txt')
  write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', outfile = 'episodic_stress-ns-sub.txt')
  write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', outfile = 'chronic_stress-ns-sub.txt')
  write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', outfile = 'negative_affect-ns-sub.txt')
  write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR', outfile = 'sleep-ns-sub.txt')
  
  range_text <- function(v, caption, v_range_code = "%0.2f", v_iqr_code = "%0.2f-%0.2f"){
    v_range <- paste(sprintf(v_range_code, range(v)), collapse = ', ') 
    v_iqr <- sprintf(v_iqr_code, quantile(v, .25), quantile(v, .75))
    N <- length(v)
    s <- sprintf(paste0(caption, " [%s] (IQR %s; N = %i)."), v_range, v_iqr, N)
    return(s)
  }
  
  #prop var sess-id not-smoothed
  cat(paste(
    c(
      range_text(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & roi_num %in% c(401, 403, 412, 414), Estimate],
                 "Range of ICC_within-person across ventricle, and white matter:"),
      range_text(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & is.na(label_text) & roi_num %in% 401:421, Estimate],
                 "Range of ICC_within-person across subcortial (no brain-stem, ventricle, or white matter):"),
      range_text(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & roi_num %in% 1:400, Estimate],
                 "Range of ICC_within-person across Sch400 cortical regions:")#,
      # range_text(oldicc_ests_meta_cont_cope_dt[roi_tag == 'HO' & ROI %in% HO_sub_cortex, Estimate],
      #            "Range of ICC_between-person across HO cortex:")
    ),
    collapse = '\n'))
  
  
  #prop var sess-id smoothed
  cat(paste(
    c(
      range_text(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '' & roi_num %in% c(401, 403, 412, 414), Estimate],
                 "Range of ICC_within-person across ventricle, and white matter:"),
      range_text(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '' & roi_num %in% c(404:407, 409:411, 415:421), Estimate],
                 "Range of ICC_within-person across subcortial (no brain-stem, ventricle, or white matter):"),
      range_text(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '' & roi_num %in% 1:400, Estimate],
                 "Range of ICC_within-person across Sch400 cortical regions:")#,
      # range_text(oldicc_ests_meta_cont_cope_dt[roi_tag == 'HO' & ROI %in% HO_sub_cortex, Estimate],
      #            "Range of ICC_between-person across HO cortex:")
    ),
    collapse = '\n'))
  
  
  
  #alpha, 15 voxels not-smoothed
  cat(paste(
    c(
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '-sub' & ns == '-ns' & roi_num %in% c(401, 403, 408, 412, 414), Estimate],
                 "Range of Alpha across brain-stem, ventricle, and white matter:"),
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '-sub' & ns == '-ns' & roi_num %in% c(404:407, 409:411, 415:421), Estimate],
                 "Range of Alpha across subcortial (no brain-stem, ventricle, or white matter):"),
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '-sub' & ns == '-ns' & roi_num %in% 1:400, Estimate],
                 "Range of Alpha across Sch400 cortical regions:")#,
      # range_text(oldicc_ests_meta_cont_cope_dt[roi_tag == 'HO' & ROI %in% HO_sub_cortex, Estimate],
      #            "Range of ICC_between-person across HO cortex:")
    ),
    collapse = '\n'))
  
  #alpha, 15 voxels smoothed
  cat(paste(
    c(
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '-sub' & ns == '' & roi_num %in% c(401, 403, 408, 412, 414), Estimate],
                 "Range of Alpha across brain-stem, ventricle, and white matter:"),
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '-sub' & ns == '' & roi_num %in% c(404:407, 409:411, 415:421), Estimate],
                 "Range of Alpha across subcortial (no brain-stem, ventricle, or white matter):"),
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '-sub' & ns == '' & roi_num %in% 1:400, Estimate],
                 "Range of Alpha across Sch400 cortical regions:")#,
      # range_text(oldicc_ests_meta_cont_cope_dt[roi_tag == 'HO' & ROI %in% HO_sub_cortex, Estimate],
      #            "Range of ICC_between-person across HO cortex:")
    ),
    collapse = '\n'))
  
  #alpha, all voxels unsmoothed
  cat(paste(
    c(
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '' & ns == '-ns' & roi_num %in% c(401, 403, 408, 412, 414), Estimate],
                 "Range of Alpha across brain-stem, ventricle, and white matter:"),
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '' & ns == '-ns' & roi_num %in% c(404:407, 409:411, 415:421), Estimate],
                 "Range of Alpha across subcortial (no brain-stem, ventricle, or white matter):"),
      range_text(reliability_dt_labels[rn == 'alpha' & sub == '' & ns == '-ns' & roi_num %in% 1:400, Estimate],
                 "Range of Alpha across Sch400 cortical regions:")#,
      # range_text(oldicc_ests_meta_cont_cope_dt[roi_tag == 'HO' & ROI %in% HO_sub_cortex, Estimate],
      #            "Range of ICC_between-person across HO cortex:")
    ),
    collapse = '\n'))
  
  #range
  #"bad" regions
  cat(paste0("Range of ICC_within-person across brainstem, ventricle, and white matter: [", 
             paste(sprintf("%0.2f", 
                           range(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & !is.na(label_text), 
                                                       Estimate])),
                   collapse = ', '),
             "]"))
  #no cortical parcels and subcortex only
  cat(paste0("Range of ICC_within-person across subcortial (no brain-stem, ventricle, or white matter) and cortical regions: [", 
             paste(sprintf("%0.2f", 
                           range(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & is.na(label_text), 
                                                       Estimate])),
                   collapse = ', '),
             "]"))
  #HO cortex
  cat(paste0("Range of ICC_within-person across HO cortex: [", 
             paste(sprintf("%0.2f", 
                           range(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & !is.na(label_text) & grepl('Cortex', label_text), 
                                                       Estimate])),
                   collapse = ', '),
             "]"))
  
  #correspondence between 15 voxel sample and all voxels
  #unsmoothed
  ggplot(dcast(reliability_dt[rn == 'prop_var' & ns == '-ns', c('Estimate', 'sub', 'roi')], roi ~ sub, value.var = 'Estimate'),
         aes(x = V1, y = `-sub`)) + 
    geom_point(alpha = .4, size = .5) + 
    geom_abline(slope = 1, intercept = 0) + 
    coord_cartesian(x = c(0, 1), y = c(0, 1)) +
    labs(title = 'Unsmoothed', x = 'All voxels', y = 'Subsampled (15 voxels)') +
    
    ggplot(dcast(reliability_dt[sub == '-sub' & ns == '-ns', c('Estimate', 'rn', 'roi')], roi ~ rn, value.var = 'Estimate'),
           aes(y = alpha, x = prop_var)) + 
    geom_point(alpha = .4, size = .5) + 
    coord_cartesian(x = c(0, 1), y = c(0, 1)) +
    labs(x = 'Proportion sub-sess variance', y = ~paste(alpha, ' (15 voxels)')) +
    
    #smoothed
    ggplot(dcast(reliability_dt[rn == 'prop_var' & ns == '' , c('Estimate', 'sub', 'roi')], roi ~ sub, value.var = 'Estimate'),
           aes(x = V1, y = `-sub`)) + 
    geom_point(alpha = .4, size = .5) + 
    geom_abline(slope = 1, intercept = 0) +
    coord_cartesian(x = c(0, 1), y = c(0, 1)) +
    labs(title = 'Smoothed', x = 'All voxels', y = 'Subsampled (15 voxels)') +
    
    ggplot(dcast(reliability_dt[sub == '-sub' & ns == '' , c('Estimate', 'rn', 'roi')], roi ~ rn, value.var = 'Estimate'),
           aes(y = alpha, x = prop_var)) + 
    geom_point(alpha = .4, size = .5) + 
    coord_cartesian(x = c(0, 1), y = c(0, 1)) +
    labs(x = 'Proportion sub-sess variance', y = ~paste(alpha, ' (15 voxels)')) + 
    patchwork::plot_layout(design = '
  AB
  CD')
  
  psych::corr.test(dcast(reliability_dt[rn == 'prop_var' & ns == '-ns', c('Estimate', 'sub', 'roi')], 
                         roi ~ sub, value.var = 'Estimate')[, -'roi'])
  
  #correspondence between prop_var and prop_var_id
  ggplot(dcast(reliability_dt[(rn == 'prop_var' | rn == 'prop_var_id') & ns == '-ns' & sub == '-sub', c('Estimate', 'rn', 'roi')], roi ~ rn, value.var = 'Estimate'),
         aes(x = prop_var, y = prop_var_id)) + 
    geom_point(alpha = .4, size = .5) + 
    geom_line(stat = 'smooth')
  
  ##Write so we can convert to NII and surface
  
  sch400_HO_dump <- fread('Sch400_HO_maskdump.tsv', sep = ' ',  col.names = c('i', 'j', 'k', 'roi_num'), header = FALSE)
  
  prop_var_sub_out_dt <- reliability_dt[rn == 'prop_var' & sub == '-sub' & ns == '-ns', ][sch400_HO_dump, on = 'roi_num']
  prop_var_sub_out_dt[is.na(Estimate), c('Estimate', 'Q2.5', 'Q97.5') := 0]
  
  readr::write_delim(prop_var_sub_out_dt[, c('i', 'j', 'k', 'Estimate')], 'prop_id-sess_var_15voxsub.tsv', delim = ' ', col_names = FALSE)
}
