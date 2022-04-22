library(data.table)
library(brms)
CPUS_DEFAULT <- 7
cpus_total <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))

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
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ggsegSchaefer)
  library(ggseg)
  
  plot_preds <- function(x, variable, outfile){
    plot <- ggplot(x[rn == variable & sub == '-sub' &
                       ns == '-ns' &
                       ! grepl('Cerebral Cortex', label)], 
                   aes(x = Estimate, xmin = Q0.006, xmax = Q99.994, y = roi_num_reverse)) +
      geom_rect(ymin = -19, ymax = 0, xmin = -1, xmax = 1, fill = 'lightgrey') + 
      geom_errorbar(alpha = .5, width = 0, size = .25, aes(color = credible)) +
      geom_point(size = .5, aes(color = credible)) + 
      scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 400, seq(320, 80, -80), 'Schaefer-400 1')) +
      theme_minimal() + 
      scale_color_manual(breaks = c(TRUE, FALSE), values = c('black', 'grey')) +
      labs(x = 'Estimate (99.988% CI)', y = 'Parcel') + 
      theme(legend.position = 'none')
    ggsave(outfile, width = 4, height = 5, dpi = 150)
    return(plot)    
  }
  
  write_cifti_text <- function(x, variable, outfile){
    #x <- copy(reliability_dt_labels)
    x <- x[rn == variable & sub == '-sub' &
             ns == '-ns' &
             !grepl('Cerebral Cortex', label), 
           c('roi_num', 'Estimate', 'credible')]
    setkey(x, roi_num)
    rois <- data.table(roi_num = 1:400)
    out_data <- x[rois]
    out_data[! credible %in% TRUE, Estimate := 0]
    write.table(out_data[, Estimate], file = outfile, row.names = FALSE, col.names = FALSE)
    return(TRUE)
  }
  
  get_selected_parcels <- function(x, variable, labels){
    # x <- copy(reliability_dt_labels)
    # variable = 'b_TIMECENTER'
    # labels = copy(ho_sub_labels)
    x <- x[rn == variable & sub == '-sub' &
             ns == '-ns' &
             ! grepl('Cerebral Cortex', label), 
           c('roi_num', 'Estimate', 'credible', 'rn', 'roi', 'sub', 'ns', 'pred')]
    setkey(x, roi_num)
    setkey(labels, roi_num)
    x <- labels[x][credible %in% TRUE]
    return(x)
  }
  
  get_parcel_var_model_df <- function(x, variable, labels){
    x <- get_selected_parcels(x, variable, labels)
    x[, mod_fn := paste0('alpha-fearGTcalm-ROI', roi, sub, ns, pred, '.rds')]
    x[, c('abs_est', 'sign') := list(abs(Estimate), sign(Estimate))]
    setkey(x, sign, abs_est)
    min_med_max_rows <- function(subx){
      n <- dim(subx)[[1]]
      idx <- c(1, round(n/2), n)
      retx <- subx[idx]
      retx[, idx := idx]
      return(retx)
    }
    x <- x[, min_med_max_rows(.SD), by = sign]
    return(x)
  }
  
  range_text <- function(v, caption, v_range_code = "%0.2f", v_iqr_code = "%0.2f-%0.2f"){
    v_range <- paste(sprintf(v_range_code, range(v)), collapse = ', ') 
    v_iqr <- sprintf(v_iqr_code, quantile(v, .25), quantile(v, .75))
    N <- length(v)
    s <- sprintf(paste0(caption, " [%s] (IQR %s; N = %i)."), v_range, v_iqr, N)
    return(s)
  }
  
  create_cifti <- function(text_file, cifti_template, cifti_out){
    module_cmd <- 'module load connectome-workbench/1.3.2-fasrc01'
    wb_cmd <- paste('wb_command', '-cifti-convert', '-from-text', text_file, cifti_template, cifti_out, sep = ' ')
    if(!file.exists(cifti_template)){
      cifti_dscalar <- gsub('pscalar', 'dscalar', cifti_template)
      cifti_dlabel <- gsub('pscalar', 'dlabel', cifti_template)
      prereqs <- c(cifti_dscalar, cifti_dlabel)
      template_prereq_exist <- file.exists(prereqs)
      if(!all(template_prereq_exist)){
        stop(paste('Not all template prerequisites exist:', 
                   paste(prereqs[!template_prereq_exist], collapse = '\n')))
      }
      mk_template_cmd <- paste('wb_command', '-cifti-parcellate', cifti_dscalar, cifti_dlabel, 'COLUMN', cifti_template, sep = ' ')
      mk_template_cmd <- paste(module_cmd, mk_template_cmd, sep = '; ')
      system(mk_template_cmd)
    }
    cmd <- paste(module_cmd, wb_cmd, sep = '; ')
    system(cmd)
    return(TRUE)
  }
  
  plot_conditional_effect <- function(x, breaks = NA, labels = NA, ylim = c(-2, 2)){
    variable <- gsub('b_', '', x$rn)
    xmx <- readRDS(file.path('fits', x$mod_fn))
    sign <- x$sign
    colors <- c('#358fc9', '#fbaa2c')
    color <- colors[sign/2 + 1.5]
    
    CEdf <- conditional_effects(xmx, effects = variable, prob = .95, re_formula = NA)[[1]]
    pointsdf <- xmx$data
    
    p <- ggplot() + 
      theme_minimal() + 
      theme(panel.grid.minor = element_blank(), axis.text.y = element_blank()) + 
      labs(x = '', y = '') +
      stat_summary(data = pointsdf, aes_string(group = sprintf('interaction(id, %s)', variable),
                                               x = variable, y = 'y'),
                   fun = mean, fun.args = list(mult = 2),
                   alpha = .2, color = color, size = .2,
                   position = position_jitter(w = .25)) +
      geom_ribbon(data = CEdf, aes_string(x = variable, y = 'estimate__', 
                                          ymin = 'lower__', ymax = 'upper__'), 
                  alpha = .5, fill = color) + 
      geom_line(data = CEdf, aes_string(x = variable, y = 'estimate__'), color = color) +
      coord_cartesian(y = ylim)
    if(!is.na(breaks)){
      p <- p + scale_x_continuous(breaks = breaks, labels = labels)
    }
    if(sign < 0){
      p <- p + theme(axis.text.x = element_blank())
    }
    return(p)
  }
  
  make_multi_plot <- function(x, xlab, file = NULL, ...){
    p_list <- lapply(1:6, function(i){
      plot_conditional_effect(x = x[i], ...)
    })
    
    p <- wrap_plots(p_list, design = "123\n456")
    gt <- patchwork::patchworkGrob(p)
    pf <- gridExtra::grid.arrange(gt, left = "Fear - Calm", bottom = xlab)
    if(!is.null(file)){
      ggsave(file, plot = pf, width = 6, height = 4)
    }
    return(pf)
  }
  
  make_parcel_icons <- function(outdir, dpi = 300){
    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
    
    region_data <- as.data.frame(schaefer7_400$data)[!is.na(schaefer7_400$data$region), c('hemi', 'side', 'region')]
    nada <- lapply(1:dim(region_data)[[1]], \(x){
      x <- region_data[x, ]
      x$selected <- 'selected'
      
      outfile <- file.path(outdir, sprintf('%s_%s.png', x$region[[1]], x$side[[1]]))
      
      outplot <- ggplot(x) + 
        geom_brain(atlas = schaefer7_400,
                   hemi = x$hemi,
                   side = x$side,
                   aes(fill = selected, size = selected, alpha = selected, color = selected))  + 
        scale_fill_manual(breaks = c('selected'), values = c('black'), na.value = 'white') + 
        scale_size_manual(breaks = c('selected'), values = c(1), na.value = .05) +
        scale_alpha_manual(breaks = c('selected'), values = c(1), na.value = 0) +
        scale_color_manual(breaks = c('selected'), values = c(alpha('black', .5)), na.value = 'black') + 
        theme_minimal() + 
        theme(panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = unit(c(0, 0, 0, 0), units = 'in'),
              legend.position = 'none') + 
        coord_sf(expand = FALSE,
                 xlim = c(0, 310),
                 ylim = c(0, 210))
      ggsave(outfile, outplot, width = 77.5, height = 52.5, units = 'px', dpi = dpi, scale = dpi/72)
      return(TRUE)
    })
    return(dir(outdir))
  }
  
  write_rez_table <- function(x, variable, tables_dir = 'tables'){
    #x <- reliability_dt_labels
    #variable <- 'b_TIMECENTER'
    Sch7_400_regions <- fread('Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order_info.txt', header = FALSE, col.names = 'region', sep = ',')
    Sch7_400_regions <- Sch7_400_regions[seq(1, 800, 2)]
    Sch7_400_regions[, roi_num := 1:400]
    setkey(Sch7_400_regions, roi_num)
    x <- x[rn == variable & sub == '-sub' &
             ns == '-ns' &
             !grepl('Cerebral Cortex', label),
           c('roi_num', 'Estimate', 'Q0.006', 'Q99.994', 'credible')]
    setkey(x, roi_num)
    
    x <- x[Sch7_400_regions]
    this_region <- x[credible == TRUE]
    
    schdat <- as.data.table(schaefer7_400$data)
    setorder(schdat, region, side)
    schdat <- schdat[!duplicated(region)]
    this_region <- schdat[this_region, on = 'region']
    this_region[, Location := pander::pandoc.image.return(file.path(parcel_icon_dir, sprintf('%s_%s.png', region, side)))]
    this_region[, 'Parcel Label' := stringr::str_replace_all(stringr::str_remove(region, '7Networks_(RH|LH)_'), '_', ' ')]
    this_region[, Hemisphere := stringr::str_to_title(hemi)]
    this_region[, 'Estimate [99.988% CI]' := sprintf('% 0.3f [% 0.3f, % 0.3f]', Estimate, Q0.006, Q99.994)]
    setorder(this_region, roi_num)
    newtable <- pander::Pandoc$new(author = 'John C. Flournoy', title = variable)
    newtable$add(
      pandoc.table.return(
        this_region[, c('Hemisphere', 'Parcel Label', 'Estimate [99.988% CI]', 'Location')],
        split.table = Inf,
        justify = c('left', 'left', 'right', 'center')))
    newtable$format <- 'docx'
    if(!dir.exists(tables_dir)) dir.create(tables_dir)
    r1 <- newtable$export(variable, open = FALSE)
    r2 <- file.rename(paste0(variable, '.md'), file.path(tables_dir, paste0(variable, '.md')))
    r3 <- file.rename(paste0(variable, '.docx'), file.path(tables_dir, paste0(variable, '.docx')))
    return(all(c(r2, r3)))
  }
  
  parcel_icon_dir <- 'plots/parcel_icons'
  if(length(dir(parcel_icon_dir)) != 400) make_parcel_icons(parcel_icon_dir, dpi = 300)
  
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
  
  # tmp <- data.table(readRDS('post_samples/alpha-fearGTcalm-ROI001-sub-ns-psamp.rds'))
  
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
  # 
  #unique(reliability_dt_labels$rn)
  
  #make sleep neg:
  reliability_dt_labels_sleepneg <- copy(
    reliability_dt_labels[rn == 'b_WCEN_SLEEPDUR'])[, c('rn', 'pred', 'Estimate', 
                                                        'Q2.5', 'Q97.5', 
                                                        'Q0.006', 'Q99.994') := 
                                                      list(sprintf('%s_NEG', rn),
                                                           sprintf('%s_neg', pred),
                                                           -Estimate, -Q2.5, -Q97.5, 
                                                           -Q0.006, -Q99.994)]
  reliability_dt_labels <- data.table::rbindlist(list(reliability_dt_labels, reliability_dt_labels_sleepneg))
  plot_preds(x = reliability_dt_labels, variable = 'b_TIMECENTER', outfile = 'time-ns-sub.svg')
  plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', outfile = 'episodic_stress-ns-sub.svg')
  plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', outfile = 'chronic_stress-ns-sub.svg')
  plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', outfile = 'negative_affect-ns-sub.svg')
  plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR_NEG', outfile = 'sleep-ns-sub.svg')
  # theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  write_cifti_text(x = reliability_dt_labels, variable = 'b_TIMECENTER', outfile = 'time-ns-sub.txt')
  write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', outfile = 'episodic_stress-ns-sub.txt')
  write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', outfile = 'chronic_stress-ns-sub.txt')
  write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', outfile = 'negative_affect-ns-sub.txt')
  write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR_NEG', outfile = 'sleep-ns-sub.txt')
  
  create_cifti(text_file = 'time-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'time-ns-sub_Sch400.pscalar.nii')
  create_cifti(text_file = 'episodic_stress-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'episodic_stress-ns-sub_Sch400.pscalar.nii')
  create_cifti(text_file = 'chronic_stress-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'chronic_stress-ns-sub_Sch400.pscalar.nii')
  create_cifti(text_file = 'negative_affect-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'negative_affect-ns-sub_Sch400.pscalar.nii')
  create_cifti(text_file = 'sleep-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'sleep-ns-sub_Sch400.pscalar.nii')
  
  write_rez_table(reliability_dt_labels, 'b_TIMECENTER')
  write_rez_table(reliability_dt_labels, 'b_WCEN_EPISODICTOT')
  write_rez_table(reliability_dt_labels, 'b_WCEN_CHRONICSEV')
  write_rez_table(reliability_dt_labels, 'b_WCEN_PANASNEG')
  write_rez_table(reliability_dt_labels, 'b_WCEN_SLEEPDUR')
  
  model_fn_list_for_plots <- lapply(list('Month' = 'b_TIMECENTER',
                                         'Stressfull Life Events' = 'b_WCEN_EPISODICTOT',
                                         'Chronic Stress' = 'b_WCEN_CHRONICSEV',
                                         'Negative Mood' = 'b_WCEN_PANASNEG',
                                         'Less Sleep' = 'b_WCEN_SLEEPDUR'), get_parcel_var_model_df, 
                                    x = reliability_dt_labels, labels = ho_sub_labels)
  
  model_fn_df_for_plots <- rbindlist(model_fn_list_for_plots)
  write.table(model_fn_df_for_plots[, mod_fn], file = 'rsync_model_list.txt', 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
  
  make_multi_plot(model_fn_list_for_plots[[1]], xlab = 'Month', ylim = c(-1.75, 1.75),
                  breaks = c(-4, -1, 2, 5), labels = c(1, 4, 7, 10))
  
  mplots <- mapply(make_multi_plot, model_fn_list_for_plots,
                   file = sprintf('ex-scatters_%s.svg', 
                                  stringr::str_to_lower(
                                    stringr::str_replace_all(names(model_fn_list_for_plots),
                                                         ' ', '-'))),
                   xlab = names(model_fn_list_for_plots),
                   breaks = c(list(c(-4, -1, 2, 5)), rep(NA, 5)), 
                   labels = c(list(c(1, 4, 7, 10)), rep(NA, 5)))
  
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
