library(data.table)
library(ggplot2)
library(ggseg)
library(ggsegSchaefer)

PLOT <- FALSE

parse_filenames <- function(x, regex, colnames){
  require(stringi)
  require(data.table)
  
  parsed_dt <- as.data.table(transpose(stri_match_all_regex(x, regex)))
  setnames(parsed_dt, c('file', colnames))
  return(parsed_dt)
}
read_rds_file <- function(x){
  require(brms)
  fits <- lapply(x, readRDS)
  cfit <- brms::combine_models(mlist = fits)
  return(cfit)
}
read_time_fac <- function(x){
  #x <- filename_list[[1]][[1]]
  require(brms)
  fit <- read_rds_file(x$file)
  newdata <- data.table(sess = sprintf('%02d', 1:10))
  sess_posterior_epred_summary <- as.data.table(
    posterior_summary(posterior_epred(fit, newdata = newdata, re_formula = ~ (1 | sess)),
                      robust = TRUE, probs = c(.025, .975)))
  sess_posterior_epred_summary <- cbind(newdata, sess_posterior_epred_summary)
  sess_posterior_epred_summary[, sess := sprintf('epred_%s', sess)]
  setnames(sess_posterior_epred_summary, 'sess', 'variable')
  sess_posterior_summary <- as.data.table(
    posterior_summary(fit, variable = '(r_sess.*|b_Intercept)', regex = TRUE, robust = TRUE, probs = c(.025, .975)),
    keep.rownames = 'variable')
  
  time_fac_sess_data <- cbind(x, rbind(sess_posterior_summary, sess_posterior_epred_summary))
  return(time_fac_sess_data)
}

parallel_read_time_fac_files <- function(x, split_cols, cpus_total){
  filename_list <- split(
    split(x,
          by = split_cols,
          drop = TRUE), 
    1:cpus_total)
  
  cl <- parallel::makeCluster(cpus_total)
  
  parallel::clusterExport(cl, c('read_time_fac', 'read_rds_file'))
  nada <- parallel::clusterEvalQ(cl, {library(data.table); library(brms); setDTthreads(1)})
  system.time({data_out_list <- parallel::parLapply(cl = cl, filename_list, function(some_fns){
    lapply(some_fns, read_time_fac)
  })})
  
  parallel::stopCluster(cl)
  data_out <- rbindlist(unlist(data_out_list, recursive = FALSE))
  return(data_out)
}
plot_random_effects <- function(x, xvar = 'session', yvar = 'Estimate', ncol = 20, title = 'Random effects', file = NULL, xlab = 'Session', scales = 'free_y', width = 12, height = 12, units = 'in', dpi = 72, label = NULL){
  require(ggplot2)
  
  if(!is.null(label)){
    lookup_ <- unique(x[, .SD, .SDcols = c('roi', label)])
    lookup <- as.character(lookup_[[2]])
    names(lookup) <- lookup_$roi
    plot_labeller <- labeller(roi = lookup)
    plot_theme <- theme(strip.text.x = element_text(), panel.grid = element_blank(), axis.text = element_blank(),
                        panel.spacing = unit(.025, 'in')) 
  } else {
    plot_labeller <- labeller()
    plot_theme <- theme(strip.text.x = element_blank(), panel.grid = element_blank(), axis.text = element_blank(),
                        panel.spacing = unit(.025, 'in')) 
  }
  
  
  plot <- ggplot(x, aes(x = !! sym(xvar), y = !! sym(yvar))) + 
    geom_hline(yintercept = 0, linewidth = .5, color = 'gray') + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = .2) + 
    geom_line(linewidth = .25) + 
    facet_wrap( ~ roi, ncol = ncol, scales = scales, labeller = plot_labeller) + 
    labs(x = xlab, title = title) + 
    theme_minimal() + plot_theme
   
  #plot
  if (! is.null(file)){
    ggsave(filename = file, plot = plot, width = width, height = height, units = units, dpi = dpi)
  }
  return(plot)
}


time_fac_data_fn <- 'time_fac_data.rds'
if(!file.exists(time_fac_data_fn)){
  fit_dir <- 'fits'
  fit_files <- parse_filenames(dir(fit_dir, pattern = '.*time.*fac.*rds', full.names = TRUE),
                               regex = 'fits/alpha-fearGT(calm|scramble)-ROI([0-9]{1,3}).*rds',
                               colnames = c('contrast', 'roi'))
  
  split_cols <- c('contrast', 'roi')
  cpus_total <- 4
  time_fac_data <- parallel_read_time_fac_files(fit_files, split_cols = split_cols, cpus_total = cpus_total)
  saveRDS(time_fac_data, file = time_fac_data_fn)
} else {
  time_fac_data <- readRDS(file = time_fac_data_fn)
}

if(PLOT){
  HO_labels <- fread('HO_sub.csv')
  time_fac_data[, c('variable', 'session') := 
                  transpose(
                    stringi::stri_match_all_regex(variable, '(b_Intercept|r_sess|epred)(?:\\[|_)*(\\d{2})*.*'))[2:3]]
  time_fac_data[, session := as.numeric(session)]
  time_fac_data[, index := ifelse(as.numeric(roi) > 400, as.numeric(roi) - 400, NA)]
  time_fac_data <- merge(time_fac_data, HO_labels, all.x = TRUE, all.y = FALSE, by = 'index')
  
  cortex =1:400
  subcortex = 400 + c(4:7, 9:11, 15:21)
  
  plot <- plot_random_effects(time_fac_data[variable == 'r_sess' & roi %in% sprintf('%03d', cortex) & contrast == 'calm'], 
                              title = 'Cortex Calm',
                              file = 'calm-cortex-r_session.png',
                              width = 20, height = 10, dpi = 150)
  plot <- plot_random_effects(time_fac_data[variable == 'r_sess' & roi %in% sprintf('%03d', cortex) & contrast == 'scramble'], 
                              title = 'Cortex Scramble',
                              file = 'scramble-cortex-r_session.png',
                              width = 20, height = 10, dpi = 150)
  plot <- plot_random_effects(time_fac_data[variable == 'r_sess' & roi %in% sprintf('%03d', subcortex) & contrast == 'calm'], 
                              title = 'Subcortex Calm',
                              file = 'calm-subcort-r_session.png',
                              label = 'label',
                              width = 10, height = 5, dpi = 150, ncol = 7)
  plot <- plot_random_effects(time_fac_data[variable == 'r_sess' & roi %in% sprintf('%03d', subcortex) & contrast == 'scramble'], 
                              title = 'Subcortex Scramble',
                              file = 'scramble-subcort-r_session.png',
                              label = 'label',
                              width = 10, height = 5, dpi = 150, ncol = 7)
  plot <- plot_random_effects(time_fac_data[variable == 'epred' & roi %in% sprintf('%03d', cortex) & contrast == 'calm'], 
                              title = 'Cortex Calm',
                              file = 'calm-cortex-epred.png',
                              width = 20, height = 10, dpi = 150, scales = 'fixed')
  plot <- plot_random_effects(time_fac_data[variable == 'epred' & roi %in% sprintf('%03d', cortex) & contrast == 'scramble'], 
                              title = 'Cortex Scramble',
                              file = 'scramble-cortex-epred.png',
                              width = 20, height = 10, dpi = 150, scales = 'fixed')
  plot <- plot_random_effects(time_fac_data[variable == 'epred' & roi %in% sprintf('%03d', subcortex) & contrast == 'calm'], 
                              title = 'Subcortex Calm',
                              file = 'calm-subcort-epred.png',
                              label = 'label',
                              width = 10, height = 5, dpi = 150, ncol = 7, scales = 'fixed')
  plot <- plot_random_effects(time_fac_data[variable == 'epred' & roi %in% sprintf('%03d', subcortex) & contrast == 'scramble'], 
                              title = 'Subcortex Scramble',
                              file = 'scramble-subcort-epred.png',
                              label = 'label',
                              width = 10, height = 5, dpi = 150, ncol = 7, scales = 'fixed')
  
  schaefer_parcel_order <- fread(file = 'Schaefer2018_400Parcels_7Networks_order.txt',
                                 header = FALSE, 
                                 col.names = c('roi', 'label'),
                                 select = 1:2)
  time_fac_parcel_data <- time_fac_data[variable == 'epred' & session %in% 1:10]
  time_fac_parcel_data[, c('index', 'label') := c(NULL, NULL)]
  time_fac_parcel_data[, roi := as.numeric(roi)]
  time_fac_parcel_data <- time_fac_parcel_data[schaefer_parcel_order, on = 'roi']
  time_fac_parcel_data[, roi := NULL]
  setnames(time_fac_parcel_data, 'label', 'region')
  
  time_fac_parcel_data <- rbind(time_fac_parcel_data, 
                                as.data.table(expand.grid(session = 1:10, 
                                                          contrast = c('calm', 'scramble'),
                                                          region = NA)), fill = TRUE)
  time_fac_parcel_data[, sig_Estimate := ifelse(sign(Q2.5) == sign(Q97.5), Estimate, NA)]
  gradient_values <- (\(x) (c(min(x), -.1, 0, .1, max(x)) - min(x)) / (max(x) - min(x)))(na.omit(time_fac_parcel_data$Estimate))
  cortical_expectations_over_time_plot <- ggplot(time_fac_parcel_data[session %in% 1:10]) + 
    geom_brain(atlas = schaefer7_400, aes(fill = Estimate)) + 
    facet_grid(session ~ contrast) + 
    scale_fill_gradientn(colors = c( 'blue', 'lightblue', 'white', 'yellow', 'red'),
                         breaks = c(-.2, -.1, 0, .1, .2),
                         values = gradient_values,
                         na.value = 'white', name = 'Expected\nContrast\nValue') + 
    theme_brain() + 
    theme(axis.text = element_blank(),
          text =element_text(family = 'sans'),
          plot.background = element_rect(fill = 'white'))
  # print(cortical_expectations_over_time_plot)
  ggsave(filename = 'cortical_expectations_over_time_plot.png', 
         plot = cortical_expectations_over_time_plot,
         width = 12, height = 10, units = 'in', dpi = 150)
  
  sig_gradient_values <- (\(x) (c(min(x), -.1, 0, .1, max(x)) - min(x)) / (max(x) - min(x)))(na.omit(time_fac_parcel_data$sig_Estimate))
  sig_cortical_expectations_over_time_plot <- ggplot(time_fac_parcel_data[session %in% 1:10]) + 
    geom_brain(atlas = schaefer7_400, aes(fill = sig_Estimate)) + 
    facet_grid(session ~ contrast) + 
    scale_fill_gradientn(colors = c( 'blue', 'lightblue', 'white', 'yellow', 'red'),
                         breaks = c(-.2, -.1, 0, .1, .2),
                         values = sig_gradient_values,
                         na.value = 'white', name = 'Expected\nContrast\nValue') + 
    theme_brain() + 
    theme(axis.text = element_blank(),
          text =element_text(family = 'sans'),
          plot.background = element_rect(fill = 'white'))
  # print(sig_cortical_expectations_over_time_plot)
  ggsave(filename = 'sig_cortical_expectations_over_time_plot.png', 
         plot = sig_cortical_expectations_over_time_plot,
         width = 12, height = 10, units = 'in', dpi = 150)
}
