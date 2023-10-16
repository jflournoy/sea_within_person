#This chunk is everything we need to start
{
source('collect_results.R')

#Necessary to get path to workbench in environment
PATH_TO_WB_COMMAND = "/n/helmod/apps/centos7/Core/connectome-workbench/1.3.2-fasrc01/bin_rh_linux64"
Sys.setenv(PATH=sprintf("%s:%s", Sys.getenv('PATH'), PATH_TO_WB_COMMAND))

library(ggplot2)
library(ggrepel)
library(patchwork)
library(ggsegSchaefer)
library(ggseg)


#Functions----

plot_preds <- function(x, variable, outfile, outcome = NULL, contrast_name = 'calm'){
  #x = reliability_dt_labels; variable = 'bsp_meWCEN_BRAINWB_SE'; outfile = 'GAD7-ns.svg'; outcome = 'GAD7'
  if(!is.null(outcome)){
    x <- x[rn == variable & contrast == contrast_name &
             ns == '-ns' & alt == '' &
             pred == sprintf('-brainpred_%s', outcome) &
             ! grepl('Cerebral Cortex', label)]
  } else {
    x <- x[rn == variable & sub == '-sub' & contrast == contrast_name &
             ns == '-ns' & alt == '-alt' &
             ! grepl('Cerebral Cortex', label)]
  }
  plot <- ggplot(x, 
                 aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
    geom_rect(ymin = -19, ymax = 0, xmin = -1, xmax = 1, fill = 'lightgrey') + 
    geom_errorbar(alpha = .5, width = 0, size = .25, aes(color = credible)) +
    geom_point(size = .5, aes(color = credible)) + 
    scale_y_continuous(breaks = c(-19, 1, seq(80, 400, 80)), labels = c('Harvard/Oxford\nSubcortex', 400, seq(320, 80, -80), 'Schaefer-400 1')) +
    theme_minimal() + 
    scale_color_manual(breaks = c(TRUE, FALSE), values = c('black', 'grey')) +
    labs(x = 'Estimate (95% CI)', y = 'Parcel') + 
    theme(legend.position = 'none')
  ggsave(outfile, width = 4, height = 5, dpi = 150)
  return(plot)    
}

write_cifti_text <- function(x, variable, outfile, outcome = NULL, unthreshed = FALSE, psign = FALSE, credible = 'credible', contrast_name = 'calm'){
  #x = copy(reliability_dt_labels); variable = 'bsp_meWCEN_BRAINWB_SE'; outfile = 'GAD7-ns.txt'; outcome = 'GAD7'
  colnames = c('roi_num', 'Estimate', credible, 'p_SSEst')
  if(!is.null(outcome)){
    x <- x[rn == variable & contrast == contrast_name & 
             ns == '-ns' &  alt == ''  &
             pred == sprintf('-brainpred_%s', outcome) &
             !grepl('Cerebral Cortex', label), 
           ..colnames]
  } else {
    x <- x[rn == variable & sub == '-sub' & contrast == contrast_name &
             ns == '-ns' &  alt == '-alt' &
             !grepl('Cerebral Cortex', label), 
           ..colnames]
  }
  setnames(x, colnames, c('roi_num', 'Estimate', 'credible', 'p_SSEst'))
  setkey(x, roi_num)
  rois <- data.table(roi_num = 1:400)
  out_data <- x[rois]
  if(!unthreshed){
    out_data[! credible %in% TRUE, Estimate := 0]
  }
  if(psign){
    out_data[, Estimate := sign(Estimate) * p_SSEst]
  }
  write.table(out_data[, Estimate], file = outfile, row.names = FALSE, col.names = FALSE)
  return(TRUE)
}

get_selected_parcels <- function(x, variable, labels, contrast_name = 'calm'){
  # x <- copy(reliability_dt_labels)
  # variable = 'b_TIMECENTER'
  # labels = copy(ho_sub_labels)
  x <- x[rn == variable & sub == '-sub' & contrast == contrast_name &
           ns == '-ns' & alt == '-alt' &
           ! grepl('Cerebral Cortex', label), 
         c('roi_num', 'Estimate', 'credible', 'rn', 'roi', 'sub', 'ns', 'pred', 'alt')]
  setkey(x, roi_num)
  setkey(labels, roi_num)
  x <- labels[x][credible %in% TRUE]
  return(x)
}

get_parcel_var_model_df <- function(x, variable, labels, contrast_name = 'calm'){
  x <- get_selected_parcels(x, variable, labels)
  x[, mod_fn := paste0('alpha-fearGT', contrast_name, '-ROI', roi, sub, ns, pred, '.rds')]
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

create_cifti <- function(text_file, cifti_template, cifti_out, txt = FALSE){
  #Relies on wb_command to be in the path.
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
  cmd <- wb_cmd
  if(txt){
    return(cmd)
  } else {
    return(system(cmd))  
  }
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
  p_list <- lapply(1:dim(x)[[1]], function(i){
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

write_rez_table <- function(x, variable, tables_dir = 'tables', outcome = NULL, contrast_name = 'calm'){
  #x <- reliability_dt_labels
  #variable <- 'b_TIMECENTER'
  #x = copy(reliability_dt_labels); variable = 'bsp_meWCEN_BRAINWB_SE'; outcome = 'GAD7'; tables_dir = 'tables'
  Sch7_400_regions <- fread('Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order_info.txt', header = FALSE, col.names = 'region', sep = ',')
  Sch7_400_regions <- Sch7_400_regions[seq(1, 800, 2)]
  Sch7_400_regions[, roi_num := 1:400]
  setkey(Sch7_400_regions, roi_num)
  if(!is.null(outcome)){
    x <- x[rn == variable &
             ns == '-ns' & alt == '' & contrast == contrast_name &
             pred == sprintf('-brainpred_%s', outcome) &
             !grepl('Cerebral Cortex', label),
           c('roi_num', 'Estimate', 'Q2.5', 'Q97.5', 'credible')]
  } else {
    x <- x[rn == variable & sub == '-sub' & contrast == contrast_name &
             ns == '-ns' & alt == '-alt' &
             !grepl('Cerebral Cortex', label),
           c('roi_num', 'Estimate', 'Q2.5', 'Q97.5', 'credible')]  
  }
  
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
  this_region[, 'Estimate [95% CI]' := sprintf('% 0.3f [% 0.3f, % 0.3f]', Estimate, Q2.5, Q97.5)]
  setorder(this_region, roi_num)
  title <- paste(c(variable, outcome), collapse = '-')
  newtable <- pander::Pandoc$new(author = 'John C. Flournoy', title = title)
  newtable$add(
    pander::pandoc.table.return(
      this_region[, c('Hemisphere', 'Parcel Label', 'Estimate [95% CI]', 'Location')],
      split.table = Inf,
      justify = c('left', 'left', 'right', 'center')))
  newtable$format <- 'docx'
  if(!dir.exists(tables_dir)) dir.create(tables_dir)
  r1 <- newtable$export(title, open = FALSE)
  r2 <- file.rename(paste0(title, '.md'), file.path(tables_dir, paste0(title, '.md')))
  r3 <- file.rename(paste0(title, '.docx'), file.path(tables_dir, paste0(title, '.docx')))
  return(all(c(r2, r3)))
}

#Code----

##Parcel Icons----
parcel_icon_dir <- 'plots/parcel_icons'
if(length(dir(parcel_icon_dir)) == 0) make_parcel_icons(parcel_icon_dir, dpi = 300)

##Prep----
gsub_regex <- 'alpha-fearGT(calm|scramble)-ROI([0-9]{1,3})(-more_\\d{2})*(-sub)*(-ns)*(-(?:brain)*pred(-alt)*_\\w+)*(-cntrlreg)*(-network)*-psamp.rds'

postsamp_files <- data.table(file = dir('post_samples', pattern = 'alpha-fearGT.*rds'))
postsamp_files[, c('contrast', 'roi', 'more', 'sub', 'ns', 'pred', 'alt', 'cntrlreg', 'network') := 
                 lapply(list('\\1', '\\2', '\\3', '\\4', '\\5', '\\6', '\\7', '\\8', '\\9'), gsub, pattern = gsub_regex, x = file)]

postsamp_filenames <- split(postsamp_files, 1:cpus_total)

prop_var_post_samps_fn <- 'prop_var_post_samps.rds'
if(!file.exists(prop_var_post_samps_fn)){
  cl <- parallel::makeCluster(cpus_total)
  
  nada <- parallel::clusterEvalQ(cl, {library(data.table); setDTthreads(1)})
  system.time({post_samps_list <- parallel::parLapply(cl = cl, postsamp_filenames, function(some_fns){
    lapply(1:dim(some_fns)[[1]], function(i){
      d <- data.table(readRDS(file.path('post_samples', some_fns[i, file])))
      d <- cbind(d, some_fns[i, -c('file')])
      return(d)
    })
  })})
  
  parallel::stopCluster(cl)
  post_samps_dt <- rbindlist(unlist(post_samps_list, recursive = FALSE), fill = TRUE)
  saveRDS(post_samps_dt, prop_var_post_samps_fn)
} else {
  post_samps_dt <- readRDS(prop_var_post_samps_fn)
}

###load MA results----

source('collect_ma_results.R')
reliability_dt <- rbindlist(list(reliability_dt, ma_dt), use.names = TRUE, fill = TRUE)

###average voxels per parcel----
reliability_dt[sub == '' & ns == '-ns' & rn == 'prop_var' & 
                 pred == '' & contrast == 'calm' & network == ''][, .(nvoxel_median = quantile(nvoxel, .5),
                                                                      nvoxel_25 = quantile(nvoxel, .25),
                                                                      nvoxel_75 = quantile(nvoxel, .75))]

###Process posterior samples----
source('get_roi_data.R')
cntrlregions <- unique(load_control_regions('tsv_nosmooth')[, c('eq_ROI', 'target_ROI', 'anatomy', 'target_N')])
setnames(post_samps_dt, 
         old = c('prop_var', 'prop_var_id'),
         new = c('sess', 'id'))
post_samps_dt[, id_sess := id + sess]

post_samps_dt[, roi_num := as.numeric(roi)]
cntrlregions[, roi_num := eq_ROI]
post_samps_dt <- merge(post_samps_dt, cntrlregions[, c('roi_num', 'anatomy')], by = 'roi_num', all.x = TRUE)

post_samps_dt[, group := NA_character_]
post_samps_dt[
  pred == '', 
  group := fifelse(roi <= 400 & cntrlreg == '' & network == '', 'cortical',
                   fifelse(roi %in% c(401, 403, 412, 414) & cntrlreg == '' & network == '', 'reference',
                           fifelse(roi %in% c(404:407, 409:411, 415:421) & cntrlreg == '' & network == '', 'subcortical',
                                   fifelse(cntrlreg == '-cntrlreg' & anatomy == 'subcort', 'SMC', 
                                           fifelse(cntrlreg == '-cntrlreg' & anatomy == 'cort', 'SMC', 
                                                   fifelse(network == '-network' & roi < 100, 'network', 
                                                           fifelse(network == '-network' & roi > 100, 'SMC_network', NA_character_)))))))]


post_samps_dt_split <- split(post_samps_dt[!is.na(group)], by = c('group', 'contrast'))

brain_group_prop_var <- rbindlist(lapply(post_samps_dt_split, function(adt){
  rbindlist(lapply(c('id', 'id_sess'), function(p){
    d <- data.table(posterior_summary(adt[, get(p)], robust = TRUE))
    d[, param:= p]
    return(d)
  }))
}), id = 'group')

brain_group_prop_var_table <- brain_group_prop_var[, list(
  'Estimate [95% CI]' = sprintf('%0.3f [%0.3f, %0.3f]', Estimate, Q2.5, Q97.5)
), by = c('group', 'param')]
brain_group_prop_var_table[, param := ifelse(param == 'id', 'Between-person', 'Within-person')]
brain_group_prop_var_table[, c('contrast', 'group') :=
                             .(gsub('.*\\.(calm|scramble)', '\\1', group),
                               gsub('(.*)\\.(calm|scramble)', '\\1', group))]
brain_group_prop_var_table <- dcast(brain_group_prop_var_table, group + contrast ~ param, value.var = 'Estimate [95% CI]')
brain_group_prop_var_table[, group := factor(group, levels = c('reference', 'SMC', 'subcortical', 'cortical', 'SMC_network', 'network'))]
setorder(brain_group_prop_var_table, contrast, group)
knitr::kable(brain_group_prop_var_table[group %in% c('SMC', 'subcortical', 'cortical'), c('contrast', 'group', 'Between-person', 'Within-person')], 
             caption = 'Proportion variance by ROI type Estimates and 95% CIs', format = 'simple')

brain_group_prop_var[, c('contrast', 'group') :=
                       .(gsub('.*\\.(calm|scramble)', '\\1', group),
                         gsub('(.*)\\.(calm|scramble)', '\\1', group))]

brain_group_prop_var_dens <- rbindlist(lapply(post_samps_dt_split, function(adt){
  rbindlist(lapply(c('id', 'id_sess'), function(p){
    d <- as.data.table(density(adt[, get(p)], n = 512)[c('x','y')])
    d[, param:= p]
    return(d)
  }))
}), id = 'group')
brain_group_prop_var[, Est.Error := NULL]
brain_group_prop_var_CI <- brain_group_prop_var[, c('group', 'Q2.5', 'Q97.5', 'param', 'contrast')]
brain_group_prop_var <-  brain_group_prop_var[, c('group', 'Estimate', 'param', 'contrast')]
brain_group_prop_var_w <- dcast(brain_group_prop_var, group + contrast ~ param, value.var = 'Estimate')
setnames(brain_group_prop_var_w, c('group', 'contrast', 'id_xmax', 'sess_xmax'))
brain_group_prop_var_w[, c('id_xmin', 'sess_xmin', 'err_xmin', 'err_xmax') :=
                         list(0, id_xmax, sess_xmax, 1)]
brain_group_prop_var_w_l <- melt(brain_group_prop_var_w, id.vars = c('group', 'contrast'))
brain_group_prop_var_w_l[, c('source', 'stat', 'variable') :=
                           list(gsub('(\\w+)_(\\w+)', '\\1', variable),
                                gsub('(\\w+)_(\\w+)', '\\2', variable),
                                NULL)]
brain_group_prop_var_forplot <- dcast(brain_group_prop_var_w_l, ... ~ stat, value.var = 'value')

region_groups <- c('cortical' = 4,
                   'reference' = 1,
                   'subcortical' = 3,
                   'SMC' = 2,
                   'network' = 5,
                   'SMC_network' = 6)
brain_group_prop_var_forplot[, group_num := region_groups[group]]
brain_group_prop_var_CI[, group_num := region_groups[group]]
brain_group_prop_var_CI[, param := factor(param, levels = c('id', 'id_sess'), labels = c('id', 'sess'))]
brain_group_prop_var_dens[, c('contrast', 'group') :=
    .(gsub('.*\\.(calm|scramble)', '\\1', group),
      gsub('(.*)\\.(calm|scramble)', '\\1', group))]
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
  rn == 'prop_var' & sub == '-sub' & ns == '-ns' & pred == '' & alt == '' & cntrlreg == '' & contrast == 'calm' & network == '', ][
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

setorder(reliability_dt, rn, sub, ns, pred, alt, p_nSSEst)
reliability_dt[grepl('-brainpred_', pred) & is.na(draws), draws := 1600]
reliability_dt[, p_SSEst := fifelse(p_SSEst == 1, draws / (draws + 1), p_SSEst)]
reliability_dt[, p_nSSEst := fifelse(p_nSSEst == 0, 1 - p_SSEst, p_nSSEst)]
reliability_dt[pred != '' & pred != '-pred_all', 
               Type_S := 1 - cumprod(p_SSEst),
               by = c('rn', 'sub', 'ns', 'pred', 'alt', 'contrast', 'cntrlreg', 'network')]
reliability_dt[pred != '' & pred != '-pred_all', 
               Type_S_cumsum := cumsum(p_nSSEst),
               by = c('rn', 'sub', 'ns', 'pred', 'alt', 'contrast', 'cntrlreg', 'network')]

reliability_dt[, credible := Type_S < .05]
reliability_dt[, credible_cumsum := Type_S_cumsum < .05]
reliability_dt[, credible_99 := (sign(Q0.006) * sign(Q99.994)) == 1]
reliability_dt[, credible_95 := (sign(Q2.5) * sign(Q97.5)) == 1]

reliability_dt_labels <- ho_sub_labels[reliability_dt, on = 'roi_num']
reliability_dt_labels[, label_text := fifelse((Estimate < .22 & rn == 'prop_var') |
                                                (Estimate < .75 & rn == 'alpha'), label, NA_character_)]

reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & pred == '' & roi > 400, c('label', 'label_text', 'Estimate')]

reliability_dt_labels[, roi_num_reverse := -1*(roi_num - 401)]

reliability_dt_labels_prop_w <- dcast(reliability_dt_labels[rn %in% c('prop_var_id', 'idsess_right') & 
                                                              ns == '-ns' & sub == '-sub' & pred == '' &
                                                              !roi_num %in% c(402, 413), #exclude cerebral cortex from HO
                                                            c('rn', 'roi_num', 'Estimate', 'roi_num_reverse', 'contrast', 'cntrlreg', 'network')],
                                      ... ~ rn, value.var = 'Estimate')

setnames(reliability_dt_labels_prop_w, 
         old = c('idsess_right', 'prop_var_id'), 
         new = c('sess_xmax', 'id_xmax'))
reliability_dt_labels_prop_w[, c('id_xmin', 'sess_xmin', 'err_xmin', 'err_xmax') :=
                               list(0, id_xmax, sess_xmax, 1)]
reliability_dt_labels_prop_l <- melt(reliability_dt_labels_prop_w, id.vars = c('roi_num', 'roi_num_reverse', 'contrast', 'cntrlreg', 'network'))
reliability_dt_labels_prop_l[, c('source', 'stat', 'variable') :=
                               list(gsub('(\\w+)_(\\w+)', '\\1', variable),
                                    gsub('(\\w+)_(\\w+)', '\\2', variable),
                                    NULL)]
reliability_dt_labels_prop_l_w <- dcast(reliability_dt_labels_prop_l, ... ~ stat, value.var = 'value')
}
##Amygdala----

x_labels <- c("b_TIMECENTER" = 'Time',
              "b_WCEN_CHRONICSEV" = 'Chronic Stress',
              "b_WCEN_EPISODICTOT" = 'Episodic Stress',
              "b_WCEN_PANASNEG" = 'Negative Affect',
              "b_WCEN_SLEEPDUR" = 'Sleep Duration')

amygdala_plot <- ggplot(reliability_dt_labels[roi %in% c(410, 420) & alt == '-alt' & ns == '-ns' & sub == '-sub' & contrast == 'calm'],
                        aes(x = rn, y = Estimate)) + 
  geom_hline(yintercept = 0) + 
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0, alpha = .8) + 
  geom_point() +
  scale_x_discrete(breaks = names(x_labels), labels = x_labels) + 
  facet_grid(label ~ ., labeller = labeller(label = label_value, contrast = c(calm = 'Neutral', scramble = 'Scramble'))) + 
  theme_minimal() + 
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 360-45, hjust = .5)) + 
  labs(x = '', y = 'Estimate and 95% credible interval')
print(amygdala_plot)
ggsave('amygdala_plot.svg', plot = amygdala_plot, width = 8, height = 5, dpi = 150)
ggsave('amygdala_plot.png', plot = amygdala_plot, width = 8, height = 5, dpi = 150)

# Variance proportions

##Plot Reliability----
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
  theme_minimal() + 
  facet_wrap( ~ contrast + cntrlreg + network)
print(proportion_by_ROI_plot)
ggsave('prop_id_sess_err_variance_by_ROI.svg', plot = proportion_by_ROI_plot, width = 8, height = 5, dpi = 150)

reliability_dt_labels_prop_l_w[, group := fifelse(roi_num <= 400 & cntrlreg == '' & network == '', 'Observed_cortical',
                                                  fifelse(roi_num %in% c(401, 403, 412, 414) & cntrlreg == '' & network == '', 'reference',
                                                          fifelse(roi_num %in% c(404:407, 409:411, 415:421) & cntrlreg == '' & network == '', 'Observed_subcortical',
                                                                  fifelse(cntrlreg == '-cntrlreg' & roi_num %in% 1:14, 'SMC', 
                                                                          fifelse(cntrlreg == '-cntrlreg' & roi_num %in% 15:24, 'SMC', 
                                                                                  fifelse(network == '-network' & roi_num < 100, 'Observed_network', 
                                                                                          fifelse(network == '-network' & roi_num > 100, 'SMC_network', NA_character_)))))))]
reliability_dt_labels_prop_l_w[, c('group', 'anatomy') := tstrsplit(group, '_')]
proportion_by_ROI_plot_redo <- ggplot(reliability_dt_labels_prop_l_w[group %in% c('Observed', 'SMC') & anatomy %in% c(NA, 'cortical', 'subcortical')], aes(x = NULL, y = factor(roi_num))) + 
  geom_ribbon(aes(xmin = xmin, xmax = xmax, fill = source, group = source), 
              alpha = 1) + 
  coord_cartesian(xlim = c(0,1)) +
  labs(x = 'Cumulative Proportion', y = 'Parcel / Region') + 
  scale_fill_manual(breaks = c('err', 'id', 'sess')[source_order], 
                    values = c('#FFF8F3', '#EDC12C', '#71A1D4')[source_order],
                    labels = c('Error', 'Participant means', 'Within-session means')[source_order],
                    name = 'Source of variance') + 
  scale_x_continuous(breaks = c(0, .5, 1)) + 
  theme_minimal() + 
  facet_wrap(~contrast + anatomy + group, scales = 'free', ncol = 3) + 
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
print(proportion_by_ROI_plot_redo)
ggsave('prop_id_sess_err_variance_by_ROI_redo.svg', plot = proportion_by_ROI_plot_redo, width = 8, height = 5, dpi = 150)

proportion_by_ROI_plot_redo_calm <- ggplot(reliability_dt_labels_prop_l_w[group %in% c('Observed', 'SMC') & contrast == 'calm' & anatomy %in% c(NA, 'cortical', 'subcortical')], aes(x = NULL, y = factor(roi_num))) + 
  geom_ribbon(aes(xmin = xmin, xmax = xmax, fill = source, group = source), 
              alpha = 1) + 
  coord_cartesian(xlim = c(0,1)) +
  labs(x = 'Cumulative Proportion', y = 'Parcel / Region') + 
  scale_fill_manual(breaks = c('err', 'id', 'sess')[source_order], 
                    values = c('#FFF8F3', '#EDC12C', '#71A1D4')[source_order],
                    labels = c('Error', 'Participant means', 'Within-session means')[source_order],
                    name = 'Source of variance') + 
  scale_x_continuous(breaks = c(0, .5, 1)) + 
  theme_minimal() + 
  facet_wrap(~group + anatomy, scales = 'free', ncol = 3) + 
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
print(proportion_by_ROI_plot_redo_calm)
ggsave('prop_id_sess_err_variance_by_ROI_redo_calm.svg', plot = proportion_by_ROI_plot_redo_calm, width = 8, height = 5, dpi = 150)


# Variance proportions by region-groups
region_group_labels <- c('cortical' = 'Cortical',
                         'reference' = 'Reference',
                         'subcortical' = 'Subcortical',
                         'SMC' = 'SMC',
                         'network' = 'Network',
                         'SMC_network' = 'SMC Network')
brain_group_prop_var_dens[, y_pct := y/max(y), by = c('group', 'param', 'contrast')]
brain_group_prop_var_dens[, x := fifelse(x < 0, NA_real_, x)]

group_nums <- c(2, 4, 3)
brain_group_prop_var_forplot_BW <- brain_group_prop_var_forplot[!is.na(group) & group_num %in% group_nums & contrast == 'calm']
brain_group_prop_var_CI_BW <- brain_group_prop_var_CI[!is.na(group) & group_num %in% group_nums & contrast == 'calm']
brain_group_prop_var_forplot_BW <- brain_group_prop_var_forplot_BW[, group := factor(group, levels = c('SMC', 'subcortical', 'cortical'), labels = c('SMC', 'Subcortical', 'Cortical'))]
brain_group_prop_var_CI_BW <- brain_group_prop_var_CI_BW[, group := factor(group, levels = c('SMC', 'subcortical', 'cortical'), labels = c('SMC', 'Subcortical', 'Cortical'))]

grouped_prop_var <- ggplot(brain_group_prop_var_forplot_BW, aes(x = NULL, y = group_num)) + 
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = group_num - .5, ymax = group_num + .5, group = source, fill = source),
            alpha = .85) + 
  # geom_rect(data = brain_group_prop_var_CI,
  #           aes(xmin = Q2.5, xmax = Q97.5, ymin = group_num - .45, ymax = group_num + .45, group = param, fill = param),
  #           alpha = .5) +
  geom_segment(data = brain_group_prop_var_CI_BW,
               aes(x = Q2.5, xend = Q97.5, y = group_num, yend = group_num, color = param),
               size = 1) +
  geom_vline(xintercept = c(0,1), color = '#666666') + 
  # geom_segment(data = brain_group_prop_var_dens,
  #              aes(x = x, xend = x, y = group_num + .5, yend = group_num - .5, alpha = y_pct,
  #                  group = param, color = param)) +
  scale_fill_manual(breaks = c('err', 'id', 'sess')[source_order],
                    values = c('#FFFFFF', '#666666', '#AAAAAA')[source_order],
                    labels = c('Error', 'Participant means', 'Within-session means')[source_order],
                    name = 'Source of variance', aesthetics = c('fill', 'color')) +
  scale_y_continuous(breaks = 2:4, labels = rev(region_group_labels[names(region_groups[region_groups %in% group_nums])])) + 
  scale_alpha_continuous(range = c(0, 1)) + 
  coord_cartesian(xlim = c(0, 1)) + 
  theme_minimal() +
  theme(axis.line = element_blank(), panel.grid = element_blank(), legend.position = 'bottom') + 
  labs(x = 'Cumulative proportion variance (95% CI)', y = '', title = 'Mean Proportion Variance Across Anatomy Types and Size-Matched Controls')
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

##ICC----

reliability_dt_labels[rn == 'icc_ma_id' & sub == '' &
                        ns == '-ns' & pred == '' &
                        network == '' & contrast == 'calm' &
                        roi_num %in% c(404:407, 409:411, 415:421), roi_num_reverse := -1:-.N]

reliability_dt_labels[rn == 'icc_ma_id' & sub == '' &
                        ns == '-ns' & pred == '' &
                        network == '' & contrast == 'calm' &
                        roi_num %in% c(1:400) & cntrlreg == '-cntrlreg', roi_num_reverse := -15:(-14 - .N)]

between_icc_plot <- ggplot(reliability_dt_labels[rn == 'icc_ma_id' & sub == '' &
                                                   ns == '-ns' & pred == '' &
                                                   network == '' & contrast == 'calm' &
                                                   roi_num %in% c(1:400, 404:407, 409:411, 415:421)], 
                           aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
  annotate('rect', ymin = -14.5, ymax = -.5, xmin = -.025, xmax = 1, color = '#666666', linewidth = .25, fill = 'white') + 
  annotate('rect', ymin = -38.5, ymax = -14.5, xmin = -.025, xmax = 1, color = '#666666', linewidth = .25, fill = 'white') + 
  geom_errorbar(alpha = .5, width = 0, linewidth = .10) +
  geom_point(size = 1.5, aes(fill = Estimate), shape = 21) + 
  scale_y_continuous(breaks = c(-26, -7, 240), labels = c('Size-Matched\nControl\n\n\n', '\n\n\nHarvard/Oxford\nSubcortex', 'Schaefer\n400')) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        axis.text.y = element_text(angle = 90, vjust = 0, hjust = 0)) +
  coord_cartesian(x = c(-.025, 1)) + 
  scale_fill_gradientn(colors = c('#ff0000', '#ff0000', '#ffff00', '#ffff00'), 
                       values = c(0, 0.10, 0.80, 1),
                       limits = c(0, 1),
                       rescaler = function(x, from = y) x) + 
  labs(x = 'Between-Person Reliability (Proportion Variance Due to Participant)', y = '')
print(between_icc_plot)
ggsave('icc_ma_id-ns.svg', width = 8, height = 5, dpi = 150)
ggsave('icc_ma_id-ns.png', width = 8, height = 5, dpi = 150)

reliability_dt_labels[rn == 'prop_var_id' & sub == '-sub' &
                        ns == '-ns' & pred == '' &
                        network == '' & contrast == 'calm' &
                        roi_num %in% c(404:407, 409:411, 415:421), roi_num_reverse := -1:-.N]

reliability_dt_labels[rn == 'prop_var_id' & sub == '-sub' &
                        ns == '-ns' & pred == '' &
                        network == '' & contrast == 'calm' &
                        roi_num %in% c(1:400) & cntrlreg == '-cntrlreg', roi_num_reverse := -15:(-14 - .N)]

between_icc_plot_all_ses <- ggplot(reliability_dt_labels[rn == 'prop_var_id' & sub == '-sub' &
                                                   ns == '-ns' & pred == '' &
                                                   network == '' & contrast == 'calm' &
                                                   roi_num %in% c(1:400, 404:407, 409:411, 415:421)], 
                           aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
  annotate('rect', ymin = -14.5, ymax = -.5, xmin = -.025, xmax = 1, color = '#666666', linewidth = .25, fill = 'white') + 
  annotate('rect', ymin = -38.5, ymax = -14.5, xmin = -.025, xmax = 1, color = '#666666', linewidth = .25, fill = 'white') + 
  geom_errorbar(alpha = .5, width = 0, linewidth = .10) +
  geom_point(size = 1.5, aes(fill = Estimate), shape = 21) + 
  scale_y_continuous(breaks = c(-26, -7, 240), labels = c('Size-Matched\nControl\n\n\n', '\n\n\nHarvard/Oxford\nSubcortex', 'Schaefer\n400')) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        axis.text.y = element_text(angle = 90, vjust = 0, hjust = 0)) +
  coord_cartesian(x = c(-.025, 1)) + 
  scale_fill_gradientn(colors = c('#ff0000', '#ff0000', '#ffff00', '#ffff00'), 
                       values = c(0, 0.10, 0.80, 1),
                       limits = c(0, 1),
                       rescaler = function(x, from = y) x) + 
  labs(x = 'Between-Person Reliability (Proportion Variance Due to Participant)', y = '')
print(between_icc_plot_all_ses)
ggsave('prop_var_id-ns.svg', width = 8, height = 5, dpi = 150)
ggsave('prop_var_id-ns.png', width = 8, height = 5, dpi = 150)

ns_icc_range <- range(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & ! grepl('Cerebral Cortex', label), Estimate])

reliability_dt_labels[rn == 'prop_var' & sub == '-sub' &
                        ns == '-ns' & pred == '' &
                        network == '' & contrast == 'calm' &
                        roi_num %in% c(404:407, 409:411, 415:421), roi_num_reverse := -1:-.N]

reliability_dt_labels[rn == 'prop_var' & sub == '-sub' &
                        ns == '-ns' & pred == '' &
                        network == '' & contrast == 'calm' &
                        roi_num %in% c(1:400) & cntrlreg == '-cntrlreg', roi_num_reverse := -15:(-14 - .N)]
#      Estimate Est.Error Q2.5 Q97.5
# [1,]     0.44      0.09 0.27  0.61
# [2,]     0.36      0.10 0.17  0.58
# [3,]     0.34      0.11 0.15  0.60

#      Estimate Est.Error Q2.5 Q97.5
# [1,]     0.44      0.00 0.43  0.45
# [2,]     0.36      0.03 0.31  0.42
# [3,]     0.35      0.03 0.30  0.41

within_icc_plot <- ggplot(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' &
                                                  ns == '-ns' & pred == '' &
                                                  network == '' & contrast == 'calm' &
                                                  roi_num %in% c(1:400, 404:407, 409:411, 415:421)], 
                          aes(x = Estimate, xmin = Q2.5, xmax = Q97.5, y = roi_num_reverse)) +
  annotate('rect', ymin = -14.5, ymax = -.5, xmin = -.025, xmax = 1, color = '#666666', linewidth = .25, fill = 'white') + 
  annotate('rect', ymin = -38.5, ymax = -14.5, xmin = -.025, xmax = 1, color = '#666666', linewidth = .25, fill = 'white') +  
  annotate('rect', ymin = -38.5, ymax = -14.5, xmin = .30, xmax = .41, fill = 'black', alpha = .2) + 
  annotate('rect', ymin = -14.5, ymax = -.5, xmin = .31, xmax = .42, fill = 'black', alpha = .2) + 
  annotate('rect', ymin = -.5, ymax = 400, xmin = .43, xmax = .45, fill = 'black', alpha = .2) + 
  geom_errorbar(alpha = .5, width = 0, linewidth = .10) +
  geom_point(size = 1.5, aes(fill = Estimate), shape = 21) + 
  scale_y_continuous(breaks = c(-26, -7, 240), labels = c('Size-Matched\nControl\n\n\n', '\n\n\nHarvard/Oxford\nSubcortex', 'Schaefer\n400')) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = 'none',
        axis.text.y = element_text(angle = 90, vjust = 0, hjust = 0)) +
  coord_cartesian(x = c(-.025, 1)) +
  scale_fill_gradientn(colors = c('#ff0000', '#ff0000', '#ffff00', '#ffff00'), 
                        values = c(0, 0.10, 0.80, 1),
                        limits = c(0, 1),
                        rescaler = function(x, from = y) x) + 
  labs(x = 'Within-Person Reliability (Proportion Variance Due to Participant-Session)', y = '')
print(within_icc_plot)
ggsave('prop_var_id-session-ns.svg', width = 8, height = 5, dpi = 150)



figure_3 <- between_icc_plot + theme(plot.tag.position = c(0,1)) + 
  within_icc_plot + theme(plot.tag.position = c(0,1)) + 
  plot_layout(ncol = 1) + 
  plot_annotation(tag_levels = 'A') 
ggsave('figure_3.svg', figure_3, width = 6, height = 7.52, dpi = 150)

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
  labs(x = 'Proportion variance due to participant-session', y = '')
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

tau_pi_delta_prop_plot <- ggplot(reliability_dt_labels[rn == 'tau_pi_delta_prop' & sub == '-sub' &
                                                         ns == '-ns' & pred == '-pred_all' & alt =='' & 
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


##Predictors----
#preds versus pred-alts
ggplot(dcast(reliability_dt_labels[grepl('-pred', pred) & pred != '-pred_all' & contrast == 'calm', c('rn', 'roi_num', 'Estimate', 'alt')], ... ~ alt, value.var = 'Estimate'),
       aes(x = V1, y = `-alt`)) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  geom_point(alpha = .5, aes(color = rn)) + 
  theme_minimal() + 
  labs(title = 'Correlation between estimates\nfor correct and incorrect specifications',
       x = 'RE by `id` (un-nested, incorrect)', y = 'RE by `id/sess` (nested, correct)', color = 'Predictor')

previously_credible <- reliability_dt_labels[grepl('-pred', pred) & contrast == 'calm' & pred != '-pred_all' & credible_99 %in% c(TRUE) & alt == '', c('rn', 'roi_num')]

ggplot(dcast(reliability_dt_labels[grepl('-pred', pred) & contrast == 'calm' & pred != '-pred_all', c('rn', 'roi_num', 'Estimate', 'alt')], ... ~ alt, value.var = 'Estimate')[
  previously_credible, on = c('rn', 'roi_num')],
  aes(x = V1, y = `-alt`)) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  geom_point(alpha = .5, aes(color = rn)) + 
  theme_minimal() + 
  labs(title = 'Correlation between estimates\nfor correct and incorrect specifications',
       subtitle = 'From previously significant parcels',
       x = 'RE by `id` (un-nested, incorrect)', y = 'RE by `id/sess` (nested, correct)', color = 'Predictor')

ggplot(dcast(reliability_dt_labels[grepl('-pred', pred) & contrast == 'calm' & pred != '-pred_all', c('rn', 'roi_num', 'Est.Error', 'alt')], ... ~ alt, value.var = 'Est.Error'),
       aes(x = V1, y = `-alt`)) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  geom_point(alpha = .5, aes(color = rn)) + 
  theme_minimal() + 
  labs(title = 'Increased uncertainty in correctly specified models',
       x = 'RE by `id` (un-nested, incorrect)', y = 'RE by `id/sess` (nested, correct)', color = 'Predictor')

ggplot(dcast(reliability_dt_labels[grepl('-pred', pred) & contrast == 'calm' & pred != '-pred_all', c('rn', 'roi_num', 'Est.Error', 'alt')], ... ~ alt, value.var = 'Est.Error')[
  previously_credible, on = c('rn', 'roi_num')],
  aes(x = V1, y = `-alt`)) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  geom_point(alpha = .5, aes(color = rn)) + 
  theme_minimal() + 
  labs(title = 'Increased uncertainty in correctly specified models',
       subtitle = 'From previously significant parcels',
       x = 'RE by `id` (un-nested, incorrect)', y = 'RE by `id/sess` (nested, correct)', color = 'Predictor')

reliability_dt_labels[grepl('-pred', pred) & contrast == 'calm' & pred != '-pred_all', c('rn', 'Estimate', 'alt', 'pred')]

plot_preds(x = reliability_dt_labels, variable = 'b_TIMECENTER', outfile = 'time-ns-sub.svg')
plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', outfile = 'episodic_stress-ns-sub.svg')
plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', outfile = 'chronic_stress-ns-sub.svg')
plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', outfile = 'negative_affect-ns-sub.svg')
plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR_NEG', outfile = 'sleep-ns-sub.svg')

plot_preds(x = reliability_dt_labels, variable = 'b_TIMECENTER', outfile = 'time-scramble-ns-sub.svg', contrast_name = 'scramble')
plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', outfile = 'episodic_stress-scramble-ns-sub.svg', contrast_name = 'scramble')
plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', outfile = 'chronic_stress-scramble-ns-sub.svg', contrast_name = 'scramble')
plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', outfile = 'negative_affect-scramble-ns-sub.svg', contrast_name = 'scramble')
plot_preds(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR_NEG', outfile = 'sleep-scramble-ns-sub.svg', contrast_name = 'scramble')

plot_preds(x = reliability_dt_labels, variable = 'bsp_meWCEN_BRAINWB_SE', outfile = 'GAD7-ns.svg', outcome = 'GAD7')
plot_preds(x = reliability_dt_labels, variable = 'bsp_meWCEN_BRAINWB_SE', outfile = 'PHQ9-ns.svg', outcome = 'PHQ9')
# theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

write_cifti_text(x = reliability_dt_labels, variable = 'b_TIMECENTER', outfile = 'time-ns-sub.txt')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', outfile = 'episodic_stress-ns-sub.txt')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', outfile = 'chronic_stress-ns-sub.txt')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', outfile = 'negative_affect-ns-sub.txt')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR_NEG', outfile = 'sleep-ns-sub.txt')

write_cifti_text(x = reliability_dt_labels, variable = 'b_TIMECENTER', outfile = 'time-scramble-ns-sub.txt', contrast_name = 'scramble')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', outfile = 'episodic_stress-scramble-ns-sub.txt', contrast_name = 'scramble')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', outfile = 'chronic_stress-scramble-ns-sub.txt', contrast_name = 'scramble')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', outfile = 'negative_affect-scramble-ns-sub.txt', contrast_name = 'scramble')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR_NEG', outfile = 'sleep-scramble-ns-sub.txt', contrast_name = 'scramble')

write_cifti_text(x = reliability_dt_labels, variable = 'bsp_meWCEN_BRAINWB_SE', outfile = 'GAD7-ns.txt', outcome = 'GAD7')
write_cifti_text(x = reliability_dt_labels, variable = 'bsp_meWCEN_BRAINWB_SE', outfile = 'PHQ9-ns.txt', outcome = 'PHQ9')

create_cifti(text_file = 'time-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'time-ns-sub_Sch400.pscalar.nii')
create_cifti(text_file = 'episodic_stress-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'episodic_stress-ns-sub_Sch400.pscalar.nii')
create_cifti(text_file = 'chronic_stress-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'chronic_stress-ns-sub_Sch400.pscalar.nii')
create_cifti(text_file = 'negative_affect-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'negative_affect-ns-sub_Sch400.pscalar.nii')
create_cifti(text_file = 'sleep-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'sleep-ns-sub_Sch400.pscalar.nii')

create_cifti(text_file = 'time-scramble-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'time-scramble-ns-sub_Sch400.pscalar.nii')
create_cifti(text_file = 'episodic_stress-scramble-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'episodic_stress-scramble-ns-sub_Sch400.pscalar.nii')
create_cifti(text_file = 'chronic_stress-scramble-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'chronic_stress-scramble-ns-sub_Sch400.pscalar.nii')
create_cifti(text_file = 'negative_affect-scramble-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'negative_affect-scramble-ns-sub_Sch400.pscalar.nii')
create_cifti(text_file = 'sleep-scramble-ns-sub.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'sleep-scramble-ns-sub_Sch400.pscalar.nii')

create_cifti(text_file = 'GAD7-ns.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', 
             cifti_out = 'GAD7-ns_Sch400.pscalar.nii')
create_cifti(text_file = 'PHQ9-ns.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', 
             cifti_out = 'PHQ9-ns_Sch400.pscalar.nii')

write_rez_table(reliability_dt_labels, 'b_TIMECENTER')
write_rez_table(reliability_dt_labels, 'b_WCEN_EPISODICTOT')
write_rez_table(reliability_dt_labels, 'b_WCEN_CHRONICSEV')
write_rez_table(reliability_dt_labels, 'b_WCEN_PANASNEG')
write_rez_table(reliability_dt_labels, 'b_WCEN_SLEEPDUR')

write_rez_table(reliability_dt_labels, 'b_TIMECENTER', contrast_name = 'scramble')
write_rez_table(reliability_dt_labels, 'b_WCEN_EPISODICTOT', contrast_name = 'scramble')
write_rez_table(reliability_dt_labels, 'b_WCEN_CHRONICSEV', contrast_name = 'scramble')
write_rez_table(reliability_dt_labels, 'b_WCEN_PANASNEG', contrast_name = 'scramble')
write_rez_table(reliability_dt_labels, 'b_WCEN_SLEEPDUR', contrast_name = 'scramble')

write_rez_table(reliability_dt_labels, 'bsp_meWCEN_BRAINWB_SE', outcome = 'GAD7')
write_rez_table(reliability_dt_labels, 'bsp_meWCEN_BRAINWB_SE', outcome = 'PHQ9')

#View(reliability_dt_labels[alt == '-alt' & roi_num > 400, c('label', 'rn', 'Estimate', 'p_nSSEst', 'Type_S')])
reliability_dt_labels[alt == '-alt' & roi_num > 400 & credible == TRUE, c('label', 'rn', 'Estimate', 'p_nSSEst', 'Type_S')]

write_cifti_text(x = reliability_dt_labels, variable = 'b_TIMECENTER', unthreshed = TRUE, psign = TRUE, outfile = 'time-ns-sub_psign_unthreshed.txt')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_EPISODICTOT', unthreshed = TRUE, psign = TRUE, outfile = 'episodic_stress-ns-sub_psign_unthreshed.txt')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_CHRONICSEV', unthreshed = TRUE, psign = TRUE, outfile = 'chronic_stress-ns-sub_psign_unthreshed.txt')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_PANASNEG', unthreshed = TRUE, psign = TRUE, outfile = 'negative_affect-ns-sub_psign_unthreshed.txt')
write_cifti_text(x = reliability_dt_labels, variable = 'b_WCEN_SLEEPDUR_NEG', unthreshed = TRUE, psign = TRUE, outfile = 'sleep-ns-sub_psign_unthreshed.txt')

create_cifti(text_file = 'time-ns-sub_psign_unthreshed.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'time-ns-sub_psign_unthreshed_Sch400.pscalar.nii')
create_cifti(text_file = 'episodic_stress-ns-sub_psign_unthreshed.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'episodic_stress-ns-sub_psign_unthreshed_Sch400.pscalar.nii')
create_cifti(text_file = 'chronic_stress-ns-sub_psign_unthreshed.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'chronic_stress-ns-sub_psign_unthreshed_Sch400.pscalar.nii')
create_cifti(text_file = 'negative_affect-ns-sub_psign_unthreshed.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'negative_affect-ns-sub_psign_unthreshed_Sch400.pscalar.nii')
create_cifti(text_file = 'sleep-ns-sub_psign_unthreshed.txt', cifti_template = 'Schaefer2018_CIFTI/Schaefer2018_400Parcels_7Networks_order.pscalar.nii', cifti_out = 'sleep-ns-sub_psign_unthreshed_Sch400.pscalar.nii')


reliability_dt_labels[credible %in% TRUE & 
                        roi_num > 400, lapply(.SD, round, 3), 
                      .SDcols = c('Estimate', 'Q2.5', 'Q97.5'), 
                      by = c('label', 'pred')]


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

##Ranges----
#prop var sess-id not-smoothed
cat(paste(
  c(
    range_text(reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'prop_var' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% 1:400 & cntrlreg == '-cntrlreg', Estimate],
               "Calm - Range of ICC_within-person across control regions:"),
    range_text(reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'prop_var' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% c(401, 403, 412, 414) & cntrlreg == '', Estimate],
               "Calm - Range of ICC_within-person across ventricle, and white matter:"),
    range_text(reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'prop_var' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% c(404:407, 409:411, 415:421) & cntrlreg == '', Estimate],
               "Calm - Range of ICC_within-person across subcortial (no brain-stem, ventricle, or white matter):"),
    range_text(reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'prop_var' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% 1:400 & cntrlreg == '' & network == '', Estimate],
               "Calm - Range of ICC_within-person across Sch400 cortical regions:"),
    range_text(reliability_dt_labels[contrast == 'scramble' & pred == '' & rn == 'prop_var' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% c(401, 403, 412, 414) & cntrlreg == '', Estimate],
               "Scramble - Range of ICC_within-person across ventricle, and white matter:"),
    range_text(reliability_dt_labels[contrast == 'scramble' & pred == '' & rn == 'prop_var' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% c(404:407, 409:411, 415:421) & cntrlreg == '', Estimate],
               "Scramble - Range of ICC_within-person across subcortial (no brain-stem, ventricle, or white matter):"),
    range_text(reliability_dt_labels[contrast == 'scramble' & pred == '' & rn == 'prop_var' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% 1:400 & cntrlreg == '', Estimate],
               "Scramble - Range of ICC_within-person across Sch400 cortical regions:")#,
    # range_text(oldicc_ests_meta_cont_cope_dt[roi_tag == 'HO' & ROI %in% HO_sub_cortex, Estimate],
    #            "Range of ICC_between-person across HO cortex:")
  ),
  collapse = '\n'))

reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'icc_ma_id' & sub == '' & ns == '-ns' & cntrlreg == '']
reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'icc_ma_id' & sub == '' & ns == '-ns' & cntrlreg == '-cntrlreg']

cat(paste(
  c(
    range_text(reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'icc_ma_id' & sub == '' & ns == '-ns' & 
                                       roi_num %in% 1:400 & cntrlreg == '-cntrlreg', Estimate],
               "Calm - Range of ICC_betwee-person across control regions:"),
    range_text(reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'icc_ma_id' & sub == '' & ns == '-ns' & 
                                       roi_num %in% c(401, 403, 412, 414) & cntrlreg == '', Estimate],
               "Calm - Range of ICC_betwee-person across ventricle, and white matter:"),
    range_text(reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'icc_ma_id' & sub == '' & ns == '-ns' & 
                                       roi_num %in% c(404:407, 409:411, 415:421) & cntrlreg == '', Estimate],
               "Calm - Range of ICC_betwee-person across subcortial (no brain-stem, ventricle, or white matter):"),
    range_text(reliability_dt_labels[contrast == 'calm' & pred == '' & rn == 'icc_ma_id' & sub == '' & ns == '-ns' & 
                                       roi_num %in% 1:400 & cntrlreg == '', Estimate],
               "Calm - Range of ICC_betwee-person across Sch400 cortical regions:"),
    range_text(reliability_dt_labels[contrast == 'scramble' & pred == '' & rn == 'prop_var_id' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% c(401, 403, 412, 414) & cntrlreg == '', Estimate],
               "Scramble - Range of ICC_betwee-person across ventricle, and white matter:"),
    range_text(reliability_dt_labels[contrast == 'scramble' & pred == '' & rn == 'prop_var_id' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% c(404:407, 409:411, 415:421) & cntrlreg == '', Estimate],
               "Scramble - Range of ICC_betwee-person across subcortial (no brain-stem, ventricle, or white matter):"),
    range_text(reliability_dt_labels[contrast == 'scramble' & pred == '' & rn == 'prop_var_id' & sub == '-sub' & ns == '-ns' & 
                                       roi_num %in% 1:400 & cntrlreg == '', Estimate],
               "Scramble - Range of ICC_betwee-person across Sch400 cortical regions:")#,
    # range_text(oldicc_ests_meta_cont_cope_dt[roi_tag == 'HO' & ROI %in% HO_sub_cortex, Estimate],
    #            "Range of ICC_between-person across HO cortex:")
  ),
  collapse = '\n'))

#prop var sess-id smoothed
cat(paste(
  c(
    range_text(reliability_dt_labels[contrast == 'calm' & rn == 'prop_var' & sub == '-sub' & ns == '' & roi_num %in% c(401, 403, 412, 414), Estimate],
               "Calm - Range of ICC_within-person across ventricle, and white matter:"),
    range_text(reliability_dt_labels[contrast == 'calm' & rn == 'prop_var' & sub == '-sub' & ns == '' & roi_num %in% c(404:407, 409:411, 415:421), Estimate],
               "Calm - Range of ICC_within-person across subcortial (no brain-stem, ventricle, or white matter):"),
    range_text(reliability_dt_labels[contrast == 'calm' & rn == 'prop_var' & sub == '-sub' & ns == '' & roi_num %in% 1:400, Estimate],
               "Calm - Range of ICC_within-person across Sch400 cortical regions:")#,
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

##Effect-size and within-person reliability


##Smoothness----
unique(reliability_dt_labels$rn)
setkey(reliability_dt_labels, roi_num, cntrlreg)
smoothness_data[, roi := as.numeric(roi)]
setnames(smoothness_data, 'roi', 'roi_num')
smoothness_data[, cntrlreg := fifelse(anatomy == 'cntrl', '-cntrlreg', '')]

smoothness_data_avg <- smoothness_data[, .(smoothness = mean(smoothness, na.rm = TRUE)), 
                                       by = c('ID', 'roi_num', 'anatomy', 'cntrlreg')][, .(smoothness = mean(smoothness, na.rm = TRUE)), 
                                                                                       by = c('roi_num', 'anatomy', 'cntrlreg')]

smoothness_rel_data <- merge(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & pred == '' & contrast == 'calm', 
                                                   c('roi_num', 'Estimate', 'cntrlreg')],
                             smoothness_data_avg, on = c('roi_num', 'cntrlreg'), all = TRUE)
smoothness_rel_data <- smoothness_rel_data[roi_num %in% c(1:400, 404:407, 409:411, 415:421)]

smoothness_rel_data[, c('rank_Estimate', 'rank_smoothness') := list(rank(Estimate, na.last = 'keep'), rank(smoothness, na.last = 'keep')), by = 'anatomy']

sch400HO <- fread(file.path('tsv_nosmooth', 'Sch400_HO_combined.tsv'), sep = ' ',  col.names = c('i', 'j', 'k', 'ROI'))
cntrlregions <- load_control_regions('tsv_nosmooth')
cntrlregions_sizes <- cntrlregions[, .(N_voxels = .N), by = c('eq_ROI')]
cntrlregions_sizes[, c('cntrlreg') := list(rep('-cntrlreg', .N))]
setnames(cntrlregions_sizes, 'eq_ROI', 'roi_num')
sch400HO_sizes <- sch400HO[, .(N_voxels = .N), by = c('ROI')]
sch400HO_sizes[, c('cntrlreg') := list(rep('', .N))]
setnames(sch400HO_sizes, 'ROI', 'roi_num')
smoothness_rel_data_full <- merge(reliability_dt_labels[rn == 'prop_var' & sub == '-sub' & ns == '-ns' & pred == '' & contrast == 'calm' & network == '', 
                                                        c('roi_num', 'Estimate', 'Est.Error', 'cntrlreg')],
                                  smoothness_data, on = c('roi_num', 'cntrlreg'), all = TRUE)
smoothness_rel_data_full <- merge(smoothness_rel_data_full, 
                                  rbindlist(list(cntrlregions_sizes, sch400HO_sizes)), 
                                  on = c('roi_num', 'cntrlreg'))
smoothness_rel_data_full <- smoothness_rel_data_full[roi_num %in% c(1:400, 404:407, 409:411, 415:421)]

range(smoothness_rel_data_full$smoothness, na.rm = TRUE)



smoothness_rel_data_full_sum <- smoothness_rel_data_full[, .(smooth_mean = mean(smoothness, na.rm = TRUE),
                                                             smooth_se = sd(smoothness, na.rm = TRUE) / sqrt(sum(!is.na(smoothness)))),
                                                         by = c('roi_num', 'Estimate', 'Est.Error', 'anatomy', 'N_voxels')]
smoothness_rel_data_full_sum[, .(N = sum(!is.na(smooth_se)), R = sum(is.na(smooth_se)), pR = sum(is.na(smooth_se)) / sum(!is.na(smooth_se)))]

library(brms)
smoothness_rel_data_full_sum[is.na(smooth_se), smooth_mean := NA]
smoothness_rel_data_full_sum[is.na(smooth_se), smooth_se := Inf]
smoothness_rel_data_full_sum[, scale_N_voxels := (N_voxels - 500) / 100]
smoothness_rel_data_full_sum[, anatomy_fac := factor(anatomy, levels = c('cort', 'subcort', 'cntrl'))]
hist(smoothness_rel_data_full_sum$scale_N_voxels)
hist(smoothness_rel_data_full_sum$Estimate)
hist(smoothness_rel_data_full_sum$smooth_mean)
smoothness_rel_data_full_sum_nomiss <- smoothness_rel_data_full_sum[!is.na(smooth_mean)]
smooth_rel_formula <- 
  bf(Estimate | mi(Est.Error) ~ 1 + anatomy_fac*scale_N_voxels + anatomy_fac*smooth_mean, 
     phi ~ 1 + (1 | anatomy_fac), 
     family = Beta())
get_prior(smooth_rel_formula, data = smoothness_rel_data_full_sum_nomiss)  
smooth_rel_prior <- c(prior('normal(0, 1)', class = 'b'))
smooth_rel_fit <- brm(formula = smooth_rel_formula, 
                      prior = smooth_rel_prior, 
                      data = smoothness_rel_data_full_sum_nomiss,
                      chains = 4, cores = 4, iter = 2000, warmup = 1000,
                      control = list(adapt_delta = .9999, max_treedepth = 20),
                      file = 'smooth_rel', file_refit = 'on_change', backend = 'rstan')

smooth_rel_s_formula <- 
  bf(Estimate | mi(Est.Error) ~ 1 + anatomy_fac + s(scale_N_voxels, by = anatomy_fac) + s( smooth_mean, by = anatomy_fac), 
     phi ~ 1 + (1 | anatomy_fac), 
     family = Beta())
get_prior(smooth_rel_s_formula, data = smoothness_rel_data_full_sum_nomiss)  
smooth_rel_s_prior <- c(prior('normal(0, 1)', class = 'b'),
                        prior('student_t(3, 0, 2)', class = 'sds', lb = 0))
smooth_rel_s_fit <- brm(formula = smooth_rel_s_formula, 
                      prior = smooth_rel_s_prior, 
                      data = smoothness_rel_data_full_sum_nomiss,
                      chains = 4, cores = 4, iter = 2000, warmup = 1000,
                      control = list(adapt_delta = .9999, max_treedepth = 20),
                      file = 'smooth_rel_s', file_refit = 'on_change', backend = 'rstan')
prior_summary(smooth_rel_s_fit)

smooth_rel_s_alt_formula <- 
  bf(Estimate | mi(Est.Error) ~ 1 + s(scale_N_voxels, anatomy_fac, bs = 'fs') + s(smooth_mean, anatomy_fac, bs = 'fs'), 
     phi ~ 1 + (1 | anatomy_fac), 
     family = Beta())
get_prior(smooth_rel_s_alt_formula, data = smoothness_rel_data_full_sum_nomiss)  
smooth_rel_s_alt_prior <- c(prior('student_t(3, 0, 2)', class = 'sds', lb = 0))
smooth_rel_s_alt_fit <- brm(formula = smooth_rel_s_alt_formula, 
                        prior = smooth_rel_s_alt_prior, 
                        data = smoothness_rel_data_full_sum_nomiss,
                        chains = 4, cores = 4, iter = 2000, warmup = 1000,
                        control = list(adapt_delta = .9999, max_treedepth = 20),
                        file = 'smooth_rel_s_alt', file_refit = 'on_change', backend = 'rstan')

smooth_rel_s_null_formula <-
  bf(Estimate | mi(Est.Error) ~ 1 + scale_N_voxels +  smooth_mean + (1 | anatomy_fac), 
     phi ~ 1 + (1 | anatomy_fac), 
     family = Beta())
get_prior(smooth_rel_s_null_formula, data = smoothness_rel_data_full_sum_nomiss)
smooth_rel_s_null_prior <- c(prior('normal(0, 1)', class = 'b'))
smooth_rel_s_null_fit <- brm(formula = smooth_rel_s_null_formula,
                        prior = smooth_rel_s_null_prior,
                        data = smoothness_rel_data_full_sum_nomiss,
                        chains = 4, cores = 4, iter = 2000, warmup = 1000,
                        control = list(adapt_delta = .9999, max_treedepth = 20),
                        file = 'smooth_rel_s_null', file_refit = 'on_change', backend = 'rstan')

smooth_rel_s_null_null_formula <-
  bf(Estimate | mi(Est.Error) ~ 1 + scale_N_voxels +  (1 | anatomy_fac), 
     phi ~ 1 + (1 | anatomy_fac), 
     family = Beta())
get_prior(smooth_rel_s_null_null_formula, data = smoothness_rel_data_full_sum_nomiss)
smooth_rel_s_null_null_prior <- c(prior('normal(0, 1)', class = 'b'))
smooth_rel_s_null_null_fit <- brm(formula = smooth_rel_s_null_null_formula,
                                  prior = smooth_rel_s_null_null_prior,
                                  data = smoothness_rel_data_full_sum_nomiss,
                                  chains = 4, cores = 4, iter = 2000, warmup = 1000,
                                  control = list(adapt_delta = .99999, max_treedepth = 20),
                                  file = 'smooth_rel_s_null_null', file_refit = 'on_change', backend = 'rstan')

smooth_rel_s_null_null_alt_formula <-
  bf(Estimate | mi(Est.Error) ~ 1 + scale_N_voxels + anatomy_fac, 
     phi ~ 1 + (1 | anatomy_fac), 
     family = Beta())
get_prior(smooth_rel_s_null_null_alt_formula, data = smoothness_rel_data_full_sum_nomiss)
smooth_rel_s_null_null_alt_prior <- c(prior('normal(0, 1)', class = 'b'))
smooth_rel_s_null_null_alt_fit <- brm(formula = smooth_rel_s_null_null_alt_formula,
                                  prior = smooth_rel_s_null_null_alt_prior,
                                  data = smoothness_rel_data_full_sum_nomiss,
                                  chains = 4, cores = 4, iter = 2000, warmup = 1000,
                                  control = list(adapt_delta = .99999, max_treedepth = 20),
                                  file = 'smooth_rel_s_null_null_alt', file_refit = 'on_change', backend = 'rstan')


smooth_rel_s_null_nullest_formula <-
  bf(Estimate | mi(Est.Error) ~ 1 + scale_N_voxels, 
     phi ~ 1 + (1 | anatomy_fac), 
     family = Beta())
get_prior(smooth_rel_s_null_nullest_formula, data = smoothness_rel_data_full_sum_nomiss)
smooth_rel_s_null_nullest_prior <- c(prior('normal(0, 1)', class = 'b'))
smooth_rel_s_null_nullest_fit <- brm(formula = smooth_rel_s_null_nullest_formula,
                                  prior = smooth_rel_s_null_nullest_prior,
                                  data = smoothness_rel_data_full_sum_nomiss,
                                  chains = 4, cores = 4, iter = 2000, warmup = 1000,
                                  control = list(adapt_delta = .99999, max_treedepth = 20),
                                  file = 'smooth_rel_s_null_nullest', file_refit = 'on_change', backend = 'rstan')


smooth_rel_s_alt_fit_CE_data <- as.data.table(smooth_rel_s_alt_fit$data)
knitr::kable(as.data.table(smooth_rel_s_alt_fit_CE_data)[, transpose(.(quantile(scale_N_voxels, probs = c(0, .25, .5, .75, 1), na.rm = TRUE))), by = c('anatomy_fac')],
             col.names = c('Anatomy', 'Q0', 'Q25', 'Q50', 'Q75', 'Q100'), format = 'simple', caption = 'N Voxels Scaled')
knitr::kable(as.data.table(smooth_rel_s_alt_fit_CE_data)[, transpose(.(quantile(scale_N_voxels * 100 + 500, probs = c(0, .25, .5, .75, 1), na.rm = TRUE))), by = c('anatomy_fac')],
             col.names = c('Anatomy', 'Q0', 'Q25', 'Q50', 'Q75', 'Q100'), format = 'simple', caption = 'N Voxels')
knitr::kable(as.data.table(smooth_rel_s_alt_fit_CE_data)[, transpose(.(quantile(smooth_mean, probs = c(0, .25, .5, .75, 1), na.rm = TRUE))), by = c('anatomy_fac')],
             col.names = c('Anatomy', 'Q0', 'Q25', 'Q50', 'Q75', 'Q100'), format = 'simple', caption = 'Smoothness')

smooth_rel_s_alt_fit_cor_dt <- as.data.table(smooth_rel_s_alt_fit$data)

sdcols <- c('Estimate', 'scale_N_voxels', 'smooth_mean')
smooth_rel_s_alt_fit_cor_dt[, paste('rank', sdcols, sep = '_') :=
                              lapply(.SD, rank), .SDcols = sdcols, by = 'anatomy_fac']

anatomy_fac_labeller <- labeller(anatomy_fac = c('cort' = 'Cortex',
                                                 subcort = 'Subcortex',
                                                 cntrl = 'Control'))
ggplot(smooth_rel_s_alt_fit_cor_dt, aes(x = rank_scale_N_voxels, y = rank_smooth_mean)) + 
  geom_point() + 
  geom_line(stat = 'smooth', method = 'lm', color = 'blue') + 
  facet_wrap(~ anatomy_fac, nrow = 1, scales = 'free', labeller = anatomy_fac_labeller) + 
  # scale_x_continuous(position = 'top') + 
  labs(x = '', y = 'Smoothness rank') + 
  theme_minimal() + 
  theme(axis.text = element_blank(), panel.grid = element_blank(), 
        panel.background = element_rect(fill = '#F9F9F9', color = '#666666')) + 
  ggplot(smooth_rel_s_alt_fit_cor_dt, aes(x = rank_scale_N_voxels, y = rank_Estimate)) + 
    geom_point() + 
    geom_line(stat = 'smooth', method = 'lm', color = 'blue') + 
    facet_wrap(~ anatomy_fac, nrow = 1, scales = 'free', labeller = anatomy_fac_labeller) + 
    labs(x = 'Size rank', y = 'Within-person reliability rank') + 
    theme_minimal() + 
    theme(axis.text = element_blank(), panel.grid = element_blank(), strip.text = element_blank(),
          panel.background = element_rect(fill = '#F9F9F9', color = '#666666')) + 
  ggplot(smooth_rel_s_alt_fit_cor_dt, aes(x = rank_smooth_mean, y = rank_Estimate)) + 
    geom_point() + 
    geom_line(stat = 'smooth', method = 'lm', color = 'blue') + 
    facet_wrap(~ anatomy_fac, nrow = 1, scales = 'free', labeller = anatomy_fac_labeller) + 
    labs(x = 'Smoothness rank', y = '') + 
    theme_minimal() + 
    theme(axis.text = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = '#F9F9F9', color = '#666666')) + 
  plot_layout(design = 'A#\nBC')
  

hist(smooth_rel_s_alt_fit_CE_data$scale_N_voxels * 100 + 500)

scale_N_voxels_values <- quantile(smooth_rel_s_alt_fit_CE_data$scale_N_voxels, probs =  (c(1/3, 2/3, 1) - c(0, 1/3, 2/3)) / 2 + c(0, 1/3, 2/3))
names(scale_N_voxels_values) <- sprintf('Parcel Size: %0.0f voxels', (scale_N_voxels_values*100+500))
smooth_rel_s_alt_fit_CE <- conditional_effects(smooth_rel_s_alt_fit, effects = c('smooth_mean:anatomy_fac'),
                                               conditions = data.frame('scale_N_voxels' = scale_N_voxels_values))
smooth_rel_s_alt_fit_CE_data[, cond__ := cut(scale_N_voxels, breaks = quantile(scale_N_voxels, probs = c(0, 1/3, 2/3, 1)), 
                                                          include.lowest = TRUE,
                                                          labels = names(scale_N_voxels_values))]
smooth_rel_s_alt_fit_CE_data[, effect2__ := anatomy_fac]
smooth_model_unconstrained_plot <- ggplot(smooth_rel_s_alt_fit_CE$`smooth_mean:anatomy_fac`,
                                          aes(x = smooth_mean, y = estimate__)) + 
  geom_errorbar(data = smooth_rel_s_alt_fit_CE_data, 
                aes(x = smooth_mean, y = Estimate, color = effect2__,
                    ymin = Estimate - 2*Est.Error, ymax = Estimate + 2*Est.Error), width = 0, alpha = .2) + 
  geom_point(data = smooth_rel_s_alt_fit_CE_data, aes(x = smooth_mean, y = Estimate, color = effect2__), size = 1) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__, group = effect2__, fill = effect2__), alpha = .4) + 
  geom_line(aes(group = effect2__, color = effect2__), alpha = 1) + 
  scale_x_continuous(trans = scales::modulus_trans(p = -.5), breaks = c(5, 10, 16)) + 
  scale_color_manual(name = 'Anatomy', aesthetics = c('color', 'fill'), values = c('#7570b3', '#d95f02', '#1b9e77'),
                     breaks = c('cort', 'subcort', 'cntrl'), labels = c('Cortex', 'Subcortex', 'Control')) + 
  facet_grid(~ cond__) + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) + 
  labs(y = 'Within-person reliability estimate (95% CI)',
       x = 'Parcel mean smoothness')
print(smooth_model_unconstrained_plot)
ggsave('smooth_model_unconstrained_plot.svg', plot = smooth_model_unconstrained_plot, width = 8, height = 5, dpi = 150)

smooth_rel_s_null_fit_conditions <- as.data.frame(expand.grid('scale_N_voxels' = scale_N_voxels_values))
smooth_rel_s_null_fit_CE <- conditional_effects(smooth_rel_s_null_null_fit, effects = 'anatomy_fac',
                                                conditions = smooth_rel_s_null_fit_conditions,
                                                re_formula = NULL)

smooth_rel_s_null_fit_conditions$facet <- names(scale_N_voxels_values)
smooth_rel_s_null_fit_CE$anatomy_fac <- as.data.table(smooth_rel_s_null_fit_CE$anatomy_fac)
smooth_rel_s_null_fit_CE$anatomy_fac <- smooth_rel_s_null_fit_CE$anatomy_fac[smooth_rel_s_null_fit_conditions, on = c('scale_N_voxels')]
smooth_rel_s_null_fit_CE$anatomy_fac[, cond__ := facet]
smooth_model_constrained_plot <- ggplot(smooth_rel_s_null_fit_CE$anatomy_fac,
       aes(x = anatomy_fac, y = estimate__)) + 
  geom_point(data = smooth_rel_s_alt_fit_CE_data, aes(x = anatomy_fac, y = Estimate), 
             alpha = .1,
             position = position_jitter(h = 0, w = .1)) + 
  geom_errorbar(aes(ymin = lower__, ymax = upper__), alpha = 1, width = .1, linewidth = 1) + 
  geom_point(aes(group = anatomy_fac), alpha = 1, size = 3) + 
  scale_x_discrete(breaks = c('cort', 'subcort', 'cntrl'), labels = c('Cortex', 'Subcortex', 'Control')) + 
  facet_grid(~ cond__) + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) + 
  labs(y = 'Within-person reliability estimate (95% CI)',
       x = '')
print(smooth_model_constrained_plot)
ggsave('smooth_model_constrained_plot.svg', plot = smooth_model_constrained_plot, width = 8, height = 5, dpi = 150)

variables(smooth_rel_s_null_fit)
(smooth_rel_s_null_null_fit_H <- hypothesis(smooth_rel_s_null_null_fit, robust = TRUE,
           hypothesis = c('r_anatomy_fac__phi[cort,Intercept] - (r_anatomy_fac__phi[subcort,Intercept] + r_anatomy_fac__phi[cntrl,Intercept])/2 > 0',
                          'r_anatomy_fac__phi[cort,Intercept] - r_anatomy_fac__phi[subcort,Intercept] > 0',
                          'r_anatomy_fac__phi[cort,Intercept] - r_anatomy_fac__phi[cntrl,Intercept] > 0'),
           class = NULL))
(smooth_rel_s_alt_fit_H <- hypothesis(smooth_rel_s_alt_fit, robust = TRUE,
           hypothesis = c('r_anatomy_fac__phi[cort,Intercept] - (r_anatomy_fac__phi[subcort,Intercept] + r_anatomy_fac__phi[cntrl,Intercept])/2 > 0',
                          'r_anatomy_fac__phi[cort,Intercept] - r_anatomy_fac__phi[subcort,Intercept] > 0',
                          'r_anatomy_fac__phi[cort,Intercept] - r_anatomy_fac__phi[cntrl,Intercept] > 0'),
           class = NULL))


round(posterior_summary(posterior_epred(smooth_rel_s_null_null_fit, 
                                        newdata = data.frame(scale_N_voxels = (312 - 500) / 100,
                                                             anatomy_fac = c('cort', 'subcort', 'cntrl'),
                                                             Est.Error = 1),
                                        re_formula = NULL), 
                        robust = TRUE), 2)

round(posterior_summary(posterior_predict(smooth_rel_s_null_null_fit, 
                                          newdata = data.frame(scale_N_voxels = (312 - 500) / 100,
                                                               anatomy_fac = c('cort', 'subcort', 'cntrl'),
                                                               Est.Error = 1),
                                          re_formula = NULL), 
                        robust = TRUE), 2)

# smooth_rel_s_alt_fit <- add_criterion(smooth_rel_s_alt_fit, criterion = c('loo', 'bayes_R2'))
# smooth_rel_s_null_fit <- add_criterion(smooth_rel_s_null_fit, criterion = c('loo', 'bayes_R2'))
# smooth_rel_s_null_null_fit <- add_criterion(smooth_rel_s_null_null_fit, criterion = c('loo', 'bayes_R2'))
# smooth_rel_s_null_nullest_fit <- add_criterion(smooth_rel_s_null_nullest_fit, criterion = c('loo', 'bayes_R2'))
# smooth_rel_s_null_null_alt_fit <- add_criterion(smooth_rel_s_null_null_alt_fit, criterion = c('loo', 'bayes_R2'))
loo_compare(smooth_rel_s_alt_fit, smooth_rel_s_null_fit, smooth_rel_s_null_null_fit, smooth_rel_s_null_nullest_fit)
loo_compare(smooth_rel_s_alt_fit, smooth_rel_s_null_null_fit, model_names = c('Unconstrained', 'Constrained'))
loo_compare(smooth_rel_s_alt_fit, smooth_rel_s_null_nullest_fit)
loo_compare(smooth_rel_s_null_null_fit, smooth_rel_s_null_nullest_fit)
loo(smooth_rel_s_fit)
loo(smooth_rel_s_null_null_fit)
loo_compare(smooth_rel_s_fit, smooth_rel_s_null_null_fit, model_names = c('Unconstrained', 'Constrained'))
loo_compare(smooth_rel_s_fit, smooth_rel_s_null_null_alt_fit, model_names = c('Unconstrained', 'Constrained'))
loo_compare(smooth_rel_s_alt_fit, smooth_rel_s_null_null_alt_fit)
loo_compare(smooth_rel_s_null_null_fit, smooth_rel_s_null_null_alt_fit, smooth_rel_s_null_nullest_fit)

pp_check(smooth_rel_s_null_null_fit) + theme_minimal()
pp_check(smooth_rel_s_null_nullest_fit) + theme_minimal()
#Roughly the same, but we prefer the one that allows us to estimate the difference between anatomy
#                               elpd_diff se_diff
# smooth_rel_s_null_null_fit      0.0       0.0  
# smooth_rel_s_null_fit          -0.6       0.7  
# smooth_rel_s_null_nullest_fit  -4.3       3.7  
# smooth_rel_s_alt_fit          -31.2      23.9 
list(bayes_R2(smooth_rel_s_alt_fit),
     bayes_R2(smooth_rel_s_null_fit),
     bayes_R2(smooth_rel_s_null_null_fit))

pp_check(smooth_rel_s_alt_fit, type = 'dens_overlay_grouped', group = 'anatomy_fac') + theme_minimal()
pp_check(smooth_rel_s_null_fit, type = 'dens_overlay_grouped', group = 'anatomy_fac') + theme_minimal()


##Priors-----

test_retest_ex_model <- readRDS('~/3dICC/bicc/fits/fit-ROI1-01_02-icc.rds')
test_retest_ex_model$formula
test_retest_ex_model$family
prior_summary(test_retest_ex_model)
test_retest_MA_ex_model <- readRDS('~/3dICC/bicc/fits/fit-ROI1-meta_icc.rds')
test_retest_MA_ex_model$formula
test_retest_MA_ex_model$family
prior_summary(test_retest_MA_ex_model)
interncons_ex_model <- readRDS('fits/alpha-fearGTcalm-ROI001-ns.rds')
interncons_ex_model$formula
interncons_ex_model$family
prior_summary(interncons_ex_model)

smooth_rel_s_alt_fit$formula
smooth_rel_s_alt_fit$family
prior_summary(smooth_rel_s_alt_fit)
smooth_rel_s_null_null_fit$formula
smooth_rel_s_null_null_fit$family
prior_summary(smooth_rel_s_null_null_fit)


##Prediction correspondence----
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
ctrl_regions <- load_control_regions(tsv_dir = 'tsv_nosmooth')
readr::write_delim(ctrl_regions[, c('i', 'j', 'k', 'eq_ROI')], 'ctrl_regions.tsv', delim = ' ', col_names = FALSE)

#prop_var_sub_out_dt <- reliability_dt[rn == 'prop_var' & sub == '-sub' & ns == '-ns', ][sch400_HO_dump, on = 'roi_num']
#prop_var_sub_out_dt[is.na(Estimate), c('Estimate', 'Q2.5', 'Q97.5') := 0]

#readr::write_delim(prop_var_sub_out_dt[, c('i', 'j', 'k', 'Estimate')], 'prop_id-sess_var_15voxsub.tsv', delim = ' ', col_names = FALSE)


##Reliability and effect sizes correlations----

rough_cor_pred_vals <- c('-pred-alt_time', '-pred-alt_chsev', '-pred-alt_eptot',
                         '-pred-alt_panasneg', '-pred-alt_sleep')

reliability_dt_labels_rough_cor_data_pred <- reliability_dt_labels[contrast == 'calm' & 
                                                                     pred %in% rough_cor_pred_vals & 
                                                                     sub == '-sub' & ns == '-ns' & 
                                                                     roi_num %in% c(1:400, 404:407, 409:411, 415:421) & 
                                                                     cntrlreg == '' &
                                                                     network == '',
                                                                   c('rn', 'roi_num', 'Estimate')]
reliability_dt_labels_rough_cor_data_pred[, Estimate := abs(Estimate)]
reliability_dt_labels_rough_cor_data_pred_w <- dcast(reliability_dt_labels_rough_cor_data_pred, 
                                                     roi_num ~ rn, value.var = 'Estimate')

reliability_dt_labels_rough_cor_data_rel <- reliability_dt_labels[contrast == 'calm' & 
                                                                    rn == 'prop_var' & 
                                                                    pred == '' & 
                                                                    sub == '-sub' & ns == '-ns' & 
                                                                    roi_num %in% c(1:400, 404:407, 409:411, 415:421) & 
                                                                    cntrlreg == '' & 
                                                                    network == '', 
                                                                  c('rn', 'roi_num', 'Estimate')]
setnames(reliability_dt_labels_rough_cor_data_rel, 'Estimate', 'int_cons')
reliability_dt_labels_rough_cor_data <- merge(reliability_dt_labels_rough_cor_data_rel, 
                                              reliability_dt_labels_rough_cor_data_pred_w,
                                              on = c('roi_num'), all = TRUE)

rel_pred_cors <- psych::corr.test(reliability_dt_labels_rough_cor_data[,3],
                                  reliability_dt_labels_rough_cor_data[,-(1:3)], 
                                  method = 'spearman')
paste(gsub('0\\.', '.', 
           sprintf('%s: r = %0.2f (p = %0.3f)',
                   dimnames(rel_pred_cors$r)[[2]],
                   rel_pred_cors$r,
                   rel_pred_cors$p.adj)),
      collapse = '; ')

variable_names <- c('Time', 'Chronic Stress', 'Stressful Life Events', 'Negative Affect', 'Sleep') 
names(variable_names) <- dimnames(rel_pred_cors$r)[[2]]
reliability_dt_labels_rough_cor_plot_data <- merge(reliability_dt_labels_rough_cor_data_rel[, -1], 
                                              reliability_dt_labels_rough_cor_data_pred,
                                              on = c('roi_num'), all = TRUE)
reliability_dt_labels_rough_cor_plot_data[, rank_int_cons := rank(int_cons), by = 'rn']
reliability_dt_labels_rough_cor_plot_data[, rank_Estimate := rank(Estimate), by = 'rn']

rel_pred_cor_plot <- ggplot(reliability_dt_labels_rough_cor_plot_data, 
                            aes(x = rank_int_cons, y = rank_Estimate)) + 
  geom_point() + 
  geom_line(stat = 'smooth', method = 'gam') + 
  facet_wrap(~ rn, labeller = as_labeller(variable_names)) + 
  theme_minimal()
print(rel_pred_cor_plot)

rel_pred_cor_raw_plot <- ggplot(reliability_dt_labels_rough_cor_plot_data, 
                            aes(x = int_cons, y = Estimate)) + 
  geom_point() + 
  geom_line(stat = 'smooth', method = 'gam') + 
  facet_wrap(~ rn, labeller = as_labeller(variable_names), scales = 'free') + 
  theme_minimal()
print(rel_pred_cor_raw_plot)

###Z scored----

reliability_dt_labels_rough_cor_data_predZ <- reliability_dt_labels[contrast == 'calm' & 
                                                                      pred %in% rough_cor_pred_vals & 
                                                                      sub == '-sub' & ns == '-ns' & 
                                                                      roi_num %in% c(1:400, 404:407, 409:411, 415:421) & 
                                                                      cntrlreg == '' &
                                                                      network == '',
                                                                    c('rn', 'roi_num', 'Estimate', 'Est.Error')]
reliability_dt_labels_rough_cor_data_predZ[, Zprime := abs(Estimate / Est.Error)]
reliability_dt_labels_rough_cor_data_predZ_w <- dcast(reliability_dt_labels_rough_cor_data_predZ, 
                                                      roi_num ~ rn, value.var = 'Zprime')

reliability_dt_labels_rough_cor_Z_data <- merge(reliability_dt_labels_rough_cor_data_rel, 
                                              reliability_dt_labels_rough_cor_data_predZ_w,
                                              on = c('roi_num'), all = TRUE)

rel_pred_Z_cors <- psych::corr.test(reliability_dt_labels_rough_cor_Z_data[,3],
                                  reliability_dt_labels_rough_cor_Z_data[,-(1:3)], 
                                  method = 'spearman')
paste(gsub('0\\.', '.', 
           sprintf('%s: r = %0.2f (p = %0.3f)',
                   dimnames(rel_pred_Z_cors$r)[[2]],
                   rel_pred_Z_cors$r,
                   rel_pred_Z_cors$p.adj)),
      collapse = '; ')

reliability_dt_labels_rough_cor_Z_plot_data <- merge(reliability_dt_labels_rough_cor_data_rel[, -1], 
                                                     reliability_dt_labels_rough_cor_data_predZ,
                                                     on = c('roi_num'), all = TRUE)
reliability_dt_labels_rough_cor_Z_plot_data[, rank_int_cons := rank(int_cons), by = 'rn']
reliability_dt_labels_rough_cor_Z_plot_data[, rank_Zprime := rank(Zprime), by = 'rn']

rel_pred_cor_Z_plot <- ggplot(reliability_dt_labels_rough_cor_Z_plot_data, 
                            aes(x = rank_int_cons, y = rank_Zprime)) + 
  geom_point() + 
  geom_line(stat = 'smooth', method = 'gam') + 
  facet_wrap(~ rn, labeller = as_labeller(variable_names)) + 
  theme_minimal()
print(rel_pred_cor_Z_plot)

rel_pred_cor_Z_raw_plot <- ggplot(reliability_dt_labels_rough_cor_Z_plot_data, 
                                aes(x = int_cons, y = Zprime)) + 
  geom_point() + 
  geom_line(stat = 'smooth', method = 'gam') + 
  facet_wrap(~ rn, labeller = as_labeller(variable_names), scales = 'free') + 
  theme_minimal()
print(rel_pred_cor_Z_raw_plot)

###Error----
sch400HO <- fread(file.path('tsv_nosmooth', 'Sch400_HO_combined.tsv'), sep = ' ',  col.names = c('i', 'j', 'k', 'ROI'))
sch400HO_sizes <- sch400HO[, .(N_voxels = .N), by = c('ROI')]
setnames(sch400HO_sizes, 'ROI', 'roi_num')

reliability_dt_labels_rough_cor_data_prederr <- reliability_dt_labels[contrast == 'calm' & 
                                                                      pred %in% rough_cor_pred_vals & 
                                                                      sub == '-sub' & ns == '-ns' & 
                                                                      roi_num %in% c(1:400, 404:407, 409:411, 415:421) & 
                                                                      cntrlreg == '' &
                                                                      network == '',
                                                                    c('rn', 'roi_num', 'Est.Error')]
reliability_dt_labels_rough_cor_data_prederr_w <- dcast(reliability_dt_labels_rough_cor_data_prederr, 
                                                      roi_num ~ rn, value.var = 'Est.Error')

reliability_dt_labels_rough_cor_err_data <- merge(reliability_dt_labels_rough_cor_data_rel, 
                                                reliability_dt_labels_rough_cor_data_prederr_w,
                                                on = c('roi_num'), all = TRUE)
reliability_dt_labels_rough_cor_err_data <- merge(sch400HO_sizes, 
                                                  reliability_dt_labels_rough_cor_err_data, 
                                                  on = 'roi_num')
rel_pred_err_cors <- psych::corr.test(reliability_dt_labels_rough_cor_err_data[,c(2,4)],
                                    reliability_dt_labels_rough_cor_err_data[,-(1:3)], 
                                    method = 'spearman')

paste(gsub('0\\.', '.', 
           sprintf('%s: r = %0.2f (p = %0.3f)',
                   paste(dimnames(rel_pred_err_cors$r)[[1]], 
                         rep(dimnames(rel_pred_err_cors$r)[[2]], each = 2)),
                   rel_pred_err_cors$r,
                   rel_pred_err_cors$p.adj)),
      collapse = '; ')

reliability_dt_labels_rough_cor_err_plot_data <- merge(reliability_dt_labels_rough_cor_data_rel[, -1], 
                                                     reliability_dt_labels_rough_cor_data_prederr,
                                                     on = c('roi_num'), all = TRUE)
reliability_dt_labels_rough_cor_err_plot_data[, rank_int_cons := rank(int_cons), by = 'rn']
reliability_dt_labels_rough_cor_err_plot_data[, rank_Est.Error := rank(Est.Error), by = 'rn']

rel_pred_cor_err_plot <- ggplot(reliability_dt_labels_rough_cor_err_plot_data, 
                              aes(x = rank_int_cons, y = rank_Est.Error)) + 
  geom_point() + 
  geom_line(stat = 'smooth', method = 'gam') + 
  facet_wrap(~ rn, labeller = as_labeller(variable_names)) + 
  theme_minimal()
print(rel_pred_cor_err_plot)

rel_pred_cor_err_raw_plot <- ggplot(reliability_dt_labels_rough_cor_err_plot_data, 
                                  aes(x = int_cons, y = Est.Error)) + 
  geom_point() + 
  geom_line(stat = 'smooth', method = 'gam') + 
  facet_wrap(~ rn, labeller = as_labeller(variable_names), scales = 'free') + 
  theme_minimal()
print(rel_pred_cor_err_raw_plot)

reliability_dt_labels_rough_cor_err_size_plot_data <- merge(sch400HO_sizes, 
                                                            reliability_dt_labels_rough_cor_err_plot_data,
                                                            on = c('roi_num'), all = TRUE)

rel_pred_cor_err_size_raw_plot <- ggplot(reliability_dt_labels_rough_cor_err_size_plot_data[!is.na(Est.Error)], 
                                    aes(x = N_voxels, y = Est.Error)) + 
  geom_point() + 
  geom_line(stat = 'smooth', method = 'gam') + 
  facet_wrap(~ rn, labeller = as_labeller(variable_names), scales = 'free') + 
  theme_minimal()
print(rel_pred_cor_err_size_raw_plot)

reliability_dt_labels_rough_cor_err_size_plot_data[, err_var := Est.Error^2]
saveRDS(reliability_dt_labels_rough_cor_err_size_plot_data, 'uberwootdata.rds')

library(brms)
intcon_err_formula <- bf(err_var ~ s(N_voxels, k = 10) + s(int_cons) + (1 | roi_num) + (1 | rn),
                         family = lognormal())
get_prior(intcon_err_formula, data = reliability_dt_labels_rough_cor_err_size_plot_data)
intcon_err_prior <- c(prior('normal(0,1)', class = 'b'), 
                      prior('student_t(3, 0, 1)', class = 'sd'))
intcon_err_fit <- brm(intcon_err_formula, 
                      data = reliability_dt_labels_rough_cor_err_size_plot_data,
                      prior = intcon_err_prior,
                      chains = 4, cores = 4, 
                      control = list(adapt_delta = .999, max_treedepth = 20),
                      file = 'intcon_err', file_refit = 'on_change')
