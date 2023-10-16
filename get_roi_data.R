get_roi_data <- function(roi = 1, cope = 'fearGTcalm', ncpus = 1, tsv_dir = 'tsv', cntrl_regions = FALSE){
  require(data.table)
  #tsv_dir = 'tsv_nosmooth'
  #ncpus = 4
  f <- dir(tsv_dir, pattern = sprintf('\\d{4}.*%s.tsv', cope), full.names = TRUE)
  message(sprintf('Number of files found: %d', length(f)))
  fregex <- '.*/(\\d{4})_(\\d{2})_(\\w+).tsv'
  
  if(cntrl_regions){
    sch400HO <- load_control_regions(tsv_dir = tsv_dir)[, c('i', 'j', 'k', 'eq_ROI')]
    setnames(sch400HO, 'eq_ROI', 'ROI')
    message(sprintf('Selecting from control regions (%d-%d)', min(sch400HO$ROI), max(sch400HO$ROI)))
  } else {
    sch400HO <- fread(file.path(tsv_dir, 'Sch400_HO_combined.tsv'), sep = ' ',  col.names = c('i', 'j', 'k', 'ROI'))
  }
  if(length(roi > 1)){
    sch400HO <- sch400HO[ROI %in% roi, c('i', 'j', 'k', 'ROI')]
    sch400HO[, v := 1:.N, by = 'ROI']
    outcols <- c('v','y','ROI')
  } else {
    sch400HO <- sch400HO[ROI == roi, c('i', 'j', 'k')]
    sch400HO[, v := 1:.N]  
    outcols <- c('v','y')
  }
  
  cl <- parallel::makePSOCKcluster(ncpus)
  parallel::clusterExport(cl = cl, varlist = c('fregex', 'sch400HO'), envir = rlang::env(fregex = fregex, sch400HO = sch400HO))
  nada <- parallel::clusterEvalQ(cl, {library(data.table); setDTthreads(1)})
  f_split <- split(f, 1:ncpus)
  sea_data_list <- parallel::parLapply(cl = cl, X = f_split, fun = function(f_chunk){
    #f_chunk=f_split[[1]]
    data_list <- lapply(1:length(f_chunk), function(i){
      #i=1
      d <- fread(f_chunk[[i]], sep = ' ',  col.names = c('i', 'j', 'k', 'y'))[sch400HO, on = c('i', 'j', 'k')][, ..outcols]
      id <- gsub(fregex, '\\1', f_chunk[[i]])
      sess <- gsub(fregex, '\\2', f_chunk[[i]])
      cope <- gsub(fregex, '\\3', f_chunk[[i]])
      d[, c('id', 'sess', 'cope') := list(id, sess, cope)]
      return(d)
    })
    
    r <- rbindlist(data_list)
    return(r)
  })
  parallel::stopCluster(cl)
  
  sea_data_dt <- rbindlist(sea_data_list)
}

make_size_matched_control_regions <- function(tsv_dir, max.iter = 50, tol = .95){
  require(data.table)
  radius_for_vol <- Vectorize(function(vol){
    r = (vol * 3 / (4*pi))^(1/3)
    return(r)
  })
  equivalent_region <- function(x, map, max.iter = 50, tol = .95){
    #x <- subcort_sizes[1,]
    #map <- sch400HO_control_regions
    iter = 0
    N = 0
    target_N <- x[['N']]
    done <- FALSE
    
    x[, radius := radius_for_vol(N)]
    random_index <- sample(1:dim(map)[[1]], max.iter)
    
    while (!done){
      random_coord <- map[random_index[[iter + 1]]][, c('i', 'j', 'k')]
      map[, c('dist_i', 'dist_j', 'dist_k') := .(i - random_coord[['i']],
                                                 j - random_coord[['j']],
                                                 k - random_coord[['k']])]
      map[, dist := sqrt(dist_i^2 + dist_j^2 + dist_k^2)]
      map[, in_sphere := dist < x[['radius']]]
      region <- map[in_sphere == TRUE]
      if(region[, .N] < tol * target_N) {
        region <- NULL
        iter <- iter + 1
      }
      if(iter >= max.iter | (!is.null(region) && region[, .N] >= tol * target_N)){
        done <- TRUE
      }
    }
    return(region)
  }
  
  control_rois <-  c(401, 403, 412, 414) #LR White matter and Ventricle
  sch400HO <- fread(file.path(tsv_dir, 'Sch400_HO_combined.tsv'), sep = ' ',  col.names = c('i', 'j', 'k', 'ROI'))
  setkey(sch400HO, 'ROI')
  sch400HO_sizes <- sch400HO[, .N, by = 'ROI']
  control_avg_voxels <- sch400HO_sizes[ROI %in% control_rois, list(mean = mean(N), sd = sd(N))]
  subcort_sizes <- sch400HO_sizes[ROI %in% c(404:407, 409:411, 415:421)]
  cort_avg_voxels <- sch400HO_sizes[ROI <= 400, list(mean = mean(N), sd = sd(N))]
  cort_size_sample <- sch400HO_sizes[ROI %in% sample(1:400, size = 10)]
  sch400HO_control_regions <- sch400HO[ROI %in% control_rois]
  eq_region_target <- rbindlist(list(subcort = subcort_sizes, cort = cort_size_sample), idcol = 'anatomy')
  eq_region_list <- lapply(1:dim(eq_region_target)[[1]], \(idx){
    target_region <- eq_region_target[idx,]
    eq_region <- equivalent_region(x = target_region, map = sch400HO_control_regions, max.iter = max.iter, tol = tol)
    setnames(target_region, c('ROI', 'N'), c('target_ROI', 'target_N'))
    return(cbind(eq_region, target_region))
  })
  eq_regions <- rbindlist(eq_region_list, idcol = 'eq_ROI', fill = TRUE)
  return(eq_regions)
}

load_control_regions <- function(tsv_dir, filename = 'control_region_map.rds', max.iter = 1000, tol = .95, replace = FALSE){
  if(!file.exists(file.path(tsv_dir, filename)) | replace == TRUE){
    message('File does not exist, creating new file: ', file.path(tsv_dir, filename))
    control_regions <- make_size_matched_control_regions(tsv_dir = tsv_dir, max.iter = max.iter, tol = tol)
    saveRDS(control_regions, file.path(tsv_dir, filename))
  } else {
    message('Loading: ', file.path(tsv_dir, filename))
    control_regions <- readRDS(file.path(tsv_dir, filename))
  }
  return(control_regions)
}

# cr <- load_control_regions(tsv_dir = 'tsv_nosmooth', max.iter = 3000, tol = .90, replace = FALSE)


