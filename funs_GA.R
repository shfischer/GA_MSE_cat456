### ------------------------------------------------------------------------ ###
### objective function for multi species run ####
### ------------------------------------------------------------------------ ###
mp_fitness <- function(params, inp_file, path, check_file = FALSE,
                       summarise_runs = FALSE,
                       MP, 
                       return_res = FALSE,
                       collapse_correction = TRUE,
                       obj_SSB = FALSE, obj_C = FALSE, obj_F = FALSE,
                       obj_risk = FALSE, obj_ICV = FALSE, obj_ICES_PA = FALSE,
                       obj_ICES_PA2 = FALSE, obj_ICES_MSYPA = FALSE,
                       stat_yrs = "all",
                       risk_threshold = 0.05,
                       ...) {
  
  ### housekeeping
  invisible(gc())
  if (exists("res_mp")) {
    rm(res_mp)
    invisible(gc())
  }
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  
  ### rounding of arguments
  if (identical(MP, "rfb")) {
    params[1:4] <- round(params[1:4])
    params[5:7] <- round(params[5:7], 1)
    params[8] <- round(params[8])
    params[9] <- round(params[9], 2)
    params[10:11] <- round(params[10:11], 2)
    ### fix NaN for upper_constraint
    if (is.nan(params[10])) params[10] <- Inf
  } else if (identical(MP, "hr")) {
    ### idxB_lag, idxB_range_3, interval [years]
    params[c(1, 2, 5)] <- round(params[c(1, 2, 5)])
    ### exp_b, comp_b_multiplier
    params[c(3, 4)] <- round(params[c(3, 4)], 1)
    ### multiplier, upper_constraint, lower_constraint
    params[c(6, 7, 8)] <- round(params[c(6, 7, 8)], 2)
    ### fix NaN for upper_constraint
    if (is.nan(params[7])) params[7] <- Inf
  } else if (identical(MP, "CC_f")) {
    params[1] <- round(params[1], 2) ### Lref_mult
    params[2] <- round(params[2], 2) ### pa_size
    params[3] <- round(params[3]) ### interval
    params[4] <- round(params[4]) ### multiplier
    params[c(5, 6)] <- round(params[c(5, 6)], 2) ### upper/lower constr
    ### fix NaN for upper_constraint
    if (is.nan(params[5])) params[5] <- Inf
  } else if (identical(MP, "CL")) {
    params[1] <- round(params[1], 0) ### interval
    params[c(2, 3)] <- round(params[c(2, 3)], 2) ### lambda lower/upper
    params[c(4, 5)] <- round(params[c(4, 5)], 2) ### lambda lower/upper
    params[c(6, 7)] <- round(params[c(6, 7)], 2) ### r/l treshold
    params[8] <- round(params[8], 2) ### Lref_mult
    params[9] <- round(params[9], 2) ### multiplier
    params[10:11] <- round(params[10:11], 2) ### catch limits
  }
  
  ### check for files?
  run_mp <- TRUE ### initialise variable
  if (isTRUE(check_file)) {
    ### current run
    run_i <- paste0(params, collapse = "_")
    ### get current stock(s)
    stock_i <- strsplit(x = tail(strsplit(x = path, split = "/")[[1]], 1), 
                        split = "_")[[1]]
    base_path <- paste0(paste0(head(strsplit(x = path, split = "/")[[1]], -1), 
                               collapse = "/"), "/")
    ### check if path exists
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    
    ### check if summary file is used
    if (isTRUE(summarise_runs) & 
        isTRUE(file.exists(paste0(path, "runs.rds")))) {
      ### load all available runs
      stats_all <- readRDS(paste0(path, "runs.rds"))
      ### check if current exists and use results
      if (isTRUE(run_i %in% names(stats_all))) {
        stats <- stats_all[[run_i]]$stats
        run_mp <- FALSE
      } else {
        run_mp <- TRUE
      }
    }
    
    ### check if run already exists in individual file (legacy approach)
    if (isTRUE(file.exists(paste0(path, run_i, ".rds"))) &
        isTRUE(run_mp)) {
      ### load stats
      stats <- readRDS(paste0(path, run_i, ".rds"))
      ### set flag for running MP
      run_mp <- FALSE
      ### use different period to calculate stats?
      if (!any(stat_yrs %in% c("all", "more"))) {
        if (!any(grepl(x = rownames(stats), pattern = stat_yrs))) run_mp <- TRUE
      }
    } else if (isTRUE(run_mp)) {
      ### check if run exist in larger group
      dir_i <- paste0(stock_i, collapse = "_")
      dirs_i <- setdiff(x = dir(path = base_path, pattern = dir_i),
                        y = dir_i)
      if (isTRUE(length(dirs_i) > 0)) {
        dirs_i <- dirs_i[which(sapply(dirs_i, function(x) {
          tmp <- strsplit(x = x, split = "_")[[1]]
          ifelse(isFALSE(dir_i %in% tmp), FALSE, TRUE)
        }))]
        files_tmp <- lapply(dirs_i, function(x) {
          #browser()
          path_tmp <- paste0(base_path, x, "/", run_i, ".rds")
          if (isTRUE(file.exists(path = path_tmp))) {
            return(path_tmp)
          } else {
            return(NA)
          }
        })
        files_tmp[is.na(files_tmp)] <- NULL
        if (isTRUE(length(files_tmp) > 0)) {
          ### load stats from larger group
          stats <- readRDS(files_tmp[[1]])
          ### subset to current group
          stats <- stats[, stock_i]
          ### do not run MP
          run_mp <- FALSE
        } 
      } 
    }
  }
  
  if (isTRUE(run_mp)) {
    
    ### load input file from disk
    input <- readRDS(inp_file)
    
    ### insert arguments into input object for mp
    ### rfb rule
    if (identical(MP, "rfb")) {
      input <- lapply(input, function(x) {
        x$ctrl$est@args$idxB_lag     <- params[1]
        x$ctrl$est@args$idxB_range_1 <- params[2]
        x$ctrl$est@args$idxB_range_2 <- params[3]
        x$ctrl$est@args$catch_range  <- params[4]
        x$ctrl$est@args$comp_m <- params[9]
        x$ctrl$phcr@args$exp_r <- params[5]
        x$ctrl$phcr@args$exp_f <- params[6]
        x$ctrl$phcr@args$exp_b <- params[7]
        x$ctrl$hcr@args$interval <- params[8]
        x$ctrl$isys@args$interval <- params[8]
        x$ctrl$isys@args$upper_constraint <- params[10]
        x$ctrl$isys@args$lower_constraint <- params[11]
        return(x)
      })
    ### harvest rates
    } else if (identical(MP, "hr")) {
      input <- lapply(input, function(x) {
        ### biomass index 
        x$ctrl$est@args$idxB_lag <- params[1]
        x$ctrl$est@args$idxB_range_3 <- params[2]
        ### biomass safeguard
        x$ctrl$phcr@args$exp_b <- params[3]
        ### change Itrigger? (default: Itrigger=1.4*Iloss)
        if (isFALSE(params[4] == 1.4)) {
          x$ctrl$est@args$I_trigger <- x$ctrl$est@args$I_trigger/1.4*params[4]
        }
        ### multiplier
        x$ctrl$est@args$comp_m <- params[6]
        ### catch interval (default: 1)
        if (is.numeric(params[5])) {
          x$ctrl$hcr@args$interval <- params[5]
          x$ctrl$isys@args$interval <- params[5]
        }
        ### catch constraint
        x$ctrl$isys@args$upper_constraint <- params[7]
        x$ctrl$isys@args$lower_constraint <- params[8]
        #x$ctrl$isys@args$cap_below_b <- params[]
        return(x)
      })
    ### constant catch with conditional PA buffer
    } else if (identical(MP, "CC_f")) {
      input <- lapply(input, function(x) {
        
        ### Lref multiplier
        x$ctrl$est@args$Lref_mult <- params[1]
        ### PA buffer size & duration
        x$ctrl$est@args$pa_size <- params[2]
        x$ctrl$est@args$pa_duration <- params[3]
        ### multiplier
        x$ctrl$est@args$comp_m <- params[4]
        ### catch interval
        if (is.numeric(params[3])) {
          x$ctrl$hcr@args$interval <- params[3]
          x$ctrl$isys@args$interval <- params[3]
        }
        ### catch constraint
        x$ctrl$isys@args$upper_constraint <- params[5]
        x$ctrl$isys@args$lower_constraint <- params[6]
        return(x)
      })
    } else if (identical(MP, "CL")) {
      input <- lapply(input, function(x) {
        ### interval
        x$ctrl$est@args$interval <- params[1]
        x$ctrl$phcr@args$interval <- params[1]
        x$ctrl$hcr@args$interval <- params[1]
        x$ctrl$isys@args$interval <- params[1]
        x$ctrl$est@args$first_catch  <- params[10]
        x$ctrl$est@args$catch_limit  <- params[11]
        ### lambda
        x$ctrl$phcr@args$lambda_lower <- params[2]
        x$ctrl$phcr@args$lambda_upper <- params[3]
        ### gamma
        x$ctrl$phcr@args$gamma_lower <- params[4]
        x$ctrl$phcr@args$gamma_upper <- params[5]
        ### r/l thresholds
        x$ctrl$phcr@args$r_threshold <- params[6]
        x$ctrl$phcr@args$l_threshold <- params[7]
        ### Lref_mult
        x$ctrl$phcr@args$Lref_mult <- params[8]
        ### multiplier
        x$ctrl$phcr@args$multiplier <- params[9]
        ### catch limit (absolute relative to past)
        x$ctrl$isys@args$catch_limit  <- params[11]
        return(x)
      })
    }
    
    ### if group of stocks, check if results for individual stocks exist
    group <- ifelse(isTRUE(length(input) > 1) & isTRUE(check_file), TRUE, FALSE)
    if (group) {
      ### get paths
      group_stocks <- names(input)
      path_base <- gsub(x = path, 
                        pattern = paste0(paste0(group_stocks, collapse = "_"), 
                                         "/"),
                        replacement = "")
      path_stocks <- paste0(path_base, group_stocks, "/")
      ### check for files
      run_exists <- file.exists(paste0(path_stocks, run_i, ".rds"))
      group <- ifelse(any(run_exists), TRUE, FALSE)
      
      ### do some results exist?
      if (group) {
        ### load results
        files_exist <- paste0(path_stocks, run_i, ".rds")[run_exists]
        stats_group <- lapply(files_exist, readRDS)
        names(stats_group) <- group_stocks[run_exists]
        ### get stocks which require simulation
        run_stocks <- group_stocks[!run_exists]
        ### subset input
        input <- input[run_stocks]
        
      }
      
    }
    
    ### run MP for each list element
    res_mp <- lapply(input, function(x) {
      if (getDoParWorkers() > 1)
        . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
      do.call(mp, x)
    })
    
    if (isTRUE(return_res)) {
      return(res_mp)
    }
    
    ### calculate stats
    stat_yrs_calc <- "more"
    stats <- mp_stats(input = input, res_mp = res_mp, stat_yrs = stat_yrs_calc,
                      collapse_correction = collapse_correction)
    
    ### add existing results for stock groups
    if (group) {
      
      ### split old stats into list
      if (isTRUE(length(stats) > 0)) {
        stats <- asplit(stats, MARGIN = 2)
      }
      ### stats_group is already a list
      ### combine new and existing stats
      stats <- c(stats_group, stats)
      ### sort and coerce into matrix
      stats <- stats[group_stocks]
      stats <- do.call(cbind, stats)
      
    }
    
    ### save result in file
    if (isTRUE(check_file)) {
      saveRDS(stats, paste0(path, run_i, ".rds"))
    }
    
  }
  
  ### prepare stats for objective function
  if (identical(stat_yrs, "all") | identical(stat_yrs, "more")) {
    SSB_rel <- stats["SSB_rel", ]
    Catch_rel <- stats["Catch_rel", ]
    Fbar_rel <- stats["Fbar_rel", ]
    risk_Blim <- stats["risk_Blim", ]
    ICV <- stats["ICV", ]
  } else if (identical(stat_yrs, "last10")) {
    SSB_rel <- stats["SSB_rel_last10", ]
    Catch_rel <- stats["Catch_rel_last10", ]
    Fbar_rel <- stats["Fbar_rel_last10", ]
    risk_Blim <- stats["risk_Blim_last10", ]
    ICV <- stats["ICV_last10", ]
  }
  ### objective function
  obj <- 0
  ### MSY objectives: target MSY reference values
  if (isTRUE(obj_SSB)) obj <- obj - sum(abs(unlist(SSB_rel) - 1))
  if (isTRUE(obj_C)) obj <- obj - sum(abs(unlist(Catch_rel) - 1))
  if (isTRUE(obj_F)) obj <- obj - sum(abs(unlist(Fbar_rel) - 1))
  ### reduce risk & ICV
  if (isTRUE(obj_risk)) obj <- obj - sum(unlist(risk_Blim))
  if (isTRUE(obj_ICV)) obj <- obj - sum(unlist(ICV))
  ### ICES approach: maximise catch while keeping risk <5%
  if (isTRUE(obj_ICES_PA)) {
    obj <- obj + sum(unlist(Catch_rel))
    ### penalise risk above 5%
    obj <- obj - sum(ifelse(test = unlist(risk_Blim) <= 0.05,
                            yes = 0,
                            no = 10)) 
  }
  if (isTRUE(obj_ICES_PA2)) {
    obj <- obj + sum(unlist(Catch_rel))
    ### penalise risk above 5% - gradual
    obj <- obj - sum(penalty(x = unlist(risk_Blim), 
                             negative = FALSE, max = 1, inflection = 0.06, 
                             steepness = 0.5e+3))
  }
  ### MSY target but replace risk with PA objective
  if (isTRUE(obj_ICES_MSYPA)) {
    obj <- obj - sum(abs(unlist(SSB_rel) - 1)) -
      sum(abs(unlist(Catch_rel) - 1)) -
      sum(unlist(ICV)) -
      sum(penalty(x = unlist(risk_Blim), 
                             negative = FALSE, max = 5, 
                             inflection = risk_threshold + 0.01, 
                             steepness = 0.5e+3))
      ### max penalty: 5
      ### for pollack zero catch has fitness of -4.7
  }
  
  ### housekeeping
  suppressWarnings(rm(res_mp, input))
  invisible(gc())
  if (getDoParWorkers() > 1)
    . <- foreach(i = 1:getDoParWorkers()) %dopar% {invisible(gc())}
  
  ### return objective function (fitness) value
  return(obj)
  
}

### ------------------------------------------------------------------------ ###
### stats from MSE run(s) ####
### ------------------------------------------------------------------------ ###

### function for calculating stats
mp_stats <- function(input, res_mp, stat_yrs = "all", 
                     collapse_correction = TRUE) {
  
  mapply(function(input_i, res_mp_i) {
    
    ### stock metrics
    SSBs <- FLCore::window(ssb(res_mp_i@om@stock), start = 101)
    Fs <- FLCore::window(fbar(res_mp_i@om@stock), start = 101)
    Cs <- FLCore::window(catch(res_mp_i@om@stock), start = 101)
    yrs <- dim(SSBs)[2]
    its <- dim(SSBs)[6]
    ### collapse correction
    if (isTRUE(collapse_correction)) {
      ### find collapses
      cd <- sapply(seq(its), function(x) {
        min_yr <- min(which(SSBs[,,,,, x] < 1))
        if (is.finite(min_yr)) {
          all_yrs <- min_yr:yrs
        } else {
          all_yrs <- NA
        }
        all_yrs + (x - 1)*yrs
      })
      cd <- unlist(cd)
      cd <- cd[which(!is.na(cd))]
      ### remove values
      SSBs@.Data[cd] <- 0
      Cs@.Data[cd] <- 0
      Fs@.Data[cd] <- 0
    }
    ### extend Catch to include ICV calculation from last historical year
    Cs_long <- FLCore::window(Cs, start = 100)
    Cs_long[, ac(100)] <- catch(res_mp_i@om@stock)[, ac(100)]
    ### refpts
    Bmsy <- c(attr(input_i, "refpts")["msy", "ssb"])
    Fmsy <- c(attr(input_i, "refpts")["msy", "harvest"])
    Cmsy <- c(attr(input_i, "refpts")["msy", "yield"])
    Blim <- attr(input_i, "Blim")
    ### TAC interval
    TAC_intvl <- input_i$ctrl$hcr@args$interval
    
    ### some stats
    stats_list <- function(SSBs, Cs, Fs, Cs_long, Blim, Bmsy, Fmsy, Cmsy,
                           TAC_intvl) {
      list(
        risk_Blim = mean(c(SSBs < Blim), na.rm = TRUE),
        risk_Bmsy = mean(c(SSBs < Bmsy), na.rm = TRUE),
        risk_halfBmsy = mean(c(SSBs < Bmsy/2), na.rm = TRUE),
        risk_collapse = mean(c(SSBs < 1), na.rm = TRUE),
        SSB = median(c(SSBs), na.rm = TRUE), Fbar = median(c(Fs), na.rm = TRUE),
        Catch = median(c(Cs), na.rm = TRUE),
        SSB_rel = median(c(SSBs/Bmsy), na.rm = TRUE),
        Fbar_rel = median(c(Fs/Fmsy), na.rm = TRUE),
        Catch_rel = median(c(Cs/Cmsy), na.rm = TRUE),
        ICV = iav(Cs_long, from = 100, period = TAC_intvl,
                  summary_all = median)
      )
    }
    stats_i <- stats_list(SSBs = SSBs, Cs = Cs, Fs = Fs, 
                          Cs_long = Cs_long, 
                          Blim = Blim, Bmsy = Bmsy, Fmsy = Fmsy, Cmsy = Cmsy,
                          TAC_intvl = TAC_intvl)
    ### additional time period?
    if (identical(stat_yrs, "last10")) {
      yrs10 <- tail(dimnames(SSBs)$year, 10)
      yrs10p1 <- tail(dimnames(SSBs)$year, 11)
      stats_i_last10 <- c(stats_list(SSBs = SSBs[, yrs10], Cs = Cs[, yrs10],
                                     Fs = Fs[, yrs10], Cs_long = Cs[, yrs10p1], 
                                     Blim = Blim, Bmsy = Bmsy, Fmsy = Fmsy, 
                                     Cmsy = Cmsy, TAC_intvl = TAC_intvl))
      names(stats_i_last10) <- paste0(names(stats_i_last10), "_last10")
      stats_i <- c(stats_i, stats_i_last10)
    } else if (identical(stat_yrs, "more")) {
      yrs_for_stats <- c("first10", "41to50", "last10", "firsthalf",
                         "lastfhalf", "11to50")
      stats_add <- lapply(yrs_for_stats, function(x) {
        ### define years for summary statistics
        yrs_tmp <- switch(x,
                         "first10" = head(dimnames(SSBs)$year, 10), 
                          "41to50" = ac(141:150), 
                          "last10" = tail(dimnames(SSBs)$year, 10), 
                          "firsthalf" = head(dimnames(SSBs)$year, 
                                             length(dimnames(SSBs)$year)/2),
                          "lastfhalf" = tail(dimnames(SSBs)$year, 
                                             length(dimnames(SSBs)$year)/2), 
                          "11to50" = ac(111:150))
        if (!any(yrs_tmp %in% dimnames(SSBs)$year)) return(NULL)
        yrs_tmpp1 <- ac(seq(from = min(as.numeric(yrs_tmp)) - 1, 
                            to = max(as.numeric(yrs_tmp))))
        stats_tmp <- c(stats_list(SSBs = SSBs[, yrs_tmp], Cs = Cs[, yrs_tmp],
                                  Fs = Fs[, yrs_tmp], Cs_long = Cs_long[, yrs_tmpp1],
                                  Blim = Blim, Bmsy = Bmsy, Fmsy = Fmsy, 
                                  Cmsy = Cmsy, TAC_intvl = TAC_intvl))
        names(stats_tmp) <- paste0(names(stats_tmp), "_", x)
        return(stats_tmp)
      })
      stats_i <- c(stats_i, unlist(stats_add))
    }
    
    return(stats_i)
  }, input, res_mp)
  
}

### ------------------------------------------------------------------------ ###
### penalty function ####
### ------------------------------------------------------------------------ ###

penalty <- function(x, negative = FALSE, max = 1,
                    inflection = 0.06, steepness = 0.5e+3) {
  y <- max / (1 + exp(-(x - inflection)*steepness))
  if (isTRUE(negative)) y <- -y
  return(y)
}

### ------------------------------------------------------------------------ ###
### postFitness - collate results and delete files from individual runs  ####
### ------------------------------------------------------------------------ ###
mp_postFitness <- function(x, path, check_file, MP, summarise_runs = FALSE,
                           ...) {
  
    if (isTRUE(summarise_runs)) {
      
    ### population - MP parameters
    pop <- as.data.frame(x@population)
    
    ### round parameters
    if (identical(MP, "hr")) {
      colnames(pop) <- c("idxB_lag", "idxB_range_3", "exp_b", "comp_b_multiplier",
                         "interval", "multiplier",
                         "upper_constraint", "lower_constraint")
      ### idxB_lag, idxB_range_3, interval [years]
      pop$idxB_lag <- round(pop$idxB_lag)
      pop$idxB_range_3 <- round(pop$idxB_range_3)
      pop$interval <- round(pop$interval)
      ### exp_b, comp_b_multiplier
      pop$exp_b <- round(pop$exp_b, 1)
      pop$comp_b_multiplier <- round(pop$comp_b_multiplier, 1)
      ### multiplier, upper_constraint, lower_constraint
      pop$multiplier <- round(pop$multiplier, 2)
      pop$upper_constraint <- round(pop$upper_constraint, 2)
      pop$lower_constraint <- round(pop$lower_constraint, 2)
      ### fix NaN for upper_constraint
      if (any(is.nan(pop$upper_constraint))) {
        pop$upper_constraint[is.nan(pop$upper_constraint)] <- Inf
      }
    } else if (identical(MP, "CL")) {
      colnames(pop) <- c("interval", "lambda_lower", "lambda_upper", 
                         "gamma_lower", "gamma_upper", "r_threshold",
                         "l_threshold", "Lref_mult", "multiplier")
      pop$interval <- round(pop$interval, 0) ### interval
      pop$lambda_lower <- round(pop$lambda_lower, 2) ### lambda lower/upper
      pop$lambda_upper <- round(pop$lambda_upper, 2) ### lambda lower/upper
      pop$gamma_lower <- round(pop$gamma_lower, 2) ### gamma lower/upper
      pop$gamma_upper <- round(pop$gamma_upper, 2) ### gamma lower/upper
      pop$r_threshold <- round(pop$r_threshold, 2) ### r_threshold
      pop$l_threshold <- round(pop$l_threshold, 2) ### l_threshold
      pop$Lref_mult <- round(pop$Lref_mult, 2) ### Lref_mult
      pop$multiplier <- round(pop$multiplier, 2) ### multiplier
    }
    ### file names
    pop$name <- apply(pop, 1, paste0, collapse = "_")
    ### remove duplicates
    pop <- unique(pop)
    
    ### load summary statistics for each run
    stats_new <- lapply(seq(nrow(pop)), function(y) {
      ### get MP parameters
      pars_y <- unlist(pop[y, -ncol(pop)])
      ### check if file exists
      if (file.exists(paste0(path, pop$name[y], ".rds"))) {
        ### load summary statistics
        stats_y <- readRDS(paste0(path, pop$name[y], ".rds"))
        ### combine and return
        list(pars = pars_y, stats = stats_y)
      } else {
        return(NULL)
      }
    })
    names(stats_new) <- pop$name
    
    ### check if summary file exists
    if (!file.exists(paste0(path, "runs.rds"))) {
      ### save file, in first generation of GA search
      saveRDS(stats_new, file = paste0(path, "runs.rds"))
    } else {
      ### load file
      stats <- readRDS(paste0(path, "runs.rds"))
      ### add new runs
      stats_add_names <- setdiff(names(stats_new), names(stats))
      if (isTRUE(length(stats_add_names) > 0)) {
        stats <- append(stats,
                        stats_new[names(stats_new) %in% stats_add_names])
        stats <- stats[unique(names(stats))]
        ### save updated file
        saveRDS(stats, file = paste0(path, "runs.rds"))
      }
  
    }
    
    ### delete individual files
    . <- suppressWarnings(file.remove(paste0(path, pop$name, ".rds")))
    
  }
  
  return(x)
  
}

