
### ------------------------------------------------------------------------ ###
### run MSE ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### arguments ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
if (exists(x = "args_local")) args <- append(args, args_local)
print("arguments passed on to this script:")
print(args)

### evaluate arguments, if they are passed to R:
if (length(args) > 0) {
  
  ### extract arguments
for (i in seq_along(args)) eval(parse(text = args[[i]]))
  ### set default arguments
  ### parallelization
  if (!exists("use_MPI")) use_MPI <- FALSE
  if (!exists("n_blocks")) n_blocks <- 1
  if (!exists("n_workers")) n_workers <- 0
  
  ### split OM into blocks?
  if (!exists("n_parts")) n_parts <- 1
  if (!exists("part")) part <- 1
  
  ### projection details
  if (!exists("n_iter")) n_iter <- 500
  if (!exists("n_yrs")) n_yrs <- 50
  if (!exists("fhist")) fhist <- "one-way"
  
  ### MP parameters
  if (!exists("MP")) MP <- "CL"
  ### MP - constant catch
  if (identical(MP, "constant_catch")) {
    if (!exists("multiplier")) multiplier <- 1
    if (!exists("comp_r")) comp_r <- FALSE
    if (!exists("comp_f")) comp_f <- FALSE
    if (!exists("comp_b")) comp_b <- FALSE
    if (!exists("pa_buffer")) pa_buffer <- TRUE
    if (!exists("pa_size")) pa_size <- 0.8
    if (!exists("pa_duration")) pa_duration <- 3
    if (!exists("interval")) interval <- 3
    if (!exists("idxB_lag")) idxB_lag <- 1
    if (!exists("upper_constraint")) upper_constraint <- Inf
    if (!exists("lower_constraint")) lower_constraint <- 0
  }
  ### MP - CC_f
  
  ### MP - CL
  if (identical(MP, "CL")) {
    if (!exists("interval")) interval <- 3
    if (!exists("lambda_upper")) lambda_upper <- 0.1
    if (!exists("lambda_lower")) lambda_lower <- 0.2
    if (!exists("gamma_lower")) gamma_lower <- 0.2
    if (!exists("gamma_upper")) gamma_upper <- 0.1
    if (!exists("r_threshold")) r_threshold <- 0.05
    if (!exists("l_threshold")) l_threshold <- 0.1
    if (!exists("Lref_mult")) Lref_mult <- 1
    if (!exists("multiplier")) multiplier <- 1
  }
  
  if (!exists("stat_yrs")) stat_yrs <- "all"
  if (!exists("scenario")) scenario <- "baseline"
  
  ### OM specifications
  ### observation uncertainty
  if (!exists("sigmaL")) sigmaL <- 0.1
  if (!exists("sigmaL_rho")) sigmaL_rho <- 0
  if (!exists("sigmaC")) sigmaC <- 0.1
  ### implementation error
  if (!exists("sigmaIEM")) sigmaIEM <- 0.1
  ### recruitment variability
  if (!exists("sigmaR")) sigmaR <- 0.6
  if (!exists("sigmaR_rho")) sigmaR_rho <- 0.0
  ### recruitment steepness
  if (!exists("steepness")) steepness <- 0.75
  
  ### what to save
  if (!exists("check_file")) check_file <- TRUE
  if (!exists("saveMP")) saveMP <- TRUE
  if (!exists("stats")) stats <- TRUE
  if (!exists("collate")) collate <- FALSE
  
  ### GA search
  if (!exists("ga_search")) ga_search <- FALSE
  if (isTRUE(ga_search)) {
    if (!exists("popSize")) stop("popSize missing")
    if (!exists("maxiter")) stop("maxiter missing")
    if (!exists("stock_id")) stop("stock_id missing")
    if (!exists("run")) run <- maxiter
    if (!exists("collate")) collate <- FALSE
    ### objective function elements
    if (!exists("obj_SSB")) obj_SSB <- FALSE
    if (!exists("obj_F")) obj_F <- FALSE
    if (!exists("obj_C")) obj_C <- FALSE
    if (!exists("obj_risk")) obj_risk <- FALSE
    if (!exists("obj_ICV")) obj_ICV <- FALSE
    if (!exists("obj_ICES_PA")) obj_ICES_PA <- FALSE
    if (!exists("obj_ICES_PA2")) obj_ICES_PA2 <- FALSE
    if (!exists("obj_ICES_MSYPA")) obj_ICES_MSYPA <- TRUE
    if (!exists("risk_threshold")) risk_threshold <- 0.05
    ### GA
    if (!exists("add_suggestions")) add_suggestions <- FALSE
    if (!exists("stat_yrs")) stat_yrs <- "all"
  }

} else {
  
  stop("no argument passed to R")

}

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
### GA fork from GitHub remotes::install_github("shfischer/GA")
### use mse fork from shfischer/mse, branch mseDL2.0 
### remotes::install_github("shfischer/mse", ref = "mseDL2.0)
req_pckgs <- c("mse", "tidyr", "dplyr", "doParallel", "GA", "doRNG")
for (i in req_pckgs) 
  suppressPackageStartupMessages(library(package = i, character.only = TRUE))

### load additional functions
source("funs.R")
source("funs_GA.R")
source("funs_OM.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###

### hybrid MPI
if (isTRUE(use_MPI)) {
  ### 1st: doMPI cluster with 1 worker per node
  message("starting doMPI")
  library(doMPI)
  cl1 <- startMPIcluster()
  message("startMPIcluster() succeeded")
  print(cl1)
  registerDoMPI(cl1)
  cl_length_1 <- cl1$workerCount
  cl_length_1
  
  ### 2nd: doParallel workers inside doMPI workers
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    ### load packages and functions into MPI workers
    for (i in req_pckgs) 
      suppressPackageStartupMessages(library(package = i, character.only = TRUE,
                                             warn.conflicts = FALSE, 
                                             verbose = FALSE, quietly = TRUE))
  }
  message("MPI package loading succeeded")
  . <- foreach(i = seq(cl_length_1)) %dopar% {
    source("funs.R", echo = FALSE)
    source("funs_GA.R", echo = FALSE)
  }
  message("MPI script loading succeeded")
  ### start doParallel inside MPI processes
  if (isTRUE(n_workers > 1)) {
    . <- foreach(i = seq(cl_length_1)) %dopar% {
      cl2 <- makeCluster(n_workers)
      registerDoParallel(cl2)
      cl_length_2 <- length(cl2)
      ### load packages and functions into parallel workers
      . <- foreach(i = seq(cl_length_2)) %dopar% {
        for (i in req_pckgs) 
          suppressPackageStartupMessages(library(package = i, 
                                                 character.only = TRUE,
                                                 warn.conflicts = FALSE, 
                                                 verbose = FALSE, 
                                                 quietly = TRUE))
        source("funs.R", echo = FALSE)
        source("funs_GA.R", echo = FALSE)
      }
    }
  }
  message("setting up doParallel inside MPI succeeded")
} else {
  if (isTRUE(n_workers > 1)) {
    ### start doParallel cluster
    cl1 <- makeCluster(n_workers)
    registerDoParallel(cl1)
    cl_length_1 <- length(cl1)
    ### load packages and functions into parallel workers
    . <- foreach(i = seq(cl_length_1)) %dopar% {
      for (i in req_pckgs) 
        suppressPackageStartupMessages(library(package = i, 
                                               character.only = TRUE,
                                               warn.conflicts = FALSE, 
                                               verbose = FALSE, quietly = TRUE))
      source("funs.R", echo = FALSE)
      source("funs_GA.R", echo = FALSE)
    }
  } else {
    cl1 <- FALSE
  }
}

### ------------------------------------------------------------------------ ###
### MP parameters ####
### ------------------------------------------------------------------------ ###

### HR rule parameters & uncertainty
if (isFALSE(exists("stock_id"))) stop("'stock_id' is missing!")
hr_params <- c(
  ### OM
  "stock_id", "fhist", "n_iter", "n_yrs", "MP", "scenario", "n_blocks",
  ### uncertainty
  "sigmaL", "sigmaL_rho", "sigmaC", "sigmaR", "sigmaR_rho", "steepness",
  "sigmaIEM",
  ### MP parameters
  "multiplier", "interval", "pa_buffer", "pa_size", "pa_duration",
  "upper_constraint", "lower_constraint",
  "lambda_upper", "lambda_lower", "gamma_lower", "gamma_upper", "r_threshold",
  "l_threshold", "Lref_mult"
)
hr_params <- as.data.frame(mget(hr_params, ifnotfound = NA))
names(hr_params)[1] <- "stocks"

### ------------------------------------------------------------------------ ###
### manual runs ####
### ------------------------------------------------------------------------ ###
if (isFALSE(ga_search)) {

  ### ---------------------------------------------------------------------- ###
  ### go through runs ####
  ### ---------------------------------------------------------------------- ###
  
  if (isTRUE(n_workers > 1 & n_blocks == 1)) {
    `%do_tmp%` <- `%dopar%`
  } else {
    `%do_tmp%` <- `%do%`
  }
  
  . <- foreach(hr_i = seq(nrow(hr_params))) %do_tmp% {
    
    par_i <- hr_params[hr_i, ]

    ### -------------------------------------------------------------------- ###
    ### generate MP input ####
    ### -------------------------------------------------------------------- ###
    
    input_i <- do.call(input_mp, as.list(par_i))
    
    ### -------------------------------------------------------------------- ###
    ### paths ####
    ### -------------------------------------------------------------------- ###
    ### generate file name
    pars_OM <- c(par_i$sigmaL, par_i$sigmaL_rho, par_i$sigmaC, par_i$sigmaR, 
                 par_i$sigmaR_rho, par_i$steepness, par_i$sigmaIEM)
    if (identical(MP, "CL")) 
      pars_MP <- c(par_i$interval, par_i$lambda_upper, par_i$lambda_lower, 
                   par_i$gamma_lower, par_i$gamma_upper, par_i$r_threshold, 
                   par_i$l_threshold, par_i$Lref_mult)
    file_out <- paste0(c(pars_OM, "", pars_MP), collapse = "_")
    path_out <- paste0("output/", MP, "/", n_iter, "_", n_yrs, "/", scenario, "/",
                       fhist, "/", paste0(names(input_i), collapse = "_"), "/")
    dir.create(path_out, recursive = TRUE)
    ### skip if run already exists
    if (file.exists(paste0(path_out, "stats_", file_out, ".rds"))) return(NULL)
    
    ### -------------------------------------------------------------------- ###
    ### run  ####
    ### -------------------------------------------------------------------- ###
    
    if (isTRUE(length(par_i$stock) > 1))
      stop("Individual MP runs only possible for one stock at a time!")
    res <- do.call(mp, input_i[[1]])

    ### -------------------------------------------------------------------- ###
    ### save ####
    ### -------------------------------------------------------------------- ###
    
    if (isTRUE(saveMP))
      saveRDS(object = res, file = paste0(path_out, "mp_", file_out, ".rds"))
    
    ### -------------------------------------------------------------------- ###
    ### stats ####
    ### -------------------------------------------------------------------- ###
    
    if (isTRUE(stats)) {
      res_stats <- mp_stats(input = input_i, res_mp = list(res), 
                            collapse_correction = TRUE,
                            stat_yrs = stat_yrs)
      res_stats <- cbind(stock = names(input_i), par_i, t(res_stats))
      saveRDS(object = res_stats, 
              file = paste0(path_out, "stats_", file_out, ".rds"))
    }
  
  }
  
  ### ---------------------------------------------------------------------- ###
  ### collate stats ####
  ### ---------------------------------------------------------------------- ###
  
  # if (isTRUE(stats) & isTRUE(collate) & isTRUE(nrow(hr_params) > 1)) {
  #   files <- paste0("stats_", hr, "_", 
  #                   sapply(seq(nrow(hr_params)), 
  #                          function(x) paste0(hr_params[x,], collapse = "_")),
  #                   ".rds")
  #   files <- paste0("output/hr/", n_iter, "_", n_yrs, "/", scenario, "/",
  #                   fhist, "/", paste0(stock, collapse = "_"), "/",
  #                   files)
  #   stats_all <- lapply(files, readRDS)
  #   stats_all <- do.call(rbind, stats_all)
  #   
  #   saveRDS(stats_all, file = paste0(
  #     "output/hr/", n_iter, "_", n_yrs, "/", scenario, "/", fhist, "/", 
  #     paste0(stock, collapse = "_"), "/",
  #     "collated_stats_", hr, "_", 
  #     paste0(apply(hr_params, 2, function(x) {
  #       ifelse(isTRUE(length(unique(x)) > 1), paste0(range(x), collapse = "-"),
  #              x[1])
  #     }), collapse = "_"), ".rds"))
  #   
  # }

### ------------------------------------------------------------------------ ###
### GA search ####
### ------------------------------------------------------------------------ ###

} else {
  
  ### ---------------------------------------------------------------------- ###
  ### prepare OM ####
  ### ---------------------------------------------------------------------- ###

  if (isTRUE(nrow(hr_params) > 1))
    stop("GA search only possible for one parameter (set)!")
  input <- do.call(input_mp, as.list(hr_params))

  ### ---------------------------------------------------------------------- ###
  ### GA set-up ####
  ### ---------------------------------------------------------------------- ###
  
  ### GA arguments
  ### MP: CC_f
  if (identical(MP, "CC_f")) {
    ga_names <- c("Lref_mult", "pa_size",
                  "interval", "multiplier",
                  "upper_constraint", "lower_constraint")
    ga_default <- c(1, 0.8, 3, 1, Inf, 0)
    ga_lower <-   c(0,   0, 1, 0,   0, 0)
    ga_upper <-   c(2,   1, 5, 2, Inf, 1)
    ga_suggestions <- rbind(c(1, 0.8, 3, 1, Inf, 0), ### default
                            c(0, 0.8, 3, 1, Inf, 0), ### zero catch
                            expand.grid(0:2, c(0.5, 0.7, 0.8, 0.9, 1), 1:5,
                                        0:2, c(1.2, Inf), c(0, 0.8)))
  } else if (identical(MP, "CL")) {
    ga_names <- c("interval",     ### 0 decimal digits; 1-5
                  "lambda_lower", ### 2 decimal digits; 0.00-0.50
                  "lambda_upper", ### 2 decimal digits; 0.00-0.50
                  "gamma_lower",  ### 2 decimal digits; 0.00-0.50
                  "gamma_upper",  ### 2 decimal digits; 0.00-0.50
                  "r_threshold",  ### 2 decimal digits; 0.00-0.50
                  "l_threshold",  ### 2 decimal digits; 0.00-0.50
                  "Lref_mult",    ### 2 decimal digits; 0.00-2.00
                  "multiplier"    ### 2 decimal digits; 0.00-2.00
                  )
    ga_default <- c(3,  0.2,  0.1,  0.2,  0.1, 0.05,  0.1, 1, 1) ### default values
    ga_lower <-   c(1,    0,    0,    0,    0,    0,    0, 0, 0) ### minima
    ga_upper <-   c(5,  0.5,  0.5,  0.5,  0.5,  0.5,  0.5, 2, 2) ### maxima
    ### ga() samples uniform real (double) values from lower ga_lower to 
    ### ga_upper and these are then rounded to the significant digits
    ### -> adjust ga_lower/upper so that minima/maxima have same probability
    ga_step  <-   c(1, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
    ga_lower <- ga_lower - (ga_step/2 - .Machine$double.eps)
    ga_upper <- ga_upper + (ga_step/2 - .Machine$double.eps)
    ### add some suggested parameterisations
    ga_suggestions <- rbind(### default
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            ### zero catch
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 0), 
                            ### interval
                            c(1, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            c(2, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            c(4, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            c(5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            ### lambda/gamma
                            c(3, 0.1, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            c(3, 0.3, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            c(3, 0.2, 0.05, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            c(3, 0.2, 0.2, 0.2, 0.1, 0.05, 0.1, 1, 1), 
                            c(3, 0.2, 0.1, 0.1, 0.1, 0.05, 0.1, 1, 1), 
                            c(3, 0.2, 0.1, 0.3, 0.1, 0.05, 0.1, 1, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.05, 0.05, 0.1, 1, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.3, 0.05, 0.1, 1, 1), 
                            ### thresholds
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.01, 0.1, 1, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.10, 0.1, 1, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.20, 0.1, 1, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.05, 1, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.2, 1, 1), 
                            ### multipliers
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 0.8, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 0.9, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1.1, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1.2, 1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 0.8), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 0.9), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1.1), 
                            c(3, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1, 1, 1.2)
                            )
  }
  
  ### turn of parameters not requested, i.e. limit to default value
  pos_default <- which(sapply(mget(ga_names, ifnotfound = FALSE), isFALSE))
  ga_lower[pos_default] <- ga_default[pos_default]
  ga_upper[pos_default] <- ga_default[pos_default]
  ### fix parameters?
  pos_fixed <- which(sapply(mget(ga_names, ifnotfound = FALSE), is.numeric))
  par_fixed <- names(pos_fixed)
  val_fixed <- as.vector(unlist(mget(ga_names, ifnotfound = FALSE)[pos_fixed]))
  ga_lower[pos_fixed] <- val_fixed
  ga_upper[pos_fixed] <- val_fixed
  ### remove not requested parameters from suggestions
  ga_suggestions[, pos_default] <- rep(ga_default[pos_default],
                                       each = nrow(ga_suggestions))
  ga_suggestions[, pos_fixed] <- rep(val_fixed,
                                     each = nrow(ga_suggestions))
  ga_suggestions <- unique(ga_suggestions)
  names(ga_suggestions) <- ga_names
  
  
  

  ### multiplier only: run all possible values
  if (isTRUE(multiplier) &
      !any(sapply(mget(setdiff(ga_names, "multiplier"), ifnotfound = FALSE),
                    isTRUE))) {
    m_vals <- seq(from = ga_lower[6], to = ga_upper[6], by = 0.01)
    ga_suggestions[1, ] <- ga_lower
    ga_suggestions <- ga_suggestions[rep(1, length(m_vals)), ]
    ga_suggestions$multiplier <- m_vals
    ### adapt GA dimensions
    maxiter <- run <- 1
    popSize <- length(m_vals)
    run_all <- TRUE
  } else {
    run_all <- FALSE
  }


  ### ---------------------------------------------------------------------- ###
  ### paths ####
  ### ---------------------------------------------------------------------- ###

  ### output path
  ### set name depending on which GA parameters are used
  scn_pars <- ga_names[setdiff(seq_along(ga_names), pos_default)]
  ### add fixed parameters
  scn_pars[which(scn_pars %in% par_fixed)] <- paste0(
    scn_pars[which(scn_pars %in% par_fixed)], val_fixed)
  scn_pars_c <- paste0(scn_pars, collapse = "-")
  path_out <- paste0("output/", MP, "/", n_iter, "_", n_yrs, "/",
                     scenario, "/", fhist, "/",
                     paste0(names(input), collapse = "_"), "/")
  dir.create(path_out, recursive = TRUE)

  ### objective function elements
  obj_fun <- c("SSB", "F", "C", "risk", "ICV", "ICES_PA", "ICES_PA2",
               "ICES_MSYPA")
  obj_fun_use <- mget(x = paste0("obj_", obj_fun),
                      ifnotfound = FALSE)
  for (i in seq_along(obj_fun)) {
    assign(x = paste0("obj_", obj_fun[i]), obj_fun_use[[i]])
  }
  obj_desc <- obj_fun[unlist(obj_fun_use)]
  obj_desc <- paste0("obj_", paste0(obj_desc, collapse = "_"), collapse = "")

  ### store input data in temp file
  inp_file <- tempfile()
  saveRDS(object = input, file = inp_file, compress = FALSE)
  rm(input)
  gc()

  ### ---------------------------------------------------------------------- ###
  ### check if previous solutions can be used as suggestions ####
  ### ---------------------------------------------------------------------- ###

  ### years for summary statistics
  file_ext <- ifelse(stat_yrs == "all", "_res",
                     paste0("_res_", stat_yrs))
  ### suffix if different risk limit used
  file_ext <- ifelse(isTRUE(!identical(risk_threshold, 0.05) &
                              isTRUE(obj_ICES_MSYPA)),
                     paste0(file_ext, "_", risk_threshold),
                     file_ext)
  file_ext <- paste0(file_ext, ".rds")

  if (isTRUE(add_suggestions)) {
    ### find files
    avail <- list.files(path_out, pattern = paste0("--", obj_desc, file_ext))
    avail <- gsub(x = avail, pattern = paste0("--", obj_desc, file_ext),
                  replacement = "")
    avail <- strsplit(x = avail, split = "-")
    ### need to have fewer parameters
    avail <- avail[which(sapply(avail, length) < length(scn_pars))]
    ### if some parameters fixed, remove suggestions without them
    if (isTRUE(length(avail) > 0)) {
      avail <- avail[which(sapply(avail, function(x)
        all(paste0(par_fixed, val_fixed) %in% x)))]
      ### skip parameters not used
      if (isTRUE(length(avail) > 0)) {
        avail <- avail[which(sapply(avail, function(x) all(x %in% scn_pars)))]
        if (isTRUE(length(avail) > 0)) {
          ### load results
          res_add <- lapply(avail, function(x) {
            tmp <- readRDS(file =
              paste0(path_out, paste0(x, collapse = "-"), "--", obj_desc,
                     "_res",
                     ifelse(identical(stat_yrs, "all"), "",
                            paste0("_", stat_yrs)),
                     ".rds"))
            tmp <- tmp@solution[1, ]
            if (is.na(tmp[which("upper_constraint" == names(tmp))])) {
              tmp[which("upper_constraint" == names(tmp))] <- Inf
            }
            return(tmp)
          })
          res_add <- do.call(rbind, res_add)
          if (isTRUE(nrow(res_add) > 1)) {
            res_add <- data.frame(res_add, stringsAsFactors = FALSE)
          } else {
            res_add <- data.frame(res_add, stringsAsFactors = FALSE)
          }
          cat("adding GA suggestions:\n")
          print(res_add)
          ### add to GA suggestions
          ga_suggestions <- rbind(ga_suggestions, res_add)
          ga_suggestions <- unique(ga_suggestions)
        }
      }
    }
  }

  ### ---------------------------------------------------------------------- ###
  ### run MSE with GA ####
  ### ---------------------------------------------------------------------- ###

  ### set random seed for reproducibility
  registerDoRNG(123)
  set.seed(1)

  ### run GA
  system.time({
    res <- ga(type = "real-valued", fitness = mp_fitness, inp_file = inp_file,
              obj_SSB = obj_SSB, obj_F = obj_F, obj_C = obj_C,
              obj_risk = obj_risk, obj_ICV = obj_ICV, obj_ICES_PA = obj_ICES_PA,
              obj_ICES_PA2 = obj_ICES_PA2, obj_ICES_MSYPA = obj_ICES_MSYPA,
              stat_yrs = stat_yrs, risk_threshold = risk_threshold,
              path = path_out, check_file = check_file,
              MP = MP,
              suggestions = ga_suggestions, lower = ga_lower, upper = ga_upper,
              names = ga_names,
              maxiter = maxiter, popSize = popSize, run = run,
              summarise_runs = TRUE,
              postFitness = mp_postFitness,
              keepBest = TRUE, parallel = cl1, seed = 1)
  })

  ### save result
  saveRDS(object = res, file = paste0(path_out, scn_pars_c,
                                      "--", obj_desc, file_ext))

  ### ---------------------------------------------------------------------- ###
  ### collate runs ####
  ### ---------------------------------------------------------------------- ###

  if (isTRUE(collate)) {
    files <- list.files(path = path_out, pattern = "[0-9]*[0-9].rds",
                        full.names = FALSE)
    files <- files[grep(x = files, pattern = "--", invert = TRUE)]
    names(files) <- sapply(files, function(x) {
      sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE)
    })
    scns <- lapply(files, function(x) {
      pars <- an(strsplit(sub(x = x, pattern = ".rds", replacement = "", fixed = TRUE),
                          split = "_")[[1]])
      names(pars) <- ga_names
      ### only keep scenarios where requested parameters are changed
      if (!all(ga_default[pos_default] == pars[pos_default])) return(NULL)
      if (!isTRUE(run_all)) {
        if (!all(val_fixed == pars[pos_fixed])) return(NULL)
      }
      stats <- readRDS(paste0(path_out, x))
      list(pars = pars, stats = stats)
    })
    scns[sapply(scns, is.null)] <- NULL
    #scns <- scns[order(sapply(scns, "[[", "obj"), decreasing = TRUE)]
    saveRDS(scns,
            file = paste0(path_out, scn_pars_c, "--", obj_desc, "_runs",
                          ifelse(identical(stat_yrs, "last10"), "_last10", ""),
                          ".rds"))
  }
  
}
  
### ------------------------------------------------------------------------ ###
### quit ####
### ------------------------------------------------------------------------ ###

quit(save = "no")
