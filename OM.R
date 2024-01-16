### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

args <- commandArgs(TRUE)
if (exists(x = "args_local")) args <- append(args, args_local)
print("arguments passed on to this script:")
print(args)

### evaluate arguments passed to R
for (i in seq_along(args)) eval(parse(text = args[[i]]))

### ------------------------------------------------------------------------ ###
### set up environment ####
### ------------------------------------------------------------------------ ###

### load packages
### use mse fork from shfischer/mse, branch mseDL2.0 
### remotes::install_github("shfischer/mse", ref = "mseDL2.0)
req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLife", 
               "tidyr", "dplyr", "foreach", "doParallel")
for (i in req_pckgs) library(package = i, character.only = TRUE)

### load additional functions
source("funs.R")

### ------------------------------------------------------------------------ ###
### setup parallel environment ####
### ------------------------------------------------------------------------ ###

if (isTRUE(n_workers > 1)) {
  ### start doParallel cluster
  cl <- makeCluster(n_workers)
  registerDoParallel(cl)
  cl_length <- length(cl)
  ### load packages and functions into parallel workers
  . <- foreach(i = seq(cl_length)) %dopar% {
    for (i in req_pckgs) library(package = i, character.only = TRUE,
                                 warn.conflicts = FALSE, verbose = FALSE,
                                 quietly = TRUE)
    source("funs.R", echo = FALSE)
  }
} else {
  cl <- FALSE
}

### ------------------------------------------------------------------------ ###
### fishing history dimensions ####
### ------------------------------------------------------------------------ ###

# n_iter <- 500
# yrs_hist <- 100
# yrs_proj <- 50

set.seed(2)

### ------------------------------------------------------------------------ ###
### with uniform distribution and random F trajectories ####
### ------------------------------------------------------------------------ ###
# fhist <- "random"#"one-way"
if (identical(fhist, "random")) {
  start <- rep(0, n_iter)
  middle <- runif(n = n_iter, min = 0, max = 1)
  end <- runif(n = n_iter, min = 0, max = 1)
  df <- t(sapply(seq(n_iter), 
    function(x) {
      c(approx(x = c(1, yrs_hist/2), 
               y = c(start[x], middle[x]), 
               n = yrs_hist/2)$y,
        approx(x = c(yrs_hist/2, yrs_hist + 1), 
               y = c(middle[x], end[x]), 
               n = (yrs_hist/2) + 1)$y[-1])
    }))
  ### convert into FLQuant for FLasher
  f_quant <- FLQuant(NA, dimnames = list(year = seq(yrs_hist),
                                         iter = seq(n_iter)))
  f_quant[] <- t(df)
}

### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

### get lhist for stocks
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)

### BRPs from Fischer et al. (2020)
brps <- readRDS("input/brps.rds")

### create FLStocks
stocks_subset <- stocks$stock[stock_id]#"bll"

if (exists("OM")) { 
  
  if (isTRUE(OM)) {
    
    stks_hist <- foreach(stock = stocks_subset, .errorhandling = "pass", 
                         .packages = c("FLCore", "FLash", "FLBRP")) %dopar% {
      stk <- as(brps[[stock]], "FLStock")
      refpts <- refpts(brps[[stock]])
      stk <- qapply(stk, function(x) {#browser()
        dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
      })
      stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
      stk <- propagate(stk, n_iter)
      ### create stock recruitment model
      stk_sr <- FLSR(params = params(brps[[stock]]), model = model(brps[[stock]]))
      ### create residuals for (historical) projection
      set.seed(0)
      residuals(stk_sr) <- rlnoise(dim(stk)[6], rec(stk) %=% 0, 
                                   sd = 0.6, b = 0)
      ### replicate residuals from GA paper
      set.seed(0)
      residuals(stk_sr)[, ac(0:150)] <- rlnoise(dim(stk)[6], 
                                                rec(stk)[, ac(0:150)] %=% 0,
                                                sd = 0.6, b = 0)
      ### replicate residuals from catch rule paper for historical period
      set.seed(0)
      residuals <- rlnoise(dim(stk)[6], (rec(stk) %=% 0)[, ac(1:100)], 
                           sd = 0.6, b = 0)
      residuals(stk_sr)[, ac(1:100)] <- residuals[, ac(1:100)]
      
      ### fishing history from previous paper
      if (isTRUE(fhist == "one-way")) {
        
        ### 0.5Fmsy until year 75, then increase to 0.8Fcrash
        fs <- rep(c(refpts["msy", "harvest"]) * 0.5, 74)
        f0 <- c(refpts["msy", "harvest"]) * 0.5
        fmax <- c(refpts["crash", "harvest"]) * 0.8
        rate <- exp((log(fmax) - log(f0)) / (25))
        fs <- c(fs, rate ^ (1:25) * f0)
        
        ### control object
        ctrl <- fwdControl(data.frame(year = 2:100, quant = "f", value = fs))
        
      ### roller-coaster
      } else if (isTRUE(fhist == "roller-coaster")) {
        
        ### 0.5Fmsy until year 75, 
        ### increase to 0.8Fcrash in 10 years
        ### keep at 0.8Fcrash for 5 years
        ### reduce to Fmsy in last 5 years
        fs <- rep(c(refpts["msy", "harvest"]) * 0.5, 75)
        f0_up <- c(refpts["msy", "harvest"]) * 0.5
        fmax_up <- c(refpts["crash", "harvest"]) * 0.8
        yrs_up <- 15
        rate_up <- exp((log(fmax_up) - log(f0_up)) / yrs_up)
        yrs_down <- 6
        f0_down <- c(refpts["msy", "harvest"])
        rate_down <- exp((log(fmax_up) - log(f0_down)) / yrs_down)
        fs <- c(fs, rate_up ^ seq(yrs_up) * f0_up, rep(fmax_up, 3),
                rev(rate_down ^ seq(yrs_down) * f0_down))
        
        ### control object
        ctrl <- fwdControl(data.frame(year = 2:100, quant = "f", value = fs))
        
      ### random F trajectories
      } else if (isTRUE(fhist == "random")) {
        
        ### target * Fcrash
        ctrl <- fwdControl(target = f_quant * c(refpts["crash", "harvest"]), 
                           quant = "f")

      }
      
      ### project fishing history
      stk_stf <- fwd(stk, control = ctrl, sr = stk_sr, 
                     residuals = residuals(stk_sr), maxF = 5)
      
      ### save OM files
      name(stk_stf) <- stock
      path <- paste0("input/", n_iter, "_", yrs_proj, "/OM/", fhist, 
                     "/", stock, "/")
      dir.create(path, recursive = TRUE)
      
      ### stock & sr (full history)
      saveRDS(stk_stf, file = paste0(path, "stk.rds"))
      saveRDS(stk_sr, file = paste0(path, "sr.rds"))
      
      return(NULL)
      #return(list(stk = stk_stf, sr = stk_sr))
    }
  }
}

# source("funs_OM.R")
# debugonce(input_mp)
# input <- input_mp(stock = "pol", fhist = "one-way", n_iter = 500, n_yrs = 50,
#                   MP = "hr")
# res <- do.call(mp, input)
# 
# input2 <- input_mp(stock = "pol", fhist = "one-way", n_iter = 500, n_yrs = 50,
#                   MP = "hr", interval = 2)
# res2 <- do.call(mp, input2)

# debugonce(input_mp)
# input_list <- input_mp(stocks = "pol", fhist = "one-way", n_iter = 500, n_yrs = 50,
#                        MP = "hr")
# res_ <- do.call(mp, input_list[[1]])

### ------------------------------------------------------------------------ ###
### gc() ####
### ------------------------------------------------------------------------ ###

gc()
if (!isFALSE(cl)) clusterEvalQ(cl, {gc()})

