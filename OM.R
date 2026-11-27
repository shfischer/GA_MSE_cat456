### ------------------------------------------------------------------------ ###
### create OMs ####
### ------------------------------------------------------------------------ ###

### some default parameters
n_workers <- 10
fhist_list <- c("one-way", "roller-coaster", "random")
# fhist <- "one-way"
stock_id <- 1:29
# stock_id <- 12 # pollack
OM <- TRUE
n_iter <- 500
yrs_hist <- 100
yrs_proj <- 100
selectivity_list <- c("default", "dome")
# selectivity <- "default"#"dome"



# args <- commandArgs(TRUE)
# if (exists(x = "args_local")) args <- append(args, args_local)
# print("arguments passed on to this script:")
# print(args)

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

### recruitment bias correction
### -> ensure that reference points from FLBRP match OM dynamics
### -> adjust recruitment residuals
rec_bias_correction <- TRUE

### ------------------------------------------------------------------------ ###
### with uniform distribution and random F trajectories ####
### ------------------------------------------------------------------------ ###
# fhist <- "random"#"one-way"
if (isTRUE("random" %in% fhist_list)) {
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
### -> use updated version
brps <- readRDS("input/brps_len.rds")
brps_dome <- readRDS("input/brps_dome_len.rds")

### create FLStocks
stocks_subset <- stocks$stock[stock_id]#"bll"

if (exists("OM")) { 
  
  if (isTRUE(OM)) {
    
    stks_hist <- foreach(fhist = fhist_list) %:%
      foreach(selectivity = selectivity_list) %:% 
      foreach(stock = stocks_subset, .errorhandling = "pass", 
              .packages = c("FLCore", "FLasher", "FLBRP")) %dopar% {
      
      if (identical(selectivity, "default")) {
        brps_i <- brps
      } else if (identical(selectivity, "dome")) {
        brps_i <- brps_dome
      }
      stk <- as(brps_i[[stock]], "FLStock")
      refpts <- refpts(brps_i[[stock]])
      stk <- qapply(stk, function(x) {#browser()
        dimnames(x)$year <- as.numeric(dimnames(x)$year) - 1; return(x)
      })
      stk <- stf(stk, yrs_hist + yrs_proj - dims(stk)$year + 1)
      stk <- propagate(stk, n_iter)
      ### create stock recruitment model
      stk_sr <- FLSR(params = params(brps_i[[stock]]), 
                     model = model(brps_i[[stock]]))
      ### create residuals for (historical) projection
      
      ### recruitment bias correction
      if (isTRUE(rec_bias_correction)) {
        # m <- 1
        # cv <- 0.6
        # v <- (cv*m)^2 ### variance
        # mu <- log(m^2/sqrt(m^2+v))
        # sR <- sqrt(log(v/m^2+1))
        sR <- 0.6
        mu <- 0 - (sR^2)/2
      } else {
        sR <- 0.6
        mu <- 0
      }
      set.seed(0)
      residuals(stk_sr) <- rlnoise(dim(stk)[6], rec(stk) %=% mu, 
                                   sd = sR, b = 0)
      # ### replicate residuals from GA paper
      # set.seed(0)
      # residuals(stk_sr)[, ac(0:150)] <- rlnoise(dim(stk)[6], 
      #                                           rec(stk)[, ac(0:150)] %=% mu,
      #                                           sd = sR, b = 0)
      # ### replicate residuals from catch rule paper for historical period
      # set.seed(0)
      # residuals <- rlnoise(dim(stk)[6], (rec(stk) %=% mu)[, ac(1:100)], 
      #                      sd = sR, b = 0)
      # residuals(stk_sr)[, ac(1:100)] <- residuals[, ac(1:100)]
      
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
        ### increase to 0.8Fcrash in 15 years
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
      path <- paste0("input/", n_iter, "_", yrs_proj, "/",
                     switch(selectivity, 
                            "default" = "OM",
                            "dome" = "OM_dome"), 
                     "_len",
                     "/", fhist, "/", stock, "/")
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

# gc()
# if (!isFALSE(cl)) clusterEvalQ(cl, {gc()})
stopCluster(cl)
gc()


### ------------------------------------------------------------------------ ###
### tmp - explore length stuff ####
### ------------------------------------------------------------------------ ###
req_pckgs <- c("FLCore", "FLasher", "FLBRP", "mse", "FLife", 
               "tidyr", "dplyr", "foreach", "doParallel")
for (i in req_pckgs) library(package = i, character.only = TRUE)
library(FLRef)

source("funs.R")
source("funs_lcc.R")

brp <- readRDS("input/brps_len.rds")[["pol"]]
stk_stf <- readRDS("input/500_100/OM_len/one-way/pol/stk.rds")
stk_sr <- readRDS("input/500_100/OM_len/one-way/pol/sr.rds")

# sweep(harvest(stk_stf), 2, apply(harvest(stk_stf), 2, max, na.rm = TRUE), "/")[, ac(50)]


gear <- attr(brp, "gear")
set.seed(1)
len_sim <- lfd.sim(catch.n(stk_stf[, ac(100)]), gear = gear[[1]], scale = FALSE)

as.data.frame(FLCore::iter(len_sim, 1)) %>%
  ggplot(aes(x = len, y = data)) +
  geom_line()
as.data.frame(FLCore::iter(len_sim, 1:5)) %>%
  ggplot(aes(x = len, y = data, colour = iter)) +
  geom_line()
as.data.frame(FLCore::iter(len_sim, 1:500)) %>%
  ggplot(aes(x = len, y = data, group = iter)) +
  geom_line(alpha = 0.05)


res_length <- rUnif_lfd(gear, years = 0:200, nsim = dim(stk_stf)[6], seed = 123)

cn <- catch.n(stk_stf[, ac(1:100)])
cn[cn <= 0] <- sqrt(.Machine$double.eps)


len_sim_full <- lfd.sim(object = cn[, ac(50:100)], 
                        gear = gear[[1]], 
                        ess_age = gear[[1]]$ess_age,
                        ess_len = gear[[1]]$ess_len,
                        scale = FALSE,
                        u_age = res_length$age$dome,
                        u_len = res_length$len$dome
)
saveRDS(object = len_sim_full, 
        file = "input/500_100/OM_len/one-way/pol/lfd.rds")

len_sim <- lfd.sim(object = cn[, ac(99:100)], 
                   gear = gear[[1]], 
                   ess_age = gear[[1]]$ess_age,
                   ess_len = gear[[1]]$ess_len,
                   scale = FALSE,
                   u_age = res_length$age$dome,
                   u_len = res_length$len$dome
                   )
oem_obs_len <- window(len_sim, end = 200)


### operating model ####
stk <- stk_stf
sr <- stk_sr
stk <- window(stk, start = 50)
sr@residuals <- window(sr@residuals, start = 50)
om <- FLom(stock = stk, ### stock 
           sr = sr, ### stock recruitment and precompiled residuals
           fleetBehaviour = mseCtrl(),
           projection = mseCtrl(method = fwd_attr))


oem_observations <- list(
  stk = stk, 
  idx = FLIndices(LFD = FLIndex(index = oem_obs_len,
                                plusgroup = NA)))
validObject(oem_observations)
validObject(oem_observations$idx)

oem <- FLoem(observations = oem_observations, 
             deviances = list(lfd = res_length),
             method = oem_lfd_sampling,
             args = list(gear = gear,
                         ess_age = 50,
                         ess_len = 500))

ctrl <- mpCtrl(list(
  est = mseCtrl(method = est_lfd_lcc,
                args = list(Linf = c(gear$dome$lhpar["linf"]))),
  phcr = mseCtrl(method = phcr_lfd_lcc,
                 args = list(m = 1, 
                             s_target = 2.2, 
                             catch_limit = 50)),
  hcr = mseCtrl(method = hcr_lcc_lfd,
                args = list(interval = 2,
                            upper_constraint = 1.2,
                            lower_constraint = 0.7))
))

args <- list(fy = 200, y0 = 75, iy = 100)


debug(mp)
undebug(mp)
debug(goFish, signature = "FLom")
undebug(goFish, signature = "FLom")
#debugonce(oem@method)
debugonce(ctrl$est@method)
debugonce(ctrl$phcr@method)
debug(ctrl$phcr@method)
undebug(ctrl$phcr@method)
debugonce(ctrl$hcr@method)
res <- mp(om = om, oem = oem, args = args, ctrl = ctrl,
          parallel = FALSE)


debugonce(ctrl.oem$method)
o.out <- ctrl.oem$method(stk = ctrl.oem$stk, observations = ctrl.oem$observations, deviances = ctrl.oem$deviances, args = ctrl.oem$args, tracking = ctrl.oem$tracking, gears = ctrl.oem$gear, ess_age = ctrl.oem$ess_age, ess_len = ctrl.oem$ess_len)

debugonce(ctrl.est$method)
out.assess <- ctrl.est$method(stk = ctrl.est$stk, idx = ctrl.est$idx, args = ctrl.est$args, tracking = ctrl.est$tracking, Linf = ctrl.est$Linf)


source("funs_OM.R")
debug(input_mp_lcc)
undebug(input_mp_lcc)

input <- input_mp_lcc(stock = "pol", fhist = "one-way", 
                      interval = 2, multiplier = 1, 
                      upper_constraint = 1.2, lower_constraint = 0.7,
                      catch_limit = 0.5)
res <- do.call(mp, input)
