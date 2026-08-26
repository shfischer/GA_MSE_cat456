### ------------------------------------------------------------------------ ###
### update BRPs ####
### ------------------------------------------------------------------------ ###
### original brps created with early development version of FLife
### update with newer version of FLife


library(FLife)
library(FLasher)
library(foreach)

### stocks and life-history parameters
stocks_lh <- read.csv("input/stocks.csv")

### set steepness to 0.75
stocks_lh$s <- 0.75

### default t0 = -0.1
stocks_lh$t0[is.na(stocks_lh$t0)] <- -0.1

### original brps
brps_original <- readRDS("input/brps.rds")

### fbar ranges
stocks_lh$minfbar <- sapply(brps_original, function(x) range(x)[["minfbar"]])
stocks_lh$maxfbar <- sapply(brps_original, function(x) range(x)[["maxfbar"]])

brps <- foreach(selectivity = c("standard", "dome")) %:% 
  foreach(stk = stocks_lh$stock) %do% { # stk="pol"
  lh_i <- stocks_lh[stocks_lh$stock == stk, ]
  lh_pars <- FLPar(lh_i[c("linf", "k", "t0", "a", "b", "l50", "a50", "s")])
  lh_pars <- lhPar(lh_pars)
  
  ### dome shaped selectivity - if needed
  if (identical(selectivity, "dome")) {
    ### dome shaped selectivity
    ### sel1 - maximum selectivity at 1st age at 100% maturity
    ### (a50+ato95; default)
    ### sel2 - sd for left limb: 1 (default, steep increase)
    ### sel3 - sd for right limb: 5 (default for dome, slower decrease)
    lh_pars["sel2"] <- 1
    lh_pars["sel3"] <- 5
  }
  
  ### Max age: age at l = 0.95 * linf
  max_age <- ceiling(log(0.05)/(-c(lh_pars["k"])) + c(lh_pars["t0"]))
  
  ### fbar age range
  #fbar_age <- c(lh_i$minfbar, lh_i$maxfbar)
  ### don't use -> use all ages
  
  ### create brp
  brp <- lhEql(lh_pars, range = c(min = 1, max = max_age, 
                                  #minfbar = fbar_age[1], maxfbar = fbar_age[2], 
                                  minfbar = 1, maxfbar = max_age,
                                  plusgroup = max_age))
  
  ### changes to original brp:
  ### m.spwn & harvest.spwn = 0 
  ### before: difference from a50 to closest integer below...
  ### stock length at age
  ### - before: age + m.spwn
  ### - new: age (beginning of year/age)
  ### catch length at age - same
  ### weights at age
  ### - before: catch weights incorrectly set to stock weights
  ### - new: fixed, mid-year weights
  ### - stock weights now at beginning of year/age (m.spwn = 0)
  ### fisheries selectivity - double normal
  ### - before: 1st age with full selectivity = a50
  ### - new: a50+1, less steep
  ###        -> less fishing of young fish
  ### fbar ages -> use all ages for simplicity
  
  ### replicate WKLIFE VII fisheries selectivity
  if (FALSE) {
    lh_pars2 <- lh_pars
    lh_pars2["sel1"] <- lh_pars2["a50"]
    lh_pars2["sel2"] <- 2
    brp2 <- lhEql(lh_pars2, range = c(min = 1, max = max_age, 
                                      minfbar = fbar_age[1], maxfbar = fbar_age[2], 
                                      plusgroup = max_age))
  }
  
  return(brp)
  
}

names(brps) <- c("default", "dome")
brps <- lapply(brps, function(x) {
  names(x) <- stocks_lh$stock
  return(x)
})

# lapply(brps, range)
# lapply(brps, refpts)

### Blim
### SSB responding to 30% recruitment impairment
### identical for all stocks because all have same 
### - B0 (1000)
### - recruitment steepness h (0.75)
R_diff <- function(RR0, brp, ssb) {
  R <- c(params(brp)["a"]) * ssb / (c(params(brp)["b"]) + ssb)
  R0 <- c(refpts(brp)["virgin", "rec"])
  return((R - R0 * RR0)^2)
}
res <- optimise(f = R_diff, 
                RR0 = 0.7, brp = brps$default$pol,
                lower = 0, upper = 1000)
Blim <- round(res$minimum, 2)

### add as attribute to brps
brps <- lapply(brps, function(x) {
  lapply(x, function(y) {
    attr(y, "Blim") <- Blim
    return(y)
  })
})


saveRDS(brps$default, file = "input/brps_new.rds")
saveRDS(brps$dome, file = "input/brps_dome_new.rds")


### some explorations with recruitment residuals...
# if (FALSE) {
#   
# brp_ <- brps$default$pol
# refpts(brp_)["fmax", ] <- NA
# brp_old <- brps_original$pol
# brp_dome <- brps$dome$pol
# refpts(brp_dome)["fmax", ] <- NA
# plot(brp_old)
# plot(brp_)
# plot(brp_dome)
#   
#   
#   stk <- as(brp, "FLStock")
#   stk[, 1] <- stk[, 2]
#   stk <- propagate(stk, 500)
#   stk_sr <- FLSR(params = params(brp), model = model(brp))
#   
#   set.seed(0)
#   residuals <- rlnoise(500, rec(stk) %=% 0, 
#                        sd = 0.6, b = 0)
#   
#   ctrl <- fwdControl(data.frame(year = 2:100, quant = "f", 
#                                 value = c(refpts(brp)["msy", "harvest"])))
#   stk_fwd <- fwd(stk, sr = stk_sr, 
#                  control = ctrl,
#                  residuals = residuals)
#   plot(stk_fwd)
#   fbar(stk_fwd)
#   refpts(brp)["msy", "harvest"]
#   ssb(stk_fwd)
#   refpts(brp)["msy", "ssb"]
#   
#   
#   stk_0 <- fwd(stk, sr = stk_sr, 
#                control = fwdControl(data.frame(year = 2:100, quant = "f", 
#                                                value = 0)),
#                residuals = residuals)
#   plot(stk_0)
#   ssb(stk_0)
#   
#   ### simplified example
#   lh_pars <- lhPar(FLPar(linf = 100))
#   brp <- lhEql(lh_pars)
#   stk <- as(brp, "FLStock")[, -101]
#   #stk[, 2:100] <- stk[, 1]
#   n_its <- 10000
#   stk10k <- propagate(stk, 10000)
#   stk500 <- propagate(stk, 500)
#   stk_sr <- FLSR(params = params(brp), model = model(brp))
#   set.seed(0)
#   residuals <- rlnoise(n_its, rec(stk) %=% 0, sd = 0.6)
#   ctrl <- fwdControl(data.frame(year = 2:100, quant = "f", value = 0))
#   stk0 <- fwd(stk, sr = stk_sr, 
#                  control = ctrl,
#                  residuals = residuals)
#   plot(stk0)
#   
#   
#   residuals2 <- residuals/mean(c(residuals))
#   stk0_2 <- fwd(stk, sr = stk_sr, 
#               control = ctrl,
#               residuals = residuals2)
#   plot(stk0_2)
#   
#   ctrl_fmsy <- fwdControl(data.frame(year = 2:100, quant = "f", 
#                                      value = c(refpts(brp)["msy", "harvest"])))
#   # ctrl_fmsy <- fwdControl(data.frame(year = 2, quant = "f", 
#   #                                    value = 0.1))
#   stk_fmsy <- fwd(stk, sr = stk_sr,
#                   control = ctrl_fmsy,
#                   residuals = residuals2)
#   plot(stk_fmsy)
#   ssb(stk_fmsy)
#   refpts(brp)
#   
#   
#   set.seed(1)
#   residuals500 <- rlnorm(500, rec(stk) %=% 0, 0.6)
#   
#   
#   
#   ### https://stackoverflow.com/questions/56821688/sample-a-lognormal-distribution-to-an-exact-mean-and-sd
#   
#   stk <- as(brp, "FLStock")
#   stk[, 1] <- stk[, 2]
#   stk <- propagate(stk, 500)
#   stk_sr <- FLSR(params = params(brp), model = model(brp))
#   
#   set.seed(0)
#   residuals <- rec(stk) %=% NA_real_
#   m <- 1
#   s <- 0.6
#   residuals[] <- rlnorm(n=length(residuals), 
#                      meanlog=log(m^2 / sqrt(s^2 + m^2)), 
#                      sdlog=sqrt(log(1 + (s^2 / m^2))))
#   residuals[] <- rlnorm(n=length(residuals), 
#                         meanlog=log(m^2 / sqrt(s^2 + m^2)), 
#                         sdlog=s)
#   residuals[] <- rlnorm(n=length(residuals), 
#                         meanlog=0 - (s^2)/2, 
#                         sdlog=s)
#   summary(residuals)
#   
#   
#   ctrl <- fwdControl(data.frame(year = 2:100, quant = "f",
#                                 value = c(refpts(brp)["msy", "harvest"])))
#   ctrl <- fwdControl(data.frame(year = 2:100, quant = "f", 
#                                 value = 0))
#   stk_fwd <- fwd(stk, sr = stk_sr, 
#                  control = ctrl,
#                  residuals = residuals)
#   ssb(stk_fwd)
#   summary(ssb(stk_fwd)[, ac(91:100)])
#   
#   stk1_fwd <- fwd(iter(stk, 1), sr = stk_sr, 
#                   control = ctrl)
#   ssb(stk1_fwd)
#   
#   set.seed(1)
#   m <- 1
#   s <- 0.6
#   data_set <- rlnorm(n=1e+5, 
#                      meanlog=log(m^2 / sqrt(s^2 + m^2)), 
#                      sdlog=sqrt(log(1 + (s^2 / m^2))))
#   mean(data_set)
#   median(data_set)
#   sd(data_set)
#   sd(data_set)/mean(data_set)
#   
#   hist(data_set, breaks = 50)
#   
#   set.seed(1)
#   res <- rnorm(1e+6, 0, 0.6)
#   summary(res)
#   
#   summary(exp(res))
#   summary(exp(res - (0.6^2)/2))
#   
#   sd(exp(res - (0.6^2)/2))
# }
