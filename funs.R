### ------------------------------------------------------------------------ ###
### observations ####
### ------------------------------------------------------------------------ ###
obs_generic <- function(stk, observations, deviances, args, tracking,
                        idxB = TRUE, ### stock/biomass index?
                        ssb_idx = FALSE, tsb_idx = FALSE, ### use SSB idx
                        idx_timing = FALSE, ### consider survey timing?
                        idx_dev = FALSE,
                        lngth = FALSE, ### catch length data?
                        lngth_dev = FALSE, 
                        lngth_par,
                        PA_status = FALSE,
                        PA_status_dev = FALSE,
                        PA_Bmsy = FALSE, PA_Fmsy = FALSE,
                        catch_dev = FALSE,
                        ...) {

  #ay <- args$ay
  ### update observations
  observations$stk <- stk
  
  ### biomass index
  if (isTRUE(idxB)) {
    ### use SSB as index?
    if (isTRUE(ssb_idx)) {
      index(observations$idx$idxB) <- ssb(observations$stk)
    ### TSB?
    } else  if (isTRUE(tsb_idx)) {
      index(observations$idx$idxB) <- tsb(observations$stk)
    ### otherwise calculate biomass index
    } else {
      sn <- stk@stock.n
      ### reduce by F and M?
      if (isTRUE(idx_timing)) {
        sn <- sn * exp(-(harvest(stk) * harvest.spwn(stk) +
                           m(stk) * m.spwn(stk)))
      }
      index(observations$idx$idxB) <- quantSums(sn * stk@stock.wt * 
                                        index(observations$idx$sel))
    }
  } else {
    ### insert dummy index
    index(observations$idx$idxB) <- ssb(stk) %=% NA_real_
    
  }
  ### use mean length in catch?
  if (isTRUE(lngth)) {
    index(observations$idx$idxL) <- lmean(stk = stk, params = lngth_par)
  }
  ### stock status for PA buffer?
  if (isTRUE(PA_status)) {
    index(observations$idx$PA_status)[] <- 
      ssb(observations$stk) > 0.5*PA_Bmsy & fbar(observations$stk) < PA_Fmsy
  }
  
  ### observation model
  stk0 <- observations$stk
  idx0 <- observations$idx
  
  ### add deviances to catch?
  if (isTRUE(catch_dev)) {
    ### add uncertainty to total catch
    catch(stk0) <- catch(observations$stk) * deviances$stk$catch
  }
  
  ### add deviances to index?
  if (isTRUE(idxB) & isTRUE(idx_dev)) {
    if (isTRUE(ssb_idx) | isTRUE(tsb_idx)) {
      index(idx0$idxB) <- index(observations$idx$idxB) * index(deviances$idx$idxB)
    } else {
      index(idx0$idxB) <- quantSums(stk@stock.n * stk@stock.wt * 
                            index(observations$idx$sel) * index(deviances$idx$sel))
      if (isTRUE("idxB" %in% names(deviances$idx)) & 
          all.equal(dim(index(deviances$idx$idxB)), dim(index(idx0$idxB))))
        index(idx0$idxB) <- index(idx0$idxB) * index(deviances$idx$idxB)
    }
  }
  ### uncertainty for catch length
  if (isTRUE(lngth) & isTRUE(lngth_dev)) {
    index(idx0$idxL) <- index(observations$idx$idxL) * index(deviances$idx$idxL)
  }
  ### uncertainty for stock status for PA buffer
  if (isTRUE(PA_status) & isTRUE(PA_status_dev)) {
    index(idx0$PA_status) <- ifelse(index(observations$idx$PA_status) == TRUE, 
                                index(deviances$idx$PA_status)["positive", ],
                                index(deviances$idx$PA_status)["negative", ])
  }
  
  return(list(stk = stk0, idx = idx0, observations = observations,
              tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### estimator ####
### ------------------------------------------------------------------------ ###

### category 4-6 - catch and length data
est_CL <- function(stk, idx, tracking, args,
                   n_catch = 5, n_length_1 = 5, n_length_2 = n_length_1,
                   r_catch = TRUE, r_length = TRUE, length_average = TRUE,
                   lag_catch = 1, lag_length = 1,
                   interval = 3,
                   ...) {
  
  ay <- args$ay
  iy <- args$iy
  
  ### only run "model" when new advice is required
  if ((ay - iy) %% interval == 0) {
  
    ### catch/advice
    ### first year - use catch
    if (identical(ay, iy)) {
      ### use current year (ay)
      advice_current <- catch(stk)[, ac(ay)]
    } else {
      ### other years - use advice
      advice_current <- tracking[[1]]["isys", ac(ay)]
    }
    
    ### linear model of (mean standardised) length indicator
    if (isTRUE(r_length)) {
      r_length <- apply(index(idx$idxL)[, ac(seq(to = ay - lag_length, 
                                          length.out = n_length_1))], 
                        6, function(x) {
        tmp_data <- as.data.frame(FLQuant(x))
        tmp_data$data <- tmp_data$data/mean(tmp_data$data, na.rm = TRUE)
        out <- try(lm(data ~ year, data = tmp_data)$coefficients[["year"]], 
                   silent = TRUE)
        if (is(out, "try-error")) {
          return(0) ### return 0, i.e. no trend detected
        } else {
          return(out)
        }
      })
    } else {
      r_length <- 0
    }
    
    ### linear model of (mean standardised) catch
    if (isTRUE(r_catch)) {
      r_catch <- apply(catch(stk)[, ac(seq(to = ay - lag_catch, 
                                           length.out = n_catch))], 
                       6, function(x) {
        tmp_data <- as.data.frame(FLQuant(x))
        tmp_data$data <- tmp_data$data/mean(tmp_data$data, na.rm = TRUE)
        out <- try(lm(data ~ year, data = tmp_data)$coefficients[["year"]], 
                   silent = TRUE)
        if (is(out, "try-error")) {
          return(0) ### return 0, i.e. no trend detected
        } else {
          return(out)
        }
      })
    } else {
      r_catch <- 0
    }
    
    ### average length
    if (isTRUE(length_average)) {
      length_average <- apply(index(idx$idxL)[, ac(seq(to = ay - lag_length, 
                                                length.out = n_length_2))],
                              6, mean, na.rm = TRUE)
    } else {
      length_average <- NA
    }
  
  } else {
    
    ### dummy values when no new advice calculated
    r_length <- r_catch <- length_average <- advice_current <- NA
    
  }
  
  ### save results
  tracking[[1]]["r_length", ac(ay)] <- r_length
  tracking[[1]]["r_catch", ac(ay)] <- r_catch
  tracking[[1]]["length_average", ac(ay)] <- length_average
  tracking[[1]]["A_last", ac(ay)] <- advice_current
  
  return(list(stk = stk, tracking = tracking))
  
}

est_comps <- function(stk, idx, tracking, args,
                      comp_r = FALSE, comp_f = FALSE, comp_b = FALSE,
                      comp_i = FALSE, comp_c = FALSE, comp_A = TRUE,
                      comp_m = FALSE, comp_hr = FALSE,
                      idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                      idxB_range_3 = 1,
                      catch_lag = 1, catch_range = 1,
                      Lref, I_trigger, Lref_mult = 1,
                      idxL_lag = 1, idxL_range = 1,
                      pa_buffer = FALSE, pa_size = 0.8, pa_duration = 3,
                      pa_buffer_conditional = FALSE,
                      Bmsy = NA,
                      ...) {
  
  ay <- args$ay
  iy <- args$iy
  
  ### component r: index trend
  if (isTRUE(comp_r)) {
    r_res <- est_r(idx = index(idx$idxB), ay = ay,
                   idxB_lag = idxB_lag, idxB_range_1 = idxB_range_1, 
                   idxB_range_2 = idxB_range_2)
  } else {
    r_res <- 1
  }
  tracking[[1]]["comp_r", ac(ay)] <- r_res
  
  ### component f: length data
  if (isTRUE(comp_f)) {
    f_res <- est_f(idx = index(idx$idxL), ay = ay,
                   Lref = Lref, idxL_range = idxL_range, idxL_lag = idxL_lag)
  } else {
    f_res <- 1
  }
  tracking[[1]]["comp_f", ac(ay)] <- f_res
  
  ### component b: biomass safeguard
  if (isTRUE(comp_b)) {
    b_res <- est_b(idx = index(idx$idxB), ay = ay,
                   I_trigger = I_trigger, idxB_lag = idxB_lag, 
                   idxB_range_3 = idxB_range_3)
  } else {
    b_res <- 1
  }
  
  ### PA buffer
  if (isTRUE(pa_buffer)) {
    b_res <- est_pa(idx = index(idx$PA_status), ay = ay, 
                    tracking = tracking, idxB_lag = idxB_lag,
                    pa_size = pa_size, pa_duration = pa_duration)
  }
  
  ### conditional PA buffer - use length
  if (isTRUE(pa_buffer_conditional)) {
    b_res <- est_pa_conditional(ay = ay, 
                                tracking = tracking, 
                                pa_size = pa_size, 
                                pa_duration = pa_duration,
                                idx = index(idx$idxL),
                                Lref = Lref, idxL_range = idxL_range, 
                                Lref_mult = Lref_mult,
                                idxL_lag = idxL_lag)
  }
  tracking[[1]]["comp_b", ac(ay)] <- b_res
  
  ### component i: index value
  if (isTRUE(comp_i)) {
    i_res <- est_i(idx = index(idx$idxB), ay = ay,
                   idxB_lag = idxB_lag, idxB_range_3 = idxB_range_3)
  } else {
    i_res <- 1
  }
  tracking[[1]]["comp_i", ac(ay)] <- i_res
  
  ### current catch/advice
  if (isTRUE(comp_c)) {
    c_res <- est_c(ay = ay, catch = catch(stk), catch_lag = catch_lag, 
                   catch_range = catch_range)
  } else if (isTRUE(comp_A)) {
    c_res <- est_A(ay = ay, iy = iy, catch = catch(stk), catch_lag = catch_lag,
                   catch_range = catch_range, tracking = tracking)
  } else {
    c_res <- 1
  }
  tracking[[1]]["comp_c", ac(ay)] <- c_res
  
  ### component m: multiplier
  if (!isFALSE(comp_m)) {
    m_res <- comp_m
    ### subset to iteration when simultion is split into blocks
    if (isTRUE(length(comp_m) > dims(stk)$iter)) {
      m_res <- comp_m[as.numeric(dimnames(stk)$iter)]
    }
  } else {
    m_res <- 1
  }
  tracking[[1]]["multiplier", ac(ay)] <- m_res
  
  ### component hr: harvest rate (catch/idx)
  if (!isFALSE(comp_hr)) {
    hr_res <- comp_hr
    ### subset to iteration when simultion is split into blocks
    if (isTRUE(length(comp_hr) > dims(stk)$iter)) {
      hr_res <- comp_hr[as.numeric(dimnames(stk)$iter)]
    }
  } else {
    hr_res <- 1
  }
  tracking[[1]]["comp_hr", ac(ay)] <- hr_res
  
  return(list(stk = stk, tracking = tracking))
  
}

### biomass index trend
est_r <- function(idx, ay,
                  idxB_lag, idxB_range_1, idxB_range_2,
                  ...) {
  
  ### index ratio
  yrs_a <- seq(to = c(ay - idxB_lag), length.out = idxB_range_1)
  yrs_b <- seq(to = min(yrs_a) - 1, length.out = idxB_range_2)
  idx_a <- yearMeans(idx[, ac(yrs_a)])
  idx_b <- yearMeans(idx[, ac(yrs_b)])
  idx_ratio <- c(idx_a / idx_b)
  
  return(idx_ratio)
  
}

### length data
est_f <- function(idx, ay, 
                  Lref, idxL_range, idxL_lag,
                  ...) {
  
  ### if fewer iterations provided expand
  if (isTRUE(length(Lref) < dims(idx)$iter)) {
    Lref <- rep(Lref, dims(idx)$iter)
    ### if more iterations provided, subset
  } else if (isTRUE(length(Lref) > dims(idx)$iter)) {
    Lref <- Lref[an(dimnames(idx)$iter)]
  }
  
  ### get mean length in catch
  idx_yrs <- seq(to = ay - idxL_range, length.out = idxL_lag)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### length relative to reference
  idx_ratio <- c(idx_mean / Lref)
  ### avoid negative values
  idx_ratio <- ifelse(idx_ratio > 0, idx_ratio, 0)
  ### avoid NAs, happens if catch = 0
  idx_ratio <- ifelse(is.na(idx_ratio), 1, idx_ratio)
  return(idx_ratio)
}

### biomass index trend
est_b <- function(idx, ay, 
                  I_trigger, idxB_lag, idxB_range_3,
                  ...) {
  
  ### if fewer iterations provided expand
  if (isTRUE(length(I_trigger) < dims(idx)$iter)) {
    I_trigger <- rep(I_trigger, dims(idx)$iter)
  ### if more iterations provided, subset
  } else if (isTRUE(length(I_trigger) > dims(idx)$iter)) {
    I_trigger <- I_trigger[an(dimnames(idx)$iter)]
  }
  
  ### calculate index mean
  idx_yrs <- seq(to = ay - idxB_lag, length.out = idxB_range_3)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### ratio
  idx_ratio <- c(idx_mean / I_trigger)
  ### b is 1 or smaller
  idx_ratio <- ifelse(idx_ratio < 1, idx_ratio, 1)
  
  return(idx_ratio)
  
}

### biomass index trend
est_pa <- function(idx, ay, tracking, pa_size, pa_duration, idxB_lag,
                   ...) {
  
  ### find last year in which buffer was applied
  last <- apply(tracking[[1]]["comp_b",,, drop = FALSE], 6, FUN = function(x) {#browser()
    ### positions (years) where buffer was applied
    yr <- dimnames(x)$year[which(x < 1)]
    ### return -Inf if buffer was never applied
    ifelse(length(yr) > 0, as.numeric(yr), -Inf)
  })
  ### find iterations to check 
  pos_check <- which(last <= (ay - pa_duration))
  ### find negative stock status (SSB<0.5Bmsy or F>Fmsy)
  pos_negative <- which(idx[, ac(ay - idxB_lag)] == 0)
  ### apply only if buffer applications need to be checked and status is negative
  pos_apply <- intersect(pos_check, pos_negative)
  
  return(ifelse(seq(dims(last)$iter) %in% pos_apply, pa_size, 1))
  
}

est_pa_conditional <- function(Lref, ### reference length (LF=M)
                               Lref_mult = 1, ### multiplier for Lref
                               idx, ### length index
                               ay, ### current (intermediate) year
                               idxL_range = 1, ### number of years for length
                               idxL_lag = 1, ### time lag relative to ay
                               tracking, ### track metrics
                               pa_duration, ### frequency of buffer application
                               pa_size, 
                               ...) {
  
  ### 1st: compare mean catch length to reference -> PA status
  
  ### if fewer iterations provided expand
  if (isTRUE(length(Lref) < dims(idx)$iter)) {
    Lref <- rep(Lref, dims(idx)$iter)
    ### if more iterations provided, subset
  } else if (isTRUE(length(Lref) > dims(idx)$iter)) {
    Lref <- Lref[an(dimnames(idx)$iter)]
  }
  
  ### get mean length in catch
  idx_yrs <- seq(to = ay - idxL_range, length.out = idxL_lag)
  idx_mean <- yearMeans(idx[, ac(idx_yrs)])
  ### length relative to reference
  idx_ratio <- c(idx_mean / (Lref * Lref_mult))
  ### avoid negative values
  idx_ratio <- ifelse(idx_ratio > 0, idx_ratio, 0)
  ### avoid NAs, happens if catch = 0
  idx_ratio <- ifelse(is.na(idx_ratio), 1, idx_ratio)
  
  ### buffer
  
  ### find last year in which buffer was applied
  last <- apply(tracking[[1]]["comp_b",,, drop = FALSE], 6, FUN = function(x) {#browser()
    ### positions (years) where buffer was applied
    yr <- dimnames(x)$year[which(x < 1)]
    ### return -Inf if buffer was never applied
    ifelse(length(yr) > 0, as.numeric(yr), -Inf)
  })
  ### find iterations to check 
  pos_check <- which(last <= (ay - pa_duration))
  ### when to apply buffer - if Lmean < LF=M == idx_ratio <= 0
  pos_negative <- which(idx_ratio <= 1)
  ### apply only if buffer applications need to be checked and status is negative
  pos_apply <- intersect(pos_check, pos_negative)
  
  ### initialise buffer value
  res <- rep(1, dim(tracking[[1]])[6])
  res[pos_apply] <- pa_size
  
  return(res)
  
}

### index value
est_i <- function(idx, ay,
                  idxB_lag, idxB_range_3,
                  ...) {
  
  ### index ratio
  yrs_r <- seq(to = c(ay - idxB_lag), length.out = idxB_range_3)
  idx_i <- yearMeans(idx[, ac(yrs_r)])
  
  return(idx_i)
  
}

### recent catch
est_c <- function(catch, ay,
                  catch_lag, catch_range,
                  ...) {
  
  catch_yrs <- seq(to = ay - catch_lag, length.out = catch_range)
  catch_current <- yearMeans(catch[, ac(catch_yrs)])
  return(catch_current)
  
}
### recent advice
est_A <- function(catch, ay, iy,
                  catch_lag, catch_range, tracking,
                  ...) {
  
  catch_yrs <- seq(to = ay - catch_lag, length.out = catch_range)
  ### first year - use catch
  if (identical(ay, iy)) {
    catch_current <- yearMeans(catch[, ac(catch_yrs)])
  ### other years - use advice
  } else {
    ### use catch_yrs + 1 because advice value is stored in corresponding
    ### years in tracking object
    catch_current <- yearMeans(tracking[[1]]["isys", ac(catch_yrs + 1)])
  }
  
  return(catch_current)
  
}

### ------------------------------------------------------------------------ ###
### phcr ####
### ------------------------------------------------------------------------ ###
### parametrization of HCR

phcr_CL <- function(tracking, args,
                    multiplier = 1,
                    lambda_upper = 0.1, lambda_lower = 0.2,
                    gamma_lower = 0.2, gamma_upper = 0.1,
                    r_threshold = 0.05, l_threshold = 0.1,
                    Lref = NA, Lref_mult = 1,
                    ...) {
  
  ay <- args$ay
  
  ### get values from tracking
  hcrpars <- tracking[[1]][c("A_last",
                        "r_length", "r_catch", "length_average",
                        "A_last", "A_last", "A_last", "A_last", ### dummy values
                        "A_last", "A_last", "A_last", "A_last",
                        "A_last"), ac(ay)]
  dimnames(hcrpars)$metric[5:13] <- c("lambda_upper", "lambda_lower",
                                      "gamma_upper", "gamma_lower",
                                      "r_threshold", "l_threshold",
                                      "Lref", "Lref_mult", "multiplier")

  ### insert control rule parameters
  hcrpars["lambda_upper", ] <- lambda_upper
  hcrpars["lambda_lower", ] <- lambda_lower
  hcrpars["gamma_upper", ] <- gamma_upper
  hcrpars["gamma_lower", ] <- gamma_lower
  hcrpars["r_threshold", ] <- r_threshold
  hcrpars["l_threshold", ] <- l_threshold
  hcrpars["Lref", ] <- Lref
  hcrpars["Lref_mult", ] <- Lref_mult
  hcrpars["multiplier", ] <- multiplier
  
  ### if threshold is zero, add tiny number to avoid computational issues
  if (identical(r_threshold, 0)) 
    hcrpars["r_threshold", ] <- r_threshold + .Machine$double.eps
  if (identical(l_threshold, 0)) 
    hcrpars["l_threshold", ] <- l_threshold + .Machine$double.eps
  
  ### convert hcrpars into FLPar to make "goFish" happy...
  hcrpars <- FLPar(hcrpars)

  return(list(tracking = tracking, hcrpars = hcrpars))
  
}


phcr_comps <- function(tracking, args, 
                       exp_r = 1, exp_f = 1, exp_b = 1,
                       ...){
  
  ay <- args$ay
  
  hcrpars <- tracking[[1]][c("comp_r", "comp_f", "comp_b", "comp_i", 
                        "comp_hr", "comp_c", "multiplier",
                        "exp_r", "exp_f", "exp_b"), ac(ay)]
  hcrpars["exp_r", ] <- exp_r
  hcrpars["exp_f", ] <- exp_f
  hcrpars["exp_b", ] <- exp_b
  
  if (exp_r != 1) tracking[[1]]["exp_r", ] <- exp_r
  if (exp_f != 1) tracking[[1]]["exp_f", ] <- exp_f
  if (exp_b != 1) tracking[[1]]["exp_b", ] <- exp_b
  
  ### convert hcrpars into FLPar to make "goFish" happy...
  hcrpars <- FLPar(hcrpars)
  
  ### return results
  return(list(tracking = tracking, hcrpars = hcrpars))
  
}

### ------------------------------------------------------------------------ ###
### hcr ####
### ------------------------------------------------------------------------ ###
### apply catch rule

hcr_CL <- function(hcrpars, args, tracking, interval = 2, 
                   combine_alpha_beta = FALSE,
                   ...) {
  
  ay <- args$ay ### current year
  iy <- args$iy ### first simulation year
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
    
    ### alpha - trends from catch and length data
    
    ### split into groups:
    ### -1: value < -0.01
    ###  0: -0.01 < value < 0.01
    ### +1: value > 0.01
    r_length <- cut(c(hcrpars["r_length"]), 
                    breaks = c(-Inf, -unique(c(hcrpars["r_threshold"])), 
                               unique(c(hcrpars["r_threshold"])), Inf), 
                    labels = c(-1, 0, 1))
    if (any(is.na(r_length))) r_length[is.na(r_length)] <- 0
    r_length <- as.numeric(as.character(r_length))
    r_catch <- cut(c(hcrpars["r_catch"]), 
                   breaks = c(-Inf, -unique(c(hcrpars["r_threshold"])),
                              unique(c(hcrpars["r_threshold"])), Inf), 
                   labels = c(-1, 0, 1))
    r_catch <- as.numeric(as.character(r_catch))
    if (any(is.na(r_catch))) r_catch[is.na(r_catch)] <- 0
    
    ### status: combine catch and length trend
    r_status <- r_length + r_catch
    
    ### assign values:
    ### -2: 1 - 2*lambda_lower
    ### -1: 1 - 1*lambda_lower
    ###  0: 1
    ###  1: 1 + 1*lambda_upper
    ###  2: 1 + 2*lambda_upper
    alpha <- sapply(as.character(r_status), function(x) {
      switch(x,
             "-2" = 1 - 2*(c(hcrpars["lambda_lower"])[1]),
             "-1" = 1 - 1*(c(hcrpars["lambda_lower"])[1]),
             "0" = 1,
             "1" = 1 + 1*(c(hcrpars["lambda_upper"])[1]),
             "2" = 1 + 2*(c(hcrpars["lambda_upper"])[1])
      )
    })
    
    ### beta - status evaluation from average catch length
    ### (relative to reference length)
    length_status <- cut(c(hcrpars["length_average"]/
                             (hcrpars["Lref"]*hcrpars["Lref_mult"])), 
                         breaks = c(-Inf, 1 - unique(c(hcrpars["l_threshold"])), 
                                    1 + unique(c(hcrpars["l_threshold"])), Inf), 
                         labels = c("negative", "neutral", "positive"))
    
    if (any(is.na(length_status))) 
      length_status[is.na(length_status)] <- "neutral"
    ### assign values:
    ### -2: 1 - 2*lambda_lower
    ### -1: 1 - 1*lambda_lower
    ###  0: 1
    ###  1: 1 + 1*lambda_upper
    ###  2: 1 + 2*lambda_upper
    beta <- sapply(as.character(length_status), function(x) {
      switch(x,
             "negative" = 1 - 1*(c(hcrpars["gamma_lower"])[1]),
             "neutral" = 1,
             "positive" = 1 + 1*(c(hcrpars["gamma_upper"])[1])
      )
    })
    
    ### calculate advice

    if (!isTRUE(combine_alpha_beta)) {
      ### only apply alpha if beta != 1  
      alpha[which(beta != 1)] <- 1
    }
    advice <- hcrpars["A_last"] * alpha * beta * hcrpars["multiplier"]
    
  } else {
    
    ### use last year's advice
    advice <- tracking[[1]]["hcr", ac(ay - 1)]
    
  }
  
  ctrl <- fwdControl(FLQuant(advice, dimnames = list(year = ay + 1)), 
                     quant = "catch")
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}


hcr_comps <- function(hcrpars, args, tracking, interval = 2, 
                  ...) {
  
  ay <- args$ay ### current year
  iy <- args$iy ### first simulation year
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
  
    ### calculate advice
    advice <- hcrpars["comp_c", ] *
                (hcrpars["comp_r", ]^hcrpars["exp_r", ]) *
                (hcrpars["comp_f", ]^hcrpars["exp_f", ]) *
                (hcrpars["comp_b", ]^hcrpars["exp_b", ]) *
                 hcrpars["comp_i", ] *
                 hcrpars["comp_hr", ] *
                 hcrpars["multiplier", ] 
    #advice <- apply(X = hcrpars, MARGIN = 6, prod, na.rm = TRUE)
    
  } else {
    
    ### use last year's advice
    advice <- tracking[[1]]["hcr", ac(ay - 1)]
    
  }

  ctrl <- fwdControl(FLQuant(advice, dimnames = list(year = ay + 1)), 
                     quant = "catch")
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### implementation ####
### ------------------------------------------------------------------------ ###
### no need to convert, already catch in tonnes
### apply TAC constraint, if required

is_comps <- function(ctrl, args, tracking, interval = 2, 
                     upper_constraint = Inf, lower_constraint = 0, 
                     cap_below_b = TRUE, ...) {
  
  ay <- args$ay ### current year
  iy <- args$iy ### first simulation year
  
  advice <- as.vector(ctrl@iters[, "value", ])
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
  
    ### apply TAC constraint, if requested
    if (!is.infinite(upper_constraint) | lower_constraint != 0) {
      
      ### get last advice
      if (isTRUE(ay == iy)) {
        ### use OM value in first year of projection
        adv_last <- tracking[[1]]["C.om", ac(iy - 1)]
      } else {
        adv_last <- tracking[[1]]["isys", ac(ay)]
      }
      ### ratio of new advice/last advice
      adv_ratio <- advice/adv_last
      
      ### upper constraint
      if (!is.infinite(upper_constraint)) {
        ### find positions
        pos_upper <- which(adv_ratio > upper_constraint)
        ### turn of constraint when index below Itrigger?
        if (isFALSE(cap_below_b)) {
          pos_upper <- setdiff(pos_upper, 
                               which(c(tracking[[1]][, ac(ay)]["comp_b", ]) < 1))
        }
        ### limit advice
        if (length(pos_upper) > 0) {
          advice[pos_upper] <- adv_last[,,,,, pos_upper] * upper_constraint
        }
        ### lower constraint
      }
      if (lower_constraint != 0) {
        ### find positions
        pos_lower <- which(adv_ratio < lower_constraint)
        ### turn of constraint when index below Itrigger?
        if (isFALSE(cap_below_b)) {
          pos_lower <- setdiff(pos_lower, 
                               which(c(tracking[[1]][, ac(ay)]["comp_b", ]) < 1))
        }
        ### limit advice
        if (length(pos_lower) > 0) {
          advice[pos_lower] <- adv_last[,,,,, pos_lower] * lower_constraint
        }
      }
    }
    
  ### otherwise do nothing here and recycle last year's advice
  } else {
    
    advice <- tracking[[1]]["isys", ac(ay)]
    
  }
  
  ### update advice values in control object
  ctrl@iters[, "value", ] <- advice
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### implementation error ####
### ------------------------------------------------------------------------ ###

iem_comps <- function(ctrl, args, tracking, 
                      iem_dev = FALSE, use_dev, ...) {
  
  ay <- args$ay
  
  ### only do something if requested
  if (isTRUE(use_dev)) {
    
    ### get advice
    advice <- ctrl@iters[, "value", ]
    ### get deviation
    dev <- c(iem_dev[, ac(ay)])
    ### implement deviation
    advice <- advice * dev
    ### insert into ctrl object
    ctrl@iters[, "value", ] <- advice
    
  }
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### projection ####
### ------------------------------------------------------------------------ ###
fwd_attr <- function(om, ### includes stock, recruitment model
                     deviances, 
                     ctrl,
                     maxF = 5, ### maximum allowed Fbar
                     ...) {
  
  ### project forward with FLasher::fwd
  om@stock[] <- fwd(object = om@stock, control = ctrl, sr = om@sr, 
                    residuals = deviances, maxF = maxF)
  
  ### return stock
  return(list(om = om))
  
}


### ------------------------------------------------------------------------ ###
### iter subset  ####
### ------------------------------------------------------------------------ ###

iter_attr <- function(object, iters, subset_attributes = TRUE) {
  
  ### subset object to iter
  res <- FLCore::iter(object, iters)
  
  if (isTRUE(subset_attributes)) {
    
    ### get default attributes of object class
    attr_def <- names(attributes(new(Class = class(object))))
    
    ### get additional attributes
    attr_new <- setdiff(names(attributes(object)), attr_def)
    
    ### subset attributes
    for (attr_i in attr_new) {
      attr(res, attr_i) <- FLCore::iter(attr(res, attr_i), iters)
    }
    
  }
  
  return(res)
  
}

### ------------------------------------------------------------------------ ###
### estimtate steepness based on l50/linf ratio ####
### according to Wiff et al. 2018
### ------------------------------------------------------------------------ ###
h_Wiff <- function(l50, linf) {
  l50linf <- l50/linf
  ### linear model
  lin <- 2.706 - 3.698*l50linf
  ### logit
  h <- (0.2 + exp(lin)) / (1 + exp(lin))
  return(h)
}

### ------------------------------------------------------------------------ ###
### mean length in catch ####
### ------------------------------------------------------------------------ ###
lmean <- function(stk, params) {
  
  ### calculate length from age with a & b
  weights <- c(catch.wt(stk)[, 1,,,, 1])
  lengths <- (weights / c(params["a"]))^(1 / c(params["b"]))
  catch.n <- catch.n(stk)
  dimnames(catch.n)$age <- lengths
  ### subset to lengths > Lc
  catch.n <- catch.n[lengths > c(params["Lc"]),]
  
  ### calculate mean length
  lmean <- apply(X = catch.n, MARGIN = c(2, 6), FUN = function(x) {
    ### calculate
    res <- weighted.mean(x = an(dimnames(x)$age), 
                         w = ifelse(is.na(x), 0, x), na.rm = TRUE)
    ### check if result obtained
    ### if all catch at all lengths = 0, return 0 as mean length
    # if (is.nan(res)) {
    #   if (all(ifelse(is.na(x), 0, x) == 0)) {
    #     res[] <- 0
    #   }
    # }
    return(res)
  })
  return(lmean)
}

### ------------------------------------------------------------------------ ###
### length at first capture ####
### ------------------------------------------------------------------------ ###
calc_lc <- function(stk, a, b) {
  ### find position in age vector
  Ac <- apply(catch.n(stk), MARGIN = c(2, 6), function(x) {
    head(which(x >= (max(x, na.rm = TRUE)/2)), 1)
  })
  Ac <- an(median(Ac))
  ### calculate lengths
  weights <- c(catch.wt(stk)[, 1,,,, 1])
  lengths <- (weights / a)^(1 / b)
  ### length at Ac
  Lc <- floor(lengths[Ac]*10)/10
  return(Lc)
}

### ------------------------------------------------------------------------ ###
### inter-annual variability ####
### ------------------------------------------------------------------------ ###
#' calculate inter-annual variability of FLQuant
#'
#' This function calculates survey indices from the numbers at age of an 
#' FLStock object
#'
#' @param object Object of class \linkS4class{FLQuant} with values.
#' @param period Select every n-th year, e.g. biennial (optional).
#' @param from,to Optional year range for analysis.
#' @param summary_per_iter Function for summarising per iter. Defaults to mean.
#' @param summary Function for summarising over iter. Defaults to mean.
#' @return An object of class \code{FLQuant} with inter-annual variability.
#'
#' @export
#' 
setGeneric("iav", function(object, period, from, to, summary_per_iter, 
                           summary_year, summary_all) {
  standardGeneric("iav")
})

### object = FLQuant
#' @rdname iav
setMethod(f = "iav",
  signature = signature(object = "FLQuant"),
  definition = function(object, 
                        period, ### periodicity, e.g. use every 2nd value 
                        from, to,### year range
                        summary_per_iter, ### summarise values per iteration
                        summary_year,
                        summary_all) {
            
  ### subset years
  if (!missing(from)) object <- FLCore::window(object, start = from)
  if (!missing(to)) object <- FLCore::window(object, end = from)
  
  ### get years in object
  yrs <- dimnames(object)$year
  
  ### select every n-th value, if requested
  if (!missing(period)) {
    yrs <- yrs[seq(from = 1, to = length(yrs), by = period)]
  }
  
  ### reference years
  yrs_ref <- yrs[-length(yrs)]
  ### years to compare
  yrs_comp <- yrs[-1]
  
  ### calculate variation (absolute values, ignore pos/neg)
  res <- abs(1 - object[, yrs_comp] / object[, yrs_ref])
  
  ### replace Inf with NA (compared to 0 catch)
  res <- ifelse(is.finite(res), res, NA)
  
  ### summarise per iteration
  if (!missing(summary_per_iter)) {
    res <- apply(res, 6, summary_per_iter, na.rm = TRUE)
  }
  
  ### summarise per year
  if (!missing(summary_year)) {
    res <- apply(res, 1:5, summary_year, na.rm = TRUE)
  }
  
  ### summarise over everything
  if (!missing(summary_all)) {
    
    res <- summary_all(c(res), na.rm = TRUE)
    
  }
  
  return(res)
  
})


### ------------------------------------------------------------------------ ###
### "correct" collapses ####
### ------------------------------------------------------------------------ ###

collapse_correction <- function(stk, quants = c("catch", "ssb", "fbar"),
                                threshold = 1, yrs) {
  names(quants) <- quants
  qnt_list <- lapply(quants, function(x) get(x)(stk))
  qnt_list <- lapply(qnt_list, function(x) x[, ac(yrs)])
  
  n_yrs <- dim(qnt_list[[1]])[2]
  n_its <- dim(qnt_list[[1]])[6]
  
  ### find collapses
  cd <- sapply(seq(n_its), function(x) {
    min_yr <- min(which(qnt_list$ssb[,,,,, x] < 1))
    if (is.finite(min_yr)) {
      all_yrs <- min_yr:n_yrs
    } else {
      all_yrs <- NA
    }
    all_yrs + (x - 1)*n_yrs
  })
  cd <- unlist(cd)
  cd <- cd[which(!is.na(cd))]
  ### remove values
  qnt_list <- lapply(qnt_list, function(x) {
    x@.Data[cd] <- 0
    return(x)
  })
  return(qnt_list)
}


