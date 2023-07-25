### ------------------------------------------------------------------------ ###
### observations ####
### ------------------------------------------------------------------------ ###
obs_generic <- function(stk, observations, deviances, args, tracking,
                        ssb_idx = FALSE, tsb_idx = FALSE, ### use SSB idx
                        idx_timing = FALSE, ### consider survey timing?
                        idx_dev = FALSE,
                        lngth = FALSE, ### catch length data?
                        lngth_dev = FALSE, 
                        lngth_par,
                        PA_status = FALSE,
                        PA_status_dev = FALSE,
                        PA_Bmsy = FALSE, PA_Fmsy = FALSE,
                        ...) {

  #ay <- args$ay
  ### update observations
  observations$stk <- stk
  ### use SSB as index?
  if (isTRUE(ssb_idx)) {
    observations$idx$idxB <- ssb(observations$stk)
  ### TSB?
  } else  if (isTRUE(tsb_idx)) {
      observations$idx$idxB <- tsb(observations$stk)
  ### otherwise calculate biomass index
  } else {
    sn <- stk@stock.n
    ### reduce by F and M?
    if (isTRUE(idx_timing)) {
      sn <- sn * exp(-(harvest(stk) * harvest.spwn(stk) +
                         m(stk) * m.spwn(stk)))
    }
    observations$idx$idxB <- quantSums(sn * stk@stock.wt * 
                                       observations$idx$sel)
  }
  ### use mean length in catch?
  if (isTRUE(lngth)) {
    observations$idx$idxL <- lmean(stk = stk, params = lngth_par)
  }
  ### stock status for PA buffer?
  if (isTRUE(PA_status)) {
    observations$idx$PA_status[] <- ssb(observations$stk) > 0.5*PA_Bmsy & 
                                    fbar(observations$stk) < PA_Fmsy
  }
  
  ### observation model
  stk0 <- observations$stk
  idx0 <- observations$idx
  ### add deviances to index?
  if (isTRUE(idx_dev)) {
    if (isTRUE(ssb_idx) | isTRUE(tsb_idx)) {
      idx0$idxB <- observations$idx$idxB * deviances$idx$idxB
    } else {
      idx0$idxB <- quantSums(stk@stock.n * stk@stock.wt * 
                             observations$idx$sel * deviances$idx$sel)
      if (isTRUE("idxB" %in% names(deviances$idx)) & 
          all.equal(dim(deviances$idx$idxB), dim(idx0$idxB)))
        idx0$idxB <- idx0$idxB * deviances$idx$idxB
    }
  }
  ### uncertainty for catch length
  if (isTRUE(lngth) & isTRUE(lngth_dev)) {
    idx0$idxL <- observations$idx$idxL * deviances$idx$idxL
  }
  ### uncertainty for stock status for PA buffer
  if (isTRUE(PA_status) & isTRUE(PA_status_dev)) {
    idx0$PA_status <- ifelse(observations$idx$PA_status == TRUE, 
                             deviances$idx$PA_status["positive", ],
                             deviances$idx$PA_status["negative", ])
  }
  
  
  return(list(stk = stk0, idx = idx0, observations = observations,
              tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### estimator ####
### ------------------------------------------------------------------------ ###

est_comps <- function(stk, idx, tracking, args,
                      comp_r = FALSE, comp_f = FALSE, comp_b = FALSE,
                      comp_i = FALSE, comp_c = FALSE, comp_A = TRUE,
                      comp_m = FALSE, comp_hr = FALSE,
                      idxB_lag = 1, idxB_range_1 = 2, idxB_range_2 = 3,
                      idxB_range_3 = 1,
                      catch_lag = 1, catch_range = 1,
                      Lref, I_trigger,
                      idxL_lag = 1, idxL_range = 1,
                      pa_buffer = FALSE, pa_size = 0.8, pa_duration = 3,
                      Bmsy = NA,
                      ...) {
  
  ay <- args$ay
  iy <- args$iy
  
  ### component r: index trend
  if (isTRUE(comp_r)) {
    r_res <- est_r(idx = idx$idxB, ay = ay,
                   idxB_lag = idxB_lag, idxB_range_1 = idxB_range_1, 
                   idxB_range_2 = idxB_range_2)
  } else {
    r_res <- 1
  }
  tracking["comp_r", ac(ay)] <- r_res
  
  ### component f: length data
  if (isTRUE(comp_f)) {
    f_res <- est_f(idx = idx$idxL, ay = ay,
                   Lref = Lref, idxL_range = idxL_range, idxL_lag = idxL_lag)
  } else {
    f_res <- 1
  }
  tracking["comp_f", ac(ay)] <- f_res
  
  ### component b: biomass safeguard
  if (isTRUE(comp_b)) {
    b_res <- est_b(idx = idx$idxB, ay = ay,
                   I_trigger = I_trigger, idxB_lag = idxB_lag, 
                   idxB_range_3 = idxB_range_3)
  } else {
    b_res <- 1
  }
  
  ### PA buffer
  if (isTRUE(pa_buffer)) {
    b_res <- est_pa(idx = idx$PA_status, ay = ay, 
                    tracking = tracking, idxB_lag = idxB_lag,
                    pa_size = pa_size, pa_duration = pa_duration)
  }
  tracking["comp_b", ac(ay)] <- b_res
  
  ### component i: index value
  if (isTRUE(comp_i)) {
    i_res <- est_i(idx = idx$idxB, ay = ay,
                   idxB_lag = idxB_lag, idxB_range_3 = idxB_range_3)
  } else {
    i_res <- 1
  }
  tracking["comp_i", ac(ay)] <- i_res
  
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
  tracking["comp_c", ac(ay)] <- c_res
  
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
  tracking["multiplier", ac(ay)] <- m_res
  
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
  tracking["comp_hr", ac(ay)] <- hr_res
  
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
  last <- apply(tracking["comp_b",,, drop = FALSE], 6, FUN = function(x) {#browser()
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
    ### yrs - 1 because advice for ay+1 is saved in ay in tracking object
    catch_current <- yearMeans(tracking["metric.is", ac(catch_yrs - 1)])
  }
  
  return(catch_current)
  
}

### ------------------------------------------------------------------------ ###
### phcr ####
### ------------------------------------------------------------------------ ###
### parametrization of HCR

phcr_comps <- function(tracking, args, 
                       exp_r = 1, exp_f = 1, exp_b = 1,
                       ...){
  
  ay <- args$ay
  
  hcrpars <- tracking[c("comp_r", "comp_f", "comp_b", "comp_i", 
                        "comp_hr", "comp_c", "multiplier",
                        "exp_r", "exp_f", "exp_b"), ac(ay)]
  hcrpars["exp_r", ] <- exp_r
  hcrpars["exp_f", ] <- exp_f
  hcrpars["exp_b", ] <- exp_b
  
  if (exp_r != 1) tracking["exp_r", ] <- exp_r
  if (exp_f != 1) tracking["exp_f", ] <- exp_f
  if (exp_b != 1) tracking["exp_b", ] <- exp_b
  
  ### return results
  return(list(tracking = tracking, hcrpars = hcrpars))
  
}

### ------------------------------------------------------------------------ ###
### hcr ####
### ------------------------------------------------------------------------ ###
### apply catch rule

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
                 hcrpars["comp_i"] *
                 hcrpars["comp_hr"] *
                 hcrpars["multiplier", ] 
    #advice <- apply(X = hcrpars, MARGIN = 6, prod, na.rm = TRUE)
    
  } else {
    
    ### use last year's advice
    advice <- tracking["metric.hcr", ac(ay - 1)]
    
  }

  ctrl <- getCtrl(values = c(advice), quantity = "catch", years = ay + 1, 
                  it = dim(advice)[6])
  
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
  
  advice <- ctrl@trgtArray[ac(ay + args$management_lag), "val", ]
  
  ### check if new advice requested
  if ((ay - iy) %% interval == 0) {
  
    ### apply TAC constraint, if requested
    if (!is.infinite(upper_constraint) | lower_constraint != 0) {
      
      ### get last advice
      if (isTRUE(ay == iy)) {
        ### use OM value in first year of projection
        adv_last <- tracking["C.om", ac(iy)]
      } else {
        adv_last <- tracking["metric.is", ac(ay - 1)]
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
                               which(c(tracking[, ac(ay)]["comp_b", ]) < 1))
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
                               which(c(tracking[, ac(ay)]["comp_b", ]) < 1))
        }
        ### limit advice
        if (length(pos_lower) > 0) {
          advice[pos_lower] <- adv_last[,,,,, pos_lower] * lower_constraint
        }
      }
    }
    
  ### otherwise do nothing here and recycle last year's advice
  } else {
    
    advice <- tracking["metric.is", ac(ay - 1)]
    
  }
  ctrl@trgtArray[ac(ay + args$management_lag),"val",] <- advice
  
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
    advice <- ctrl@trgtArray[ac(ay + args$management_lag), "val", ]
    ### get deviation
    dev <- c(iem_dev[, ac(ay)])
    ### implement deviation
    advice <- advice * dev
    ### insert into ctrl object
    ctrl@trgtArray[ac(ay + args$management_lag),"val",] <- advice
    
  }
  
  return(list(ctrl = ctrl, tracking = tracking))
  
}

### ------------------------------------------------------------------------ ###
### projection ####
### ------------------------------------------------------------------------ ###
fwd_attr <- function(stk, ctrl,
                     sr, ### stock recruitment model
                     sr.residuals, ### recruitment residuals
                     sr.residuals.mult = TRUE, ### are res multiplicative?
                     maxF = 5, ### maximum allowed Fbar
                     dupl_trgt = FALSE,
                     ...) {
  
  ### avoid the issue that the catch is higher than the targeted catch
  ### can happen due to bug in FLash if >1 iteration provided
  ### sometimes, FLash struggles to get estimates and then uses F estimate from
  ### previous iteration
  ### workaround: target same value several times and force FLash to try again
  if (isTRUE(dupl_trgt)) {

    ### duplicate target
    ctrl@target <- rbind(ctrl@target, ctrl@target, ctrl@target)
    ### replace catch in second row with landings
    ctrl@target$quantity[1] <- "landings"
    ctrl@target$quantity[3] <- "catch"

    ### extract target values
    val_temp <- ctrl@trgtArray[, "val", ]

    ### extend trgtArray
    ### extract dim and dimnames
    dim_temp <- dim(ctrl@trgtArray)
    dimnames_temp <- dimnames(ctrl@trgtArray)
    ### duplicate years
    dim_temp[["year"]] <- dim_temp[["year"]] * 3
    dimnames_temp$year <- rep(dimnames_temp$year, 3)

    ### create new empty array
    trgtArray <- array(data = NA, dim = dim_temp, dimnames = dimnames_temp)

    ### fill with values
    ### first as target
    trgtArray[1, "val", ] <- val_temp
    ### then again, but as max
    trgtArray[2, "max", ] <- val_temp
    ### min F
    trgtArray[3, "max", ] <- val_temp

    ### insert into ctrl object
    ctrl@trgtArray <- trgtArray
  }
  
  ### project forward with FLash::fwd
  stk[] <- fwd(object = stk, control = ctrl, sr = sr, 
               sr.residuals = sr.residuals, 
               sr.residuals.mult = sr.residuals.mult,
               maxF = maxF)
  
  ### return stock
  return(list(object = stk))
  
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


