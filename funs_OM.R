
### ------------------------------------------------------------------------ ###
### function for preparing OM input for MP ####
### ------------------------------------------------------------------------ ###

input_mp <- function(stocks, 
                     fhist, ### fishing history, e.g. "one-way"
                     n_iter = 500, ### number of iterations
                     n_yrs = 50, ### number of years for MP
                     MP, ### e.g. "hr", "const_catch", "CC_f"
                     hist_yr_min = 50, ### first year
                     scenario = "",
                     idx_sel = "tsb", ### index selectivity
                     n_blocks = 1, 
                     ### index uncertainty & auto-correlation
                     sigmaB = 0.2,
                     sigmaL = 0.2,
                     sigmaB_rho = 0,
                     sigmaL_rho = 0,
                     ### recruitment variability
                     sigmaR = 0.6,
                     sigmaR_rho = 0,
                     ### recruitment steepness
                     steepness = 0.75,
                     ### HR definition
                     hr_value = "length",
                     ### est
                     multiplier = 1,
                     comp_b_multiplier = 1.4,
                     idxB_lag = 1,
                     idxB_range_1 = 2, idxB_range_2 = 3, ### 2 over 3 ratio
                     idxB_range_3 = 1, ### for b & i
                     catch_lag = 1, catch_range = 1,
                     idxL_lag = 1, idxL_range = 1,
                     pa_buffer = FALSE, pa_size = 0.8, pa_duration = interval,
                     Lref_mult = 1,
                     ### phcr
                     exp_r = 1, exp_f = 1, exp_b = 1, ### exponents (rfb only)
                     ### hcr
                     interval = 1, ### TAC interval
                     ### is
                     upper_constraint = 1.2, ### uncertainty cap
                     lower_constraint = 0.7, 
                     cap_below_b = FALSE ### conditional cap?
                     ### iem
                     
                     ) {
  
  ### load life-history parameters
  brps <- readRDS("input/brps.rds")
  stocks_lh <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
  
  ### stock ID
  if (is.numeric(stocks)) stocks <- stocks_lh$stock[stocks]
  
  ### go through list of stocks
  input_list <- lapply(stocks, function(stock) {  
    
    lhist <- stocks_lh[stocks_lh$stock == stock, ]
    ### get reference points
    refpts <- refpts(brps[[stock]])
    Blim <- attr(brps[[stock]], "Blim")
    
    ### load stock & sr
    stk <- readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM/", fhist, "/",
                          stock, "/stk.rds"))
    sr <- readRDS(paste0("input/", n_iter, "_", n_yrs, "/OM/", fhist, "/",
                          stock, "/sr.rds"))
    
    ### adapt recruitment variability if needed
    if (isTRUE(!identical(sigmaR, 0.6) | !identical(sigmaR_rho, 0))) {
      ### create recruitment residuals for projection period
      set.seed(0)
      residuals(sr) <- rlnoise(n_iter, rec(stk) %=% 0, 
                               sd = sigmaR, b = sigmaR_rho)
      ### replicate residuals from GA paper
      set.seed(0)
      residuals(sr)[, ac(0:150)] <- rlnoise(n_iter, rec(stk)[, ac(0:150)] %=% 0,
                                            sd = sigmaR, b = sigmaR_rho)
      ### replicate residuals from catch rule paper for historical period
      set.seed(0)
      residuals(sr)[, ac(1:100)] <- rlnoise(n_iter, (rec(stk) %=% 0)[, ac(1:100)], 
                                            sd = sigmaR, b = sigmaR_rho)
      
    }
    
    ### cut off history
    stk <- window(stk, start = hist_yr_min)
    sr@residuals <- window(sr@residuals, start = hist_yr_min)
    
    ### adapt recruitment steepness if needed
    if (!identical(steepness, 0.75)) {
      ### calculate new recruitment model parameters with new steepness
      alpha <- (4*steepness*c(refpts["virgin", "rec"])) /
        (5*steepness - 1)
      beta <- (c(refpts["virgin", "ssb"]) * (1 - steepness)) /
        (5*steepness - 1)
      ### insert values
      params(sr)[] <- c(alpha, beta)
    }
    
    ### operating model
    om <- FLom(stock = stk, ### stock 
               sr = sr, ### stock recruitment and precompiled residuals
               fleetBehaviour = mseCtrl(),
               projection = mseCtrl(method = fwd_attr,
                                    args = list(dupl_trgt = TRUE)))
  
    
    ### length data
    if (MP %in% c("rfb", "hr", "CC_f")) {
      pars_l <- FLPar(a = lhist$a,
                      b = lhist$b,
                      Lc = calc_lc(stk = stk[, ac(75:100)], 
                                   a = lhist$a, b = lhist$b))
    }
    
    ### observation (error) model
    idx <- FLQuants(
      sel = stk@mat %=% NA_real_,
      idxB = ssb(stk) %=% NA_real_,
      idxL = ssb(stk) %=% NA_real_,
      PA_status = ssb(stk) %=% NA_integer_)
    if (identical(MP, "hr")) {
      oem <- FLoem(method = obs_generic,
                   observations = list(stk = stk, idx = idx), 
                   deviances = list(stk = FLQuant(), idx = idx),
                   args = list(idx_dev = TRUE, ssb_idx = FALSE, tsb_idx = FALSE,
                               lngth = FALSE, lngth_dev = FALSE,
                               lngth_par = pars_l,
                               PA_status = FALSE, PA_status_dev = FALSE,
                               PA_Bmsy = c(refpts(brps[[stock]])["msy", "ssb"]), 
                               PA_Fmsy = c(refpts(brps[[stock]])["msy", "harvest"])))
    } else if (identical(MP, "const_catch")) {
      oem <- FLoem(method = obs_generic,
                   observations = list(stk = stk, idx = idx), 
                   deviances = list(stk = FLQuant(), idx = idx),
                   args = list(idx_dev = TRUE, ssb_idx = FALSE, tsb_idx = FALSE,
                               lngth = FALSE, lngth_dev = FALSE,
                               #lngth_par = FLPar(),
                               PA_status = TRUE, PA_status_dev = FALSE,
                               PA_Bmsy = Inf, ### i.e. B always below Bmsy/2
                               PA_Fmsy = 0 ### i.e. F always above Fmsy
                               ### -> always apply buffer
                               ))
    } else if (identical(MP, "CC_f")) {
        oem <- FLoem(method = obs_generic,
                     observations = list(stk = stk, idx = idx), 
                     deviances = list(stk = FLQuant(), idx = idx),
                     args = list(idx_dev = FALSE, ssb_idx = FALSE, tsb_idx = FALSE,
                                 lngth = TRUE, lngth_dev = TRUE,
                                 lngth_par = pars_l,
                                 PA_status = TRUE, PA_status_dev = FALSE,
                                 PA_Bmsy = Inf, ### i.e. B always below Bmsy/2
                                 PA_Fmsy = 0 ### i.e. F always above Fmsy
                                 ### -> always apply buffer
                     ))
      
    }
    ### create biomass index
    ### update index if not total stock biomass
    if (identical(idx_sel, "ssb")) {
      oem@args$ssb_idx <- TRUE
      oem@args$tsb_idx <- FALSE
      oem@args$idx_timing <- TRUE
      oem@observations$idx$sel <- mat(stk)
    } else if (identical(idx_sel, "tsb")) {
      oem@args$ssb_idx <- FALSE
      oem@args$tsb_idx <- TRUE
      oem@args$idx_timing <- TRUE
      oem@observations$idx$sel <- mat(stk) %=% 1
    } else if (identical(idx_sel, "catch")) {
      oem@args$ssb_idx <- FALSE
      oem@args$tsb_idx <- FALSE
      oem@args$idx_timing <- TRUE
      ### estimate selectivity of catch (i.e. catch numbers/stock numbers)
      cn <- oem@observations$stk@catch.n
      sn <- oem@observations$stk@stock.n
      csel <- cn/sn
      csel_max <- csel
      for (age in seq(dim(csel)[1]))
        csel_max[age, ] <- csel[dim(csel)[1], ]
      csel <- csel/csel_max
      ### standardise for all years
      csel <- yearMeans(csel)
      oem@observations$idx$sel[] <- csel
    }  else if (identical(idx_sel, "dome_shaped")) {
      oem@args$ssb_idx <- FALSE
      oem@args$tsb_idx <- FALSE
      oem@args$idx_timing <- TRUE
      ### get life-history parameters
      if (is.na(lhist$t0)) lhist$t0 <- -0.1
      if (is.na(lhist$a50))
        lhist$a50 <- -log(1 - lhist$l50/lhist$linf)/lhist$k + lhist$t0
      ages <- as.numeric(dimnames(stk)$age)
      ### define selectivity (follow FLife's double normal function)
      sel_dn <- function(t, t1, sl, sr) {
        ifelse(t < t1, 2^(-((t - t1)/sl)^2), 2^(-((t - t1)/sr)^2))
      }
      sel <- sel_dn(t = ages, t1 = lhist$a50, sl = 1, sr = 10)
      oem@observations$idx$sel[] <- sel
    } else {
      ### standard scientific survey
      oem@args$ssb_idx <- FALSE
      oem@args$tsb_idx <- FALSE
      oem@args$idx_timing <- TRUE
      sel <- 1/(1 + exp(-1*(an(dimnames(stk)$age) - dims(stk)$max/10)))
      oem@observations$idx$sel[] <- sel
    }
    
    ### calculate biomass index with selectivity
    ### get stock numbers and reduce by F and M
    sn <- stock.n(stk) * exp(-(harvest(stk) * harvest.spwn(stk) + 
                                 m(stk) * m.spwn(stk)))
    ### calculate index
    idxB <- quantSums(sn * stock.wt(stk) * oem@observations$idx$sel)
    oem@observations$idx$idxB <- idxB
    
    
      
    ### length index
    ### only include when required for MP
    if (isTRUE(MP %in% c("rfb", "CC_f"))) {
      oem@args$lngth <- TRUE
      oem@args$lngth_dev <- TRUE
      oem@args$lngth_par <- pars_l
      ### calculate mean length
      oem@observations$idx$idxL <- lmean(stk = stk, params = pars_l)
    }
  
    ### index deviations
    ### PA buffer deviations (based on SPiCT performance)
    if (isTRUE(MP %in% c("2over3"))) {
      PA_status_dev <- FLQuant(NA, dimnames = list(age = c("positive", "negative"), 
                                                   year = dimnames(stk)$year, 
                                                   iter = dimnames(stk)$iter))
      set.seed(1)
      PA_status_dev["positive"] <- rbinom(n = PA_status_dev["positive"], 
                                          size = 1, prob = 0.9886215)
      set.seed(2)
      PA_status_dev["negative"] <- rbinom(n = PA_status_dev["negative"], 
                                          size = 1, prob = 1 - 0.4216946)
      oem@deviances$idx$PA_status <- PA_status_dev
    } else if (isTRUE(MP %in% c("const_catch"))) {
      PA_status_dev <- FLQuant(NA, dimnames = list(age = c("positive", "negative"), 
                                                   year = dimnames(stk)$year, 
                                                   iter = dimnames(stk)$iter))
      set.seed(1)
      PA_status_dev["positive"] <- 1
      set.seed(2)
      PA_status_dev["negative"] <- 1
      oem@deviances$idx$PA_status <- PA_status_dev
    } else if (identical(MP, "CC_f")) {
      PA_status_dev <- FLQuant(NA, dimnames = list(age = c("positive", "negative"), 
                                                   year = dimnames(stk)$year, 
                                                   iter = dimnames(stk)$iter))
      set.seed(1)
      PA_status_dev["positive"] <- 1
      set.seed(2)
      PA_status_dev["negative"] <- 1
      oem@deviances$idx$PA_status <- PA_status_dev
      
      set.seed(696)
      oem@deviances$idx$idxL <- rlnoise(n = dims(idx$idxL)$iter, idx$idxL %=% 0, 
                                        sd = 0.2, b = 0)
        
    }
    
    oem@deviances$idx$sel <- oem@deviances$idx$sel %=% 1
    
    ### biomass index deviations
    set.seed(695)
    if (isTRUE(MP %in% c("rfb", "2over3", "hr"))) {
      oem@deviances$idx$idxB <- rlnoise(n = dims(idx$idxB)$iter, idx$idxB %=% 0, 
                                        sd = sigmaB, b = sigmaB_rho)
    }
    if (isTRUE(MP %in% c("rfb", "hr"))) {
      oem@deviances$idx$idxL <- rlnoise(n = dims(idx$idxL)$iter, idx$idxL %=% 0, 
                                        sd = sigmaL, b = sigmaL_rho)
    }
    ### replicate previous deviates from GA paper
    set.seed(696)
    if (isTRUE(MP %in% c("rfb", "2over3", "hr"))) {
      oem@deviances$idx$idxB[, ac(50:150)] <- 
        rlnoise(n = dims(oem@deviances$idx$idxB)$iter,
                window(oem@deviances$idx$idxB, end = 150) %=% 0,
                sd = sigmaB, b = sigmaB_rho)
    }
    if (isTRUE(MP %in% c("rfb", "hr"))) {
      oem@deviances$idx$idxL[, ac(50:150)] <- 
        rlnoise(n = dims(oem@deviances$idx$idxL)$iter, 
                window(oem@deviances$idx$idxB, end = 150) %=% 0,
                sd = sigmaL, b = sigmaL_rho)
    }
    
    ### biomass index reference points
    ### lowest observed index in last 50 years
    I_loss <- apply(window(oem@observations$idx$idxB * 
                             oem@deviances$idx$idxB, end = 100),
                    6, min)
    I_trigger <- I_loss * comp_b_multiplier
  
    ### iem deviation
    set.seed(205)
    iem_dev <- FLQuant(rlnoise(n = dims(stk)$iter,  catch(stk) %=% 0,
                               sd = 0.1, b = 0))
    ### iem object
    iem <- FLiem(method = iem_comps,
                 args = list(use_dev = FALSE, iem_dev = iem_dev))
    
    
    ### harvest rate (chr)
    if (identical(MP, "hr")) {
      ### load some reference harvest rates 
      hr_refs <- readRDS("input/catch_rates.rds")[stock]
      ### define target harvest rate(s)
      if (identical(hr_value, "uniform")) {
        set.seed(33)
        hr_val <- runif(n = n_iter, min = 0, max = 1)
      } else if (identical(hr_value, "Fmsy")) {
        hr_val <- switch(idx_sel, 
                         "tsb" = hr_refs[[stock]]$Fmsy$tsb,
                         "ssb" = hr_refs[[stock]]$Fmsy$ssb,
                         "survey" = hr_refs[[stock]]$Fmsy$idx,
                         hr_refs[[stock]]$Fmsy$tsb)
      } else if (identical(hr_value, "LFeM")) {
        hr_val <- switch(idx_sel, 
                         "tsb" = hr_refs[[stock]]$LFeM$tsb,
                         "ssb" = hr_refs[[stock]]$LFeM$ssb,
                         "survey" = hr_refs[[stock]]$LFeM$idx,
                         hr_refs[[stock]]$LFeM$tsb)
      } else if (is.numeric(hr_value)) {
        hr_val <- hr_value
      } else if (identical(hr_value, "length")) {
        ### determine harvest rate based on historical mean catch length
        # stk <- FLCore::window(input$om@stock, end = 100)
        Lc <- calc_lc(stk = window(stk, end = 100), a = lhist$a, b = lhist$b)
        ### reference length
        LFeM <- (lhist$linf + 2*1.5*c(Lc)) / (1 + 2*1.5)
        ### mean catch length index (including noise)
        Lmean <- lmean(stk = stk, params = pars_l)
        idxL <- Lmean * oem@deviances$idx$idxL
        idxL <- window(idxL, end = 100)
        ### relative to reference length
        idxL <- idxL / LFeM
        ### biomass index
        idxB <- oem@observations$idx$idxB * oem@deviances$idx$idxB
        idxB <- window(idxB, end = 100)
        ### historical harvest rate
        CI <- window(catch(stk), end = 100)/idxB
        ### average of harvest rates where catch length is above reference
        hr_val <- sapply(dimnames(stk)$iter, function(i){
          mean(CI[, which(idxL[,,,,, i] >= 1),,,, i])
        })
      }
      
      ### set up MP ctrl object
      ctrl <- ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_comps,
                      args = list(comp_r = FALSE, 
                                  comp_f = FALSE, 
                                  comp_c = FALSE,
                                  comp_b = TRUE,
                                  comp_m = multiplier,
                                  pa_buffer = FALSE, 
                                  comp_i = TRUE, 
                                  I_trigger = I_trigger,
                                  idxB_lag = idxB_lag,
                                  idxB_range_3 = idxB_range_3,
                                  comp_hr = hr_val
                                  )),
        phcr = mseCtrl(method = phcr_comps,
                       args = list()),
        hcr = mseCtrl(method = hcr_comps,
                      args = list(interval = interval)),
        isys = mseCtrl(method = is_comps,
                       args = list(interval = interval, 
                                   upper_constraint = upper_constraint, 
                                   lower_constraint = lower_constraint, 
                                   cap_below_b = cap_below_b))
      ))
    
    } else if (identical(MP, "const_catch")) {
 
      ### set up MP ctrl object
      ctrl <- ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_comps,
                      args = list(comp_r = FALSE, 
                                  comp_f = FALSE, 
                                  comp_c = FALSE,
                                  comp_A = TRUE,
                                  comp_b = FALSE,
                                  comp_m = FALSE,
                                  pa_buffer = TRUE, 
                                  comp_i = FALSE, 
                                  idxB_lag = idxB_lag,
                                  pa_size = pa_size, 
                                  pa_duration = interval,
                                  catch_lag = 0
                      )),
        phcr = mseCtrl(method = phcr_comps,
                       args = list()),
        hcr = mseCtrl(method = hcr_comps,
                      args = list(interval = interval)),
        isys = mseCtrl(method = is_comps,
                       args = list(interval = interval))
      ))
      
    } else if (identical(MP, "CC_f")) {
      
      ### calculate reference length FL=M
      Lref <- c(0.75 * pars_l["Lc"] + 0.25 * lhist$linf)
      ### set up MP ctrl object
      ctrl <- ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_comps,
                      args = list(comp_r = FALSE, 
                                  comp_f = FALSE, 
                                  comp_c = FALSE,
                                  comp_A = TRUE,
                                  comp_b = FALSE,
                                  comp_m = FALSE,
                                  pa_buffer = FALSE, 
                                  pa_buffer_conditional = TRUE,
                                  comp_i = FALSE, 
                                  idxB_lag = idxB_lag,
                                  pa_size = pa_size, 
                                  pa_duration = interval,
                                  catch_lag = 0,
                                  Lref = Lref,
                                  Lref_mult = Lref_mult
                      )),
        phcr = mseCtrl(method = phcr_comps,
                       args = list()),
        hcr = mseCtrl(method = hcr_comps,
                      args = list(interval = interval)),
        isys = mseCtrl(method = is_comps,
                       args = list(interval = interval))
      ))
      
    }
    
    ### tracking
    if (isTRUE(MP %in% c("rfb", "hr", "2over3", "const_catch", "CC_f"))) {
      tracking <- c("comp_c", "comp_i", "comp_r", "comp_f", "comp_b",
                    "multiplier", "comp_hr", "exp_r", "exp_f", "exp_b")
      
      
    }
    
    ### args
    args <- list(fy = dims(stk)$maxyear, ### final simulation year
                 y0 = dims(stk)$minyear, ### first data year
                 iy = 100, ### first simulation (intermediate) year
                 nsqy = 3, ### not used, but has to provided
                 nblocks = n_blocks, ### block for parallel processing
                 seed = 1, ### random number seed before starting MSE
                 seed_part = FALSE
    )
    
    ### list with input object for mp()
    input <- list(om = om, oem = oem, iem = iem, ctrl = ctrl, 
                  args = args,
                  scenario = scenario, 
                  tracking = tracking, 
                  verbose = TRUE,
                  refpts = refpts, Blim = Blim)
    
    return(input)
    
  })
  names(input_list) <- stocks
  
  #if (identical(length(input_list), 1L)) input_list <- input_list[[1]]
  
  return(input_list)
  
}


