
### ------------------------------------------------------------------------ ###
### function for preparing OM input for MP ####
### ------------------------------------------------------------------------ ###

input_mp <- function(stocks, 
                     fhist, ### fishing history, e.g. "one-way"
                     selectivity = "default", ### selectivity default/dome
                     n_iter = 500, ### number of iterations
                     n_yrs = 50, ### number of years for MP
                     MP, ### e.g. "hr", "constant_catch", "CC_f"
                     hist_yr_min = 50, ### first year
                     scenario = "",
                     idx_sel = ifelse(isTRUE(MP %in% c("hr")), "tsb", "standard"), ### index selectivity
                     n_blocks = 1, 
                     ### index uncertainty & auto-correlation
                     sigmaB = 0.2,
                     sigmaL = 0.2,
                     sigmaCC = 0.1,
                     sigmaB_rho = 0,
                     sigmaL_rho = 0,
                     sigmaCC_rho = 0,
                     ### catch observation uncertainty
                     sigmaC = 0,
                     ### Catch implementation uncertainty
                     sigmaIEM = 0,
                     ### recruitment variability
                     rec_bias_correction = FALSE,
                     sigmaR = 0.6,
                     sigmaR_rho = 0,
                     ### recruitment steepness
                     steepness = 0.75,
                     ### HR definition
                     hr_value = "length",
                     ### est
                     multiplier = 1, ### also for CL
                     comp_b_multiplier = 1.4,
                     idxB_lag = 1,
                     idxB_range_1 = 2, idxB_range_2 = 3, ### 2 over 3 ratio
                     idxB_range_3 = 1, ### for b & i
                     catch_lag = 1, catch_range = 1,
                     idxL_lag = 1, idxL_range = 1,
                     pa_buffer = FALSE, pa_size = 0.8, pa_duration = interval,
                     Lref_mult = 1, ### also for CL
                     ### phcr
                     exp_r = 1, exp_f = 1, exp_b = 1, ### exponents (rfb only)
                     ### hcr
                     interval = ifelse(isTRUE(MP %in% c("rfb", "rb")), 2, 1), ### TAC interval
                     ### is
                     upper_constraint = 1.2, ### uncertainty cap
                     lower_constraint = 0.7, 
                     cap_below_b = FALSE, ### conditional cap?
                     ### MP: CL
                     n_catch = 3, 
                     n_length_1 = 3, 
                     n_length_2 = n_length_1,
                     lambda_upper = 0.1, 
                     lambda_lower = 0.2,
                     gamma_lower = 0.2, 
                     gamma_upper = 0.2,
                     r_threshold = 0.05, 
                     l_threshold = 0.01,
                     f_threshold = 0.1,
                     cc_threshold = 0.01,
                     #Lref_mult = 1, # defined above
                     first_catch = "advice",
                     catch_limit = 0,
                     ...
                     ) {
  
  ### load life-history parameters
  if (identical(selectivity, "default")) {
    brps <- readRDS("input/brps_new.rds")
  } else if (identical(selectivity, "dome")) {
    brps <- readRDS("input/brps_dome_new.rds")
  }
  
  stocks_lh <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
  
  ### stock ID
  if (is.numeric(stocks)) stocks <- stocks_lh$stock[stocks]
  
  ### go through list of stocks ####
  input_list <- lapply(stocks, function(stock) {#browser()
    
    lhist <- stocks_lh[stocks_lh$stock == stock, ]
    ### get reference points
    refpts <- refpts(brps[[stock]])
    Blim <- attr(brps[[stock]], "Blim")
    
    ### load stock & sr
    stk <- readRDS(paste0("input/", n_iter, "_", n_yrs, 
                          switch(selectivity, 
                                 "default" = "/OM/",
                                 "dome" = "/OM_dome/"),
                          fhist, "/",
                          stock, "/stk.rds"))
    sr <- readRDS(paste0("input/", n_iter, "_", n_yrs, 
                         switch(selectivity, 
                                "default" = "/OM/",
                                "dome" = "/OM_dome/"),
                         fhist, "/", stock, "/sr.rds"))
    
    ### adapt recruitment variability if needed #####
    if (isTRUE(!identical(sigmaR, 0.6) | !identical(sigmaR_rho, 0))) {
      ### bias correction
      if (isTRUE(rec_bias_correction)) {
        # m <- 1
        # cv <- sigmaR
        # v <- (cv*m)^2 ### variance
        # mu <- log(m^2/sqrt(m^2+v))
        # sigma <- sqrt(log(v/m^2+1))
        mu <- 0 - (sigmaR^2)/2
      } else {
        sigma <- sigmaR
        mu <- 0
      }
      ### create recruitment residuals for projection period
      set.seed(0)
      residuals(sr) <- rlnoise(n_iter, rec(stk) %=% mu, 
                               sd = sigma, b = sigmaR_rho)
      ### replicate residuals from GA paper
      set.seed(0)
      residuals(sr)[, ac(0:150)] <- rlnoise(n_iter, rec(stk)[, ac(0:150)] %=% mu,
                                            sd = sigma, b = sigmaR_rho)
      ### replicate residuals from catch rule paper for historical period
      set.seed(0)
      residuals(sr)[, ac(1:100)] <- rlnoise(n_iter, (rec(stk) %=% mu)[, ac(1:100)], 
                                            sd = sigma, b = sigmaR_rho)
    }
    
    ### cut off history
    stk <- window(stk, start = hist_yr_min)
    sr@residuals <- window(sr@residuals, start = hist_yr_min)
    
    ### adapt recruitment steepness if needed #####
    if (!identical(steepness, 0.75)) {
      ### calculate new recruitment model parameters with new steepness
      alpha <- (4*steepness*c(refpts["virgin", "rec"])) /
        (5*steepness - 1)
      beta <- (c(refpts["virgin", "ssb"]) * (1 - steepness)) /
        (5*steepness - 1)
      ### insert values
      params(sr)[] <- c(alpha, beta)
    }
    
    ### operating model ####
    om <- FLom(stock = stk, ### stock 
               sr = sr, ### stock recruitment and precompiled residuals
               fleetBehaviour = mseCtrl(),
               projection = mseCtrl(method = fwd_attr,
                                    args = list(dupl_trgt = TRUE)))
  
    ### length data ####
    if (MP %in% c("rfb", "hr", "CC_f", "CL", "CL_cc", "fcc")) {
      ### allometric a & b
      ### length at first capture
      pars_l <- FLPar(a = lhist$a,
                      b = lhist$b,
                      Lc = calc_lc(stk = stk[, ac(75:100)], 
                                   a = lhist$a, b = lhist$b),
                      Linf = lhist$linf)
      ### mean catch length
      idx_Lmean <- lmean(stk = stk, params = pars_l)
      
      ### calculate reference length FL=M
      Lref <- c(0.75 * pars_l["Lc"] + 0.25 * lhist$linf)
      
      ### length catch curve
      if (isTRUE(MP %in% c("CL_cc", "fcc"))) {
        idx_cc <- lcc(ay = 100,
                      idx = idx_Lmean %=% NA_real_,
                      stk = stk, params = pars_l)
      }
      
    }
    
    ### index template ####
    idx <- FLIndices(
      sel       = FLIndex(index = stk@mat %=% NA_real_),
      idxB      = FLIndex(index = ssb(stk) %=% NA_real_),
      idxL      = FLIndex(index = ssb(stk) %=% NA_real_),
      idxCC     = FLIndex(index = ssb(stk) %=% NA_real_),
      PA_status = FLIndex(index = ssb(stk) %=% NA_integer_)
    )
    
    ### observation (error) model ####
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
    } else if (identical(MP, "rfb")) {
      oem <- FLoem(method = obs_generic,
                   observations = list(stk = stk, idx = idx), 
                   deviances = list(stk = FLQuant(), idx = idx),
                   args = list(idxB = TRUE, idx_dev = TRUE, ssb_idx = FALSE, tsb_idx = FALSE,
                               lngth = TRUE, lngth_dev = TRUE,
                               lngth_par = pars_l,
                               PA_status = FALSE, PA_status_dev = FALSE,
                               PA_Bmsy = c(refpts(brps[[stock]])["msy", "ssb"]), 
                               PA_Fmsy = c(refpts(brps[[stock]])["msy", "harvest"])))
    } else if (identical(MP, "constant_catch")) {
      oem <- FLoem(method = obs_generic,
                   observations = list(stk = stk, idx = idx), 
                   deviances = list(stk = FLQuant(), idx = idx),
                   args = list(idxB = FALSE, 
                               lngth = FALSE, lngth_dev = FALSE,
                               PA_status = TRUE, PA_status_dev = FALSE,
                               PA_Bmsy = Inf, ### i.e. B always below Bmsy/2
                               PA_Fmsy = 0 ### i.e. F always above Fmsy
                               ### -> always apply buffer
                               ))
    } else if (identical(MP, "CC_f")) {
        oem <- FLoem(method = obs_generic,
                     observations = list(stk = stk, idx = idx), 
                     deviances = list(stk = FLQuant(), idx = idx),
                     args = list(idxB = FALSE,
                                 lngth = TRUE, lngth_dev = TRUE,
                                 lngth_par = pars_l,
                                 PA_status = TRUE, PA_status_dev = FALSE,
                                 PA_Bmsy = Inf, ### i.e. B always below Bmsy/2
                                 PA_Fmsy = 0 ### i.e. F always above Fmsy
                                 ### -> always apply buffer
                     ))
    ### Catch + length data
    } else if (identical(MP, "CL")) {
      oem <- FLoem(method = obs_generic,
                   observations = list(stk = stk, idx = idx), 
                   deviances = list(stk = FLQuant(), idx = idx),
                   args = list(idxB = FALSE,
                               lngth = TRUE, lngth_dev = TRUE,
                               lngth_par = pars_l
                   ))
    ### Catch + length data CC
    } else if (identical(MP, "CL_cc")) {
      oem <- FLoem(method = obs_generic,
                   observations = list(stk = stk, idx = idx), 
                   deviances = list(stk = FLQuant(), idx = idx),
                   args = list(idxB = FALSE,
                               lngth = TRUE, lngth_dev = TRUE,
                               lngth_cc = TRUE, lngth_cc_dev = TRUE,
                               lngth_par = pars_l
                   ))
    } else if (identical(MP, "fcc")) {
      oem <- FLoem(method = obs_generic,
                   observations = list(stk = stk, idx = idx), 
                   deviances = list(stk = FLQuant(), idx = idx),
                   args = list(idxB = FALSE,
                               lngth = TRUE, lngth_dev = TRUE,
                               lngth_cc = TRUE, lngth_cc_dev = TRUE,
                               lngth_par = pars_l
                   ))
    }
    
    ### create biomass index - if needed ####
    if (MP %in% c("rfb", "hr")) {
      ### update index if not total stock biomass
      if (identical(idx_sel, "ssb")) {
        oem@args$ssb_idx <- TRUE
        oem@args$tsb_idx <- FALSE
        oem@args$idx_timing <- TRUE
        index(oem@observations$idx$sel) <- mat(stk)
      } else if (identical(idx_sel, "tsb")) {
        oem@args$ssb_idx <- FALSE
        oem@args$tsb_idx <- TRUE
        oem@args$idx_timing <- TRUE
        index(oem@observations$idx$sel) <- mat(stk) %=% 1
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
        index(oem@observations$idx$sel)[] <- csel
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
        index(oem@observations$idx$sel)[] <- sel
      } else {
        ### standard scientific survey
        oem@args$ssb_idx <- FALSE
        oem@args$tsb_idx <- FALSE
        oem@args$idx_timing <- TRUE
        sel <- 1/(1 + exp(-1*(an(dimnames(stk)$age) - dims(stk)$max/10)))
        index(oem@observations$idx$sel)[] <- sel
      }
      
      ### calculate biomass index with selectivity
      ### get stock numbers and reduce by F and M
      sn <- stock.n(stk) * exp(-(harvest(stk) * harvest.spwn(stk) + 
                                   m(stk) * m.spwn(stk)))
      ### calculate index
      idxB <- quantSums(sn * stock.wt(stk) * index(oem@observations$idx$sel))
      index(oem@observations$idx$idxB) <- idxB
    } 

    ### length index
    ### only include when required for MP
    if (isTRUE(MP %in% c("CC_f", "CL", "CL_cc", "rfb", "fcc"))) {
      oem@args$lngth <- TRUE
      oem@args$lngth_dev <- TRUE
      oem@args$lngth_par <- pars_l
      ### calculate mean length
      index(oem@observations$idx$idxL) <- idx_Lmean
      ### catch curve
      if (isTRUE(MP %in% c("CL_cc", "fcc"))) {
        index(oem@observations$idx$idxCC) <- idx_cc
      }
    }
  
    ### index deviations ####
    
    ### PA buffer deviations ####
    if (isTRUE(MP %in% c("constant_catch", "CC_f"))) {
      PA_status_dev <- FLQuant(NA, dimnames = list(age = c("positive", "negative"), 
                                                   year = dimnames(stk)$year, 
                                                   iter = dimnames(stk)$iter))
      ### always apply PA buffer
      set.seed(1)
      PA_status_dev["positive"] <- 1
      set.seed(2)
      PA_status_dev["negative"] <- 1
      index(oem@deviances$idx$PA_status) <- PA_status_dev
    } 
    
    index(oem@deviances$idx$sel) <- index(oem@deviances$idx$sel) %=% 1
    
    ### biomass index deviations ####
    set.seed(695)
    if (isTRUE(MP %in% c("hr", "rfb"))) {
      index(oem@deviances$idx$idxB) <- 
        rlnoise(n = dims(idx$idxB)$iter, index(idx$idxB) %=% 0, 
                sd = sigmaB, b = sigmaB_rho)
    }
    if (isTRUE(MP %in% c("hr", "rfb", "CC_f", "CL", "CL_cc", "fcc"))) {
      index(oem@deviances$idx$idxL) <- 
        rlnoise(n = dims(idx$idxL)$iter, index(idx$idxL) %=% 0, 
        sd = sigmaL, b = sigmaL_rho)
      if (isTRUE(MP %in% c("CL_cc", "fcc"))) {
        index(oem@deviances$idx$idxCC) <- 
          rlnoise(n = dims(idx$idxCC)$iter, index(idx$idxCC) %=% 0, 
                  sd = sigmaCC, b = sigmaCC_rho)
      }
      
    }
    ### replicate previous deviates from GA paper
    set.seed(696)
    if (isTRUE(MP %in% c("hr", "rfb"))) {
      index(oem@deviances$idx$idxB)[, ac(50:150)] <- 
        rlnoise(n = dims(oem@deviances$idx$idxB)$iter,
                window(index(oem@deviances$idx$idxB), end = 150) %=% 0,
                sd = sigmaB, b = sigmaB_rho)
    }
    if (isTRUE(MP %in% c("hr", "rfb", "CC_f", "CL", "CL_cc", "fcc"))) {
      index(oem@deviances$idx$idxL)[, ac(50:150)] <- 
        rlnoise(n = dims(oem@deviances$idx$idxL)$iter, 
                window(index(oem@deviances$idx$idxB), end = 150) %=% 0,
                sd = sigmaL, b = sigmaL_rho)
    }
    
    if (MP %in% c("hr", "rfb")) {
      ### biomass index reference points
      ### lowest observed index in last 50 years
      I_loss <- apply(window(index(oem@observations$idx$idxB) * 
                                     index(oem@deviances$idx$idxB), end = 100),
                      6, min)
      I_loss * comp_b_multiplier
    } else {
      I_loss <- NA
    }
    
    ### catch deviation ####
    if (MP %in% c("CL", "fcc")) {
      ### deviation on total catch
      set.seed(699)
      oem@deviances$stk <- FLQuants(
        catch = rlnoise(n = dims(stk)$iter, catch(stk) %=% 0, 
                        sd = sigmaC, b = 0))
      oem@args$catch_dev <- TRUE
    }
  
    ### Implementation error model (IEM) ####
    set.seed(205)
    iem_dev <- FLQuant(rlnoise(n = dims(stk)$iter,  catch(stk) %=% 0,
                               sd = sigmaIEM, b = 0))
    ### iem object
    iem <- FLiem(method = iem_comps,
                 args = list(use_dev = ifelse(isTRUE(sigmaIEM > 0),
                                              TRUE, FALSE), 
                             iem_dev = iem_dev))
    
    ### MP control object ####
    
    ### absolute catch limit
    if (isTRUE(MP %in% c("fcc", "CL", "CL_cc"))) {
      if (isTRUE(catch_limit > 0)) {
        catch_oem <- catch(oem@observations$stk) * oem@deviances$stk$catch
        catch_limit_val <- c(apply(window(catch_oem, end = 100), 6, quantile, 
                                   probs = catch_limit))
      } else {
        catch_limit_val <- Inf
      }
    }
    
    ### harvest rate (chr) ####
    if (identical(MP, "hr")) {
      ### define target harvest rate(s)
      if (identical(hr_value, "uniform")) {
        set.seed(33)
        hr_val <- runif(n = n_iter, min = 0, max = 1)
      } else if (identical(hr_value, "Fmsy")) {
        ### load some reference harvest rates 
        hr_refs <- readRDS("input/catch_rates.rds")[stock]
        hr_val <- switch(idx_sel, 
                         "tsb" = hr_refs[[stock]]$Fmsy$tsb,
                         "ssb" = hr_refs[[stock]]$Fmsy$ssb,
                         "survey" = hr_refs[[stock]]$Fmsy$idx,
                         hr_refs[[stock]]$Fmsy$tsb)
      } else if (identical(hr_value, "LFeM")) {
        ### load some reference harvest rates 
        hr_refs <- readRDS("input/catch_rates.rds")[stock]
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
        idxL <- Lmean * index(oem@deviances$idx$idxL)
        idxL <- window(idxL, end = 100)
        ### relative to reference length
        idxL <- idxL / LFeM
        ### biomass index
        idxB <- index(oem@observations$idx$idxB) * index(oem@deviances$idx$idxB)
        idxB <- window(idxB, end = 100)
        ### historical harvest rate
        CI <- window(catch(stk), end = 100)/idxB
        ### average of harvest rates where catch length is above reference
        hr_val <- sapply(dimnames(stk)$iter, function(i){
          mean(CI[, which(idxL[,,,,, i] >= 1),,,, i])
        })
      }
      
      ### set up MP ctrl object
      ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_comps,
                      args = list(comp_r = FALSE, 
                                  comp_f = FALSE, 
                                  comp_c = FALSE, comp_A = FALSE,
                                  comp_b = TRUE,
                                  comp_m = multiplier,
                                  pa_buffer = FALSE, 
                                  comp_i = TRUE, 
                                  I_loss = I_loss,
                                  comp_b_multiplier = 1.4,
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
    
    } else if (identical(MP, "rfb")) {
      
      ### set up MP ctrl object
      ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_comps,
                      args = list(comp_r = TRUE, 
                                  comp_f = TRUE, 
                                  comp_c = TRUE, comp_A = FALSE,
                                  comp_b = TRUE,
                                  comp_m = multiplier,
                                  pa_buffer = FALSE, 
                                  comp_i = FALSE, 
                                  idxB_lag = 1, 
                                  idxB_range_1 = 2, idxB_range_2 = 3,
                                  catch_lag = 1, catch_range = 1,
                                  Lref = Lref, 
                                  I_loss = I_loss,
                                  comp_b_multiplier = 1.4,
                                  idxL_lag = 1, idxL_range = 1
                      )),
        phcr = mseCtrl(method = phcr_comps,
                       args = list(exp_r = 1, exp_f = 1, exp_b = 1)),
        hcr = mseCtrl(method = hcr_comps,
                      args = list(interval = interval)),
        isys = mseCtrl(method = is_comps,
                       args = list(interval = interval, 
                                   upper_constraint = upper_constraint, 
                                   lower_constraint = lower_constraint, 
                                   cap_below_b = cap_below_b))
      ))
      
    } else if (identical(MP, "constant_catch")) {
      ### constant_catch ####
      ### set up MP ctrl object
      ctrl <- mpCtrl(list(
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
      ### CC_f ####
      ### set up MP ctrl object
      ctrl <- mpCtrl(list(
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
      
    } else if (identical(MP, "CL")) {
      ### CL ####
      ### set up MP ctrl object
      ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_CL,
                      args = list(interval = interval,
                                  n_catch = n_catch, 
                                  n_length_1 = n_length_1, 
                                  n_length_2 = n_length_2,
                                  first_catch = first_catch)),
        phcr = mseCtrl(method = phcr_CL,
                       args = list(interval = interval,
                                   lambda_upper = lambda_upper,
                                   lambda_lower = lambda_lower,
                                   gamma_lower = gamma_lower, 
                                   gamma_upper = gamma_upper,
                                   r_threshold = r_threshold, 
                                   l_threshold = l_threshold,
                                   f_threshold = f_threshold,
                                   Lref = Lref,
                                   Lref_mult = Lref_mult,
                                   multiplier = multiplier)),
        hcr = mseCtrl(method = hcr_CL,
                      args = list(interval = interval,
                                  combine_alpha_beta = FALSE)),
        isys = mseCtrl(method = is_comps,
                       args = list(interval = interval,
                                   catch_limit = catch_limit_val))
      ))
      
    } else if (identical(MP, "CL_cc")) {
      ### CL_cc ####
      ### set up MP ctrl object
      ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_CL,
                      args = list(interval = interval,
                                  n_catch = n_catch, 
                                  n_length_1 = n_length_1, 
                                  n_length_2 = n_length_2,
                                  r_catch = TRUE, r_length = FALSE, r_cc = TRUE,
                                  first_catch = first_catch)),
        phcr = mseCtrl(method = phcr_CL,
                       args = list(interval = interval,
                                   lambda_upper = lambda_upper,
                                   lambda_lower = lambda_lower,
                                   gamma_lower = gamma_lower, 
                                   gamma_upper = gamma_upper,
                                   r_threshold = r_threshold, 
                                   l_threshold = l_threshold,
                                   f_threshold = f_threshold,
                                   Lref = Lref,
                                   Lref_mult = Lref_mult,
                                   multiplier = multiplier)),
        hcr = mseCtrl(method = hcr_CL,
                      args = list(interval = interval,
                                  combine_alpha_beta = FALSE,
                                  multiple = TRUE)),
        isys = mseCtrl(method = is_comps,
                       args = list(interval = interval,
                                   catch_limit = catch_limit_val))
      ))
      
    } else if (identical(MP, "fcc")) {
      ### find reference slope from length catch curve
      ### same principle as for chr rule
      Lc <- calc_lc(stk = window(stk, end = 100), a = lhist$a, b = lhist$b)
      LFeM <- (lhist$linf + 2*1.5*c(Lc)) / (1 + 2*1.5)
      Lmean <- window(index(oem@observations$idx$idxL) * 
                        index(oem@deviances$idx$idxL),
                      end = 100)
      Lstatus <- Lmean/LFeM
      Lcc <- window(index(oem@observations$idx$idxCC) * 
                      index(oem@deviances$idx$idxCC),
                    end = 100)
      Lcc_use <- Lcc
      Lcc_use[Lstatus < 1] <- NA
      Lcc_use <- c(yearMeans(Lcc_use))
      
      ### set up MP ctrl object
      ctrl <- mpCtrl(list(
        est = mseCtrl(method = est_comps,
                      args = list(comp_fcc = TRUE, 
                                  comp_m = multiplier,
                                  catch_lag = 1, catch_range = 1,
                                  Lref = Lcc_use, 
                                  idxL_lag = idxL_lag, idxL_range = idxL_range
                      )),
        phcr = mseCtrl(method = phcr_comps,
                       args = list()),
        hcr = mseCtrl(method = hcr_comps,
                      args = list(interval = interval)),
        isys = mseCtrl(method = is_comps,
                       args = list(interval = interval, 
                                   upper_constraint = upper_constraint, 
                                   lower_constraint = lower_constraint, 
                                   cap_below_b = TRUE,
                                   catch_limit = catch_limit_val))
      ))
      
    }
    
    ### tracking ####
    if (isTRUE(MP %in% c("rfb", "hr", "2over3", "constant_catch", "CC_f", "fcc"))) {
      tracking <- c("comp_c", "comp_i", "comp_r", "comp_f", "comp_b",
                    "multiplier", "comp_hr", "exp_r", "exp_f", "exp_b")
    } else if (isTRUE(MP %in% c("CL", "CL_cc"))) {
      tracking <- c("A_last", "r_length", "r_catch", "length_average",
                    "catch_limit")
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
    
    ### list with input object for mp() ####
    input <- list(om = om, oem = oem, iem = iem, ctrl = ctrl, 
                  args = args,
                  scenario = scenario, 
                  tracking = tracking, 
                  verbose = TRUE,
                  parallel = n_blocks)
    attr(input, "refpts") <- refpts
    attr(input, "Blim") <- Blim
    
    return(input)
    
  })
  names(input_list) <- stocks
  
  #if (identical(length(input_list), 1L)) input_list <- input_list[[1]]
  
  return(input_list)
  
}


### ------------------------------------------------------------------------ ###
### function for preparing OM input for MP - lcc ####
### ------------------------------------------------------------------------ ###

input_mp_lcc <- function(
  stock, 
  fhist, ### fishing history, e.g. "one-way"
  selectivity = "default", ### selectivity default/dome
  n_iter = 500, ### number of iterations
  n_yrs = 50, ### number of years for MP
  MP = "lcc", ### e.g. "hr", "constant_catch", "CC_f"
  hist_yr_min = 50, ### first year
  scenario = "",
  n_blocks = 1, 
  ### length sampling (effective sample size)
  ess_age = 50, ### 1st step for sampling from ages
  ess_len = 500, ### 2nd step for sampling lengths
  ### catch observation uncertainty
  sigmaC = 0,
  ### Catch implementation uncertainty
  sigmaIEM = 0,
  ### recruitment variability
  rec_bias_correction = TRUE,
  sigmaR = 0.6,
  sigmaR_rho = 0,
  ### recruitment steepness
  steepness = 0.75,
  ### MP
  multiplier = 1, 
  interval = 2, ### TAC interval
  upper_constraint = 1.2, ### uncertainty cap
  lower_constraint = 0.7, 
  catch_limit = 0.5, ### percentile of historical catches
  ...
) {
  
  ### load life-history parameters
  if (identical(selectivity, "default")) {
    brps <- readRDS("input/brps_len.rds")
  } else if (identical(selectivity, "dome")) {
    brps <- readRDS("input/brps_dome_len.rds")
  }
  
  stocks_lh <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
  
  ### stock ID
  if (is.numeric(stock)) stock <- stocks_lh$stock[stock]
  
  lhist <- stocks_lh[stocks_lh$stock == stock, ]
  ### get reference points
  refpts <- refpts(brps[[stock]])
  Blim <- attr(brps[[stock]], "Blim")
  
  ### load stock & sr
  stk <- readRDS(paste0("input/", n_iter, "_", 100, 
                        switch(selectivity, 
                               "default" = "/OM_len/",
                               "dome" = "/OM_len_dome/"),
                        fhist, "/",
                        stock, "/stk.rds"))
  sr <- readRDS(paste0("input/", n_iter, "_", 100, 
                       switch(selectivity, 
                              "default" = "/OM_len/",
                              "dome" = "/OM_len_dome/"),
                       fhist, "/", stock, "/sr.rds"))
  
  ### adapt recruitment variability if needed #####
  if (isTRUE(!identical(sigmaR, 0.6) | !identical(sigmaR_rho, 0))) {
    ### bias correction
    if (isTRUE(rec_bias_correction)) {
      mu <- 0 - (sigmaR^2)/2
    } else {
      sigma <- sigmaR
      mu <- 0
    }
    ### create recruitment residuals for projection period
    set.seed(0)
    residuals(sr) <- rlnoise(n_iter, rec(stk) %=% mu, 
                             sd = sigma, b = sigmaR_rho)
  }
  
  ### cut off history
  stk <- window(stk, start = hist_yr_min)
  sr@residuals <- window(sr@residuals, start = hist_yr_min)
  
  ### years for projection
  stk <- window(stk, end = 100 + n_yrs)
  sr@residuals <- window(sr@residuals, end = 100 + n_yrs)
  
  ### adapt recruitment steepness if needed #####
  if (!identical(steepness, 0.75)) {
    ### calculate new recruitment model parameters with new steepness
    alpha <- (4*steepness*c(refpts["virgin", "rec"])) /
      (5*steepness - 1)
    beta <- (c(refpts["virgin", "ssb"]) * (1 - steepness)) /
      (5*steepness - 1)
    ### insert values
    params(sr)[] <- c(alpha, beta)
  }
  
  ### operating model ####
  om <- FLom(stock = stk, ### stock 
             sr = sr, ### stock recruitment and precompiled residuals
             fleetBehaviour = mseCtrl(),
             projection = mseCtrl(method = fwd_attr))
  
  ### generate historical length data
  ### object with gear definition
  gear <- attr(brp, "gear")
  ### random number seeds for length sampling
  ### 1st create full year range and later subset to those used
  res_length <- rUnif_lfd(gear, 
                          years = seq(from = range(stk)[["minyear"]],
                                      to = 200), 
                          nsim = dim(stk)[6], 
                          seed = 123)
  res_length <- lapply(res_length, function(x) {
    x[[1]][, dimnames(stk)$year, ]
  })
  ### catch numbers at age
  cn <- window(catch.n(stk), end = 100)
  cn[cn <= 0] <- sqrt(.Machine$double.eps)
  ### simulate length data (reproducible)
  lfd_file <- paste0("input/", n_iter, "_", 100, 
                     switch(selectivity, 
                            "default" = "/OM_len/",
                            "dome" = "/OM_len_dome/"),
                     fhist, "/",
                     stock, "/lfd.rds")
  if (file.exists(lfd_file)) {
    len_sim <- readRDS(lfd_file)
  } else {
    len_sim <- lfd.sim(object = cn, 
                       gear = gear[[1]], 
                       ess_age = ess_age,
                       ess_len = ess_len,
                       scale = FALSE,
                       u_age = res_length$age,
                       u_len = res_length$len
    )
  }
  
  ### calculate catch curve slope for all years and iterations
  idx_s <- lcc_fun(len_data = len_sim, Linf = lhist$linf)

  ### calculate length at first capture
  Lc <- apply(len_sim, MARGIN = c(2:6), function(x) {
    #browser()
    ### find position of first length class at or above half mode abundance
    lc_pos <- head(which(x >= (max(x, na.rm = TRUE)/2)), 1)
    ### Lc
    lc <- an(names(x)[lc_pos])
    if (is.na(lc)) {
      return(an(names(x)[1]))
    } else {
      return(lc)
    }
  })
  Lc <- median(Lc)
  
  ### calculate mean length and find reference years with L >= LF=M
  Lmean <- apply(X = len_sim[an(dimnames(len_sim)$len) >= Lc, ], 
                 MARGIN = c(2, 6), FUN = function(x) {#browser()
    ### calculate
    res <- weighted.mean(x = an(dimnames(x)$len), 
                         w = ifelse(is.na(x), 0, x), na.rm = TRUE)
    return(res)
  })
  ### calculate reference length FL=M
  Lref <- c(0.75 * Lc + 0.25 * lhist$linf)

  ### calculate target slope
  ### average of slope values where catch length is above reference
  s_target <- sapply(dimnames(stk)$iter, function(i){
    mean(idx_s[, which((Lmean/Lref)[,,,,, i] >= 1),,,, i])
  })
  
  ### OEM
  oem_observations <- list(
    stk = stk, 
    idx = FLIndices(
      LFD = FLIndex(index = window(len_sim, end = range(stk)[["maxyear"]]),
                    plusgroup = NA)))
  oem <- FLoem(observations = oem_observations, 
               deviances = list(lfd = res_length),
               method = oem_lfd_sampling,
               args = list(gear = gear,
                           ess_age = ess_age,
                           ess_len = ess_len))
  
  ### calculate absolute catch limit
  ### based on percentile of historical catches
  if (isTRUE(catch_limit > 0)) {
    catch_limit_val <- c(apply(window(catch(stk), end = 100), 6, quantile, 
                               probs = catch_limit))
  } else {
    catch_limit_val <- Inf
  }

  ### MP
  ctrl <- mpCtrl(list(
    est = mseCtrl(method = est_lfd_lcc,
                  args = list(Linf = c(lhist$linf))),
    phcr = mseCtrl(method = phcr_lfd_lcc,
                   args = list(m = multiplier, 
                               s_target = s_target, 
                               catch_limit = catch_limit_val)),
    hcr = mseCtrl(method = hcr_lcc_lfd,
                  args = list(interval = interval,
                              upper_constraint = upper_constraint,
                              lower_constraint = lower_constraint))
  ))
  
  ### Implementation error model (IEM) ####
  set.seed(205)
  iem_dev <- FLQuant(rlnoise(n = dims(stk)$iter,  catch(stk) %=% 0,
                             sd = sigmaIEM, b = 0))
  ### iem object
  iem <- FLiem(method = iem_comps,
               args = list(use_dev = ifelse(isTRUE(sigmaIEM > 0),
                                            TRUE, FALSE), 
                           iem_dev = iem_dev))
  
  ### args
  args <- list(fy = dims(stk)$maxyear, ### final simulation year
               y0 = dims(stk)$minyear, ### first data year
               iy = 100, ### first simulation (intermediate) year
               nsqy = 3, ### not used, but has to provided
               nblocks = n_blocks, ### block for parallel processing
               seed = 1, ### random number seed before starting MSE
               seed_part = FALSE
  )
  
  ### list with input object for mp() ####
  input <- list(om = om, 
                oem = oem, 
                iem = iem, 
                ctrl = ctrl, 
                args = args,
                scenario = scenario, 
                verbose = TRUE,
                parallel = n_blocks)
  attr(input, "refpts") <- refpts
  attr(input, "Blim") <- Blim
  
  return(input)

}


