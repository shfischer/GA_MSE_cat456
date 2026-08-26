### ------------------------------------------------------------------------ ###
### analysis of category 4-6 harvest control rules ####
### ------------------------------------------------------------------------ ###

library(mse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(doParallel)
library(scales)
library(patchwork)
library(RColorBrewer)
library(foreach)
source("funs.R")

### stock list
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
### brps
brps <- readRDS("input/brps_new.rds")

### ------------------------------------------------------------------------ ###
### "constant_catch" - current ICES approach ####
### ------------------------------------------------------------------------ ###
### constant catch with PA buffer
### - reduce triennially by 20%

### ------------------------------------------------------------------------ ###
### constant_catch - risk & projection time - all stocks ####
### ------------------------------------------------------------------------ ###
### 500 iters, 100 years, all stocks, one-way & random

stats_risk <- foreach(stock = stocks$stock, k = stocks$k, 
                      .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows)  %do% {
    #browser()
    res <- readRDS(paste0("output/constant_catch/500_100/baseline/", fhist,"/",
                          stock, "/mp_0.1_0_0.1_0.6_0_0.75_0_.rds"))
    ### collapse correction
    res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
    
    ### get fishing history
    stk <- readRDS(paste0("input/500_100/OM/", fhist, "/", stock, "/stk.rds"))
    ssb_ <- ssb(stk)
    catch_ <- catch(stk)
    ### add simulated values
    ssb_[, ac(101:200)] <- res_corrected$ssb
    catch_[, ac(101:200)] <- res_corrected$catch
    
    ### metrics (relative)
    ssb_rel <- ssb_/c(refpts(brps[[stock]])["msy", "ssb"])
    catch_rel <- catch_/c(refpts(brps[[stock]])["msy", "yield"])
    risk <- apply(ssb_ < attr(brps[[stock]], "Blim"), 2, mean)
    risk[, ac(1:100)] <- NA
    
    ### cumulative values (year 101 to year x)
    ssb_rel_cum <- iterMeans(ssb_rel) %=% NA_real_
    ssb_rel_cum[, ac(101:200)] <- sapply(1:100, function(x) {
      median(c(ssb_rel[, ac(seq(from = 101, length.out = x))]), 
             na.rm = TRUE)
    })
    catch_rel_cum <- iterMeans(catch_rel) %=% NA_real_
    catch_rel_cum[, ac(101:200)] <- sapply(1:100, function(x) {
      median(c(catch_rel[, ac(seq(from = 101, length.out = x))]), 
             na.rm = TRUE)
    })
    risk_cum <- iterMeans(ssb_rel) %=% NA_real_
    risk_cum[, ac(101:200)] <- sapply(1:100, function(x) {
      mean(c(ssb_[, ac(seq(from = 101, length.out = x))] < attr(brps[[stock]], "Blim")), 
           na.rm = TRUE)
    })
    
    ### annual values (time series)
    ssb_rel <- iterMedians(ssb_rel)
    catch_rel <- iterMedians(catch_rel)
    
    ### combine data into data.frame
    df_tmp <- as.data.frame(FLQuants(ssb_annual = ssb_rel, 
                                     ssb_cum = ssb_rel_cum,
                                     catch_annual = catch_rel, 
                                     catch_cum = catch_rel_cum,
                                     risk_annual = risk, 
                                     risk_cum = risk_cum))
    df_tmp <- df_tmp %>% 
      select(year, data, qname) %>%
      separate(col = qname, sep = "_", into = c("metric", "calculation"))
    df_tmp$fhist <- fhist
    df_tmp$stock <- stock
    df_tmp$k <- k
    return(df_tmp)
}

stats_risk <- stats_risk %>%
  mutate(metric = factor(metric,
                         levels = c("ssb", "catch", "risk"),
                         labels = c("SSB/B[MSY]",
                                    "Catch/MSY",
                                    "B[lim]~risk")),
         fhist = factor(fhist,
                        levels = c("one-way", "random")),
         calculation = factor(calculation,
                              levels = c("annual", "cum"),
                              labels = c("Annual", "Cumulative")))
### calculate median over all stocks
stats_risk <- bind_rows(
  stats_risk %>%
    mutate(source = "stocks"),
  stats_risk %>%
    group_by(year, metric, calculation, fhist) %>%
    summarise(data = median(data)) %>%
    mutate(source = "median")) %>%
  mutate(source = factor(source, 
                         levels = c("stocks", "median")))


res_def_colours <- c("one-way" = RColorBrewer::brewer.pal(n = 4, name = "Set1")[1], 
                     "random" = RColorBrewer::brewer.pal(n = 4, name = "Set1")[2])



#stats_risk %>%
#filter(stock == "pol") %>%
ggplot() +
  geom_line(data = stats_risk %>% filter(source == "stocks"),
            aes(x = year - 100, y = data,
                group = interaction(stock, fhist, source), colour = fhist),
            linewidth = 0.1, linetype = "dashed") +
  geom_line(data = stats_risk %>% filter(source == "median"),
            aes(x = year - 100, y = data, 
                colour = fhist),
            linewidth = 0.5, linetype = "solid") +
  facet_grid(metric ~ calculation, labeller = "label_parsed", switch = "y", 
             scales = "free_y") +
  scale_colour_manual("fishing history", values = res_def_colours) +
  theme_bw(base_size = 8) +
  #xlim(c(-2, NA)) +
  coord_cartesian(xlim = c(-2, 102), ylim = c(0, NA), expand = FALSE) +
  labs(x = "Year") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.6),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"))

ggsave(filename = "output/plots/constant_catch/CC_all_timeseries.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)

### annual only
ggplot() +
  geom_line(data = stats_risk %>% 
              filter(source == "stocks" & calculation == "Annual"),
            aes(x = year - 100, y = data,
                group = interaction(stock, fhist, source), colour = fhist),
            linewidth = 0.1, linetype = "dashed") +
  geom_line(data = stats_risk %>% 
              filter(source == "median" & calculation == "Annual"),
            aes(x = year - 100, y = data, 
                colour = fhist),
            linewidth = 0.5, linetype = "solid") +
  facet_wrap(~ metric, labeller = "label_parsed", strip.position = "left", 
             scales = "free_y", ncol = 1) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  theme_bw(base_size = 8) +
  coord_cartesian(xlim = c(-2, 102), ylim = c(0, NA), expand = FALSE) +
  labs(x = "Year") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.55),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"))

ggsave(filename = "output/plots/constant_catch/CC_all_annual_timeseries.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### constant_catch - risk vs depletion ####
### ------------------------------------------------------------------------ ###
### random fhist, 10,000 iterations, 100 years
### examples: pollack, turbot, herring, anglerfish

res_df <- foreach(stk = c("pol", "tur", "ang3", "her"), 
                  .combine = bind_rows) %do% {
  #browser()
  
  ### reference points
  Blim <- attr(brps[[stk]], "Blim")
  Bmsy <- brps[[stk]]@refpts["msy", "ssb"]
  MSY <- brps[[stk]]@refpts["msy", "yield"]
  
  ### MP output
  res <- readRDS(paste0("output/constant_catch/10000_100/baseline/random/",
                        stk, "/mp_0.1_0_0.1_0.6_0_0.75_0_.rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
  
  ### starting condition
  SSBs0 <- ssb(res@om@stock)[, ac(100)]
  SSBs0 <- SSBs0/c(Bmsy)
  SSBs0 <- c(SSBs0)
  max_SSB <- max(SSBs0) #3
  SSB_breaks <- seq(from = 0, to = max_SSB, by = 0.1)
  SSB_groups <- cut(SSBs0, breaks = SSB_breaks)
  SSB_levels <- unique(as.character(SSB_groups))
  ### number of replicates per group
  group_n <- sapply(SSB_levels, function(x) {
    length(which(SSB_groups %in% x))
  })
  group_n[sort(names(group_n))]
  # sum(group_n)
  ### Blim risk per group
  ### SSB is on absolute scale 
  risk_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]
    mean(tmp < Blim)
  })
  ### SSB (long-term median) per group
  SSB_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]/c(Bmsy)
    median(tmp)
  })
  ### Catch (long-term median) per group
  Catch_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$catch[,,,,, which(SSB_groups %in% x)]/c(MSY)
    median(tmp)
  })
  ### get starting conditions
  SSB_levels <- sapply(SSB_levels, function(x) {
    x <- gsub(x = x, pattern = "\\(|\\]", replacement = "")
    x <- unlist(strsplit(x, split = ","))
    mean(as.numeric(x))
  })
  pos_remove <- which(is.na(SSB_levels))
  stats_risk_depletion <- data.frame(
    risk_Blim = unlist(risk_group)[-pos_remove],
    SSB_rel = unlist(SSB_group)[-pos_remove],
    Catch_rel = unlist(Catch_group)[-pos_remove],
    SSB0_rel = unlist(SSB_levels)[-pos_remove],
    n_iter_part = unlist(group_n)[-pos_remove])
  row.names(stats_risk_depletion) <- NULL
  stats_risk_depletion <- stats_risk_depletion[order(stats_risk_depletion$SSB0_rel), ]
  stats_risk_depletion <- stats_risk_depletion %>%
    mutate(stock = stk)
  
  return(stats_risk_depletion)
  
}

df_plot <- res_df %>%
  pivot_longer(c(risk_Blim, SSB_rel, Catch_rel)) %>%
  mutate(name = factor(name, levels = c("SSB_rel", "risk_Blim", "Catch_rel"),
                       labels = c("SSB/B[MSY]", "B[lim]~risk", "Catch/MSY"))) %>%
  # left_join(stocks %>%
  #             select(stock, k)) %>%
  mutate(stock_name = factor(stock,
                             levels = c("ang3", "pol", "tur", "her"),
                             labels = c("anglerfish", "pollack",
                                        "turbot", "herring")))
  

p <- df_plot %>%
  ggplot(aes(x = SSB0_rel, y = value, 
             colour = stock_name, fill = stock_name, linetype = stock_name)) +
  stat_smooth(n = 100, span = 0.25, se = FALSE, geom = "line", linewidth = 0.4) +
  geom_point(size = 0.15, stroke = 0, shape = 21) +
  scale_colour_brewer("", palette = "Set1") +
  scale_fill_brewer("", palette = "Set1") +
  scale_linetype("") +
  facet_wrap(~ name, scales = "free_y", labeller = "label_parsed",
             strip.position = "left") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, NA)) +
  labs(x = expression(SSB[y == 0]/B[MSY])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines")
  )
p

ggsave(filename = "output/plots/constant_catch/10k_depletion.png",
       type = "cairo-png", plot = p,
       width = 16, height = 4.5, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### constant_catch - plot time series example ####
### ------------------------------------------------------------------------ ###
### pollack, random, 10k iterations

brp <- brps$pol
### quantiles
stk <- readRDS("input/10000_100/OM/random/pol/stk.rds")
res_mp <- readRDS(paste0("output/constant_catch/10000_100/baseline/random/pol/",
                         "mp_0.1_0_0.1_0.6_0_0.75_0_.rds"))
stk[, ac(101:200)] <- res_mp@om@stock[, ac(101:200)]

### collapse correction
stk_qnt <- collapse_correction(stk = stk, yrs = 101:200,
                               quants = c("ssb", "catch"))
stk_qnt$ssb <- window(stk_qnt$ssb, start = 1)
stk_qnt$catch <- window(stk_qnt$catch, start = 1)
stk_qnt$ssb[, ac(1:100)] <- ssb(stk)[, ac(1:100)]
stk_qnt$catch[, ac(1:100)] <- catch(stk)[, ac(1:100)]
### relative to MSY
stk_qnt_rel <- stk_qnt
stk_qnt_rel$catch <- stk_qnt$catch/c(refpts(brp)["msy", "yield"])
stk_qnt_rel$ssb <- stk_qnt$ssb/c(refpts(brp)["msy", "ssb"])
### get quants
qnts <- list(ssb = stk_qnt_rel$ssb,
             catch = stk_qnt_rel$catch,
             risk = apply(stk_qnt$ssb < attr(brp, "Blim"), 2, mean))
### get quantiles
qnts <- lapply(qnts, function(x) {
  quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95))
})
qnts <- as(qnts, "FLQuants")
df <- as.data.frame(qnts)


### some iterations
qnts_iter <- lapply(stk_qnt_rel, function(x) {
  iter(x, 1:5)
})
qnts_iter <- as(qnts_iter, "FLQuants")
df_iter <- as.data.frame(qnts_iter)

### prepare for plotting
df_plot_qnts <- df %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "catch", "risk"),
                        labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")),
         year = year - 100)
df_plot_iter <- df_iter %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "catch", "risk"),
                        labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")),
         year = year - 100)

df_plot_qnts %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.25) + 
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.5) +
  geom_line(aes(x = year, y = `50%`), show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "grey", linewidth = 0.3, 
             linetype = "dotted") +
  geom_line(data = df_plot_iter,
            aes(x = year, y = data, colour = iter),
            show.legend = FALSE, size = 0.2, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  facet_wrap( ~ qname, scales = "free_y", labeller = "label_parsed", 
              strip.position = "left", ncol = 1) +
  coord_cartesian(xlim = c(0, 100)) + 
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank())

ggsave(filename = "output/plots/constant_catch/pol_timeseries.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### CL - new HCR ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### CL - risk & projection time - all stocks ####
### ------------------------------------------------------------------------ ###
### 500 iters, 100 years, all stocks, one-way & random

stats_risk <- foreach(stock = stocks$stock, k = stocks$k, 
                      .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows)  %do% {
    #browser()
    res <- readRDS(paste0("output/CL/500_100/default/", fhist,"/",
                          stock, "/mp_0.1_0_0.1_0.6_0_0.75_0_",
                          "_2_0.1_0.2_0.2_0.1_0.05_0.01_0.1_1_1_0.4_0.rds"))
    ### collapse correction
    res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
    
    ### get fishing history
    stk <- readRDS(paste0("input/500_100/OM/", fhist, "/", stock, "/stk.rds"))
    ssb_ <- ssb(stk)
    catch_ <- catch(stk)
    ### add simulated values
    ssb_[, ac(101:200)] <- res_corrected$ssb
    catch_[, ac(101:200)] <- res_corrected$catch
    
    ### metrics (relative)
    ssb_rel <- ssb_/c(refpts(brps[[stock]])["msy", "ssb"])
    catch_rel <- catch_/c(refpts(brps[[stock]])["msy", "yield"])
    risk <- apply(ssb_ < attr(brps[[stock]], "Blim"), 2, mean)
    risk[, ac(1:100)] <- NA
    
    ### cumulative values (year 101 to year x)
    ssb_rel_cum <- iterMeans(ssb_rel) %=% NA_real_
    ssb_rel_cum[, ac(101:200)] <- sapply(1:100, function(x) {
      median(c(ssb_rel[, ac(seq(from = 101, length.out = x))]), 
             na.rm = TRUE)
    })
    catch_rel_cum <- iterMeans(catch_rel) %=% NA_real_
    catch_rel_cum[, ac(101:200)] <- sapply(1:100, function(x) {
      median(c(catch_rel[, ac(seq(from = 101, length.out = x))]), 
             na.rm = TRUE)
    })
    risk_cum <- iterMeans(ssb_rel) %=% NA_real_
    risk_cum[, ac(101:200)] <- sapply(1:100, function(x) {
      mean(c(ssb_[, ac(seq(from = 101, length.out = x))] < attr(brps[[stock]], "Blim")), 
           na.rm = TRUE)
    })
    
    ### annual values (time series)
    ssb_rel <- iterMedians(ssb_rel)
    catch_rel <- iterMedians(catch_rel)
    
    ### combine data into data.frame
    df_tmp <- as.data.frame(FLQuants(ssb_annual = ssb_rel, 
                                     ssb_cum = ssb_rel_cum,
                                     catch_annual = catch_rel, 
                                     catch_cum = catch_rel_cum,
                                     risk_annual = risk, 
                                     risk_cum = risk_cum))
    df_tmp <- df_tmp %>% 
      select(year, data, qname) %>%
      separate(col = qname, sep = "_", into = c("metric", "calculation"))
    df_tmp$fhist <- fhist
    df_tmp$stock <- stock
    df_tmp$k <- k
    return(df_tmp)
}

stats_risk <- stats_risk %>%
  mutate(metric = factor(metric,
                         levels = c("ssb", "catch", "risk"),
                         labels = c("SSB/B[MSY]",
                                    "Catch/MSY",
                                    "B[lim]~risk")),
         fhist = factor(fhist,
                        levels = c("one-way", "random")),
         calculation = factor(calculation,
                              levels = c("annual", "cum"),
                              labels = c("Annual", "Cumulative")))
### calculate median over all stocks
stats_risk <- bind_rows(
  stats_risk %>%
    mutate(source = "stocks"),
  stats_risk %>%
    group_by(year, metric, calculation, fhist) %>%
    summarise(data = median(data)) %>%
    mutate(source = "median")) %>%
  mutate(source = factor(source, 
                         levels = c("stocks", "median")))


res_def_colours <- 
  c("one-way" = RColorBrewer::brewer.pal(n = 4, name = "Set1")[1], 
    "random" = RColorBrewer::brewer.pal(n = 4, name = "Set1")[2])

ggplot() +
  geom_line(data = stats_risk %>% filter(source == "stocks"),
            aes(x = year - 100, y = data,
                group = interaction(stock, fhist, source), colour = fhist),
            linewidth = 0.1, linetype = "dashed") +
  geom_line(data = stats_risk %>% filter(source == "median"),
            aes(x = year - 100, y = data, 
                colour = fhist),
            linewidth = 0.5, linetype = "solid") +
  facet_grid(metric ~ calculation, labeller = "label_parsed", switch = "y", 
             scales = "free_y") +
  scale_colour_manual("fishing history", values = res_def_colours) +
  theme_bw(base_size = 8) +
  #xlim(c(-2, NA)) +
  coord_cartesian(xlim = c(-2, 102), ylim = c(0, NA), expand = FALSE) +
  labs(x = "Year") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.6),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"))

### annual only
ggplot() +
  geom_line(data = stats_risk %>% 
              filter(source == "stocks" & calculation == "Annual"),
            aes(x = year - 100, y = data,
                group = interaction(stock, fhist, source), colour = fhist),
            linewidth = 0.1, linetype = "dashed") +
  geom_line(data = stats_risk %>% 
              filter(source == "median" & calculation == "Annual"),
            aes(x = year - 100, y = data, 
                colour = fhist),
            linewidth = 0.5, linetype = "solid") +
  facet_wrap(~ metric, labeller = "label_parsed", strip.position = "left", 
             scales = "free_y", ncol = 1) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  theme_bw(base_size = 8) +
  coord_cartesian(xlim = c(-2, 102), ylim = c(0, NA), expand = FALSE) +
  labs(x = "Year") +
  theme(strip.placement.y = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        strip.text.x = element_text(size = 8),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.55),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"))

ggsave(filename = "output/plots/CL/CL_all_annual_timeseries.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### CL - plot time series example ####
### ------------------------------------------------------------------------ ###
### pollack, random, 10k iterations

brp <- brps$pol
### quantiles
stk <- readRDS("input/10000_100/OM/random/pol/stk.rds")
res_mp <- readRDS(paste0("output/CL/10000_100/baseline/random/pol/",
                         "mp_0.1_0_0.1_0.6_0_0.75_0_",
                         "_2_0.1_0.2_0.2_0.1_0.05_0.01_0.1_1_1_0.4_0.rds"))
stk[, ac(101:200)] <- res_mp@om@stock[, ac(101:200)]

### collapse correction
stk_qnt <- collapse_correction(stk = stk, yrs = 101:200,
                               quants = c("ssb", "catch"))
stk_qnt$ssb <- window(stk_qnt$ssb, start = 1)
stk_qnt$catch <- window(stk_qnt$catch, start = 1)
stk_qnt$ssb[, ac(1:100)] <- ssb(stk)[, ac(1:100)]
stk_qnt$catch[, ac(1:100)] <- catch(stk)[, ac(1:100)]
### relative to MSY
stk_qnt_rel <- stk_qnt
stk_qnt_rel$catch <- stk_qnt$catch/c(refpts(brp)["msy", "yield"])
stk_qnt_rel$ssb <- stk_qnt$ssb/c(refpts(brp)["msy", "ssb"])
### get quants
qnts <- list(ssb = stk_qnt_rel$ssb,
             catch = stk_qnt_rel$catch,
             risk = apply(stk_qnt$ssb < attr(brp, "Blim"), 2, mean))
### get quantiles
qnts <- lapply(qnts, function(x) {
  quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95))
})
qnts <- as(qnts, "FLQuants")
df <- as.data.frame(qnts)


### some iterations
qnts_iter <- lapply(stk_qnt_rel, function(x) {
  iter(x, 1:5)
})
qnts_iter <- as(qnts_iter, "FLQuants")
df_iter <- as.data.frame(qnts_iter)

### prepare for plotting
df_plot_qnts <- df %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "catch", "risk"),
                        labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")),
         year = year - 100)
df_plot_iter <- df_iter %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "catch", "risk"),
                        labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")),
         year = year - 100)

df_plot_qnts %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.25) + 
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.5) +
  geom_line(aes(x = year, y = `50%`), show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "grey", linewidth = 0.3, 
             linetype = "dotted") +
  geom_line(data = df_plot_iter,
            aes(x = year, y = data, colour = iter),
            show.legend = FALSE, linewidth = 0.2, alpha = 0.5) +
  scale_colour_brewer(palette = "Dark2") +
  facet_wrap( ~ qname, scales = "free_y", labeller = "label_parsed", 
              strip.position = "left", ncol = 1) +
  coord_cartesian(xlim = c(0, 100)) + 
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank())

ggsave(filename = "output/plots/CL/pol_timeseries.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### CL - risk vs depletion ####
### ------------------------------------------------------------------------ ###
### random fhist, 10,000 iterations, 100 years
### examples: pollack, turbot, herring, anglerfish

res_df <- foreach(stk = c("pol", "tur", "ang3", "her"), 
                  .combine = bind_rows) %do% {
  #browser()
  
  ### reference points
  Blim <- attr(brps[[stk]], "Blim")
  Bmsy <- brps[[stk]]@refpts["msy", "ssb"]
  MSY <- brps[[stk]]@refpts["msy", "yield"]
  
  ### MP output
  res <- readRDS(paste0("output/CL/10000_100/baseline/random/",
                        stk, "/mp_0.1_0_0.1_0.6_0_0.75_0_",
                        "_2_0.1_0.2_0.2_0.1_0.05_0.01_0.1_1_1_0.4_0.rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
  
  ### starting condition
  SSBs0 <- ssb(res@om@stock)[, ac(100)]
  SSBs0 <- SSBs0/c(Bmsy)
  SSBs0 <- c(SSBs0)
  max_SSB <- max(SSBs0) #3
  SSB_breaks <- seq(from = 0, to = max_SSB, by = 0.1)
  SSB_groups <- cut(SSBs0, breaks = SSB_breaks)
  SSB_levels <- unique(as.character(SSB_groups))
  ### number of replicates per group
  group_n <- sapply(SSB_levels, function(x) {
    length(which(SSB_groups %in% x))
  })
  group_n[sort(names(group_n))]
  # sum(group_n)
  ### Blim risk per group
  ### SSB is on absolute scale 
  risk_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]
    mean(tmp < Blim)
  })
  ### SSB (long-term median) per group
  SSB_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]/c(Bmsy)
    median(tmp)
  })
  ### Catch (long-term median) per group
  Catch_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$catch[,,,,, which(SSB_groups %in% x)]/c(MSY)
    median(tmp)
  })
  ### get starting conditions
  SSB_levels <- sapply(SSB_levels, function(x) {
    x <- gsub(x = x, pattern = "\\(|\\]", replacement = "")
    x <- unlist(strsplit(x, split = ","))
    mean(as.numeric(x))
  })
  pos_remove <- which(is.na(SSB_levels))
  stats_risk_depletion <- data.frame(
    risk_Blim = unlist(risk_group)[-pos_remove],
    SSB_rel = unlist(SSB_group)[-pos_remove],
    Catch_rel = unlist(Catch_group)[-pos_remove],
    SSB0_rel = unlist(SSB_levels)[-pos_remove],
    n_iter_part = unlist(group_n)[-pos_remove])
  row.names(stats_risk_depletion) <- NULL
  stats_risk_depletion <- stats_risk_depletion[order(stats_risk_depletion$SSB0_rel), ]
  stats_risk_depletion <- stats_risk_depletion %>%
    mutate(stock = stk)
  
  return(stats_risk_depletion)
  
}

df_plot <- res_df %>%
  pivot_longer(c(risk_Blim, SSB_rel, Catch_rel)) %>%
  mutate(name = factor(name, levels = c("SSB_rel", "risk_Blim", "Catch_rel"),
                       labels = c("SSB/B[MSY]", "B[lim]~risk", "Catch/MSY"))) %>%
  # left_join(stocks %>%
  #             select(stock, k)) %>%
  mutate(stock_name = factor(stock,
                             levels = c("ang3", "pol", "tur", "her"),
                             labels = c("anglerfish", "pollack",
                                        "turbot", "herring")))


p <- df_plot %>%
  ggplot(aes(x = SSB0_rel, y = value, 
             colour = stock_name, fill = stock_name, linetype = stock_name)) +
  stat_smooth(n = 100, span = 0.25, se = FALSE, geom = "line", linewidth = 0.4) +
  geom_point(size = 0.15, stroke = 0, shape = 21) +
  scale_colour_brewer("", palette = "Set1") +
  scale_fill_brewer("", palette = "Set1") +
  scale_linetype("") +
  facet_wrap(~ name, scales = "free_y", labeller = "label_parsed",
             strip.position = "left") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, NA)) +
  labs(x = expression(SSB[y == 0]/B[MSY])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines")
  )
p

ggsave(filename = "output/plots/CL/10k_depletion.png",
       type = "cairo-png", plot = p,
       width = 16, height = 4.5, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### CL - example GA optimisation ####
### ------------------------------------------------------------------------ ###
### pollack, one-way, 500 iterations
### risk limit: 10%

### stats
stats_mp <- readRDS(paste0("output/CL/500_100/all_par/one-way/pol/",
                         "1_0.1_0.48_0.48_",
                         "0.48_0.35_0.25_0.38_1.53_1.16_0.59_0.3.rds"))

brp <- brps$pol
### optimised MP
res_mp <- readRDS(paste0("output/CL/500_100/all_par/one-way/pol/",
                         "mp_0.1_0_0.1_0.6_0_0.75_0_", "_1_0.1_0.48_0.48_",
                         "0.48_0.35_0.25_0.38_1.53_1.16_0.59_0.3.rds"))


### quantiles
stk <- readRDS("input/500_100/OM/one-way/pol/stk.rds")
stk[, ac(101:200)] <- res_mp@om@stock[, ac(101:200)]

### collapse correction
stk_qnt <- collapse_correction(stk = stk, yrs = 101:200,
                               quants = c("ssb", "catch"))
stk_qnt$ssb <- window(stk_qnt$ssb, start = 1)
stk_qnt$catch <- window(stk_qnt$catch, start = 1)
stk_qnt$ssb[, ac(1:100)] <- ssb(stk)[, ac(1:100)]
stk_qnt$catch[, ac(1:100)] <- catch(stk)[, ac(1:100)]
### relative to MSY
stk_qnt_rel <- stk_qnt
stk_qnt_rel$catch <- stk_qnt$catch/c(refpts(brp)["msy", "yield"])
stk_qnt_rel$ssb <- stk_qnt$ssb/c(refpts(brp)["msy", "ssb"])
### get quants
qnts <- list(ssb = stk_qnt_rel$ssb,
             catch = stk_qnt_rel$catch,
             risk = apply(stk_qnt$ssb < attr(brp, "Blim"), 2, mean))
### get quantiles
qnts <- lapply(qnts, function(x) {
  quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95))
})
qnts <- as(qnts, "FLQuants")
df <- as.data.frame(qnts)

### some iterations
qnts_iter <- lapply(stk_qnt_rel, function(x) {
  iter(x, 1:5)
})
qnts_iter <- as(qnts_iter, "FLQuants")
df_iter <- as.data.frame(qnts_iter)

### prepare for plotting
df_plot_qnts <- df %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "catch", "risk"),
                        labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")),
         year = year - 100)
df_plot_iter <- df_iter %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "catch", "risk"),
                        labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")),
         year = year - 100)
df_plot_risk <- df_plot_qnts %>%
  select(data = `50%`, qname) %>%
  mutate(data = NA) %>%
  unique() %>%
  mutate(data = ifelse(qname == "B[lim]~risk", 0.05, NA))

df_plot_qnts %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.25) + 
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.5) +
  geom_line(aes(x = year, y = `50%`), show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "grey", linewidth = 0.3, 
             linetype = "dotted") +
  geom_line(data = df_plot_iter,
            aes(x = year, y = data, colour = iter),
            show.legend = FALSE, linewidth = 0.2, alpha = 0.5) +
  geom_hline(data = df_plot_risk, aes(yintercept = data),
             colour = "red", linetype = "1111",
             linewidth = 0.3) +
  scale_colour_brewer(palette = "Dark2") +
  facet_wrap( ~ qname, scales = "free_y", labeller = "label_parsed", 
              strip.position = "left", ncol = 1) +
  coord_cartesian(xlim = c(0, 100)) + 
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank())

ggsave(filename = "output/plots/CL/pol_ow_opt_timeseries.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### CL - pollack - exploration: 1 parameter at a time ####
### ------------------------------------------------------------------------ ###

### collate runs
runs <- foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {
  #browser()
  runs_i <- readRDS(paste0("output/CL/500_100/1par/", fhist, "/pol/runs.rds"))
  df_i <- lapply(runs_i, function(x) {
    bind_cols(data.frame(t(x$pars)), data.frame(t(x$stats)))
  })
  df_i <- do.call(bind_rows, df_i)
  df_i$fhist <- fhist
  return(df_i)
}
par_names <- c("interval", "lambda_lower", "lambda_upper", "gamma_lower",
               "gamma_upper", "r_threshold", "l_threshold", "f_threshold",
               "Lref_mult", "multiplier", "first_catch", "catch_limit")
names(runs)[1:12] <- par_names
par_default <- c(2, 0.2, 0.1, 0.2, 0.1, 0.05, 0.01, 0.1, 1, 1, 0.4, 0)

runs_groups <- foreach(par = par_names, .combine = rbind) %:%
  foreach(fhist = c("one-way", "random"), .combine = rbind) %do% {
    #browser()
    
    ### subset to fhist
    runs_i <- runs %>%
      filter(fhist == !!fhist)
    
    ### find parameters with default values
    pos <- which(par_names != par)
    ### filter these
    for (i in pos) {
      runs_i <- runs_i[runs_i[, i] == par_default[i], ]
      runs_i[, i] <- NA
    }
    
    runs_i$par_opt <- par
    
    return(unnest(runs_i, cols = everything()))
    
  }

df_plot <- runs_groups %>%
  #filter(par_opt == "interval" & fhist == "one-way") %>%
  select(fhist, par_opt, 1:12, SSB_rel, risk_Blim, Catch_rel) %>%
  mutate(par_opt = factor(par_opt, levels = par_names)) %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random"))) %>%
  pivot_longer(cols = c(SSB_rel, risk_Blim, Catch_rel),
               names_to = "metric", values_to = "value") %>%
  pivot_longer(cols = interval:catch_limit,
               names_to = "par_name", values_to = "par_value") %>%
  mutate(par_name = factor(par_name, levels = par_names)) %>%
  mutate(metric = factor(metric,
                         levels = c("SSB_rel", "risk_Blim", "Catch_rel"),
                         labels = c("SSB/B[MSY]", "B[lim]~risk", "Catch/MSY")))

df_refs <- df_plot %>%
  select(-value, -par_value) %>%
  unique() %>%
  full_join(data.frame(par_name = par_names, value_default = par_default)) %>%
  mutate(par_name = factor(par_name, levels = par_names))

p <- df_plot %>%
  #filter(par_name == "interval") %>%
  ggplot() +
  geom_vline(data = df_refs,
             aes(xintercept = value_default),
             colour = "darkgrey", linetype = "1111", linewidth = 0.3) +
  geom_line(data = df_plot,
            aes(x = par_value, y = value, colour = fhist),
            linewidth = 0.4) +
  scale_colour_brewer("Fishing history", palette = "Set1") +
  facet_grid(metric ~ par_name, scales = "free",
             labeller = label_parsed, switch = "y") +
  #scale_x_continuous(breaks = scales::pretty_breaks()) +
  labs(x = "Parameter value") +
  theme_bw(base_size = 7) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines"),
        legend.position = "bottom",
        panel.spacing.x = unit(0.5, "lines"))
p
ggsave(filename = "output/plots/CL/CL_pol_1par_stats.png",
       type = "cairo",
       width = 30, height = 12, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### more code - not used ####
### ------------------------------------------------------------------------ ###


# ### ------------------------------------------------------------------------ ###
# ### conditional PA buffer - pollack example ####
# ### ------------------------------------------------------------------------ ###
# ### apply buffer when Lmean <= LF=M
# ### every 3 years, 20%
# 
# stats_cond <- foreach(rule = c("default", "conditional"), 
#                       .combine = bind_rows) %do% {
#     #browser()
#     if (identical(rule, "default"))
#       res <- readRDS(paste0("output/const_catch/500_100/baseline/one-way/pol/",
#                             "/mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_",
#                             "0.6_0_0.75.rds"))
#     if (identical(rule, "conditional"))
#       res <- readRDS(paste0("output/const_catch/500_50/CC_f/one-way/pol/",
#                             "mp_CC_f_3_3_0.8_Inf_0_0.2_0.2_0_0_0.6_0_0.75.rds"))
#     ### collapse correction
#     res_corrected <- collapse_correction(stk = res@stock, yrs = 101:150)
#     
#     ### get fishing history
#     stk <- readRDS(paste0("input/500_50/OM/one-way/pol/stk.rds"))
#     ### combine with results
#     res_corrected <- lapply(res_corrected, window, start = 0)
#     res_corrected$catch[, ac(0:100)] <- catch(stk)[, ac(0:100)]
#     res_corrected$ssb[, ac(0:100)] <- ssb(stk)[, ac(0:100)]
#     res_corrected$fbar[, ac(0:100)] <- fbar(stk)[, ac(0:100)]
#     
#     ### relative metrics
#     res_corrected$catch <- res_corrected$catch/c(refpts(brps$pol)["msy", "yield"])
#     res_corrected$ssb <- res_corrected$ssb/c(refpts(brps$pol)["msy", "ssb"])
#     res_corrected$fbar <- res_corrected$fbar/c(refpts(brps$pol)["msy", "harvest"])
#     
#     ### quantiles
#     res_corrected <- lapply(res_corrected, quantile, 
#                             c(0.05, 0.25, 0.5, 0.75, 0.95))
#     res_corrected <- as(res_corrected, "FLQuants")
#     
#     df <- as.data.frame(res_corrected)
#     df$rule <- rule
#     
#     return(df)
# }
# 
# stats_cond <- stats_cond %>%
#   pivot_wider(names_from = iter, values_from = data) %>%
#   mutate(qname = factor(qname, 
#                         levels = c("ssb", "fbar", "catch"),
#                         labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY")),
#          rule = factor(rule, levels = c("default", "conditional"))) %>%
#   select(-age, - unit, -season, -area)
# 
# stats_cond %>%
#   mutate(year = year - 100) %>%
#   ggplot(aes(x = year, fill = rule)) +
#   geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.125) + 
#   geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.25) +
#   geom_line(aes(y = `50%`, colour = rule), show.legend = FALSE) +
#   geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
#   scale_colour_manual("PA buffer", 
#     values = c("default" = brewer.pal(3, name = "Dark2")[1],
#                "conditional" = brewer.pal(3, name = "Dark2")[2])) +
#   scale_fill_manual("PA buffer", 
#                     values = c("default" = brewer.pal(3, name = "Dark2")[1],
#                                "conditional" = brewer.pal(3, name = "Dark2")[2])) +
#   facet_wrap(~ qname, ncol = 1,
#              scales = "free", labeller = "label_parsed", 
#              switch = "y") +
#   labs(x = "Year") +
#   coord_cartesian(xlim = c(0, 50)) + 
#   theme_bw(base_size = 8) +
#   theme(strip.placement = "outside",
#         strip.background.y = element_blank(),
#         strip.text.y = element_text(size = 8),
#         axis.title.y = element_blank())
# ggsave(filename = "output/plots/constant_catch/CC_CC_f.png",
#        type = "cairo", 
#        width = 17, height = 8, units = "cm", dpi = 600)
# 
# ### ------------------------------------------------------------------------ ###
# ### conditional PA buffer - pollack example - optimised ####
# ### ------------------------------------------------------------------------ ###
# ### apply buffer when Lmean <= LF=M
# ### every 3 years, 20%
# ### + optimised
# 
# ### get optimised parameters
# ga_res <- readRDS(paste0("output/CC_f/500_50/CC_f_opt/one-way/pol/",
#                          "Lref_mult-pa_size-interval--obj_ICES_MSYPA_res.rds"))
# ga_res@solution
# 
# stats_cond <- foreach(rule = c("default", "conditional", "optimised"), 
#                       .combine = bind_rows) %do% {
#   #browser()
#   if (identical(rule, "default"))
#     res <- readRDS(paste0("output/const_catch/500_100/baseline/one-way/pol/",
#                           "/mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_",
#                           "0.6_0_0.75.rds"))
#   if (identical(rule, "conditional"))
#     res <- readRDS(paste0("output/const_catch/500_50/CC_f/one-way/pol/",
#                           "mp_CC_f_3_3_0.8_Inf_0_0.2_0.2_0_0_0.6_0_0.75.rds"))
#   if (identical(rule, "optimised"))
#     res <- readRDS(paste0("output/CC_f/500_50/CC_f_opt/one-way/pol/",
#                           "mp_1.9_0.49_5_1_Inf_0.rds"))
#   ### collapse correction
#   res_corrected <- collapse_correction(stk = res@stock, yrs = 101:150)
#   
#   ### get fishing history
#   stk <- readRDS(paste0("input/500_50/OM/one-way/pol/stk.rds"))
#   ### combine with results
#   res_corrected <- lapply(res_corrected, window, start = 0)
#   res_corrected$catch[, ac(0:100)] <- catch(stk)[, ac(0:100)]
#   res_corrected$ssb[, ac(0:100)] <- ssb(stk)[, ac(0:100)]
#   res_corrected$fbar[, ac(0:100)] <- fbar(stk)[, ac(0:100)]
#   
#   ### relative metrics
#   res_corrected$catch <- res_corrected$catch/c(refpts(brps$pol)["msy", "yield"])
#   res_corrected$ssb <- res_corrected$ssb/c(refpts(brps$pol)["msy", "ssb"])
#   res_corrected$fbar <- res_corrected$fbar/c(refpts(brps$pol)["msy", "harvest"])
#   
#   ### quantiles
#   res_corrected <- lapply(res_corrected, quantile, 
#                           c(0.05, 0.25, 0.5, 0.75, 0.95))
#   res_corrected <- as(res_corrected, "FLQuants")
#   
#   df <- as.data.frame(res_corrected)
#   df$rule <- rule
#   
#   return(df)
# }
# 
# stats_cond <- stats_cond %>%
#   pivot_wider(names_from = iter, values_from = data) %>%
#   mutate(qname = factor(qname, 
#                         levels = c("ssb", "fbar", "catch"),
#                         labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY")),
#          rule = factor(rule, 
#                        levels = c("default", "conditional", "optimised"))) %>%
#   select(-age, - unit, -season, -area)
# 
# stats_cond %>%
#   mutate(year = year - 100) %>%
#   ggplot(aes(x = year, fill = rule)) +
#   geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.125) + 
#   geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.25) +
#   geom_line(aes(y = `50%`, colour = rule), show.legend = FALSE) +
#   geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
#   scale_colour_manual("PA buffer", 
#     values = c("default" = brewer.pal(3, name = "Dark2")[1],
#                "conditional" = brewer.pal(3, name = "Dark2")[2],
#                "optimised" = brewer.pal(3, name = "Dark2")[3])) +
#   scale_fill_manual("PA buffer", 
#     values = c("default" = brewer.pal(3, name = "Dark2")[1],
#                "conditional" = brewer.pal(3, name = "Dark2")[2],
#                "optimised" = brewer.pal(3, name = "Dark2")[3])) +
#   facet_wrap(~ qname, ncol = 1,
#              scales = "free", labeller = "label_parsed", 
#              switch = "y") +
#   labs(x = "Year") +
#   coord_cartesian(xlim = c(0, 50)) + 
#   theme_bw(base_size = 8) +
#   theme(strip.placement = "outside",
#         strip.background.y = element_blank(),
#         strip.text.y = element_text(size = 8),
#         axis.title.y = element_blank())
# ggsave(filename = "output/plots/constant_catch/CC_CC_f_opt.png",
#        type = "cairo", 
#        width = 17, height = 8, units = "cm", dpi = 600)
# 
# 
# ### ------------------------------------------------------------------------ ###
# ### CL - pollack - default - without beta ####
# ### ------------------------------------------------------------------------ ###
# 
# res_rnd <- readRDS("../GA_MSE_cat456/output/CL/500_50/baseline/random/pol/mp_0.1_0_0.1_0.6_0_0.75_0.1__3_0.1_0.2_0.2_0.1_0.05_999_1.rds")
# res_ow <- readRDS("../GA_MSE_cat456/output/CL/500_50/baseline/one-way/pol/mp_0.1_0_0.1_0.6_0_0.75_0.1__3_0.1_0.2_0.2_0.1_0.05_999_1.rds")
# 
# plot(FLQuants(SSB = ssb(res_rnd@stock), 
#               Catch = catch(res_rnd@stock), 
#               Fbar = fbar(res_rnd@stock)), iter = 1:5) + 
#   theme_bw()
# ggsave(filename = "output/plots/CL/pol_default_rnd_without_beta.png",
#        type = "cairo", 
#        width = 12, height = 12, units = "cm", dpi = 600)
# plot(FLQuants(SSB = ssb(res_ow@stock), 
#               Catch = catch(res_ow@stock),
#               Fbar = fbar(res_ow@stock)), iter = 1:5) + 
#   theme_bw()
# ggsave(filename = "output/plots/CL/pol_default_ow_without_beta.png",
#        type = "cairo", 
#        width = 12, height = 12, units = "cm", dpi = 600)
# 
# ### ------------------------------------------------------------------------ ###
# ### CL - pollack - default - with beta ####
# ### ------------------------------------------------------------------------ ###
# 
# res_rnd <- readRDS("../GA_MSE_cat456/output/CL/500_50/baseline/random/pol/mp_0.1_0_0.1_0.6_0_0.75_0.1__3_0.1_0.2_0.2_0.1_0.05_0.1_1.rds")
# res_ow <- readRDS("../GA_MSE_cat456/output/CL/500_50/baseline/one-way/pol/mp_0.1_0_0.1_0.6_0_0.75_0.1__3_0.1_0.2_0.2_0.1_0.05_0.1_1.rds")
# 
# plot(FLQuants(SSB = ssb(res_rnd@stock), 
#               Catch = catch(res_rnd@stock), 
#               Fbar = fbar(res_rnd@stock)), iter = 1:5) + 
#   theme_bw()
# ggsave(filename = "output/plots/CL/pol_default_rnd_with_beta.png",
#        type = "cairo", 
#        width = 12, height = 12, units = "cm", dpi = 600)
# plot(FLQuants(SSB = ssb(res_ow@stock), 
#               Catch = catch(res_ow@stock),
#               Fbar = fbar(res_ow@stock)), iter = 1:5) + 
#   theme_bw()
# ggsave(filename = "output/plots/CL/pol_default_ow_with_beta.png",
#        type = "cairo", 
#        width = 12, height = 12, units = "cm", dpi = 600)
# 
# 
# ### optimised
# res_ow <- readRDS("../GA_MSE_cat456/output/CL/500_50/baseline/one-way/pol/mp_0.1_0_0.1_0.6_0_0.75_0.1__5_0.21_0.08_0.3_0.5_0_0_1.08.rds")
# plot(FLQuants(SSB = ssb(res_ow@stock), 
#               Catch = catch(res_ow@stock),
#               Fbar = fbar(res_ow@stock)), iter = 1:5) + 
#   theme_bw()
# ggsave(filename = "output/plots/CL/pol_default_ow_with_beta_opt.png",
#        type = "cairo", 
#        width = 12, height = 12, units = "cm", dpi = 600)

# 
# ### ------------------------------------------------------------------------ ###
# ###  ####
# ### ------------------------------------------------------------------------ ###
# 
# runs <- readRDS("output/CL/500_100/baseline/one-way/pol/runs.rds")
# runs <- runs[sapply(runs, function(x) {!is.null(x)})]
# 
# runs_fitness <- lapply(runs, function(x) {
#   #browser()
#   obj <- 0 
#   obj <- obj + sum(unlist(x$stats["Catch_rel", ]))
#   ### penalise risk above 5% - gradual
#   obj <- obj - sum(penalty(x = unlist(x$stats["risk_Blim", ]), 
#                            negative = FALSE, max = 1, inflection = 0.06, 
#                            steepness = 0.5e+3))
#   return(obj)
# })
# pos <- which.max(runs_fitness)
# pos
# runs_fitness[pos]
# runs[pos]
# 
# 


### ------------------------------------------------------------------------ ###
### old (CC) vs new (CL) - compare risk for all stocks  ####
### ------------------------------------------------------------------------ ###

stats_risk <- foreach(stock = stocks$stock, k = stocks$k, 
                      .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  foreach(hcr = c("constant_catch", "CL"), .combine = bind_rows) %do% {
    #browser()
    res <- readRDS(paste0(
      "output/", hcr, "/500_100/", 
      switch(hcr, 
             constant_catch = "baseline",
             CL = "default"),
      "/", fhist,"/", stock, "/",
      switch(hcr, 
             constant_catch = "mp_0.1_0_0.1_0.6_0_0.75_0_",
             CL = "mp_0.1_0_0.1_0.6_0_0.75_0__2_0.1_0.2_0.2_0.1_0.05_0.01_0.1_1_1_0.4_0"),
      ".rds"))
    
    ### collapse correction
    res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
    
    ### get fishing history
    stk <- readRDS(paste0("input/500_100/OM/", fhist, "/", stock, "/stk.rds"))
    qnts <- FLQuants(ssb = ssb(stk), catch = catch(stk), 
                     risk = ssb(stk) %=% NA_real_)
    ### add simulated values
    qnts$ssb[, ac(101:200)] <- res_corrected$ssb
    qnts$catch[, ac(101:200)] <- res_corrected$catch
    
    ### risk
    qnts$risk <- apply(qnts$ssb < attr(brps[[stock]], "Blim"), 2, mean)
    
    ### median for SSB and catch
    qnts$ssb <- iterMedians(qnts$ssb)
    qnts$catch <- iterMedians(qnts$catch)

    ### combine data into data.frame
    df_tmp <- as.data.frame(qnts)[, c("year", "data", "qname")]
    df_tmp$fhist <- fhist
    df_tmp$stock <- stock
    df_tmp$k <- k
    df_tmp$hcr <- hcr
    return(df_tmp)
}

stats_risk %>%
  filter(qname == "risk") %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random")),
         hcr = factor(hcr, levels = c("constant_catch", "CL"),
                      labels = c("old", "new"))) %>%
  ggplot(aes(x = year - 100, y = data, group = stock, colour = k)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid(hcr ~ fhist) +
  coord_cartesian(xlim = c(-5, 100), ylim = c(0, 1.05), expand = FALSE) +
  labs(x = "Year", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/CL/comparison_risk_all.png",
       type = "cairo",
       width = 20, height = 12, units = "cm", dpi = 600)

stats_risk %>%
  filter(qname == "ssb") %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random")),
         hcr = factor(hcr, levels = c("constant_catch", "CL"),
                      labels = c("old", "new"))) %>%
  ggplot(aes(x = year - 100, y = data, group = stock, colour = k)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid(hcr ~ fhist) +
  coord_cartesian(xlim = c(-5, 100), ylim = c(0, NA), expand = FALSE) +
  labs(x = "Year", y = "SSB") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/CL/comparison_SSB_all.png",
       type = "cairo",
       width = 20, height = 12, units = "cm", dpi = 600)
stats_risk %>%
  filter(qname == "catch") %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random")),
         hcr = factor(hcr, levels = c("constant_catch", "CL"),
                      labels = c("old", "new"))) %>%
  ggplot(aes(x = year - 100, y = data, group = stock, colour = k)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid(hcr ~ fhist) +
  coord_cartesian(xlim = c(-5, 100), ylim = c(0, 300), expand = FALSE) +
  labs(x = "Year", y = "Catch") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/CL/comparison_catch_all.png",
       type = "cairo",
       width = 20, height = 12, units = "cm", dpi = 600)




### change in risk between HCRs
p <- stats_risk %>%
  filter(qname == "risk") %>%
  pivot_wider(names_from = hcr, values_from = data) %>%
  mutate(risk_change = (CL/constant_catch - 1)*100) %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random"))) %>%
  ggplot(aes(x = year - 100, y = risk_change, group = stock, colour = k)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid( ~ fhist) +
  coord_cartesian(xlim = c(-5, 100), ylim = c(-100, 100), expand = FALSE) +
  labs(x = "Year", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)
p
ggsave(filename = "output/plots/CL/comparison_risk_all_comparison.png",
       type = "cairo",
       width = 12, height = 5, units = "cm", dpi = 600)
p + coord_cartesian(xlim = c(-4, 21), ylim = c(-100, 100), expand = FALSE)
ggsave(filename = "output/plots/CL/comparison_risk_all_comparison_zoom.png",
       type = "cairo",
       width = 12, height = 5, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### old (CC) vs new (CL) - compare risk for all stocks - 5-year trends  ####
### ------------------------------------------------------------------------ ###

stats_risk5 <- foreach(stock = stocks$stock, k = stocks$k, 
                      .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows) %:%
  foreach(hcr = c("constant_catch", "CL"), .combine = bind_rows) %do% {
    #browser()
    res <- readRDS(paste0(
      "output/", hcr, "/500_100/", 
      switch(hcr, 
             constant_catch = "baseline",
             CL = "default"),
      "/", fhist,"/", stock, "/",
      switch(hcr, 
             constant_catch = "mp_0.1_0_0.1_0.6_0_0.75_0_",
             CL = "mp_0.1_0_0.1_0.6_0_0.75_0__2_5_5_3_0.2_0.1_0.2_0.1_0.05_0.01_0.1_1_1_0.4_0"),
      ".rds"))
    
    ### collapse correction
    res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
    
    ### get fishing history
    stk <- readRDS(paste0("input/500_100/OM/", fhist, "/", stock, "/stk.rds"))
    qnts <- FLQuants(ssb = ssb(stk), catch = catch(stk), 
                     risk = ssb(stk) %=% NA_real_)
    ### add simulated values
    qnts$ssb[, ac(101:200)] <- res_corrected$ssb
    qnts$catch[, ac(101:200)] <- res_corrected$catch
    
    ### risk
    qnts$risk <- apply(qnts$ssb < attr(brps[[stock]], "Blim"), 2, mean)
    
    ### median for SSB and catch
    qnts$ssb <- iterMedians(qnts$ssb)
    qnts$catch <- iterMedians(qnts$catch)
    
    ### combine data into data.frame
    df_tmp <- as.data.frame(qnts)[, c("year", "data", "qname")]
    df_tmp$fhist <- fhist
    df_tmp$stock <- stock
    df_tmp$k <- k
    df_tmp$hcr <- hcr
    return(df_tmp)
}

stats_risk5 %>%
  filter(qname == "risk") %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random")),
         hcr = factor(hcr, levels = c("constant_catch", "CL"),
                      labels = c("old", "new"))) %>%
  ggplot(aes(x = year - 100, y = data, group = stock, colour = k)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid(hcr ~ fhist) +
  coord_cartesian(xlim = c(-5, 100), ylim = c(0, 1.05), expand = FALSE) +
  labs(x = "Year", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/CL/comparison5_risk_all.png",
       type = "cairo",
       width = 20, height = 12, units = "cm", dpi = 600)

stats_risk5 %>%
  filter(qname == "ssb") %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random")),
         hcr = factor(hcr, levels = c("constant_catch", "CL"),
                      labels = c("old", "new"))) %>%
  ggplot(aes(x = year - 100, y = data, group = stock, colour = k)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid(hcr ~ fhist) +
  coord_cartesian(xlim = c(-5, 100), ylim = c(0, NA), expand = FALSE) +
  labs(x = "Year", y = "SSB") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/CL/comparison5_SSB_all.png",
       type = "cairo",
       width = 20, height = 12, units = "cm", dpi = 600)
stats_risk5 %>%
  filter(qname == "catch") %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random")),
         hcr = factor(hcr, levels = c("constant_catch", "CL"),
                      labels = c("old", "new"))) %>%
  ggplot(aes(x = year - 100, y = data, group = stock, colour = k)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid(hcr ~ fhist) +
  coord_cartesian(xlim = c(-5, 100), ylim = c(0, 300), expand = FALSE) +
  labs(x = "Year", y = "Catch") +
  theme_bw(base_size = 8)
ggsave(filename = "output/plots/CL/comparison5_catch_all.png",
       type = "cairo",
       width = 20, height = 12, units = "cm", dpi = 600)




### change in risk between HCRs
p <- stats_risk5 %>%
  filter(qname == "risk") %>%
  pivot_wider(names_from = hcr, values_from = data) %>%
  mutate(risk_change = (CL/constant_catch - 1)*100) %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random"))) %>%
  ggplot(aes(x = year - 100, y = risk_change, group = stock, colour = k)) +
  geom_line() +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid( ~ fhist) +
  coord_cartesian(xlim = c(-5, 100), ylim = c(-100, 100), expand = FALSE) +
  labs(x = "Year", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8)
p
ggsave(filename = "output/plots/CL/comparison5_risk_all_comparison.png",
       type = "cairo",
       width = 12, height = 5, units = "cm", dpi = 600)
p + coord_cartesian(xlim = c(-4, 21), ylim = c(-100, 100), expand = FALSE)
ggsave(filename = "output/plots/CL/comparison5_risk_all_comparison_zoom.png",
       type = "cairo",
       width = 12, height = 5, units = "cm", dpi = 600)
### ------------------------------------------------------------------------ ###
### old (CC) vs new (CL) - pollack - compare recovery vs. depletion  ####
### ------------------------------------------------------------------------ ###

res_df <- foreach(stk = c("pol", "tur", "ang3", "her"), 
                  .combine = bind_rows) %:%
  foreach(hcr = c("constant_catch", "CL"), .combine = bind_rows)  %do% {
  #browser()
  
  ### reference points
  Blim <- attr(brps[[stk]], "Blim")
  Bmsy <- brps[[stk]]@refpts["msy", "ssb"]
  MSY <- brps[[stk]]@refpts["msy", "yield"]
  
  ### MP output
  res <- readRDS(paste0(
    "output/", hcr, "/10000_100/", 
    switch(hcr, 
           constant_catch = "baseline",
           CL = "baseline"),
    "/random/", stk, "/",
    switch(hcr, 
           constant_catch = "mp_0.1_0_0.1_0.6_0_0.75_0_",
           CL = "mp_0.1_0_0.1_0.6_0_0.75_0__2_0.1_0.2_0.2_0.1_0.05_0.01_0.1_1_1_0.4_0"),
    ".rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
  
  ### starting condition
  SSBs0 <- ssb(res@om@stock)[, ac(100)]
  SSBs0 <- SSBs0/c(Bmsy)
  SSBs0 <- c(SSBs0)
  max_SSB <- max(SSBs0) #3
  SSB_breaks <- seq(from = 0, to = max_SSB, by = 0.25)
  SSB_groups <- cut(SSBs0, breaks = SSB_breaks)
  SSB_levels <- sort(unique(as.character(SSB_groups)))
  ### number of replicates per group
  group_n <- sapply(SSB_levels, function(x) {
    length(which(SSB_groups %in% x))
  })
  group_n
  # sum(group_n)
  ### Blim risk per group
  ### SSB is on absolute scale 
  df <- foreach(x = SSB_levels, .combine = bind_rows) %do% {
    tmp_ssb <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]
    tmp_catch <- res_corrected$catch[,,,,, which(SSB_groups %in% x)]
    
    df_tmp <- as.data.frame(iterMeans(tmp_ssb < Blim)) %>%
      select(year, risk = data)
    df_tmp$SSB <- c(iterMedians(tmp_ssb/c(Bmsy)))
    df_tmp$Catch <- c(iterMedians(tmp_catch/c(MSY)))
    df_tmp$depletion <- x
    return(df_tmp)
  }
  df$hcr <- hcr
  df$stk <- stk
  return(df)
}

p <- res_df %>%
  mutate(depletion = factor(depletion)) %>%
  mutate(hcr = factor(hcr, levels = c("constant_catch", "CL"),
                      labels = c("current approach", "new"))) %>%
  mutate(stk = factor(stk, levels = c("ang3", "pol", "tur", "her"))) %>%
  #filter(stk == "pol") %>%
  filter(depletion %in% c("(0,0.25]", "(0.25,0.5]", "(0.5,0.75]", "(0.75,1]")) %>%
  ggplot(aes(x = year - 100, y = risk, linetype = depletion)) +
  geom_line() +
  scale_linetype_discrete(name = expression(Initial~depletion~SSB/B[MSY])) +
  facet_grid(stk ~ hcr) +
  coord_cartesian(xlim = c(0, 100), ylim = c(-0.025, 1.05), expand = FALSE) + 
  labs(x = "Year", y = expression(B[lim]~risk)) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"))
p
ggsave(filename = "output/plots/CL/comparison_10k_risk.png",
       type = "cairo",
       width = 16, height = 10, units = "cm", dpi = 600)
p + coord_cartesian(xlim = c(0, 26), ylim = c(-0.025, 1.05), expand = FALSE)
ggsave(filename = "output/plots/CL/comparison_10k_risk_zoom.png",
       type = "cairo",
       width = 16, height = 10, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### fcc - all stocks ####
### ------------------------------------------------------------------------ ###

stats_fcc <- foreach(stock = stocks$stock, k = stocks$k, 
                      .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows)  %dopar% {
    #browser()
    res <- readRDS(paste0(
      "output/fcc/500_100/default/", fhist,"/", stock, "/",
      "mp_0.1_0_0.1_0.6_0_0.75_0__2_1_1_0.5_1.1_0.7",
      ".rds"))
    
    ### collapse correction
    res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
    
    ### get fishing history
    stk <- readRDS(paste0("input/500_100/OM/", fhist, "/", stock, "/stk.rds"))
    qnts <- FLQuants(ssb = ssb(stk), catch = catch(stk), 
                     risk = ssb(stk) %=% NA_real_)
    ### add simulated values
    qnts$ssb[, ac(101:200)] <- res_corrected$ssb
    qnts$catch[, ac(101:200)] <- res_corrected$catch
    
    ### risk
    qnts$risk <- apply(qnts$ssb < attr(brps[[stock]], "Blim"), 2, mean)
    
    ### median for SSB and catch
    qnts$ssb <- iterMedians(qnts$ssb)
    qnts$catch <- iterMedians(qnts$catch)
    
    ### combine data into data.frame
    df_tmp <- as.data.frame(qnts)[, c("year", "data", "qname")]
    df_tmp$fhist <- fhist
    df_tmp$stock <- stock
    df_tmp$k <- k
    return(df_tmp)
}

stats_plot <- stats_fcc %>%
  filter(year >= 90) %>%
  #filter(qname == "risk") %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "random")),
         qname = factor(qname, levels = c("ssb", "catch", "risk"),
                        labels = c("SSB", "Catch", "B[lim]~risk")))
stats_plot %>%
  ggplot(aes(x = year - 100, y = data, group = stock, colour = k)) +
  geom_line() +
  geom_blank(data = stats_plot %>% mutate(data = data * 1.05)) +
  geom_vline(xintercept = 0) +
  scale_colour_gradientn("k/year", 
                         colours = scales::brewer_pal(palette = "Blues", 
                                                      direction = -1)(9),
                         values = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1)) +
  facet_grid(qname ~ fhist, scales = "free_y", switch = "y", 
             labeller = "label_parsed") +
  coord_cartesian(xlim = c(-5, 100), ylim = c(0, NA), expand = FALSE) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8))
ggsave(filename = "output/plots/fcc/comparison_all.png",
       type = "cairo",
       width = 20, height = 12, units = "cm", dpi = 600)

### plot pollack with iterations
res_pol <- readRDS("output/fcc/500_100/default/one-way/pol/mp_0.1_0_0.1_0.6_0_0.75_0__2_1_1_0.5_1.1_0.7.rds")
plot(stock(res_pol), iter = 1:5)
ggsave(filename = "output/plots/fcc/pol.png",
       type = "cairo",
       width = 20, height = 12, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### fcc - pollack - multiplier ####
### ------------------------------------------------------------------------ ###

res <- foreach(multiplier = seq(0.5, 1.4, 0.1),
               .combine = bind_rows) %do% {
  res_i <- readRDS(paste0("output/fcc/500_100/default/one-way/pol/mp_0.1_0_0.1_0.6_0_0.75_0__2_", multiplier, "_1_0.5_1.1_0.7.rds"))
  Blim <- attr(brps[["pol"]], "Blim")
  Bmsy <- c(refpts(brps[["pol"]])["msy", "ssb"])
  Cmsy <- c(refpts(brps[["pol"]])["msy", "yield"])
  risk_i <- mean(iterMeans(ssb(stock(om(res_i)))[, ac(101:200)] < Blim))
  ssb_i <- median(ssb(stock(om(res_i)))[, ac(101:200)]/Bmsy)
  catch_i <- median(catch(stock(om(res_i)))[, ac(101:200)]/Cmsy)
  data.frame(multiplier = multiplier, 
             risk = risk_i,
             ssb = ssb_i,
             catch = catch_i)
}

### ------------------------------------------------------------------------ ###
### fcc - risk vs depletion ####
### ------------------------------------------------------------------------ ###
### random fhist, 10,000 iterations, 100 years
### examples: pollack, turbot, herring, anglerfish

res_df <- foreach(stk = c("ang3","pol", "tur", "her"),
                  .combine = bind_rows) %do% {
  #browser()
  
  ### reference points
  Blim <- attr(brps[[stk]], "Blim")
  Bmsy <- brps[[stk]]@refpts["msy", "ssb"]
  MSY <- brps[[stk]]@refpts["msy", "yield"]
  
  ### MP output
  res <- readRDS(paste0("output/fcc/10000_100/default/random/",
                        stk, "/mp_0.1_0_0.1_0.6_0_0.75_0_",
                        "_2_1_1_0.5_1.1_0.7.rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
  
  ### starting condition
  SSBs0 <- ssb(res@om@stock)[, ac(100)]
  SSBs0 <- SSBs0/c(Bmsy)
  SSBs0 <- c(SSBs0)
  max_SSB <- max(SSBs0) #3
  SSB_breaks <- seq(from = 0, to = max_SSB, by = 0.1)
  SSB_groups <- cut(SSBs0, breaks = SSB_breaks)
  SSB_levels <- unique(as.character(SSB_groups))
  ### number of replicates per group
  group_n <- sapply(SSB_levels, function(x) {
    length(which(SSB_groups %in% x))
  })
  group_n[sort(names(group_n))]
  # sum(group_n)
  ### Blim risk per group
  ### SSB is on absolute scale 
  risk_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]
    mean(tmp < Blim)
  })
  ### SSB (long-term median) per group
  SSB_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[,,,,, which(SSB_groups %in% x)]/c(Bmsy)
    median(tmp)
  })
  ### Catch (long-term median) per group
  Catch_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$catch[,,,,, which(SSB_groups %in% x)]/c(MSY)
    median(tmp)
  })
  ### same but for long-term (last 50 years)
  risk_group_long <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[, ac(151:200),,,, which(SSB_groups %in% x)]
    mean(tmp < Blim)
  })
  SSB_group_long <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[, ac(151:200),,,, which(SSB_groups %in% x)]/c(Bmsy)
    median(tmp)
  })
  Catch_group_long <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$catch[, ac(151:200),,,, which(SSB_groups %in% x)]/c(MSY)
    median(tmp)
  })
  ### get starting conditions
  SSB_levels <- sapply(SSB_levels, function(x) {
    x <- gsub(x = x, pattern = "\\(|\\]", replacement = "")
    x <- unlist(strsplit(x, split = ","))
    mean(as.numeric(x))
  })
  pos_remove <- which(is.na(SSB_levels))
  stats_risk_depletion <- data.frame(
    riskBlim_all = unlist(risk_group)[-pos_remove],
    SSBrel_all = unlist(SSB_group)[-pos_remove],
    Catchrel_all = unlist(Catch_group)[-pos_remove],
    riskBlim_long = unlist(risk_group_long)[-pos_remove],
    SSBrel_long = unlist(SSB_group_long)[-pos_remove],
    Catchrel_long = unlist(Catch_group_long)[-pos_remove],
    SSB0_rel = unlist(SSB_levels)[-pos_remove],
    n_iter_part = unlist(group_n)[-pos_remove])
  row.names(stats_risk_depletion) <- NULL
  stats_risk_depletion <- stats_risk_depletion[order(stats_risk_depletion$SSB0_rel), ]
  stats_risk_depletion <- stats_risk_depletion %>%
    mutate(stock = stk)
  
  return(stats_risk_depletion)
  
}

res_df <- res_df %>%
  pivot_longer(cols = riskBlim_all:Catchrel_long,
               names_to = c("name", "period"),
               names_sep = "_",
               values_to = "value") %>%
  mutate(name = factor(name, levels = c("SSBrel", "riskBlim", "Catchrel"),
                       labels = c("SSB/B[MSY]", "B[lim]~risk", "Catch/MSY"))) %>%
  mutate(stock_name = factor(stock,
                             levels = c("ang3", "pol", "tur", "her"),
                             labels = c("anglerfish", "pollack",
                                        "turbot", "herring")))

p <- res_df %>%
  filter(period == "long") %>%
  ggplot(aes(x = SSB0_rel, y = value, 
             colour = stock_name, fill = stock_name, linetype = stock_name)) +
  stat_smooth(n = 100, span = 0.2, se = FALSE, geom = "line", linewidth = 0.4) +
  geom_point(size = 0.15, stroke = 0, shape = 21) +
  scale_colour_brewer("", palette = "Set1") +
  scale_fill_brewer("", palette = "Set1") +
  scale_linetype("") +
  facet_wrap(~ name, scales = "free_y", labeller = "label_parsed",
             strip.position = "left") +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, NA)) +
  labs(x = expression(SSB[y == 0]/B[MSY])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.6, "lines")
  )
p

ggsave(filename = "output/plots/fcc/default_10k_depletion_long.png",
       type = "cairo-png", plot = p,
       width = 16, height = 6, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### fcc - multiplier tuning ####
### ------------------------------------------------------------------------ ###

if (FALSE) {
  res <- foreach(fhist = c("one-way", "roller-coaster", "random"),
                 .combine = bind_rows, .errorhandling = "remove") %:%
    foreach(stock = stocks$stock, k = stocks$k, .combine = bind_rows, 
            .errorhandling = "remove") %:% 
    foreach(multiplier = seq(0, 2, 0.05), .combine = bind_rows, 
            .errorhandling = "remove") %dopar% {
      
      #browser()
      file_i <- paste0("output/fcc/500_100/default/", fhist, "/", stock, 
                       "/mp_0.1_0_0.1_0.6_0_0.75_0__2_",
                       multiplier, "_1_0.5_1.1_0.7.rds")
      if (!file.exists(file_i)) return(NULL)
      e <- try({stk_i <- stock(om(readRDS(file_i)))}, silent = TRUE)
      if (is(e, "try-error")) return(NULL)
      qnts <- collapse_correction(stk = stk_i, yrs = 101:200)
      
      ### reference points
      Blim <- attr(brps[[stock]], "Blim")
      Bmsy <- brps[[stock]]@refpts["msy", "ssb"]
      Fmsy <- brps[[stock]]@refpts["msy", "harvest"]
      MSY <- brps[[stock]]@refpts["msy", "yield"]
      
      df_i <- data.frame(
        stock = stock, fhist = fhist, multiplier = multiplier,
        SSBrel_all = median(qnts$ssb[, ac(101:200)]/c(Bmsy)),
        SSBrel_long = median(qnts$ssb[, ac(151:200)]/c(Bmsy)),
        Frel_all = median(qnts$fbar[, ac(101:200)]/c(Bmsy)),
        Frel_long = median(qnts$fbar[, ac(151:200)]/c(Bmsy)),
        Catchrel_all = median(qnts$catch[, ac(101:200)]/c(MSY)),
        Catchrel_long = median(qnts$catch[, ac(151:200)]/c(MSY)),
        riskBlim_all = mean(qnts$ssb[, ac(101:200)] < Blim),
        riskBlim_long = mean(qnts$ssb[, ac(151:200)] < Blim)
      )
      return(df_i)
  }
  saveRDS(res, "output/fcc/mult_smry.rds")
}
res <- readRDS("output/fcc/mult_smry.rds")

### stocks: jnd/san/sar implausible - Fcrash = NA
res <- res %>%
  filter(stock %in% stocks$stock[c(1:25, 28)]) %>%
  filter(multiplier > 0)


res_plot <- res %>%
  filter(stock %in% stocks$stock[c(1:24)]) %>%
  pivot_longer(cols = SSBrel_all:riskBlim_long,
               names_to = c("name", "period"),
               names_sep = "_",
               values_to = "value") %>%
  mutate(name = factor(name, levels = c("SSBrel", "Frel", "riskBlim", "Catchrel"),
                       labels = c("SSB/B[MSY]", "F[MSY]", "B[lim]~risk", 
                                  "Catch/MSY"))) %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "roller-coaster",
                                          "random"))) %>%
  left_join(stocks %>% select(stock, k)) %>%
  mutate(stock = factor(stock),
         period = factor(period)) %>%
  mutate(stock2 = paste0(stock, "-", fhist))

res_plot <- bind_rows(
  res_plot, 
  res_plot %>%
    mutate(fhist = "all")) %>%
  mutate(fhist = factor(fhist, levels = c("one-way", "roller-coaster",
                                          "random", "all")))

res_plot <- bind_rows(
  res_plot %>%
    mutate(type = "stocks"),
  res_plot %>%
    group_by(fhist, multiplier, period, name) %>%
    summarise(value = mean(value)) %>%
    mutate(type = "average"))

p_risk <- res_plot %>%
  filter(period == "long") %>%
  filter(name == "B[lim]~risk") %>%
  ggplot() +
  geom_vline(xintercept = 1, colour = "#444", linetype = "1111") +
  geom_line(aes(x = multiplier, y = value, group = stock2, 
                colour = k, linewidth = type, linetype = type)) +
  scale_linewidth_manual("", values = c(stocks = 0.2, average = 0.5)) +
  scale_linetype_manual("", values = c(stocks = "1111", average = "solid")) +
  scale_colour_gradient(expression("K (yr"^{-1}*")"), na.value = "black") +
  geom_hline(yintercept = 0.05, colour = "red", linetype = "2222") +
  facet_wrap(~ fhist, nrow = 1) + 
  labs(y = expression(B[lim]~risk)) +
  coord_cartesian(xlim = c(0, 1.99), ylim = c(0, 0.5), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
p_catch <- res_plot %>%
  filter(period == "long") %>%
  filter(name == "Catch/MSY") %>%
  ggplot() +
  geom_vline(xintercept = 1, colour = "#444", linetype = "1111") +
  geom_line(aes(x = multiplier, y = value, group = stock2, 
                colour = k, linewidth = type, linetype = type)) +
  scale_linewidth_manual("", values = c(stocks = 0.2, average = 0.5)) +
  scale_linetype_manual("", values = c(stocks = "1111", average = "solid")) +
  scale_colour_gradient(expression("K (yr"^{-1}*")"), na.value = "black") +
  facet_wrap(~ fhist, nrow = 1) + 
  labs(y = "Catch/MSY") +
  coord_cartesian(xlim = c(0, 1.99), ylim = c(0, 1), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_blank())
p_ssb <- res_plot %>%
  filter(period == "long") %>%
  filter(name == "SSB/B[MSY]") %>%
  ggplot() +
  geom_vline(xintercept = 1, colour = "#444", linetype = "1111") +
  geom_line(aes(x = multiplier, y = value, group = stock2, 
                colour = k, linewidth = type, linetype = type)) +
  scale_linewidth_manual("", values = c(stocks = 0.2, average = 0.5)) +
  scale_linetype_manual("", values = c(stocks = "1111", average = "solid")) +
  scale_colour_gradient(expression("K (yr"^{-1}*")"), na.value = "black") +
  facet_wrap(~ fhist, nrow = 1) + 
  labs(y = expression(SSB/B[MSY]), x = "Multiplier") +
  coord_cartesian(xlim = c(0, 1.99), ylim = c(0, NA), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.6, "lines"),
        strip.text = element_blank())

(p <- p_risk + p_catch + p_ssb + 
  plot_layout(ncol = 1, guides = "collect"))
ggsave(filename = "output/plots/fcc/mult_tuning_fhist.png",
       type = "cairo-png", plot = p,
       width = 16, height = 10, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### fcc vs ICES - trajectories ####
### ------------------------------------------------------------------------ ###
### use random fhist with 10k iterations, 100 years
### use B/Bmsy~1 for comparison
### pollack

res_df <- foreach(stk = c("pol"), .combine = bind_rows) %:%
  foreach(MP = c("fcc", "constant_catch"), .combine = bind_rows) %do% {
    #browser()
    
    ### reference points
    Blim <- attr(brps[[stk]], "Blim")
    Bmsy <- brps[[stk]]@refpts["msy", "ssb"]
    MSY <- brps[[stk]]@refpts["msy", "yield"]
    
    ### MP output
    if (identical(MP, "fcc")) {
      res <- readRDS(paste0("output/fcc/10000_100/default/random/",
                            stk, "/mp_0.1_0_0.1_0.6_0_0.75_0_",
                            "_2_1_1_0.5_1.1_0.7.rds"))
    } else if (identical(MP, "constant_catch")) {
      res <- readRDS(paste0("output/constant_catch/10000_100/baseline/random/",
                            stk, "/mp_0.1_0_0.1_0.6_0_0.75_0_.rds"))
    }

    ### collapse correction
    res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
    
    ### starting condition
    SSBs0 <- ssb(res@om@stock)[, ac(100)]
    SSBs0 <- SSBs0/c(Bmsy)
    SSBs0 <- c(SSBs0)
    
    its_Bmsy <- which(SSBs0 >= 0.9 & SSBs0 <= 1.1)
    its_thirdBmsy <- which(SSBs0 >= 0.2 & SSBs0 <= 0.4)
    
    qnts <- FLQuants(
      SSB_Bmsy = res_corrected$ssb[,,,,, its_Bmsy]/c(Bmsy),
      SSB_thirdBmsy = res_corrected$ssb[,,,,, its_thirdBmsy]/c(Bmsy),
      Catch_Bmsy = res_corrected$catch[,,,,, its_Bmsy]/c(MSY),
      Catch_thirdBmsy = res_corrected$catch[,,,,, its_thirdBmsy]/c(MSY),
      risk_Bmsy = iterMeans(res_corrected$ssb[,,,,, its_Bmsy] < Blim),
      risk_thirdBmsy = iterMeans(res_corrected$ssb[,,,,, its_thirdBmsy] < Blim)
    )
    qnts <- lapply(qnts, iterMedians)
    df <- as.data.frame(qnts)[, c("year", "data", "qname")]
    df$stock <- stk
    df$MP <- MP
    return(df)
  
}

res_df_plot <- res_df %>%
  separate(qname, into = c("quant", "level"), sep = "_") %>%
  mutate(quant = factor(quant, levels = c("Catch", "SSB", "risk"),
                        labels = c("Catch/MSY", "SSB/B[MSY]", "B[lim]~risk")),
         level = factor(level, levels = c("Bmsy", "thirdBmsy"),
                        labels = c("Sustainable", "Depleted")),
         MP = factor(MP, levels = c("fcc", "constant_catch"),
                     labels = c("New", "ICES")))

res_df_plot %>%
  ggplot(aes(x = year - 100, y = data, colour = MP, linetype = level)) +
  geom_line() +
  # geom_smooth(aes(x = year - 100, y = data, colour = MP, linetype = level),
  #             se = FALSE, span = 0.4, linewidth = 0.4) + 
  scale_colour_manual("Approach", values = c("black", "red")) +
  scale_linetype_manual("Depletion", values = c("solid", "1111")) +
  #guides(linetype = "none") +
  guides(linetype = guide_legend(override.aes = list(color = "black"))) +
  facet_wrap(~quant, labeller = "label_parsed", scales = "free_y",
             strip.position = "left") +
  labs(x = "Year") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size = 8),
        axis.title.y = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        #legend.position = "inside",
        #legend.position.inside = c(0.8, 0.5),
        legend.key = element_blank(),
        legend.background = element_blank())

ggsave(filename = "output/plots/fcc/fcc_cc_traj.png",
       type = "cairo-png", 
       width = 14, height = 5, units = "cm", dpi = 600)

### only catch and depleted state
res_df_plot %>%
  filter(quant == "Catch/MSY" & level == "Depleted") %>%
  ggplot(aes(x = year - 100, y = data, colour = MP)) +
  #geom_line() +
  geom_smooth(se = FALSE, span = 0.4, show.legend = FALSE) + 
  geom_text(data = data.frame(year = 180,
                              data = c(0.55, 0.08),
                              MP = c("New", "ICES")),
            mapping = aes(label = MP),
            show.legend = FALSE) +
  scale_colour_manual("Approach", values = c("black", "red")) +
  labs(x = "Time", y = "Catch") +
  ylim(c(0, NA)) +
  theme_bw(base_size = 8) +
  theme(axis.text = element_blank())

ggsave(filename = "output/plots/fcc/fcc_cc_traj_simple.png",
       type = "cairo-png", 
       width = 5, height = 3, units = "cm", dpi = 600)
ggsave(filename = "output/plots/fcc/fcc_cc_traj_simple.pdf",
       width = 5, height = 3, units = "cm")


### ------------------------------------------------------------------------ ###
### fcc vs ICES - catch-risk plot ####
### ------------------------------------------------------------------------ ###
### use random fhist with 10k iterations, 100 years
### use B/Bmsy~1 for comparison
### use long-term (last 50 years)
### pollack

res_df <- foreach(stk = c("ang3","pol", "tur", "her"), .combine = bind_rows) %:%
  foreach(multiplier = c(0.75, 1), .combine = bind_rows) %:%
  foreach(MP = c("fcc", "constant_catch"), .combine = bind_rows) %do% {
    #browser()
    
    ### reference points
    Blim <- attr(brps[[stk]], "Blim")
    Bmsy <- brps[[stk]]@refpts["msy", "ssb"]
    MSY <- brps[[stk]]@refpts["msy", "yield"]
    
    ### MP output
    if (identical(MP, "fcc")) {
      file <- paste0("output/fcc/10000_100/default/random/",
                     stk, "/mp_0.1_0_0.1_0.6_0_0.75_0_",
                     "_2_", multiplier, "_1_0.5_1.1_0.7.rds")
    } else if (identical(MP, "constant_catch")) {
      file <- paste0("output/constant_catch/10000_100/baseline/random/",
                     stk, "/mp_0.1_0_0.1_0.6_0_0.75_0_.rds")
    }

    if (!file.exists(file)) return(NULL)
    if (isTRUE(MP == "constant_catch" & multiplier != 1)) return(NULL)
    res <- readRDS(file)
    
    ### collapse correction
    res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
    res_corrected <- lapply(res_corrected, window, start = 151, end = 200)
    
    ### starting condition
    SSBs0 <- ssb(res@om@stock)[, ac(100)]
    SSBs0 <- SSBs0/c(Bmsy)
    SSBs0 <- c(SSBs0)
    
    its_Bmsy <- which(SSBs0 >= 0.9 & SSBs0 <= 1.1)
    its_thirdBmsy <- which(SSBs0 >= 0.2 & SSBs0 <= 0.4)
    
    df <- data.frame(
      SSB_Bmsy = median(res_corrected$ssb[,,,,, its_Bmsy]/c(Bmsy)),
      SSB_thirdBmsy = median(res_corrected$ssb[,,,,, its_thirdBmsy]/c(Bmsy)),
      Catch_Bmsy = median(res_corrected$catch[,,,,, its_Bmsy]/c(MSY)),
      Catch_thirdBmsy = median(res_corrected$catch[,,,,, its_thirdBmsy]/c(MSY)),
      risk_Bmsy = mean(iterMeans(res_corrected$ssb[,,,,, its_Bmsy] < Blim)),
      risk_thirdBmsy = mean(iterMeans(res_corrected$ssb[,,,,, its_thirdBmsy] < Blim))
    )
    df$stock <- stk
    df$MP <- MP
    df$multiplier <- multiplier
    return(df)
    
}

res_df_plot <- res_df %>%
  pivot_longer(1:6) %>%
  separate(name, into = c("quant", "level"), sep = "_") %>%
  pivot_wider(names_from = quant, values_from = value) %>%
  mutate(level = factor(level, levels = c("Bmsy", "thirdBmsy"),
                        labels = c("Sustainable", "Depleted")),
         MP = factor(MP, levels = c("fcc", "constant_catch"),
                     labels = c("New", "ICES")))

res_df_plot %>%
  filter(stock == "pol" & multiplier == 1) %>%
  ggplot(aes(x = risk, y = Catch, colour = MP, 
             group = interaction(stock, MP, multiplier))) +
  geom_point(aes(shape = level)) +
  geom_line(show.legend = FALSE) +
  scale_colour_manual("Approach", values = c("black", "red")) +
  scale_shape("Depletion") +
  labs(x = "Depletion risk", y = "Catch/MSY") +
  coord_cartesian(xlim = c(0, 0.135), ylim = c(-0.01, 0.7), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        legend.position = "top", legend.direction = "vertical")
ggsave(filename = "output/plots/fcc/fcc_cc_catch_risk.png",
       type = "cairo-png", 
       width = 5.5, height = 6, units = "cm", dpi = 600)

### with tuning
res_df_plot %>%
  filter(stock == "pol") %>%
  ggplot(aes(x = risk, y = Catch, colour = MP, 
             group = interaction(stock, MP, multiplier))) +
  geom_point(aes(shape = level)) +
  geom_line(show.legend = FALSE) +
  scale_colour_manual("Approach", values = c("black", "red")) +
  scale_shape("Depletion") +
  labs(x = "Depletion risk", y = "Catch/MSY") +
  coord_cartesian(xlim = c(0, 0.135), ylim = c(-0.01, 0.7), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        legend.position = "top", legend.direction = "vertical")
ggsave(filename = "output/plots/fcc/fcc_cc_catch_risk_tuning.png",
       type = "cairo-png", 
       width = 5.5, height = 6, units = "cm", dpi = 600)

res_df_plot %>%
  #filter(stock == "pol") %>%
  ggplot(aes(x = risk, y = Catch, colour = MP, 
             group = interaction(stock, MP))) +
  geom_point(aes(shape = level)) +
  geom_line(show.legend = FALSE) +
  scale_colour_manual("Approach", values = c("black", "red")) +
  scale_shape("Depletion") +
  labs(x = "Depletion risk", y = "Catch/MSY") +
  coord_cartesian(xlim = c(0, 0.3), ylim = c(-0.01, 0.74), expand = FALSE) +
  theme_bw(base_size = 8) +
  theme(legend.key.height = unit(0.5, "lines"),
        legend.position = "top", legend.direction = "vertical")
ggsave(filename = "output/plots/fcc/fcc_cc_catch_risk_4.png",
       type = "cairo-png", 
       width = 5.5, height = 6, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### fcc vs ICES - sensitivity analysis ####
### ------------------------------------------------------------------------ ###

### ------------------------------------------------------------------------ ###
### sensitivity to simulation assumptions ####
### ------------------------------------------------------------------------ ###
### use pollack as example
### default: 50 yrs, 500 iterations

### some generic parameters
brp <- readRDS("input/brps_new.rds")[["pol"]]
Blim <- brp@Blim
Bmsy <- c(refpts(brp)["msy", "ssb"])
MSY <- c(refpts(brp)["msy", "yield"])

### stats over time
stats_sens_time <- foreach(fhist = c("random", "one-way", "roller-coaster"),
                           .combine = rbind) %do% {
  #browser()
  file <- paste0("mp_0.1_0_0.1_0.6_0_0.75_0__2_0.75_1_0.5_1.1_0.7")
  res <- readRDS(paste0("output/fcc/500_100/default/", fhist, 
                        "/pol/", file, ".rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
  ### template
  tmp <- data.frame(year = 1:100)
  ### Blim risk
  tmp$risk_average <- sapply(1:100, function(x) {
    mean(c(res_corrected$ssb[, ac(seq(from = 101, length.out = x))] < Blim), 
         na.rm = TRUE)
  })
  tmp$risk_annual <- c(apply(res_corrected$ssb < Blim, 2, mean, na.rm = TRUE))
  ### SSB
  tmp$SSB_annual <- sapply(1:100, function(x) {
    median(c(res_corrected$ssb[, x]/Bmsy), na.rm = TRUE)
  })
  tmp$SSB_average <-  sapply(1:100, function(x) {
    median(c(res_corrected$ssb[, ac(seq(from = 101, length.out = x))]/Bmsy), 
           na.rm = TRUE)
  })
  ### Catch
  tmp$Catch_annual <- sapply(1:100, function(x) {
    median(c(res_corrected$catch[, x]/MSY), na.rm = TRUE)
  })
  tmp$Catch_average <-  sapply(1:100, function(x) {
    median(c(res_corrected$catch[, ac(seq(from = 101, length.out = x))]/MSY), 
           na.rm = TRUE)
  })
  tmp <- tmp %>%
    pivot_longer(2:7, names_to = c(".value", "period"), names_sep = "_")
  ### full data.frame
  df_i <- data.frame(
    stock = "pol", 
    interval = 2, multiplier = 0.75, idxL_range = 1, catch_limit = 0.5,
    upper_constraint = 1.1, lower_constraint = 0.7,
    sigmaL = 0.1, sigmaL_rho = 0, sigmaR = 0.6, sigmaR_rho = 0, 
    steepness = 0.75,
    risk_Blim = tmp$risk,
    SSB_rel = tmp$SSB,
    Catch_rel = tmp$Catch,
    stat_metric = tmp$period,
    fhist = fhist,
    n_yrs = tmp$year,
    n_iter = 500,
    sensitivity = "period") %>%
    arrange(stat_metric, n_yrs)
    return(df_i)
}
stats_sens_time <- stats_sens_time %>%
  pivot_longer(c(risk_Blim, SSB_rel, Catch_rel)) %>%
  mutate(name = factor(name, levels = c("Catch_rel", "SSB_rel", "risk_Blim"),
                       labels = c("Catch/MSY", "SSB/B[MSY]", "B[lim]~risk")),
         fhist = factor(fhist, levels = c("one-way", "roller-coaster",
                                          "random")))
saveRDS(stats_sens_time, "output/fcc/pol_sensitivity_time.rds")
stats_sens_time <- readRDS("output/fcc/pol_sensitivity_time.rds")

### stock status
stats_sens_status <- foreach(fhist = c("random"),
  .combine = rbind) %do% {
  #browser()
  res <- readRDS(paste0("output/fcc/10000_100/default/", fhist, 
                       "/pol/mp_0.1_0_0.1_0.6_0_0.75_0__2_0.75_1_0.5_1.1_0.7.rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@om@stock, yrs = 101:200)
  ### starting condition
  SSBs0 <- ssb(res@om@stock)[, ac(100)]
  SSBs0 <- SSBs0/Bmsy
  SSBs0 <- c(SSBs0)
  SSB_breaks <- seq(from = 0, to = max(SSBs0), by = 0.1)
  SSB_groups <- cut(SSBs0, breaks = SSB_breaks)
  SSB_levels <- unique(as.character(SSB_groups))
  ### number of replicates per group
  group_n <- sapply(SSB_levels, function(x) {
    length(which(SSB_groups %in% x))
  })
  # group_n[sort(names(group_n))]
  ### Blim risk per group
  ### SSB is on absolute scale 
  risk_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[, ac(151:200),,,, which(SSB_groups %in% x)]
    mean(tmp < (Blim))
  })
  ### SSB (long-term median) per group
  SSB_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$ssb[, ac(151:200),,,, which(SSB_groups %in% x)]/Bmsy
    median(tmp)
  })
  ### Catch (long-term median) per group
  Catch_group <- sapply(SSB_levels, function(x) {
    tmp <- res_corrected$catch[, ac(151:200),,,, which(SSB_groups %in% x)]/MSY
    median(tmp)
  })
  ### get starting conditions
  SSB_levels <- sapply(SSB_levels, function(x) {
    x <- gsub(x = x, pattern = "\\(|\\]", replacement = "")
    x <- unlist(strsplit(x, split = ","))
    mean(as.numeric(x))
  })
  pos_remove <- which(is.na(SSB_levels))
  df_i <- data.frame(
    stock = "pol", 
    interval = 2, multiplier = 0.75, idxL_range = 1, catch_limit = 0.5,
    upper_constraint = 1.1, lower_constraint = 0.7,
    sigmaL = 0.1, sigmaL_rho = 0, sigmaR = 0.6, sigmaR_rho = 0, 
    steepness = 0.75,
    risk_Blim = unlist(risk_group)[-pos_remove],
    SSB_rel = unlist(SSB_group)[-pos_remove],
    Catch_rel = unlist(Catch_group)[-pos_remove],
    SSB0_rel = unlist(SSB_levels)[-pos_remove],
    n_iter_part = unlist(group_n)[-pos_remove],
    fhist = fhist,
    n_yrs = 50,
    n_iter = 10000,
    sensitivity = "stock_status")
  row.names(df_i) <- NULL
  df_i <- df_i[order(df_i$SSB0_rel), ]
  return(df_i)
}
stats_sens_status <- stats_sens_status %>%
  pivot_longer(c(risk_Blim, SSB_rel, Catch_rel)) %>%
  mutate(name = factor(name, levels = c("Catch_rel", "SSB_rel", "risk_Blim"),
                       labels = c("Catch/MSY", "SSB/B[MSY]", "B[lim]~risk")),
         fhist = factor(fhist, levels = c("one-way", "roller-coaster",
                                          "random")))
saveRDS(stats_sens_status, "output/fcc/pol_sensitivity_status.rds")
stats_sens_status <- readRDS("output/fcc/pol_sensitivity_status.rds")

### sensitivity runs - collate stats
stats_runs <- foreach(fhist = c("one-way", "roller-coaster", "random"), 
        .combine = bind_rows) %do% {#browser()
  path <- paste0("output/fcc/500_100/default/", fhist, "/pol/")
  files <- list.files(path, pattern = "stats_")
  scns <- lapply(paste0(path, files), readRDS)
  scns <- do.call(bind_rows, scns)
  scns <- scns %>% filter(multiplier == 0.75)
  saveRDS(scns, file = paste0(path, "sensitivity_stats.rds"))
  return(scns)
}
saveRDS(stats_runs, "output/fcc/pol_sensitivity.rds")
stats_sens <- readRDS("output/fcc/pol_sensitivity.rds")

### plot sensitivity analysis
### use long term (last 50 years)
stats_sens_plot <- stats_runs %>%
  select(stock, fhist, n_iter, n_yrs, sigmaL, sigmaL_rho, sigmaR, sigmaR_rho,
         steepness, multiplier, interval, upper_constraint, lower_constraint,
         catch_limit, idxL_range, 
         SSB_rel = SSB_rel_lastfhalf,
         Catch_rel = Catch_rel_lastfhalf,
         risk_Blim = risk_Blim_lastfhalf) %>%
  mutate(SSB_rel = unlist(SSB_rel),
         Catch_rel = unlist(Catch_rel),
         risk_Blim = unlist(risk_Blim)) %>%
  pivot_longer(c(risk_Blim, SSB_rel, Catch_rel)) %>%
  mutate(name = factor(name, levels = c("Catch_rel", "SSB_rel", "risk_Blim"),
                       labels = c("Catch/MSY", "SSB/B[MSY]", "B[lim]~risk")),
         fhist = factor(fhist, levels = c("one-way", "roller-coaster",
                                          "random")))
df_blank <- data.frame(name = rep(c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk"),
                                  each = 2),
                       x = c(0, 1, 0, 1, 0, 1),
                       value = c(0, 6, 0, 1, 0, 1),
                       fhist = factor("one-way", 
                                      levels = c("one-way", "roller-coaster",
                                                 "random"))) %>%
  mutate(name = factor(name, levels = c("Catch/MSY", "SSB/B[MSY]", "B[lim]~risk")))
res_def_colours <- c("one-way" = brewer.pal(n = 4, name = "Set1")[1], 
                     "roller-coaster" = brewer.pal(n = 4, name = "Set1")[4], 
                     "random" = brewer.pal(n = 4, name = "Set1")[2])
res_def_linetype <- c("one-way" = "solid", 
                      "roller-coaster" = "1212", 
                      "random" = "3232")
p_sigmaR <- stats_sens_plot %>%
  filter(sigmaL == 0.1 & sigmaL_rho == 0 & sigmaR_rho == 0 & steepness == 0.75) %>%
  ggplot(aes(x = sigmaR, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.6, size = 0.4, colour = "grey") +
  stat_smooth(n = 20, span = 0.6, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Recruitment\nvariability'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(sigma[R])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_text(size = 8),
        #strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 4), "pt"))

p_sigmaR_rho <- stats_sens_plot %>%
  filter(sigmaL == 0.1 & sigmaL_rho == 0 & sigmaR == 0.6 & steepness == 0.75) %>%
  ggplot(aes(x = sigmaR_rho, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.0, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Recruitment\nauto-correlation'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(italic(rho)[R])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        #strip.text.x = element_text(margin = margin(8, 0, 1.2, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_steepness <- stats_sens_plot %>%
  filter(sigmaL == 0.1 & sigmaL_rho == 0 & sigmaR == 0.6 & sigmaR_rho == 0) %>%
  ggplot(aes(x = steepness, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.75, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.4, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Recruitment\nsteepness'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(italic(h))) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        #strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sigmaL <- stats_sens_plot %>%
  filter(steepness == 0.75 & sigmaL_rho == 0 & sigmaR == 0.6 & sigmaR_rho == 0) %>%
  ggplot(aes(x = sigmaL, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.2, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Observation\nuncertainty'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(italic(sigma)[L])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        #strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sigmaL_rho <- stats_sens_plot %>%
  filter(steepness == 0.75 & sigmaL == 0.1 & sigmaR == 0.6 & sigmaR_rho == 0) %>%
  ggplot(aes(x = sigmaL_rho, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 0.0, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.2, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Observation\nauto-correlation'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = expression(italic(rho)[L])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        #strip.text.x = element_text(margin = margin(8, 0, 1, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sens_status <- stats_sens_status %>%
  filter(SSB0_rel <= 2) %>%
  ggplot(aes(x = SSB0_rel, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  stat_smooth(n = 50, span = 0.4, se = FALSE, geom = "line", size = 0.4,
              show.legend = FALSE) + 
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Initial\nstock status'", scales = "free",
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(limits = c(-0.05, 2.05)) +
  labs(x = expression(SSB[y == 0]/B[MSY])) +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        ### manual margins because no letter goes below base
        #strip.text.x = element_text(margin = unit(c(3.8, 0, 3.8, 0), "pt")),
        #strip.text.x = element_text(margin = margin(8, 0, 1.2, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 2, 4, 0), "pt"))
p_sens_period <- stats_sens_time %>%
  filter(stat_metric == "annual") %>%
  ggplot(aes(x = n_yrs, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  geom_vline(xintercept = 100, size = 0.4, colour = "grey") +
  stat_smooth(n = 50, span = 0.1, se = FALSE, geom = "line", size = 0.4) + 
  geom_point(size = 0.15, stroke = 0, shape = 21) +
  geom_blank(data = df_blank, aes(x = x, y = value)) +
  facet_grid(name ~ "'Implementation\nperiod'", scales = "free", 
             labeller = "label_parsed",
             switch = "y") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.text.y = element_blank(),
        #strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
        strip.background.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.switch.pad.grid = unit(0, "pt"),
        plot.margin = unit(c(2, 4, 4, 0), "pt"),
        legend.position = "inside",
        legend.position.inside = c(0.55, 0.26),
        legend.background = element_blank(),
        legend.key.height = unit(0.5, "lines"),
        legend.key.width = unit(0.6, "lines"),
        legend.title = element_blank(),
        legend.key = element_blank())

p <- p_sigmaR + p_sigmaR_rho + p_steepness +
  p_sigmaL + p_sigmaL_rho +
  p_sens_status + p_sens_period +
  plot_layout(nrow = 1)
p

ggsave(filename = "output/plots/fcc/pol_sensitivity_stats.png", 
       type = "cairo", plot = p,
       width = 17, height = 8, units = "cm", dpi = 600)
ggsave(filename = "output/plots/fcc/pol_sensitivity_stats.pdf", plot = p,
       width = 17, height = 8, units = "cm")

