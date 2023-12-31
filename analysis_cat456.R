### ------------------------------------------------------------------------ ###
### analysis of constant harvest rate rule ####
### ------------------------------------------------------------------------ ###

library(mse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(doParallel)
library(scales)
library(patchwork)
library(RColorBrewer)
source("funs.R")

### stock list
stocks <- read.csv("input/stocks.csv", stringsAsFactors = FALSE)
### brps
brps <- readRDS("input/brps.rds")


### ------------------------------------------------------------------------ ###
### risk vs depletion ####
### ------------------------------------------------------------------------ ###
### example: pollack, random fhist, 10,000 iterations, 100 years



### stock status
input <- readRDS("input/10000_100/OM/random/pol/stk.rds")
Blim <- attr(brps$pol, "Blim")
Bmsy <- brps$pol@refpts["msy", "ssb"]
MSY <- brps$pol@refpts["msy", "yield"]

res <- readRDS(paste0("output/const_catch/10000_100/baseline/random/pol/", 
                      "mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_0.6_0_0.75.rds"))
### collapse correction
res_corrected <- collapse_correction(stk = res@stock, yrs = 101:200)
### starting condition
SSBs0 <- ssb(res@stock)[, ac(100)]
SSBs0 <- SSBs0/c(brps$pol@refpts["msy", "ssb"])
SSBs0 <- c(SSBs0)
SSB_breaks <- seq(from = 0, to = max(SSBs0), by = 0.1)
SSB_groups <- cut(SSBs0, breaks = SSB_breaks)
SSB_levels <- unique(as.character(SSB_groups))
### number of replicates per group
group_n <- sapply(SSB_levels, function(x) {
 length(which(SSB_groups %in% x))
})
group_n[sort(names(group_n))]
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
 stock = "pol",
 risk_Blim = unlist(risk_group)[-pos_remove],
 SSB_rel = unlist(SSB_group)[-pos_remove],
 Catch_rel = unlist(Catch_group)[-pos_remove],
 SSB0_rel = unlist(SSB_levels)[-pos_remove],
 n_iter_part = unlist(group_n)[-pos_remove],
 fhist = "random",
 n_yrs = 100,
 n_iter = 10000,
 steepness = 0.75,
 sensitivity = "stock_status")
row.names(stats_risk_depletion) <- NULL
stats_risk_depletion <- stats_risk_depletion[order(stats_risk_depletion$SSB0_rel), ]

### ------------------------------------------------------------------------ ###
### risk vs projection period ####
### ------------------------------------------------------------------------ ###
### example: pollack, random & one-way fhist, 10,000 iterations, 100 years

stats_sens_time <- foreach(fhist = c("random", "one-way"),
  .combine = rbind) %do% {
  #browser()
  file <- "mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_0.6_0_0.75"
  res <- readRDS(paste0("output/const_catch/10000_100/baseline/", fhist, 
                       "/pol/", file, ".rds"))
  
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@stock, yrs = 101:200)
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
   risk_Blim = tmp$risk,
   SSB_rel = tmp$SSB,
   Catch_rel = tmp$Catch,
   stat_metric = tmp$period,
   fhist = fhist,
   n_yrs = tmp$year,
   n_iter = 10000,
   steepness = 0.75,
   sensitivity = "period") %>%
   arrange(stat_metric, n_yrs)
  return(df_i)
}

### ------------------------------------------------------------------------ ###
### plot sensitivity ####
### ------------------------------------------------------------------------ ###

stats_sens <- bind_rows(stats_risk_depletion, stats_sens_time)
stats_sens_plot <- stats_sens %>%
  pivot_longer(c(risk_Blim, SSB_rel, Catch_rel)) %>%
  mutate(name = factor(name, levels = c("SSB_rel", "Catch_rel", "risk_Blim"),
                       labels = c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk")),
         fhist = factor(fhist, levels = c("one-way", "random")))
df_blank <- data.frame(name = rep(c("SSB/B[MSY]", "Catch/MSY", "B[lim]~risk"),
                                  each = 2),
                       x = c(0, 1, 0, 1, 0, 1),
                       value = c(0, 1.8, 0, 1.4, 0, 1),
                       fhist = NA)
res_def_colours <- c("one-way" = brewer.pal(n = 4, name = "Set1")[1], 
                     #"roller-coaster" = brewer.pal(n = 4, name = "Set1")[4], 
                     "random" = brewer.pal(n = 4, name = "Set1")[2])
res_def_linetype <- c("one-way" = "solid", 
                      #"roller-coaster" = "1212", 
                      "random" = "3232")
### ggplot2 bug... need to plot all elements individually
# p_sens_status <- stats_sens_plot %>%
#   filter(sensitivity == "stock_status" &
#            SSB0_rel <= 2) %>%
#   ggplot(aes(x = SSB0_rel, y = value, fill = fhist, colour = fhist,
#              linetype = fhist)) +
#   stat_smooth(n = 50, span = 0.4, se = FALSE, geom = "line", linewidth = 0.4,
#               show.legend = FALSE) +
#   geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
#   geom_blank(data = df_blank, aes(x = x, y = value)) +
#   facet_grid(name ~ "'\ \ \ \ \ \ Initial\nstock status'", scales = "free",
#              labeller = "label_parsed",
#              switch = "y") +
#   scale_linetype_manual("fishing history", values = res_def_linetype) +
#   scale_colour_manual("fishing history", values = res_def_colours) +
#   scale_fill_manual("fishing history", values = res_def_colours) +
#   scale_x_continuous(limits = c(-0.05, 2.05)) +
#   labs(x = expression(SSB[y == 0]/B[MSY])) +
#   theme_bw(base_size = 8) +
#   theme(strip.placement = "outside",
#         strip.text.y = element_text(size = 8),
#         strip.text.x = element_text(margin = margin(8, 0, 1.2, 0)),
#         strip.background.y = element_blank(),
#         axis.title.y = element_blank(),
#         strip.switch.pad.grid = unit(0, "pt"),
#         plot.margin = unit(c(2, 2, 4, 4), "pt"))

# p_sens_period <- stats_sens_plot %>%
#   filter(sensitivity == "period" &
#            stat_metric == "average") %>%
#   ggplot(aes(x = n_yrs, y = value, fill = fhist, colour = fhist, 
#              linetype = fhist)) +
#   geom_vline(xintercept = 50, size = 0.4, colour = "grey") +
#   stat_smooth(n = 50, span = 0.1, se = FALSE, geom = "line", size = 0.4) + 
#   geom_point(size = 0.15, stroke = 0, shape = 21) +
#   geom_blank(data = df_blank, aes(x = x, y = value)) +
#   facet_grid(name ~ "'Implementation\n\ \ \ \ \ \ \ period'", scales = "free", 
#              labeller = "label_parsed",
#              switch = "y") +
#   scale_linetype_manual("fishing history", values = res_def_linetype) +
#   scale_colour_manual("fishing history", values = res_def_colours) +
#   scale_fill_manual("fishing history", values = res_def_colours) +
#   labs(x = "years") +
#   theme_bw(base_size = 8) +
#   theme(strip.placement = "outside",
#         strip.text.y = element_blank(),
#         strip.text.x = element_text(margin = margin(8, 0, 0, 0)),
#         strip.background.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(), 
#         strip.switch.pad.grid = unit(0, "pt"),
#         plot.margin = unit(c(2, 4, 4, 0), "pt"),
#         legend.position = c(0.55, 0.26),
#         legend.background = element_blank(),
#         legend.key.height = unit(0.5, "lines"),
#         legend.key.width = unit(0.6, "lines"),
#         legend.title = element_blank(),
#         legend.key = element_blank())

### status - Blim risk
p_sens_status_risk <- stats_sens_plot %>%
  filter(sensitivity == "stock_status" & name == "B[lim]~risk" &
           SSB0_rel <= 2) %>%
  ggplot(aes(x = SSB0_rel, y = value, fill = fhist, colour = fhist,
             linetype = fhist)) +
  stat_smooth(n = 50, span = 0.4, se = FALSE, geom = "line", linewidth = 0.4,
              show.legend = FALSE) +
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  facet_wrap(~ "'Initial stock status'", scales = "free",
              labeller = "label_parsed") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(limits = c(-0.05, 2.05)) +
  labs(x = expression(SSB[y == 0]/B[MSY]),
       y = expression(B[lim]~risk)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 2)) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
### status - catch
p_sens_status_catch <- stats_sens_plot %>%
  filter(sensitivity == "stock_status" & name == "Catch/MSY" &
           SSB0_rel <= 2) %>%
  ggplot(aes(x = SSB0_rel, y = value, fill = fhist, colour = fhist,
             linetype = fhist)) +
  stat_smooth(n = 50, span = 0.4, se = FALSE, geom = "line", linewidth = 0.4,
              show.legend = FALSE) +
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  scale_x_continuous(limits = c(-0.05, 2.05)) +
  labs(x = expression(SSB[y == 0]/B[MSY]),
       y = expression("Catch/MSY")) +
  coord_cartesian(ylim = c(0, 1.4), xlim = c(0, 2)) +
  theme_bw(base_size = 8) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())
### status - ssb
p_sens_status_ssb <- stats_sens_plot %>%
  filter(sensitivity == "stock_status" & name == "SSB/B[MSY]" &
           SSB0_rel <= 2) %>%
  ggplot(aes(x = SSB0_rel, y = value, fill = fhist, colour = fhist,
             linetype = fhist)) +
  stat_smooth(n = 50, span = 0.4, se = FALSE, geom = "line", linewidth = 0.4,
              show.legend = FALSE) +
  geom_point(size = 0.15, stroke = 0, shape = 21, show.legend = FALSE) +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  #scale_x_continuous(limits = c(-0.05, 2.05)) +
  labs(x = expression(SSB[y == 0]/B[MSY]),
       y = expression(SSB/B[MSY])) +
  coord_cartesian(ylim = c(0, 4), xlim = c(0, 2)) +
  theme_bw(base_size = 8) +
  theme()

### duration - risk
p_sens_period_risk <- stats_sens_plot %>%
  filter(sensitivity == "period" & name == "B[lim]~risk" &
           stat_metric == "average") %>%
  ggplot(aes(x = n_yrs, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  stat_smooth(n = 50, span = 0.1, se = FALSE, geom = "line", size = 0.4) + 
  geom_point(size = 0.15, stroke = 0, shape = 21) +
  facet_wrap(~ "'Implementation period'", scales = "free", 
             labeller = "label_parsed") +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +
  theme_bw(base_size = 8) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = c(0.7, 0.7),
        legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank())
### duration - catch
p_sens_period_catch <- stats_sens_plot %>%
  filter(sensitivity == "period" & name == "Catch/MSY" &
           stat_metric == "average") %>%
  ggplot(aes(x = n_yrs, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  stat_smooth(n = 50, span = 0.1, se = FALSE, geom = "line", size = 0.4) + 
  geom_point(size = 0.15, stroke = 0, shape = 21) +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  coord_cartesian(ylim = c(0, 1.4), xlim = c(0, 100)) +
  theme_bw(base_size = 8) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank())
### duration - ssb
p_sens_period_ssb <- stats_sens_plot %>%
  filter(sensitivity == "period" & name == "SSB/B[MSY]" &
           stat_metric == "average") %>%
  ggplot(aes(x = n_yrs, y = value, fill = fhist, colour = fhist, 
             linetype = fhist)) +
  stat_smooth(n = 50, span = 0.1, se = FALSE, geom = "line", size = 0.4) + 
  geom_point(size = 0.15, stroke = 0, shape = 21) +
  scale_linetype_manual("fishing history", values = res_def_linetype) +
  scale_colour_manual("fishing history", values = res_def_colours) +
  scale_fill_manual("fishing history", values = res_def_colours) +
  coord_cartesian(ylim = c(0, 4), xlim = c(0, 100)) +
  theme_bw(base_size = 8) +
  labs(x = "Years") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.5, "lines"),
        legend.background = element_blank())

p <- (p_sens_status_risk + p_sens_period_risk)/
  (p_sens_status_catch + p_sens_period_catch)/
  (p_sens_status_ssb + p_sens_period_ssb)
p

ggsave(filename = "output/plots/constant_catch/CC_sensitivity.png",
       type = "cairo", plot = p,
       width = 17, height = 8, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### plot time series ####
### ------------------------------------------------------------------------ ###
brp <- brps$pol
### quantiles
df_time <- foreach(fhist = c("random", "one-way"),
                             .combine = rbind) %do% {
  #browser()
  file <- "mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_0.6_0_0.75"
  res <- readRDS(paste0("output/const_catch/10000_100/baseline/", fhist, 
                        "/pol/", file, ".rds"))
  
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@stock, yrs = 101:200)
  res_corrected$catch[] <- res_corrected$catch/c(refpts(brp)["msy", "yield"])
  res_corrected$ssb[] <- res_corrected$ssb/c(refpts(brp)["msy", "ssb"])
  res_corrected$fbar[] <- res_corrected$fbar/c(refpts(brp)["msy", "harvest"])
  
  
  qnts <- lapply(res_corrected, function(x) {
    quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95))
  })
  qnts <- as(qnts, "FLQuants")
  df <- as.data.frame(qnts)
  df$fhist <- fhist
  return(df)
}
### some iterations
df_time_iters <- foreach(fhist = c("random", "one-way"),
                   .combine = rbind) %do% {
  #browser()
  file <- "mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_0.6_0_0.75"
  res <- readRDS(paste0("output/const_catch/10000_100/baseline/", fhist, 
                        "/pol/", file, ".rds"))
  
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@stock, yrs = 101:200)
  res_corrected$catch[] <- res_corrected$catch/c(refpts(brp)["msy", "yield"])
  res_corrected$ssb[] <- res_corrected$ssb/c(refpts(brp)["msy", "ssb"])
  res_corrected$fbar[] <- res_corrected$fbar/c(refpts(brp)["msy", "harvest"])
  qnts <- lapply(res_corrected, function(x) {
   iter(x, 1:5)
  })
  qnts <- as(qnts, "FLQuants")
  df <- as.data.frame(qnts)
  df$fhist <- fhist
  return(df)

}

df_time %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "fbar", "catch"),
                        labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY")),
         year = year - 100) %>%
  ggplot() +
  geom_ribbon(aes(x = year, ymin = `5%`, ymax = `95%`), alpha = 0.25) + 
  geom_ribbon(aes(x = year, ymin = `25%`, ymax = `75%`), alpha = 0.5) +
  geom_line(aes(x = year, y = `50%`), show.legend = FALSE) +
  geom_line(data = df_time_iters %>%
              mutate(qname = factor(qname,
                                    levels = c("ssb", "fbar", "catch"),
                                    labels = c("SSB/B[MSY]", "F/F[MSY]",
                                               "Catch/MSY")),
                     year = year - 100),
            aes(x = year, y = data, colour = iter),
            show.legend = FALSE, size = 0.2, alpha = 0.5) +
  scale_colour_brewer(palette = "Set1") +
  facet_grid(qname ~ fhist, scales = "free", labeller = "label_parsed", 
             switch = "y") +
  labs(x = "Year") +
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank())

ggsave(filename = "output/plots/constant_catch/CC_pol_timeseries.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)


### ------------------------------------------------------------------------ ###
### risk & projection time - all stocks ####
### ------------------------------------------------------------------------ ###
### 500 iters, 100 years, all stocks, one-way & random

stats_risk <- foreach(stock = stocks$stock, k = stocks$k, 
                      .combine = bind_rows) %:%
  foreach(fhist = c("one-way", "random"), .combine = bind_rows)  %do% {
    #browser()
    res <- readRDS(paste0("output/const_catch/500_100/baseline/", fhist,"/",
                          stock, "/mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_",
                          "0.6_0_0.75.rds"))
    ### collapse correction
    res_corrected <- collapse_correction(stk = res@stock, yrs = 101:200)
    
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


res_def_colours <- c("one-way" = brewer.pal(n = 4, name = "Set1")[1], 
                     "random" = brewer.pal(n = 4, name = "Set1")[2])



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
        legend.position = c(0.9, 0.6),
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
  facet_wrap(~ metric, labeller = "label_parsed", switch = "y", 
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
        legend.position = c(0.9, 0.55),
        axis.title.y = element_blank(),
        plot.margin = unit(x = c(1, 3, 0, 3), units = "pt"))

ggsave(filename = "output/plots/constant_catch/CC_all_annual_timeseries.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### conditional PA buffer - pollack example ####
### ------------------------------------------------------------------------ ###
### apply buffer when Lmean <= LF=M
### every 3 years, 20%

stats_cond <- foreach(rule = c("default", "conditional"), 
                      .combine = bind_rows) %do% {
    #browser()
    if (identical(rule, "default"))
      res <- readRDS(paste0("output/const_catch/500_100/baseline/one-way/pol/",
                            "/mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_",
                            "0.6_0_0.75.rds"))
    if (identical(rule, "conditional"))
      res <- readRDS(paste0("output/const_catch/500_50/CC_f/one-way/pol/",
                            "mp_CC_f_3_3_0.8_Inf_0_0.2_0.2_0_0_0.6_0_0.75.rds"))
    ### collapse correction
    res_corrected <- collapse_correction(stk = res@stock, yrs = 101:150)
    
    ### get fishing history
    stk <- readRDS(paste0("input/500_50/OM/one-way/pol/stk.rds"))
    ### combine with results
    res_corrected <- lapply(res_corrected, window, start = 0)
    res_corrected$catch[, ac(0:100)] <- catch(stk)[, ac(0:100)]
    res_corrected$ssb[, ac(0:100)] <- ssb(stk)[, ac(0:100)]
    res_corrected$fbar[, ac(0:100)] <- fbar(stk)[, ac(0:100)]
    
    ### relative metrics
    res_corrected$catch <- res_corrected$catch/c(refpts(brps$pol)["msy", "yield"])
    res_corrected$ssb <- res_corrected$ssb/c(refpts(brps$pol)["msy", "ssb"])
    res_corrected$fbar <- res_corrected$fbar/c(refpts(brps$pol)["msy", "harvest"])
    
    ### quantiles
    res_corrected <- lapply(res_corrected, quantile, 
                            c(0.05, 0.25, 0.5, 0.75, 0.95))
    res_corrected <- as(res_corrected, "FLQuants")
    
    df <- as.data.frame(res_corrected)
    df$rule <- rule
    
    return(df)
}

stats_cond <- stats_cond %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "fbar", "catch"),
                        labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY")),
         rule = factor(rule, levels = c("default", "conditional"))) %>%
  select(-age, - unit, -season, -area)

stats_cond %>%
  mutate(year = year - 100) %>%
  ggplot(aes(x = year, fill = rule)) +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.125) + 
  geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.25) +
  geom_line(aes(y = `50%`, colour = rule), show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  scale_colour_manual("PA buffer", 
    values = c("default" = brewer.pal(3, name = "Dark2")[1],
               "conditional" = brewer.pal(3, name = "Dark2")[2])) +
  scale_fill_manual("PA buffer", 
                    values = c("default" = brewer.pal(3, name = "Dark2")[1],
                               "conditional" = brewer.pal(3, name = "Dark2")[2])) +
  facet_wrap(~ qname, ncol = 1,
             scales = "free", labeller = "label_parsed", 
             switch = "y") +
  labs(x = "Year") +
  coord_cartesian(xlim = c(0, 50)) + 
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank())
ggsave(filename = "output/plots/constant_catch/CC_CC_f.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)

### ------------------------------------------------------------------------ ###
### conditional PA buffer - pollack example - optimised ####
### ------------------------------------------------------------------------ ###
### apply buffer when Lmean <= LF=M
### every 3 years, 20%
### + optimised

### get optimised parameters
ga_res <- readRDS(paste0("output/CC_f/500_50/CC_f_opt/one-way/pol/",
                         "Lref_mult-pa_size-interval--obj_ICES_MSYPA_res.rds"))
ga_res@solution

stats_cond <- foreach(rule = c("default", "conditional", "optimised"), 
                      .combine = bind_rows) %do% {
  #browser()
  if (identical(rule, "default"))
    res <- readRDS(paste0("output/const_catch/500_100/baseline/one-way/pol/",
                          "/mp_const_catch_3_3_0.8_Inf_0_0.2_0.2_0_0_",
                          "0.6_0_0.75.rds"))
  if (identical(rule, "conditional"))
    res <- readRDS(paste0("output/const_catch/500_50/CC_f/one-way/pol/",
                          "mp_CC_f_3_3_0.8_Inf_0_0.2_0.2_0_0_0.6_0_0.75.rds"))
  if (identical(rule, "optimised"))
    res <- readRDS(paste0("output/CC_f/500_50/CC_f_opt/one-way/pol/",
                          "mp_1.9_0.49_5_1_Inf_0.rds"))
  ### collapse correction
  res_corrected <- collapse_correction(stk = res@stock, yrs = 101:150)
  
  ### get fishing history
  stk <- readRDS(paste0("input/500_50/OM/one-way/pol/stk.rds"))
  ### combine with results
  res_corrected <- lapply(res_corrected, window, start = 0)
  res_corrected$catch[, ac(0:100)] <- catch(stk)[, ac(0:100)]
  res_corrected$ssb[, ac(0:100)] <- ssb(stk)[, ac(0:100)]
  res_corrected$fbar[, ac(0:100)] <- fbar(stk)[, ac(0:100)]
  
  ### relative metrics
  res_corrected$catch <- res_corrected$catch/c(refpts(brps$pol)["msy", "yield"])
  res_corrected$ssb <- res_corrected$ssb/c(refpts(brps$pol)["msy", "ssb"])
  res_corrected$fbar <- res_corrected$fbar/c(refpts(brps$pol)["msy", "harvest"])
  
  ### quantiles
  res_corrected <- lapply(res_corrected, quantile, 
                          c(0.05, 0.25, 0.5, 0.75, 0.95))
  res_corrected <- as(res_corrected, "FLQuants")
  
  df <- as.data.frame(res_corrected)
  df$rule <- rule
  
  return(df)
}

stats_cond <- stats_cond %>%
  pivot_wider(names_from = iter, values_from = data) %>%
  mutate(qname = factor(qname, 
                        levels = c("ssb", "fbar", "catch"),
                        labels = c("SSB/B[MSY]", "F/F[MSY]", "Catch/MSY")),
         rule = factor(rule, 
                       levels = c("default", "conditional", "optimised"))) %>%
  select(-age, - unit, -season, -area)

stats_cond %>%
  mutate(year = year - 100) %>%
  ggplot(aes(x = year, fill = rule)) +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.125) + 
  geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.25) +
  geom_line(aes(y = `50%`, colour = rule), show.legend = FALSE) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  scale_colour_manual("PA buffer", 
    values = c("default" = brewer.pal(3, name = "Dark2")[1],
               "conditional" = brewer.pal(3, name = "Dark2")[2],
               "optimised" = brewer.pal(3, name = "Dark2")[3])) +
  scale_fill_manual("PA buffer", 
    values = c("default" = brewer.pal(3, name = "Dark2")[1],
               "conditional" = brewer.pal(3, name = "Dark2")[2],
               "optimised" = brewer.pal(3, name = "Dark2")[3])) +
  facet_wrap(~ qname, ncol = 1,
             scales = "free", labeller = "label_parsed", 
             switch = "y") +
  labs(x = "Year") +
  coord_cartesian(xlim = c(0, 50)) + 
  theme_bw(base_size = 8) +
  theme(strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y = element_text(size = 8),
        axis.title.y = element_blank())
ggsave(filename = "output/plots/constant_catch/CC_CC_f_opt.png",
       type = "cairo", 
       width = 17, height = 8, units = "cm", dpi = 600)
