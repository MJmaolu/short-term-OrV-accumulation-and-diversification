# SHORT-TERM ACCUMULATION AND DIVERSITY DYNAMICS
# MS FIGURES
# MJ OLMO-UCEDA

# 2026/04/28
################################################################################
# FIG 2A
################################################################################
df_complete.normH.summary %>%
  filter(hpi > 0) %>%
  ggplot(aes(x = hpi,
             group = rna,
             color = rna)) +
  # geom_hline(yintercept = 0,
  #            linewidth = 0.2,
  #            color = "black") +
  # geom_point(data = data.frame(hpi = c(0, 2, 4, 6, 8, 12, 18, 22, 26, 32, 38, 44),
  #                              y = rep(0, 12)),
  #            aes(x = hpi, y = y),
  #            alpha = .5,
  #            inherit.aes = FALSE) +
  # SD OF MEANS OF THE MEAN DIVERSITY PER TIMEPOINT
  stat_summary(
    aes(y = mean_Sn.plus, fill = rna),
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      data.frame(y = m, ymin = m - s, ymax = m + s)
    },
    geom = "ribbon",
    alpha = .4,
    color = NA
  ) +
  stat_summary(
    aes(y = mean_Sn.minus,
        fill = rna),
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      data.frame(y = m, ymin = m - s, ymax = m + s)
    },
    geom = "ribbon",
    alpha = .1,
    color = NA
  ) +
  
  # SAMPLE MEAN +- SEM
  geom_point(aes(y = mean_Sn.plus),
             size = 2,
             alpha = .8
             ) +
  geom_errorbar(aes(ymin = mean_Sn.plus - (sd_Sn.plus / sqrt(n)),
                    ymax = mean_Sn.plus + (sd_Sn.plus / sqrt(n))),
                width = 0.05) +
  
  geom_point(aes(y = mean_Sn.minus),
             shape = 21,
             size = 2,
             alpha = .8) +
  geom_errorbar(aes(ymin = mean_Sn.minus - (sd_Sn.minus / sqrt(n)),
                    ymax = mean_Sn.minus + (sd_Sn.minus / sqrt(n))),
                width = 0.1) +
  stat_summary(fun = "mean",
               geom = "line",
               linetype = "solid",
               aes(y = mean_Sn.plus)) +
 
  stat_summary(fun = "mean",
               geom = "line",
               linetype = "dashed",
               aes(y = mean_Sn.minus)) +
  
  
  scale_color_manual(values = color.rna) +
  scale_fill_manual(values = color.rna) +
  
  ylab(paste("Mean Shannon entropy")) +
  #ylab(paste("Mean Shannon entropy", "SEM", sep = " \u00B1 ")) +
  
  xlab("Time(hpi)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))


################################################################################
# FIG 2B
################################################################################
samples.VPMV.shannon %>% 
  mutate(point_size = ifelse(hpi == 18, "18", "rest")) %>% 
  
  ggplot(aes(x = VPMB,
             y = mean_Sn,
             color = rna,
             linetype = strand)) +
  
  geom_point(aes(shape = strand,
                 size = point_size),
             alpha = 0.8) +
  
  # IC (ribbon)
  geom_ribbon(data = pred_sat_data,
              aes(x = VPMB,
                  ymin = lower,
                  ymax = upper,
                  fill = rna,
                  group = strand
              ),
              alpha = 0.15,
              inherit.aes = FALSE) +
  
  # Curva
  geom_line(data = pred_sat_data,
            aes(x = VPMB,
                y = mean_Sn,
                color = rna,
                linetype = strand),
            linewidth = 0.8) +
  
  # geom_hline(data = plateau_df,
  #          aes(yintercept = a, 
  #              color = rna),
  #          linetype = "dotted",
  #          alpha = 0.6,
  #          inherit.aes = FALSE) +
  
  facet_grid(~rna, scales = "free") +
  
  scale_color_manual(values = c("RNA1" = "#4E3BA7",
                                "RNA2" = "#90B083")) +
  
  scale_fill_manual(values = c("RNA1" = "#4E3BA7",
                               "RNA2" = "#90B083")) +
  
  scale_linetype_manual(values = c("plus" = "solid",
                                   "minus" = "dashed")) +
  
  scale_shape_manual(values = c("plus" = 19,
                                "minus" = 21)) +
  
  scale_size_manual(values = c("18" = 3,
                               "rest" = 1.5)) +
  
  ylab(expression(mu * S[n])) +
  xlab("VPMB") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))
