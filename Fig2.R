# SHORT-TERM ACCUMULATION AND DIVERSITY DYNAMICS
# MS FIGURES
# MJ OLMO-UCEDA

# 2026/11/05
################################################################################
# FIG 2A
################################################################################

df_complete.normH.summary.q20 <- df_complete.normH_wider_q20 %>%
  group_by(
    sample,
    rna,
    hpi
  ) %>%
  summarise(
    n = n(),
    
    mean_APD.plus =
      mean(`D_j/site_pairs_plus`,
           na.rm = TRUE),
    
    sd_APD.plus =
      sd(`D_j/site_pairs_plus`,
         na.rm = TRUE),
    
    mean_APD.minus =
      mean(`D_j/site_pairs_minus`,
           na.rm = TRUE),
    
    sd_APD.minus =
      sd(`D_j/site_pairs_minus`,
         na.rm = TRUE),
    
    mean_Sn.plus =
      mean(norm.H_plus,
           na.rm = TRUE),
    
    sd_Sn.plus =
      sd(norm.H_plus,
         na.rm = TRUE),
    
    mean_Sn.minus =
      mean(norm.H_minus,
           na.rm = TRUE),
    
    sd_Sn.minus =
      sd(norm.H_minus,
         na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    replica =
      as.numeric(
        strsplit(sample, "_")[[1]][3]
      )
  )

################################################################################
# FIG 2A but with APD
################################################################################
# With APD
df_complete.normH.summary.q20 %>%
  filter(hpi > 2) %>%
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
    aes(y = mean_APD.plus, 
        fill = rna),
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      data.frame(y = m, ymin = m - s, ymax = m + s)
    },
    geom = "ribbon",
    alpha = .3,
    color = NA
  ) +
  stat_summary(
    aes(y = mean_APD.minus,
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
  geom_point(aes(y = mean_APD.plus),
             size = 2,
             alpha = .8
  ) +
  geom_errorbar(aes(ymin = mean_APD.plus - (sd_APD.plus / sqrt(n)),
                    ymax = mean_APD.plus + (sd_APD.plus / sqrt(n))),
                width = 0.05) +
  
  geom_point(aes(y = mean_APD.minus),
             shape = 21,
             size = 2,
             alpha = .8) +
  geom_errorbar(aes(ymin = mean_APD.minus - (sd_APD.minus / sqrt(n)),
                    ymax = mean_APD.minus + (sd_APD.minus / sqrt(n))),
                width = 0.1) +
  stat_summary(fun = "mean",
               geom = "line",
               linetype = "solid",
               aes(y = mean_APD.plus)) +
  
  stat_summary(fun = "mean",
               geom = "line",
               linetype = "dashed",
               aes(y = mean_APD.minus)) +
  
  
  scale_color_manual(values = color.rna) +
  scale_fill_manual(values = color.rna) +
  
  ylab(paste("Mean Average Pairwise Diversity")) +
  #ylab(paste("Mean Shannon entropy", "SEM", sep = " \u00B1 ")) +
  
  xlab("Time(hpi)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

##===============================================================================
## FIG 2A
#===============================================================================
# With Shannon
df_complete.normH.summary.q20 %>%
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
    alpha = .3,
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
                width = 0.05) +
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
                                    fill = "transparent")) +

#===============================================================================
## FIG 2B
#===============================================================================
samples.VPMV.shannon.q20 %>%
  mutate(hpi = as.numeric(hpi)) %>%
  
  ggplot(
    aes(
      x = VPMB,
      y = mean_Sn,
      color = rna,
      linetype = strand
    )
  ) +
  
  geom_point(
    aes(
      shape = strand),
    alpha = 0.8,
    size = 2,
  ) +
  
  geom_ribbon(
    data = pred_sat_data.q20,
    
    aes(
      x = VPMB,
      ymin = lower,
      ymax = upper,
      fill = rna,
      group = interaction(rna, strand)
    ),
    
    alpha = 0.15,
    colour = NA,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = pred_sat_data.q20,
    aes(
      x = VPMB,
      y = mean_Sn,
      color = rna,
      linetype = strand
    ),
    
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  
  facet_grid(
    ~rna,
    scales = "free"
  ) +
  
  scale_color_manual(
    values = c(
      RNA1 = "#4E3BA7",
      RNA2 = "#90B083"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      RNA1 = "#4E3BA7",
      RNA2 = "#90B083"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      plus = "solid",
      minus = "dashed"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      plus = 19,
      minus = 21
    )
  ) +
  
  # ylab(expression(mu*S[n])) +
  ylim(c(0,0.0065)) +
  ylab("Mean Shannon entropy") +
  xlab("VPMB") +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  )

bind_rows(
  samples.VPMV.shannon %>%
    mutate(dataset = "Original"),
  
  samples.VPMV.shannon.q20 %>%
    mutate(dataset = "Q20")
) %>%
  
  ggplot(
    aes(
      VPMB,
      mean_Sn,
      colour = strand
    )
  ) +
  
  geom_point(alpha = .4) +
  
  geom_smooth(
    method = "loess",
    se = FALSE
  ) +
  
  scale_color_manual(values = color.strand) +
  
  facet_grid(
    dataset ~ rna
  )


## MODELS

df_complete.normH.summary.long.q20 <- df_complete.normH.summary.q20 %>%
  pivot_longer(
    cols = c(
      mean_APD.plus,
      sd_APD.plus,
      mean_APD.minus,
      sd_APD.minus,
      mean_Sn.plus,
      sd_Sn.plus,
      mean_Sn.minus,
      sd_Sn.minus
    ),
    names_to = c(
      "measure_name",
      "strand"
    ),
    names_pattern =
      "(.*).(plus|minus)",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = measure_name,
    values_from = value
  )


# q20
lmer(
  mean_Sn ~
    rna * hpi * strand + (1 | sample),
  data = df_complete.normH.summary.long.q20
) %>% summary()

# original
lmer(
  mean_Sn ~ rna * hpi * strand + (1 | sample),
  data = df_complete.normH.summary.long
) %>% summary()


car::Anova(
  lmer(
    mean_Sn ~ rna * hpi * strand +
      (1|sample),
    data = df_complete.normH.summary.long.q20
  ),
  type = 3
)

emmeans::emmeans(
  lmer(
    mean_Sn ~
      rna * hpi * strand + (1 | sample),
    data = df_complete.normH.summary.long.q20
  ),
  ~ strand | rna
)

emmeans::emmeans(
  lmer(
    mean_Sn ~
      rna * hpi * strand + (1 | sample),
    data = df_complete.normH.summary.long.q20
  ),
  pairwise ~ strand
)


df_complete.normH.summary.long.q20 %>%
  group_by(
    rna,
    hpi,
    strand
  ) %>%
  
  summarise(
    m = mean(mean_Sn),
    .groups="drop"
  ) %>%
  
  pivot_wider(
    names_from = strand,
    values_from = m
  ) %>%
  
  mutate(
    delta = minus - plus
  ) %>% print(n = 24)
