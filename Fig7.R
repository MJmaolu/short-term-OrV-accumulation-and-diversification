################################################################################
# FIG 7 A-B
################################################################################
selec.A <- ggplot(
  selection_sample_summary %>%
    mutate(
      protein = factor(
        protein,
        levels = c("RdRP", "CP", "delta")
      )
    ),
  aes(protein, mean_dNdS)
) +
  geom_hline(yintercept = 0) +
  
  # geom_violin(
  #   aes(fill = factor(protein)),
  #   outlier.shape = NA,
  #   width = .4,
  #   alpha = .2
  # ) +
  # geom_quasirandom(data = selection_w4s1 %>%
  #                       mutate(
  #                         dNdS = dN - dS,
  #                         protein = factor(protein,
  #                                          levels = c("RdRP","CP","delta"))) %>%
  #                    filter(dNdS != 0),
  #             aes(protein,
  #                 dNdS,
  #                 color = protein),
  #             alpha = .2,
  #             ) +
  geom_quasirandom(
    aes(size = (hpi),
        fill = protein),
    shape = 21,
    #size = 2,
    width = .15,
    alpha = .9
  ) +
  
  labs(
    x = NULL,
    y = expression("Mean " * (d[N] - d[S])),
    fill = "hpi"
  ) +
  
  stat_summary(aes(color = protein),
               size = 1,
               fun = "median",
               fill = "white",
               shape = 21) +
  
  scale_color_manual(values = color.proteins) +
  scale_fill_manual(values = color.proteins) +
  
  ylim(c(-0.00765, 0.0025)) +
  scale_size_continuous(range = c(1, 3)) +
  
  theme_classic() +
  theme(legend.position = 'none',
        panel.border = element_rect(fill = 'transparent'),
        text = element_text(color = 'black'),
        axis.text = element_text(size = 12,
                                 color = 'black'),
        strip.text = element_text(size = 14))


# PLOTTING WITH THE FIT OF THE LMM
selec.B <- selection_sample_summary %>%
  mutate(protein = factor(protein,
                          levels = c("RdRP", "CP", "delta"))) %>%
  arrange(hpi) %>%
  ggplot(aes(hpi,
             mean_dNdS,
             color = protein)) +
  
  geom_hline(
    yintercept = 0,
    color = "grey50"
  ) +
  
  # Observed sample means
  geom_jitter(
    aes(fill = protein,
        size = hpi),
    color = "black",
    shape = 21,
    #size = 2.5,
    width = 0.4,
    height = 0
  ) +
  
  # Model confidence interval
  geom_ribbon(
    data = newdat,
    aes(
      x = hpi,
      ymin = lwr,
      ymax = upr,
      fill = protein
    ),
    alpha = 0.2,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  # Model prediction
  geom_line(
    data = newdat,
    aes(
      x = hpi,
      y = pred,
      color = protein
    ),
    linewidth = 1,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(
    values = color.proteins
  ) +
  scale_fill_manual(
    values = color.proteins
  ) +
  scale_size_continuous(range = c(1, 3)) +
  
  ylim(c(-0.00765, 0.0025)) +
  
  labs(
    x = "Time (hpi)",
    y = expression("Mean " * (d[N] - d[S])),
    fill = "Protein"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(
      fill = "transparent"
    ),
    text = element_text(
      color = "black"
    ),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    strip.text = element_text(
      size = 14
    )
  )

selec.A + selec.B +
  plot_layout(widths = c(1,2))

################################################################################
# FIG 7C
################################################################################
################################################################################
# Test whether dN-dS differs from zero in each sliding window
#
# One-sided Wilcoxon signed-rank tests are performed independently for:
#   - positive selection (dN-dS > 0)
#   - purifying selection (dN-dS < 0)
#
# P-values are adjusted independently using the Benjamini-Hochberg procedure.
################################################################################

selection_hotspots_w4s1 <-
  
  selection_w4s1 %>%
  mutate(dNdS = dN - dS) %>%
  
  group_by(protein,
           window_start
  ) %>%
  
  summarise(n = sum(!is.na(dNdS)),
            median_dNdS = median(dNdS, na.rm = TRUE),
            mean_dNdS = mean(dNdS, na.rm = TRUE),
    
    p_positive = {
      x <- dNdS[!is.na(dNdS)]
      if(length(unique(x)) < 2){
        NA_real_
      } else {
        wilcox.test(
          x,
          mu = 0,
          alternative = "greater",
          exact = FALSE
        )$p.value
      }
    },
    
    p_negative = {
      x <- dNdS[!is.na(dNdS)]
      if(length(unique(x)) < 2){
        NA_real_
      } else {
        wilcox.test(
          x,
          mu = 0,
          alternative = "less",
          exact = FALSE
        )$p.value
      }
    },
    .groups = "drop"
  ) %>% 
  
  group_by(protein) %>%
  
  # MULTIPLE TEST CORRECTION: FDR BH
  mutate(
    q_positive = p.adjust(
      p_positive,
      method = "BH"
    ),
    q_negative = p.adjust(
      p_negative,
      method = "BH"
    ),
    hotspot =
      case_when(
        q_positive < 0.05 ~ "positive",
        q_negative < 0.05 ~ "negative",
        TRUE ~ "none"
      )
  ) %>%
  ungroup()

selection_hotspots_w4s1$hotspot %>% table()
# negative     none 
# 12     1684 


selection_hotspots_w4s1 %>% 
  filter(p_positive < 0.05)

selection_hotspots_w4s1 %>% 
  filter(p_negative < 0.05)

selection_hotspots_w4s1 %>% 
  filter(q_negative < 0.05)



# DF FOR LABELLING THE HOTSPOTS
selection_w4s1_4hotspot_labels <- selection_w4s1_4plot %>%
  filter(p_positive < 0.05 | p_negative < 0.05) %>%
  mutate(
    window_mid = window_start + (window_end - window_start)/2,
    label = paste0(window_start, "-", window_end),
    
    label_type = case_when(
      q_positive < 0.05 ~ "q_positive",
      q_negative < 0.05 ~ "q_negative",
      p_positive < 0.05 ~ "p_positive",
      p_negative < 0.05 ~ "p_negative",
      
    )
  ) %>%
  distinct(
    protein,
    window_start,
    .keep_all = TRUE
  )


## HIGHLIGHTING SIGNIFICANCE
{selection_w4s1_4plot %>%
    left_join(.,
              samples_info %>% dplyr::select(c("sample", "hpi", "replica")),
              by = "sample") %>%
    rowwise() %>%
    mutate(window_mid = window_start + (window_end - window_start) / 2,
           `pN/pS` = pN / pS,
           `pN - pS` = pN - pS,
           `dN - dS` = dN - dS,
           `NS-S` = nonsynonymous_mutations - synonymous_mutations,
           protein = factor(protein, levels = c("RdRP",
                                                "CP",
                                                "delta"))) %>%
    arrange(hpi) %>%
    #filter(abs(nonsynonymous_mutations - synonymous_mutations) > 0) %>%
    filter(abs(`dN - dS`) > 0,
           protein == "RdRP") %>%
    ggplot(. ,
           aes(x = window_mid,
               y = `dN - dS`,
               #y = `NS-S`,
               color = protein
           ),
    ) +
    ## Palmprint motifs ----------------------------------------------------------
    geom_rect(
      data = palmprint.orv %>%
        mutate(
          protein = factor("RdRP", levels = c("RdRP", "CP", "delta"))),
      aes(
        xmin = aa.start,
        xmax = aa.end,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      alpha = 0.10,
      linewidth = 0.05,
      colour = "black"
    ) +
    
    geom_hline(yintercept = 0,
               linewidth = .2) +
    geom_point(
      aes(
        colour = significance,
        alpha = significance,
        size = significance
      ),
      #size = .5
    ) +
    
    # P VALUES SIGNIFICANT
    geom_label_repel(
      data = selection_w4s1_4hotspot_labels %>%
        filter(protein == "RdRP",
               label_type %in% c("p_positive","p_negative")),
      
      aes(
        x = window_mid,
        y = 0,
        label = label
      ),
      
      fill = NA,
      colour = "black",
      label.size = .25,
      size = 2.8,
      box.padding = 0,
      point.padding = .2,
      segment.size = .2,
      min.segment.length = 0,
      inherit.aes = FALSE
    ) +
    
    # P-adjusted (BH) VALUES SIGNIFICANT
    geom_label_repel(
      data = selection_w4s1_4hotspot_labels %>%
        filter(protein == "RdRP",
               label_type %in% c("q_positive","q_negative")),
      
      aes(
        x = window_mid,
        y = 0,
        label = label,
        fill = label_type
      ),
      
      colour = "black",
      label.size = .25,
      size = 2.8,
      box.padding = .3,
      point.padding = .2,
      segment.size = .2,
      min.segment.length = 0,
      inherit.aes = FALSE
    ) +
    facet_grid( ~ protein,
                scales = "free_x",
                space = "free_x",
    ) +
    #scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    coord_cartesian(clip = "off") +
    scale_colour_manual(
      values = c(
        none = "black",
        p_negative = "#ff1493",
        q_negative = "#ff1493",
        p_positive = "#6b8e23",
        q_positive = "#6b8e23")) +

    scale_alpha_manual(values = c("none" = .5, 
                                  "p_negative" = .5,
                                  "q_negative" = 1,
                                  "p_positive" = .5,
                                  "q_positive" = 1)) +
    scale_size_manual(values = c("none" = .5, 
                                 "p_negative" = .5,
                                 "q_negative" = .8,
                                 "p_positive" = .5,
                                 "q_positive" = .8)) +
    scale_fill_manual(
      values = c(
        q_negative = "#ff1493",
        q_positive = "#6b8e23"
      )
    ) +
    labs(x = "") +
    
    theme_classic() +
    theme(legend.position = 'none',
          panel.border = element_rect(fill = 'transparent'),
          text = element_text(color = 'black'),
          axis.text = element_text(size = 12,
                                   color = 'black'),
          strip.text = element_text(size = 14))} /
  
  ## CP & DELTA
  ###############
  {selection_w4s1_4plot %>%
      left_join(.,
                samples_info %>% dplyr::select(c("sample", "hpi", "replica")),
                by = c("sample", "hpi", "replica")) %>%
      rowwise() %>%
      mutate(window_mid = window_start + (window_end - window_start) / 2,
             `pN/pS` = pN / pS,
             `pN - pS` = pN - pS,
             `dN - dS` = dN - dS,
             `NS-S` = nonsynonymous_mutations - synonymous_mutations,
             protein = factor(protein, levels = c("RdRP",
                                                  "CP",
                                                  "delta"))) %>%
      arrange(hpi) %>%
      #filter(abs(nonsynonymous_mutations - synonymous_mutations) > 0) %>%
      filter(abs(`pN - pS`) > 0,
             protein != "RdRP") %>%
      ggplot(. ,
             aes(x = window_mid,
                 y = `dN - dS`,
                 #y = `NS-S`,
                 color = protein
             ),
      ) +
      geom_hline(yintercept = 0,
                 linewidth = .2) +
      geom_point(
        aes(colour = significance,
            alpha = significance,
            size = significance),
        #size = .5
        ) +
      
      # P VALUES SIGNIFICANT
      geom_label_repel(
        data = selection_w4s1_4hotspot_labels %>%
          filter(protein != "RdRP",
                 label_type %in% c("p_positive","p_negative")),
        
        aes(
          x = window_mid,
          y = 0,
          label = label
        ),
        
        fill = NA,
        colour = "black",
        label.size = .25,
        size = 2.8,
        box.padding = 0,
        point.padding = .2,
        segment.size = .2,
        min.segment.length = 0,
        inherit.aes = FALSE
      ) +
      
      # P-adjusted (BH) VALUES SIGNIFICANT
      geom_label_repel(
        data = selection_w4s1_4hotspot_labels %>%
          filter(protein != "RdRP",
                 label_type %in% c("q_positive","q_negative")),
        
        aes(
          x = window_mid,
          y = 0,
          label = label,
          fill = label_type
        ),
        
        colour = "black",
        label.size = .25,
        size = 2.8,
        box.padding = .3,
        point.padding = .2,
        segment.size = .2,
        min.segment.length = 0,
        max.overlaps = Inf,
        
        inherit.aes = FALSE
      ) +
      
      facet_grid( ~ protein,
                  scales = "free_x",
                  space = "free_x",
      ) +
      #scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(clip = "off") +
      scale_colour_manual(values = c("none" = "black",
                                     "p_negative" = "#ff1493",
                                     "q_negative" = "#ff1493",
                                     "p_positive" = "#6b8e23",
                                     "q_positive" = "#6b8e23")) +
      scale_alpha_manual(values = c("none" = .5, 
                                    "p_negative" = .5,
                                    "q_negative" = 1,
                                    "p_positive" = .5,
                                    "q_positive" = 1)) +
      scale_size_manual(values = c("none" = .5, 
                                    "p_negative" = .5,
                                    "q_negative" = .8,
                                    "p_positive" = .5,
                                    "q_positive" = .8)) +
      scale_fill_manual(
        values = c(
          q_negative = "#ff1493",
          q_positive = "#6b8e23"
        )
      ) +
      
      labs(x = "Codon") +
      
      theme_classic() +
      theme(legend.position = 'none',
            panel.border = element_rect(fill = 'transparent'),
            text = element_text(color = 'black'),
            axis.text = element_text(size = 12,
                                     color = 'black'),
            strip.text = element_text(size = 14))} + 
  
  plot_layout(nrow = 2)

################################################################################
# STATS
################################################################################
#test_dN-dS.R

# 1. Global question: does the distribution of dN−dS differ between proteins?
# 2. Temporal question: within each protein, does the distribution change over 
#    infection?
# 3. Local question: are there specific windows that consistently deviate from 
#    neutrality? (genomic hotspot analysis).
################################################################################
selection_w4s1 %>%
  mutate(
    dNdS = dN - dS
  ) %>% 
  kruskal.test(
    dNdS ~ protein,
    data = .
  )

# Kruskal-Wallis rank sum test
# 
# data:  dNdS by protein
# Kruskal-Wallis chi-squared = 24.321, df = 2, p-value = 5.232e-06

selection_w4s1 %>%
  mutate(dNdS = dN - dS) %>%
  with(
    pairwise.wilcox.test(
      x = dNdS,
      g = protein,
      p.adjust.method = "BH"
    )
  )
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  dNdS and protein 
# 
# CP      delta  
# delta 0.00049 -      
#   RdRP  0.72562 1.1e-06
# 
# P value adjustment method: BH 

# 2.

selection_sample_summary <- selection_w4s1_time %>%
  mutate(dNdS = dN - dS) %>%
  group_by(
    sample,
    hpi,
    protein) %>%
  summarise(
    median_dNdS = median(dNdS, na.rm = TRUE),
    mean_dNdS = mean(dNdS, na.rm = TRUE),
    .groups = "drop"
  )


lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
)

lmer(
  dNdS ~ protein * hpi + (1 | sample),
  data = selection_w4s1_time
)

# LA INTERACCION NO MEJORA EL MODEL (P = 0.09)
anova(lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
)
,lmer(
  dNdS ~ protein * hpi + (1 | sample),
  data = selection_w4s1_time
))

library(lmerTest)

# MAIN EFFECT OF PROTEIN
anova(lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
))

anova(lmer(
  dNdS ~ protein * hpi + (1 | sample),
  data = selection_w4s1_time
))

car::Anova(lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
), type = 3)

# Protein and hpi has an effect. The effect of hpi is not different between proteins
# post-hoc

emmeans(lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
), ~ protein) %>% 
  pairs(.,adjust = "fdr")
  #cld(.,adjust = "fdr")

emmeans(lmer(
  dNdS ~ protein +(1 | sample),
  data = selection_w4s1_time
), ~ protein) %>% 
  pairs(.,adjust = "fdr")
#cld(.,adjust = "fdr")


lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
)

summary(lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
))


# GENERATING PREDICTION 
mod <- lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
)

newdat <- expand.grid(
  protein = unique(selection_w4s1_time$protein),
  hpi = sort(unique(selection_w4s1_time$hpi))
)

# Predicción marginal (sin efecto aleatorio)
newdat$pred <- predict(
  mod,
  newdata = newdat,
  re.form = NA
)

# Design matrix para efectos fijos
X <- model.matrix(
  ~ protein + hpi,
  data = newdat
)

# Varianza-covarianza de los efectos fijos
V <- vcov(mod)

# Error estándar de cada predicción
newdat$se <- sqrt(diag(X %*% V %*% t(X)))

# IC95%
newdat <- newdat %>%
  mutate(
    lwr = pred - 1.96 * se,
    upr = pred + 1.96 * se
  )


confint(lmer(
  dNdS ~ protein + hpi + (1 | sample),
  data = selection_w4s1_time
), parm = "hpi")

# 3
# TEST: Ho: dN-dS == 0
# Per window
################################################################################
selection_w4s1 %>%
  mutate(`dN-dS` = dN - dS)

################################################################################
# Test whether dN-dS differs from zero in each sliding window
#
# One-sided Wilcoxon signed-rank tests are performed independently for:
#   - positive selection (dN-dS > 0)
#   - purifying selection (dN-dS < 0)
#
# P-values are adjusted independently using the Benjamini-Hochberg procedure.
################################################################################

selection_hotspots_w4s1 <-
  
  selection_w4s1 %>%
  mutate(dNdS = dN - dS) %>%
  
  group_by(protein,
           window_start
  ) %>%
  
  summarise(n = sum(!is.na(dNdS)),
            median_dNdS = median(dNdS, na.rm = TRUE),
            mean_dNdS = mean(dNdS, na.rm = TRUE),
            
            p_positive = {
              x <- dNdS[!is.na(dNdS)]
              if(length(unique(x)) < 2){
                NA_real_
              } else {
                wilcox.test(
                  x,
                  mu = 0,
                  alternative = "greater",
                  exact = FALSE
                )$p.value
              }
            },
            
            p_negative = {
              x <- dNdS[!is.na(dNdS)]
              if(length(unique(x)) < 2){
                NA_real_
              } else {
                wilcox.test(
                  x,
                  mu = 0,
                  alternative = "less",
                  exact = FALSE
                )$p.value
              }
            },
            .groups = "drop"
  ) %>% 
  
  group_by(protein) %>%
  
  # MULTIPLE TEST CORRECTION: FDR BH
  mutate(
    q_positive = p.adjust(
      p_positive,
      method = "BH"
    ),
    q_negative = p.adjust(
      p_negative,
      method = "BH"
    ),
    hotspot =
      case_when(
        q_positive < 0.05 ~ "positive",
        q_negative < 0.05 ~ "negative",
        TRUE ~ "none"
      )
  ) %>%
  ungroup()

selection_hotspots_w4s1$hotspot %>% table()
# negative     none 
# 12     1684 

selection_hotspots_w4s1 %>% filter(hotspot != "none")
