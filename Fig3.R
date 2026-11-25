################################################################################
# A) 
################################################################################
################################################################################
# POSITION-WISE RESIDUALS
################################################################################

meanResiduals.df.q20 <- df_complete.normH_wider_q20 %>%
  filter(hpi > 2) %>%
  group_by(rna, pos, hpi) %>%
  summarise(
    meanResidual = mean(
      Sn.residuals,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  filter(!is.nan(meanResidual))

################################################################################
# SEGMENT-SPECIFIC THRESHOLDS
################################################################################

meanSnResiduals.quantiles.q20 <- meanResiduals.df.q20 %>%
  group_by(rna) %>%
  summarise(
    Q1  = quantile(meanResidual, 0.25, na.rm = TRUE),
    Q3  = quantile(meanResidual, 0.75, na.rm = TRUE),
    P5  = quantile(meanResidual, 0.05, na.rm = TRUE),
    P95 = quantile(meanResidual, 0.95, na.rm = TRUE),
    P1  = quantile(meanResidual, 0.01, na.rm = TRUE),
    P99 = quantile(meanResidual, 0.99, na.rm = TRUE),
    .groups = "drop"
  )

################################################################################
# CLASSIFY POSITIONS
################################################################################

meanResiduals.classified.q20 <- meanResiduals.df.q20 %>%
  left_join(
    meanSnResiduals.quantiles.q20,
    by = "rna"
  ) %>%
  mutate(
    class = case_when(
      
      meanResidual < P1 |
        meanResidual > P99 ~ "extreme",
      
      meanResidual < P5 |
        meanResidual > P95 ~ "strong",
      
      meanResidual < Q1 |
        meanResidual > Q3 ~ "moderate",
      
      TRUE ~ "central"
    )
  )

################################################################################
# TOP PANEL
################################################################################

p.top.q20 <- meanResiduals.classified.q20 %>%
  ggplot(
    aes(
      x = pos,
      y = meanResidual,
      color = class
    )
  ) +
  
  geom_rect(
    data = orv.proteins %>%
      filter(protein != "putative"),
    aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf,
      fill = protein
    ),
    alpha = 0.05,
    inherit.aes = FALSE
  ) +
  
  geom_rect(
    data = palmprint.orv,
    aes(
      xmin = nt.start,
      xmax = nt.end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "grey50",
    color = "black",
    alpha = 0.15,
    linewidth = 0.05,
    inherit.aes = FALSE
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.2,
    color = "black"
  ) +
  
  geom_segment(
    aes(
      xend = pos,
      y = 0,
      yend = meanResidual
    ),
    linewidth = 0.12,
    alpha = 0.25,
    color = "black"
  ) +
  
  geom_point(
    size = 0.35,
    alpha = 0.85
  ) +
  
  scale_color_manual(
    values = c(
      central  = "grey90",
      moderate = "grey55",
      strong   = "#4A4E69",
      extreme  = "#8792F1"
    )
  ) +
  
  scale_fill_manual(
    values = color.proteins
  ) +
  
  ylab(
    expression(
      paste(
        bar(S[n]),
        "(p) residuals (",
        S[italic(np)]^"+",
        " - ",
        S[italic(np)]^"-",
        ")"
      )
    )
  ) +
  
  xlab(NULL) +
  
  scale_x_continuous(expand = c(0,0)) +
  
  facet_grid(
    ~rna,
    scales = "free_x",
    space = "free_x"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    )
  )


################################################################################
# B)
################################################################################
################################################################################
# SLIDING WINDOWS
################################################################################

window.size <- 50
step.size   <- 20

sliding.df.q20 <- meanResiduals.df.q20 %>%
  group_by(rna) %>%
  arrange(pos) %>%
  group_modify(~{
    
    df <- .x
    
    max.pos <- max(
      df$pos,
      na.rm = TRUE
    )
    
    windows <- tibble(
      window.start = seq(
        1,
        max.pos - window.size,
        by = step.size
      )
    ) %>%
      mutate(
        window.end =
          window.start +
          window.size - 1
      )
    
    windows %>%
      rowwise() %>%
      mutate(
        mean.window.residual =
          mean(
            df$meanResidual[
              df$pos >= window.start &
                df$pos <= window.end
            ],
            na.rm = TRUE
          )
      ) %>%
      ungroup()
    
  }) %>%
  ungroup()

sliding.df.q20 <- meanResiduals.df.q20 %>%
  group_by(rna) %>%
  arrange(pos) %>%
  group_modify(~{
    
    df <- .x
    
    max.pos <- max(df$pos)
    
    windows <- tibble(
      window.start = seq(
        1,
        max.pos - window.size,
        by = step.size
      )
    ) %>%
      mutate(
        window.end = window.start + window.size - 1
      )
    
    windows %>%
      rowwise() %>%
      mutate(
        
        mean.window.residual =
          mean(
            df$meanResidual[
              df$pos >= window.start &
                df$pos <= window.end
            ],
            na.rm = TRUE
          ),
        
        sd.window.residual =
          sd(
            df$meanResidual[
              df$pos >= window.start &
                df$pos <= window.end
            ],
            na.rm = TRUE
          ),
        
        n =
          sum(
            !is.na(
              df$meanResidual[
                df$pos >= window.start &
                  df$pos <= window.end
              ]
            )
          ),
        
        sem.window.residual =
          sd.window.residual / sqrt(n)
        
      ) %>%
      ungroup()
    
  }) %>%
  ungroup()

################################################################################
# B)
################################################################################

p.bottom.q20 <- sliding.df.q20 %>%
  ggplot(
    aes(
      x = window.start,
      y = mean.window.residual
    )
  ) +
  
  geom_rect(
    data = orv.proteins %>%
      filter(protein != "putative"),
    aes(
      xmin = start,
      xmax = end,
      ymin = -Inf,
      ymax = Inf,
      fill = protein
    ),
    alpha = 0.05,
    inherit.aes = FALSE
  ) +
  
  geom_rect(
    data = palmprint.orv,
    aes(
      xmin = nt.start,
      xmax = nt.end,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "grey50",
    color = "black",
    alpha = 0.15,
    linewidth = 0.05,
    inherit.aes = FALSE
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.25
  ) +
  
  geom_ribbon(
    aes(
      ymin = mean.window.residual - sem.window.residual,
      ymax = mean.window.residual + sem.window.residual
    ),
    fill = "grey20",
    alpha = 0.4
  ) +
  
  geom_line(
    linewidth = 0.4,
    color = "black"
  ) +
  
  scale_fill_manual(
    values = color.proteins
  ) +
  
  scale_x_continuous(expand = c(0,0)) +
  
  ylab(
    "Sliding-window \nmean residuals"
  ) +
  
  xlab(NULL) +
  
  facet_grid(
    ~rna,
    scales = "free_x",
    space = "free_x"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    ),
    strip.text = element_blank()
  )


################################################################################
# C)
################################################################################
p.mfe <-ggplot(rnafold.df %>%
                 mutate(rna = ifelse(chrom == "RNA1",
                                     "RNA 1",
                                     "RNA 2")),
               aes(x = window_mid,
                   y = MFE_plot,
                   color = strand,
                   fill = strand)) +
  # Coding regions
  ###########################################
  geom_rect(data = orv.proteins %>%
              filter(protein != "putative"),
            aes(xmin = start,
                xmax = end,
                ymin = -Inf,
                ymax = Inf,
                fill = protein),
            alpha = 0.05,
            inherit.aes = FALSE) +
  
  
  # Palmprint motifs
  ###########################################
  geom_rect(data = palmprint.orv,
            aes(xmin = nt.start,
                xmax = nt.end,
                ymin = -Inf,
                ymax = Inf),
            color = "black",
            fill = "grey50",
            alpha = 0.15,
            linewidth = 0.05,
            inherit.aes = FALSE) +
  
  # # high asymmetry regions (by windows without promediate)
  # geom_rect(
  #   data = high.diff %>%
  #     mutate(rna = ifelse(chrom == "RNA1",
  #                         "RNA 1",
  #                         "RNA 2")),
  #   inherit.aes = FALSE,
  #   aes(
  #     xmin = window_mid - 25,
  #     xmax = window_mid + 25,
  #     ymin = 20,
  #     ymax = Inf
  #   ),
  #   fill = "gold",
  #   alpha = 0.12
  # ) +

  # # Structural-bias hotspots
  # ###########################################
  geom_rect(
    data = merged.df %>%
      filter(mfe.hotspot),

    inherit.aes = FALSE,

    aes(
      xmin = window.start,
      xmax = window.end,
      ymin = -Inf,
      ymax = Inf
    ),

    fill = "gold",
    alpha = 0.4
  ) +
  # uncertainty ribbon
  geom_ribbon(
    aes(
      ymin = MFE_plot - diversity_scaled,
      ymax = MFE_plot + diversity_scaled
    ),
    alpha = 0.7,
    linewidth = 0
  ) +
  
  geom_line(linewidth = 0.6) +
  scale_x_continuous(expand = c(0,0)) +
  
  facet_grid( ~ rna,
              scales = "free_x",
              space = "free_x") +
  
  xlab("Genomic position (nt)") +
  ylab("-MFE (kcal/mol)") +
  
  theme_classic() +
  scale_color_manual(values = color.strand) +
  scale_fill_manual(values = c(color.strand,
                               color.proteins)) +
  
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.line = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA
    ),
    #axis.title.x = element_blank(),
    strip.text = element_blank()
  )




p.top / p.bottom / p.mfe +
  plot_layout(heights = c(1.5, .5, 1))

################################################################################
# STATS
################################################################################
# STRUCTURE BIAS vs DIVERSITY BIAS
# Q20 REANALYSIS
################################################################################

# Reproduces all statistics reported in the manuscript using:
#   - sliding.df.q20 (new diversity residuals)
#   - rnafold.df      (unchanged structure predictions)
#   - gc.df           (unchanged GC windows)
#
# Objects assumed to already exist:
#   sliding.df.q20
#   rnafold.df
#   gc.df
#
# MJ OLMO-UCEDA
# 2026/11/17
################################################################################

library(tidyverse)
library(mgcv)
library(lme4)

################################################################################
# BUILD MFE BIAS DATASET
################################################################################

# MFE bias:
# positive values  -> more structure in (+)
# negative values  -> more structure in (-)

# mfe.bias.df <- rnafold.df %>%
#   
#   mutate(
#     rna = ifelse(
#       chrom == "RNA1",
#       "RNA 1",
#       "RNA 2"
#     )
#   ) %>%
#   
#   group_by(
#     rna,
#     window_mid
#   ) %>%
#   
#   summarise(
#     
#     mfe.plus =
#       mean(
#         MFE_plot[strand == "+"],
#         na.rm = TRUE
#       ),
#     
#     mfe.minus =
#       mean(
#         MFE_plot[strand == "-"],
#         na.rm = TRUE
#       ),
#     
#     mfe.bias =
#       mfe.plus - mfe.minus,
#     
#     .groups = "drop"
#     
#   )

################################################################################
# MERGE STRUCTURE + DIVERSITY (Q20)
################################################################################

merged.df.q20 <- sliding.df.q20 %>%
  rowwise() %>%
  mutate(
    mfe.bias = mean(
      mfe.bias.df$mfe.bias[
        mfe.bias.df$rna == rna &
          mfe.bias.df$window_mid >= window.start &
          mfe.bias.df$window_mid <= window.end
      ],
      na.rm = TRUE
    )
  ) %>%
  ungroup()

################################################################################
# ADD GC CONTENT
################################################################################

merged.df.q20 <- merged.df.q20 %>%
  mutate(
    GC = purrr::pmap_dbl(
      list(
        rna,
        window.start,
        window.end
      ),
      ~ mean(
        gc.df$GC_content[
          gc.df$rna == ..1 &
            gc.df$window.mid >= ..2 &
            gc.df$window.mid <= ..3],
        na.rm = TRUE)))

################################################################################
# DIVERSITY HOTSPOTS 
################################################################################

threshold.div <- quantile(
  abs(merged.df.q20$mean.window.residual),
  0.95,
  na.rm = TRUE
)

merged.df.q20 <- merged.df.q20 %>%
  mutate(
    hotspot =
      abs(mean.window.residual) >=
      threshold.div
  )

table(
  merged.df.q20$hotspot
)

################################################################################
# STRUCTURAL HOTSPOTS (segment-specific)
################################################################################

q95.mfe.rna1 <- quantile(
  abs(
    merged.df.q20$mfe.bias[
      merged.df.q20$rna == "RNA 1"
    ]
  ),
  0.95,
  na.rm = TRUE
)

q95.mfe.rna2 <- quantile(
  abs(
    merged.df.q20$mfe.bias[
      merged.df.q20$rna == "RNA 2"
    ]
  ),
  0.95,
  na.rm = TRUE
)

merged.df.q20 <- merged.df.q20 %>%
  mutate(
    mfe.hotspot =
      abs(mfe.bias) >
      ifelse(
        rna == "RNA 1",
        q95.mfe.rna1,
        q95.mfe.rna2
      )
  )

table(
  merged.df.q20$mfe.hotspot
)

################################################################################
# GLOBAL CORRELATION
#
# Reported in manuscript:
# rs = ...
# P = ...
################################################################################

cor.global.q20 <- cor.test(
  merged.df.q20$mean.window.residual,
  merged.df.q20$mfe.bias,
  method = "spearman"
)

cor.global.q20

################################################################################
# RNA1
################################################################################

cor.rna1.q20 <- cor.test(
  
  merged.df.q20 %>%
    filter(
      rna == "RNA 1"
    ) %>%
    pull(
      mean.window.residual
    ),
  
  merged.df.q20 %>%
    filter(
      rna == "RNA 1"
    ) %>%
    pull(
      mfe.bias
    ),
  
  method = "spearman"
  
)

cor.rna1.q20

################################################################################
# RNA2
################################################################################

cor.rna2.q20 <- cor.test(
  
  merged.df.q20 %>%
    filter(
      rna == "RNA 2"
    ) %>%
    pull(
      mean.window.residual
    ),
  
  merged.df.q20 %>%
    filter(
      rna == "RNA 2"
    ) %>%
    pull(
      mfe.bias
    ),
  method = "spearman"
)

cor.rna2.q20

################################################################################
# DEGREES OF FREEDOM
################################################################################

n.global.df <- nrow(merged.df.q20) - 2

n.rna1.df <- merged.df.q20 %>%
  filter(rna == "RNA 1") %>%
  nrow() - 2

n.rna2.df <- merged.df.q20 %>%
  filter(rna == "RNA 2") %>%
  nrow() - 2

n.global.df
n.rna1.df
n.rna2.df

################################################################################
# MAGNITUDE RELATIONSHIP ONLY
################################################################################

cor.abs.q20 <- cor.test(
  abs(merged.df.q20$mean.window.residual),
  abs(merged.df.q20$mfe.bias),
  method = "spearman"
)

cor.abs.q20

################################################################################
# PERMUTATION TEST
#
# Reported in manuscript
################################################################################

set.seed(123)

obs.cor <- cor(
  merged.df.q20$mean.window.residual,
  merged.df.q20$mfe.bias,
  method = "spearman"
)

perm.cor <- replicate(
  5000,
  cor(
    sample(
      merged.df.q20$mean.window.residual
    ),
    merged.df.q20$mfe.bias,
    method = "spearman"
  )
)

perm.p <- mean(
  abs(perm.cor) >= abs(obs.cor)
)
perm.p

################################################################################
# DO DIVERSITY HOTSPOTS SHOW LARGER STRUCTURAL BIAS?
#
# Reported in manuscript
################################################################################

wilcox.div.hotspot <- wilcox.test(
  abs(mfe.bias) ~ hotspot,
  alternative = "less",
  data = merged.df.q20
)
wilcox.div.hotspot

################################################################################
# DO STRUCTURAL HOTSPOTS SHOW LARGER DIVERSITY BIAS?
#
# Reported in manuscript
################################################################################

wilcox.mfe.hotspot <- wilcox.test(
  abs(mean.window.residual) ~ mfe.hotspot,
  alternative = "less",
  data = merged.df.q20
)

wilcox.mfe.hotspot

################################################################################
# DIRECTIONAL HOTSPOTS
#
# Used for Fisher tests in original analysis
################################################################################

q95.plus.div <- quantile(
  merged.df.q20$mean.window.residual[
    merged.df.q20$mean.window.residual > 0
  ],
  0.95,
  na.rm = TRUE
)

q05.minus.div <- quantile(
  merged.df.q20$mean.window.residual[
    merged.df.q20$mean.window.residual < 0
  ],
  0.05,
  na.rm = TRUE  
)

merged.df.q20 <- merged.df.q20 %>%
  mutate(
    hotspot.diversity.sign = case_when(
      mean.window.residual >= q95.plus.div
      ~ "plus_diversity",
      mean.window.residual <= q05.minus.div
      ~ "minus_diversity",
      TRUE
      ~ "none"
    )
  )

################################################################################
# DIRECTIONAL STRUCTURAL HOTSPOTS
################################################################################

q95.plus.mfe.rna1 <- quantile(
  merged.df.q20$mfe.bias[
    merged.df.q20$rna == "RNA 1" &
      merged.df.q20$mfe.bias > 0
  ],
  0.95,
  na.rm = TRUE
  
)

q05.minus.mfe.rna1 <- quantile(
  merged.df.q20$mfe.bias[
    merged.df.q20$rna == "RNA 1" &
      merged.df.q20$mfe.bias < 0
  ],
  0.05,
  na.rm = TRUE
)

q95.plus.mfe.rna2 <- quantile(
  merged.df.q20$mfe.bias[
    merged.df.q20$rna == "RNA 2" &
      merged.df.q20$mfe.bias > 0
  ],
  0.95,
  na.rm = TRUE
  
)

q05.minus.mfe.rna2 <- quantile(
  merged.df.q20$mfe.bias[
    merged.df.q20$rna == "RNA 2" &
      merged.df.q20$mfe.bias < 0
  ],
  0.05,
  na.rm = TRUE
)

merged.df.q20 <- merged.df.q20 %>%
  mutate(
    hotspot.structure.sign = case_when(
      rna == "RNA 1" &
        mfe.bias >= q95.plus.mfe.rna1
      ~ "plus_stable",
      rna == "RNA 1" &
        mfe.bias <= q05.minus.mfe.rna1
      ~ "minus_stable",
      rna == "RNA 2" &
        mfe.bias >= q95.plus.mfe.rna2
      ~ "plus_stable",
      rna == "RNA 2" &
        mfe.bias <= q05.minus.mfe.rna2
      ~ "minus_stable",
      TRUE
      ~ "none"
    )
  )

################################################################################
# FISHER TEST
#
# Reproduces directional hotspot overlap test
################################################################################

hotspots.q20.table <- table(
  merged.df.q20$hotspot.diversity.sign,
  merged.df.q20$hotspot.structure.sign
)

hotspots.q20.table

fisher.test(
  hotspots.q20.table
)

################################################################################
# LINEAR MODEL
#
# Reported in exploratory section
################################################################################

lm.bias.q20 <- lm(
  abs(mean.window.residual) ~
    abs(mfe.bias),
  data = merged.df.q20
)
summary(lm.bias.q20)

################################################################################
# GC AS CONFOUNDER
################################################################################

lm.bias.gc.q20 <- lm(
  abs(mean.window.residual) ~
    abs(mfe.bias) +
    GC,
  data = merged.df.q20 
)
summary(lm.bias.gc.q20)

################################################################################
# NON-LINEAR MODEL
################################################################################

gam.bias.q20 <- gam(
  abs(mean.window.residual) ~
    s(abs(mfe.bias)) +
    s(GC),
  data = merged.df.q20
  
)

summary(gam.bias.q20)

################################################################################
# SAVE DATASET
################################################################################

write.csv(
  merged.df.q20,
  "merged_df_structure_diversity_q20.csv",
  row.names = FALSE
)

################################################################################
# QUICK SUMMARY TABLE
################################################################################

tibble(
  analysis = c(
    "global_spearman",
    "rna1_spearman",
    "rna2_spearman",
    "permutation",
    "wilcox_div_hotspots",
    "wilcox_structure_hotspots"
  ),
  
  statistic = c(
    unname(cor.global.q20$estimate),
    unname(cor.rna1.q20$estimate),
    unname(cor.rna2.q20$estimate),
    perm.p,
    wilcox.div.hotspot$p.value,
    wilcox.mfe.hotspot$p.value
  )
)

