################################################################################
# FIGURE 6 - Q20 VERSION
#
# LoFreq:
#   --min-bq 20
#   --min-alt-bq 20
#
# Uses:
#   all.snv.byStrand.bed_wider.q20
#   both.snv.af0.01.binMat.q20
#   both.snv.af0.01.pAFMat.q20
#   both.snv.af0.01.pAFMat_inAllPlus.q20
#   both.snv.up3samples_annot.q20

# MJ OLMO-UCEDA
# 2026/11/10
################################################################################

library(tidyverse)
library(ggrepel)
library(reshape2)
library(ComplexHeatmap)
library(ggnewscale)
library(patchwork)
library(tibble)

################################################################################
# FIG 6A
# AF CORRELATION
################################################################################

all.snv.byStrand.bed_wider.q20 %>%
  
  filter(
    `AF_+` > 0.01,
    `AF_-` > 0.01
  ) %>%
  
  ggplot(
    aes(
      `AF_-`,
      `AF_+`,
      fill = region
    )
  ) +
  
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dotted",
    linewidth = 0.2
  ) +
  
  geom_smooth(
    method = "lm",
    linewidth = 0.6,
    aes(color = region)
  ) +
  
  geom_point(
    shape = 21,
    size = 2,
    alpha = 0.8
  ) +
  
  geom_text_repel(
    data = \(x)
    x %>%
      filter(
        `AF_+` > 0.2 |
          `AF_-` > 0.2
      ) %>%
      mutate(
        label =
          if_else(
            aa_change == "synonymous",
            id_snv,
            SNV
          )
      ),
    aes(
      label = label,
      color = region
    ),
    size = 3,
    max.overlaps = Inf
  ) +
  
  scale_fill_manual(
    values = c(
      color.proteins,
      "5'UTR" = "gray30",
      "3'UTR" = "white"
    )
  ) +
  
  scale_color_manual(
    values = color.proteins
  ) +
  
  facet_grid(
    rna ~ aa_change
  ) +
  
  coord_equal() +
  
  lims(
    x = c(0,1),
    y = c(0,1)
  ) +
  
  labs(
    x = "AF(-)",
    y = "AF(+)"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    panel.border = element_rect(
      colour = "black",
      fill = "transparent"
    )
  )


all.snv.byStrand.bed_wider.q20 %>%
  filter(
    `AF_+` > 0.01,
    `AF_-` > 0.01,
    region %in% c("RdRP","CP","delta P")
  ) %>%
  group_by(
    region,
    aa_change,
    id_snv_chrm
  ) %>%
  summarise(
    recurrence = n_distinct(sample),
    .groups = "drop"
  ) 


################################################################################
# FIG 6A
# n = total observations
# N = unique SNVs
################################################################################

fig6A_labels.q20 <- all.snv.byStrand.bed_wider.q20 %>%
  
  filter(
    `AF_+` > 0.01,
    `AF_-` > 0.01
  ) %>%
  
  group_by(
    rna,
    region,
    aa_change
  ) %>%
  
  summarise(
    n = n(),
    N = n_distinct(id_snv_chrm),
    .groups = "drop"
  ) %>%
  
  arrange(
    rna,
    region,
    aa_change
  )

fig6A_labels.q20 %>%
   mutate(
    label = paste0(
      "n = ", n,
      "\nN = ", N
    )
  ) %>%
  
  dplyr::select(
    rna,
    region,
    aa_change,
    label
  )

################################################################################
# FIG 6B
# DISTRIBUTION OF THE RECURRENT SNVS (ONLY BOTH)
# with labels and space adapted to the distribution
################################################################################

all.snv.byStrand.bed_wider.q20 %>%
  filter(
    `AF_+` > 0.01,
    `AF_-` > 0.01,
    region %in% c("RdRP", "CP", "delta P")
  ) %>%
  
  # 1. recurrencia por SNV
  group_by(region, 
           aa_change, 
           id_snv_chrm) %>%
  summarise(
    recurrence = n_distinct(sample),
    .groups = "drop"
  ) %>%
  
  # 2. número de SNVs por combinación real
  group_by(region, aa_change, recurrence) %>%
  summarise(n = n(), .groups = "drop") %>%
  
  # 3. COMPLETE POR REGIÓN (clave)
  group_by(region, aa_change) %>%
  group_modify(~{
    df <- .x
    tibble(
      recurrence = seq(1, max(df$recurrence, na.rm = TRUE))
    ) %>%
      left_join(df, by = "recurrence") %>%
      mutate(n = tidyr::replace_na(n, 0))
  }) %>%
  ungroup() %>%
  
  mutate(
    region = factor(
      region,
      levels = c("5'UTR", "RdRP", "CP", "delta P", "3'UTR")
    ),
    recurrence = factor(recurrence)
  ) %>%
  
  ggplot(
    aes(
      x = recurrence,
      y = n,
      group = region,
      fill = region,
      linewidth = aa_change,
      linetype = aa_change,
      color = aa_change
    )
  ) +
  
  geom_col(
    position = position_dodge2(
      width = 0.9,
      preserve = "single"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      color.proteins,
      "5'UTR" = "gray30",
      "3'UTR" = "white"
    )
  ) +
  
  new_scale_fill() +
  
  geom_label(
    data = \(x) dplyr::filter(x, n > 0),
    aes(
      x = recurrence,
      y = n,
      label = n,
      fill = aa_change
    ),
    inherit.aes = FALSE,
    position = position_dodge2(width = 0.9, preserve = "single"),
    size = 3,
    color = "black",
    label.size = 0.25,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      "synonymous" = "white",
      "non-synonymous" = "gray50"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "synonymous" = "dotted",
      "non-synonymous" = "solid"
    )
  ) +
  
  scale_linewidth_manual(
    values = c(
      "synonymous" = 0.8,
      "non-synonymous" = 0.4
    )
  ) +
  
  scale_color_manual(
    values = c(
      "synonymous" = "white",
      "non-synonymous" = "black"
    )
  ) +
  
  facet_grid(~region, 
             scales = "free_x", 
             space = "free_x") +
  
  scale_y_continuous(
    breaks = seq(0, 
                 80, 
                 by = 5),
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Samples recurrence",
    y = "Number of SNVs\n(both strands AF > 0.01)"
  ) +
  
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    panel.background = element_blank(),
    strip.text = element_text(size = 12),
    panel.border = element_rect(
      colour = "black",
      fill = "transparent"
    )
  )



################################################################################
# FIG 6C - Q20
# SNVs present in both strands (AF > 0.01)
# and detected in at least 4 samples
# The AF+ is shown also in teh case is in more 
################################################################################

reshape2::melt(
  as.matrix(
    both.snv.af0.01.pAFMat_inAllPlus.q20
  )
) %>%
  
  rowwise() %>%
  
  mutate(
    
    timepoint =
      strsplit(
        as.character(Var2),
        "_"
      )[[1]][2],
    
    hpi =
      as.numeric(
        stringr::str_replace(
          timepoint,
          "h",
          ""
        )
      ),
    
    replica =
      strsplit(
        as.character(Var2),
        "_"
      )[[1]][3],
    
    sample =
      paste0(
        timepoint,
        " r",
        replica
      ),
    
    sample =
      factor(
        sample,
        levels = c(
          "0h r1","0h r2","0h r3",
          "2h r1","2h r2","2h r3",
          "4h r1","4h r2","4h r3",
          "6h r1","6h r2","6h r3",
          "8h r1","8h r2","8h r3",
          "12h r1","12h r2","12h r3",
          "18h r1","18h r2","18h r3",
          "22h r1","22h r2","22h r3",
          "26h r1","26h r2","26h r3",
          "32h r1","32h r2","32h r3",
          "38h r1","38h r2","38h r3",
          "44h r1","44h r2","44h r3"
        )
      )
    
  ) %>%
  
  arrange(hpi) %>%
  
  merge(
    all.snv.byStrand.bed.q20 %>%
      dplyr::select(
        rna,
        id_snv,
        id_snv_chrm,
        cds_region,
        aa_change,
        aa_ref,
        aa_alt,
        SNV,
        region
      ) %>%
      distinct(),
    by.x = "Var1",
    by.y = "id_snv_chrm"
  ) %>%
  
  distinct() %>%
  
  ggplot(
    aes(
      x = sample,
      y = Var1,
      fill = value
    )
  ) +
  
  geom_tile(
    color = "black",
    linewidth = 0.1
  ) +
  ## ADDING INFO ABOUT WHERE IS FOUND IF NOT IN BOTH STRANDS
  # only (+)
  geom_tile(
    data = fig6C.categories.q20 %>%
      filter(category == "plus"),
    aes(
      x = sample_plot,
      y = id_snv_chrm
    ),
    fill = NA,
    color = "black",
    linewidth = .4
  ) +
  
  # only (-)
  geom_tile(
    data = fig6C.categories.q20 %>%
      filter(category == "minus"),
    aes(
      x = sample_plot,
      y = id_snv_chrm
    ),
    fill = NA,
    color = "#dc267fff",
    linewidth = .4
  ) +
  
  geom_text(
    data = . %>%
      filter(value >= 0.1),
    aes(
      label = sprintf("%.2f", value)
    ),
    size = 2.5,
    color = "white"
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "black",
    limits = c(
      0,
      max(
        both.snv.af0.01.pAFMat_inAllPlus.q20,
        na.rm = TRUE
      )
    )
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_tile(
    data = both.snv.up3samples_annot.q20,
    aes(
      x = "rna",
      y = id_snv_chrm,
      fill = rna
    ),
    color = "black",
    linewidth = 0.1,
    inherit.aes = FALSE
  ) +
  
  geom_tile(
    data = both.snv.up3samples_annot.q20,
    aes(
      x = "region",
      y = id_snv_chrm,
      fill = region
    ),
    color = "black",
    linewidth = 0.1,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      color.proteins,
      color.rna
    ),
    name = "Region"
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Samples",
    y = "",
    fill = "AF(+)"
  ) +
  
  coord_equal() +
  
  theme(
    legend.position = "none",
    
    axis.text.x =
      element_text(
        size = 8,
        angle = 90,
        hjust = 1,
        color = "black"
      ),
    
    axis.text.y =
      element_text(
        face =
          ifelse(
            
            both.snv.up3samples_annot.q20[
              match(
                both.snv.up3samples_annot.q20$id_snv_chrm,
                rownames(
                  both.snv.af0.01.pAFMat_inAllPlus.q20
                )
              ),
              "aa_change"
            ] == "non-synonymous",
            
            "bold",
            
            "italic"
            
          ),
        
        color = "black"
      )
  )

################################################################################
# COUNTS FOR FIG 6C (Q20)
#
# both = detected in both strands
# plus = only (+) strand above threshold
# minus = only (-) strand above threshold
################################################################################

fig6C_counts.q20 <- all.snv.byStrand.bed_wider.q20 %>%
  
  filter(
    id_snv_chrm %in%
      both.snv.up3samples_annot.q20$id_snv_chrm
  ) %>%
  
  mutate(
    
    category = case_when(
      
      `AF_+` > 0.01 &
        `AF_-` > 0.01
      ~ "both",
      
      `AF_+` > 0.01 &
        `AF_-` <= 0.01
      ~ "plus",
      
      `AF_+` <= 0.01 &
        `AF_-` > 0.01
      ~ "minus",
      
      TRUE ~ NA_character_
      
    )
    
  ) %>%
  
  filter(!is.na(category)) %>%
  
  group_by(
    id_snv_chrm,
    category
  ) %>%
  
  summarise(
    n_samples = n_distinct(sample),
    .groups = "drop"
  ) %>%
  
  tidyr::pivot_wider(
    names_from = category,
    values_from = n_samples,
    values_fill = 0
  ) %>%
  
  right_join(
    both.snv.up3samples_annot.q20 %>%
      dplyr::select(
        id_snv_chrm,
        rna,
        region,
        aa_change
      ),
    by = "id_snv_chrm"
  ) %>%
  
  dplyr::select(
    id_snv_chrm,
    rna,
    region,
    aa_change,
    both,
    plus,
    minus
  )

fig6C_counts.q20 %>% 
  arrange(desc(both))

fig6C.categories.q20 <- all.snv.byStrand.bed_wider.q20 %>%
  
  mutate(
    
    category = case_when(
      
      `AF_+` > 0.01 &
        `AF_-` > 0.01 ~ "both",
      
      `AF_+` > 0.01 &
        `AF_-` <= 0.01 ~ "plus",
      
      `AF_+` <= 0.01 &
        `AF_-` > 0.01 ~ "minus",
      
      TRUE ~ "none"
      
    )
    
  ) %>%
  
  filter(
    id_snv_chrm %in%
      rownames(
        both.snv.af0.01.pAFMat_inAllPlus.q20
      )
  ) %>%
  
  mutate(
    
    timepoint =
      strsplit(sample, "_") |>
      sapply(`[`, 2),
    
    replica =
      strsplit(sample, "_") |>
      sapply(`[`, 3),
    
    sample_plot =
      factor(
        paste0(timepoint, " r", replica),
        levels = c(
          "0h r1","0h r2","0h r3",
          "2h r1","2h r2","2h r3",
          "4h r1","4h r2","4h r3",
          "6h r1","6h r2","6h r3",
          "8h r1","8h r2","8h r3",
          "12h r1","12h r2","12h r3",
          "18h r1","18h r2","18h r3",
          "22h r1","22h r2","22h r3",
          "26h r1","26h r2","26h r3",
          "32h r1","32h r2","32h r3",
          "38h r1","38h r2","38h r3",
          "44h r1","44h r2","44h r3"
        )
      )
    
  ) %>%
  
  dplyr::select(
    id_snv_chrm,
    sample_plot,
    category
  )

################################################################################
# FIGURE 6 - STATISTICS (Q20)
################################################################################
# Shared SNVs detected on both strands
# LoFreq:
#   --min-bq 20
#   --min-alt-bq 20
#
# MJ Olmo-Uceda
# Updated Q20 reanalysis
# 2026/11/20
################################################################################

################################################################################
# DATASET
################################################################################

shared.snvs.q20 <- all.snv.byStrand.bed_wider.q20 %>%
  filter(
    `AF_+` > 0.01,
    `AF_-` > 0.01
  )

shared.snvs.cds.q20 <- shared.snvs.q20 %>%
  filter(
    region %in% c("RdRP","CP","delta P")
  )

################################################################################
# FIGURE 6A
# AF CORRELATION BETWEEN STRANDS
################################################################################

cor.af.q20 <- cor.test(
  shared.snvs.q20$`AF_+`,
  shared.snvs.q20$`AF_-`
)

cor.af.q20

################################################################################
# EFFECT OF SYNONYMOUS/NON-SYNONYMOUS STATUS
################################################################################

mod.af.effect.q20 <- lm(
  `AF_+` ~ `AF_-` * aa_change,
  data = shared.snvs.cds.q20
)

summary(mod.af.effect.q20)

anova(mod.af.effect.q20)

################################################################################
# ESTIMATED SLOPES
################################################################################

slopes.aa.q20 <- emtrends(
  mod.af.effect.q20,
  ~ aa_change,
  var = "AF_-"
)

slopes.aa.q20

pairs(slopes.aa.q20)

################################################################################
# EFFECT OF PROTEIN
################################################################################

mod.af.region.q20 <- lm(
  `AF_+` ~ `AF_-` * region + aa_change,
  data = shared.snvs.cds.q20
)

summary(mod.af.region.q20)

anova(mod.af.region.q20)

################################################################################
# PROTEIN-SPECIFIC SLOPES
################################################################################

slopes.region.q20 <- emtrends(
  mod.af.region.q20,
  ~ region,
  var = "AF_-"
)

slopes.region.q20

pairs(
  slopes.region.q20,
  adjust = "BH"
)

6################################################################################
# COUNTS SHOWN IN FIGURE 6A
################################################################################

fig6A.labels.q20 <- shared.snvs.q20 %>%
  group_by(
    rna,
    region,
    aa_change
  ) %>%
  summarise(
    n = n(),
    N = n_distinct(id_snv_chrm),
    .groups = "drop"
  )

fig6A.labels.q20

################################################################################
# SNV-LEVEL DATASET
# One row per unique SNV
################################################################################

snv.recurrence.q20 <- shared.snvs.cds.q20 %>%
  group_by(
    id_snv_chrm,
    region,
    aa_change
  ) %>%
  summarise(
    recurrence = n_distinct(sample),
    .groups = "drop"
  ) %>%
  mutate(
    recurrent_ge2 = recurrence >= 2,
    recurrent_ge3 = recurrence >= 3,
    recurrent_ge4 = recurrence >= 4
  )

################################################################################
# OVERALL RECURRENCE
################################################################################

snv.recurrence.q20 %>%
  summarise(
    total = n(),
    recurrent = sum(recurrent_ge2),
    proportion = 100 * recurrent / total
  )

################################################################################
# RECURRENCE TABLE
################################################################################

snv.recurrence.q20 %>%
  count(
    region,
    aa_change,
    recurrence
  )

################################################################################
# ARE SYNONYMOUS VARIANTS MORE RECURRENT?
# >= 2 samples
################################################################################

mod.recurrent.ge2.q20 <- glm(
  recurrent_ge2 ~ aa_change + region,
  family = binomial,
  data = snv.recurrence.q20
)

summary(mod.recurrent.ge2.q20)

drop1(
  mod.recurrent.ge2.q20,
  test = "Chisq"
)

anova(
  mod.recurrent.ge2.q20,
  test = "Chisq"
)

################################################################################
# ODDS RATIOS
################################################################################

OR.ge2.q20 <- exp(
  cbind(
    OR = coef(mod.recurrent.ge2.q20),
    confint(mod.recurrent.ge2.q20)
  )
)

OR.ge2.q20

################################################################################
# >= 3 samples
################################################################################

mod.recurrent.ge3.q20 <- glm(
  recurrent_ge3 ~ aa_change + region,
  family = binomial,
  data = snv.recurrence.q20
)

summary(mod.recurrent.ge3.q20)

drop1(
  mod.recurrent.ge3.q20,
  test = "Chisq"
)

anova(
  mod.recurrent.ge3.q20,
  test = "Chisq"
)

################################################################################
# ODDS RATIOS
################################################################################

OR.ge3.q20 <- exp(
  cbind(
    OR = coef(mod.recurrent.ge3.q20),
    confint(mod.recurrent.ge3.q20)
  )
)

OR.ge3.q20

################################################################################
# NEGATIVE BINOMIAL MODEL
# RECURRENCE INTENSITY
################################################################################

mod.nb.q20 <- glm.nb(
  recurrence ~ aa_change + region,
  data = snv.recurrence.q20
)

summary(mod.nb.q20)

anova(
  mod.nb.q20,
  test = "Chisq"
)

################################################################################
# ARE SOME PROTEINS ENRICHED FOR NON-SYNONYMOUS VARIANTS?
################################################################################

mod.aa.region.q20 <- glm(
  aa_change ~ region,
  family = binomial,
  data = snv.recurrence.q20 %>%
    mutate(
      aa_change = relevel(
        factor(aa_change),
        ref = "non-synonymous"
      )
    )
)

summary(mod.aa.region.q20)

drop1(
  mod.aa.region.q20,
  test = "Chisq"
)

anova(
  mod.aa.region.q20,
  test = "Chisq"
)

################################################################################
# PAIRWISE CONTRASTS
################################################################################

pairs(
  emmeans(
    mod.aa.region.q20,
    ~ region
  ),
  adjust = "BH"
)

################################################################################
# RNA1 vs RNA2
# NON-SYNONYMOUS PROPORTION
################################################################################

rna.effect.q20 <- shared.snvs.q20 %>%
  filter(
    aa_change %in%
      c("synonymous","non-synonymous")
  ) %>%
  distinct(
    id_snv_chrm,
    aa_change,
    rna
  )

tab.rna.q20 <- table(
  rna.effect.q20$rna,
  rna.effect.q20$aa_change
)

tab.rna.q20

fisher.rna.q20 <- fisher.test(
  tab.rna.q20
)

fisher.rna.q20

################################################################################
# RECURRENCE PERCENTAGES USED IN TEXT
################################################################################

recurrence.summary.q20 <- snv.recurrence.q20 %>%
  group_by(
    region,
    aa_change
  ) %>%
  summarise(
    total = n(),
    recurrent = sum(recurrent_ge2),
    percentage = 100 * recurrent / total,
    .groups = "drop"
  )

recurrence.summary.q20

################################################################################
# FIGURE 6C
# MOST RECURRENT SHARED SNVs
################################################################################

shared.top.q20 <- shared.snvs.cds.q20 %>%
  group_by(
    id_snv_chrm,
    SNV,
    region,
    aa_change
  ) %>%
  summarise(
    recurrence = n_distinct(sample),
    max_AF_plus = max(`AF_+`),
    max_AF_minus = max(`AF_-`),
    mean_AF_plus = mean(`AF_+`),
    mean_AF_minus = mean(`AF_-`),
    .groups = "drop"
  ) %>%
  arrange(desc(recurrence))

shared.top.q20

################################################################################
# TABLE OF SHARED SNVs PRESENT IN >=4 SAMPLES
################################################################################

shared.top.q20 %>%
  filter(recurrence >= 4) %>%
  arrange(desc(recurrence))

################################################################################
# WILCOXON ANALYSIS
# RECURRENT SNVs VS BACKGROUND
################################################################################

################################################################################
# RECURRENT >1 SAMPLE
################################################################################

both.snv.up1sample.q20 <- shared.snvs.q20 %>%
  group_by(id_snv_chrm) %>%
  summarise(
    recurrence = n_distinct(sample),
    .groups = "drop"
  ) %>%
  filter(recurrence > 1) %>%
  pull(id_snv_chrm)

################################################################################
# RECURRENT >2 SAMPLES
################################################################################

both.snv.up2sample.q20 <- shared.snvs.q20 %>%
  group_by(id_snv_chrm) %>%
  summarise(
    recurrence = n_distinct(sample),
    .groups = "drop"
  ) %>%
  filter(recurrence > 2) %>%
  pull(id_snv_chrm)

################################################################################
# RECURRENT >3 SAMPLES
################################################################################

both.snv.up3sample.q20 <- shared.snvs.q20 %>%
  group_by(id_snv_chrm) %>%
  summarise(
    recurrence = n_distinct(sample),
    .groups = "drop"
  ) %>%
  filter(recurrence > 3) %>%
  pull(id_snv_chrm)

################################################################################
# FUNCTION
################################################################################

wilcox_recurrent.q20 <- function(recurrent_snvs){
  
  all.snv.byStrand.bed_wider.q20 %>%
    
    mutate(
      is_recurrent = id_snv_chrm %in% recurrent_snvs,
      AF_mean = (`AF_+` + `AF_-`) / 2
    ) %>%
    
    group_by(sample) %>%
    
    summarise(
      
      hpi = unique(hpi),
      
      n_recurrent = sum(is_recurrent),
      
      n_background = sum(!is_recurrent),
      
      p = if(
        n_recurrent > 0 &
        n_background > 0
      ){
        
        wilcox.test(
          AF_mean[is_recurrent],
          AF_mean[!is_recurrent],
          alternative = "greater",
          exact = FALSE
        )$p.value
        
      } else {
        
        NA_real_
        
      },
      
      median_recurrent =
        if(n_recurrent > 0)
          median(AF_mean[is_recurrent], na.rm = TRUE)
      else
        NA_real_,
      
      median_background =
        if(n_background > 0)
          median(AF_mean[!is_recurrent], na.rm = TRUE)
      else
        NA_real_,
      
      effect_size =
        median_recurrent -
        median_background,
      
      .groups = "drop"
      
    ) %>%
    
    mutate(
      FDR = p.adjust(
        p,
        method = "BH"
      )
    ) %>%
    
    arrange(hpi)
  
}

################################################################################
# CALCULATE ALL
################################################################################

results.wilcox.q20 <- list(
  
  recurrent_gt1 =
    wilcox_recurrent.q20(
      both.snv.up1sample.q20
    ),
  
  recurrent_gt2 =
    wilcox_recurrent.q20(
      both.snv.up2sample.q20
    ),
  
  recurrent_gt3 =
    wilcox_recurrent.q20(
      both.snv.up3sample.q20
    )
  
)

results.wilcox.q20

################################################################################
# EXPORT
################################################################################

writexl::write_xlsx(
  results.wilcox.q20,
  "REVIEW_GENETICS/wilcox.AF.recurrents_vs_background.Q20.xlsx"
)

################################################################################
# END
################################################################################
