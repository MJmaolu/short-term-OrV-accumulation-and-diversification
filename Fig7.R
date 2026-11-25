################################################################################
# FIGURE 7 (Q20 VERSION)
#
# Recalculation of dN-dS sliding-window analyses after Phred >=20 filtering
#
# INPUT REQUIRED:
#   selection_table.q20
#
# OUTPUT:
#   selection_w4s1.q20
#   selection_w4s2.q20
#   selection_sample_summary.q20
#   selection_hotspots_w4s1.q20
#   mod.q20
#   newdat.q20
#
# MJ OLMO-UCEDA
# 2026/11/10
################################################################################
################################################################################
# CALCULATE pN, pS, dN, dS AND OMEGA IN SLIDING WINDOWS
################################################################################

calculate_selection_windows <- function(data,
                                        window_size = 4,
                                        step_size = 1){
  
  data <- data %>%
    arrange(protein, codon_position)
  
  window_results <- list()
  result_index <- 1
  
  for(current_protein in unique(data$protein)){
    
    protein_data <- data %>%
      filter(protein == current_protein)
    
    last_window_start <-
      max(protein_data$codon_position) -
      window_size +
      1
    
    if(last_window_start < 1){
      next
    }
    
    for(window_start in seq(
      from = 1,
      to = last_window_start,
      by = step_size
    )){
      
      window_end <- window_start + window_size - 1
      
      window_data <- protein_data %>%
        filter(
          codon_position >= window_start,
          codon_position <= window_end
        )
      
      synonymous_sites <-
        sum(window_data$synonymous_sites)
      
      nonsynonymous_sites <-
        sum(window_data$nonsynonymous_sites)
      
      synonymous_mutations <-
        sum(window_data$synonymous_mutations)
      
      nonsynonymous_mutations <-
        sum(window_data$nonsynonymous_mutations)
      
      pS <- if(synonymous_sites > 0){
        synonymous_mutations /
          synonymous_sites
      } else {
        NA_real_
      }
      
      pN <- if(nonsynonymous_sites > 0){
        nonsynonymous_mutations /
          nonsynonymous_sites
      } else {
        NA_real_
      }
      
      dS <- if(
        !is.na(pS) &&
        (1 - 4 * pS / 3) > 0
      ){
        -(3/4) * log(1 - 4 * pS / 3)
      } else {
        NA_real_
      }
      
      dN <- if(
        !is.na(pN) &&
        (1 - 4 * pN / 3) > 0
      ){
        -(3/4) * log(1 - 4 * pN / 3)
      } else {
        NA_real_
      }
      
      omega <- if(
        !is.na(dS) &&
        dS > 0
      ){
        dN / dS
      } else {
        NA_real_
      }
      
      window_results[[result_index]] <-
        
        tibble(
          
          sample =
            unique(window_data$sample),
          
          protein =
            current_protein,
          
          window_start =
            window_start,
          
          window_end =
            window_end,
          
          window_size =
            window_size,
          
          synonymous_sites =
            synonymous_sites,
          
          nonsynonymous_sites =
            nonsynonymous_sites,
          
          synonymous_mutations =
            synonymous_mutations,
          
          nonsynonymous_mutations =
            nonsynonymous_mutations,
          
          pS = pS,
          pN = pN,
          
          dS = dS,
          dN = dN,
          
          omega = omega
          
        )
      
      result_index <- result_index + 1
    }
  }
  
  bind_rows(window_results)
  
}

################################################################################
# Q20 VERSION
################################################################################

snvs_4selection_annotated.q20 <- all.snv.byStrand.bed.q20 %>%
  filter(
    strand == "+",
    AF > 0.01,
    DP > 20,
    cds_region == "coding"
  ) %>%
  left_join(
    genomic_position_annotation,
    by = c("rna", "POS")
  )

################################################################################
# OBSERVED MUTATIONS
################################################################################

observed_mutations.q20 <- snvs_4selection_annotated.q20 %>%
  group_by(
    sample,
    protein,
    rna,
    codon_position
  ) %>%
  summarise(
    synonymous_mutations =
      sum(AF[aa_change == "synonymous"]),
    
    nonsynonymous_mutations =
      sum(AF[aa_change == "non-synonymous"]),
    
    synonymous_variants =
      sum(aa_change == "synonymous"),
    
    nonsynonymous_variants =
      sum(aa_change == "non-synonymous"),
    
    .groups = "drop"
  )

################################################################################
# MERGE OBSERVED + EXPECTED
################################################################################

selection_table.q20 <- sample_codon_table %>%
  
  left_join(
    observed_mutations.q20,
    by = c(
      "sample",
      "protein",
      "rna",
      "codon_position"
    )
  ) %>%
  
  mutate(
    synonymous_mutations =
      replace_na(synonymous_mutations, 0),
    
    nonsynonymous_mutations =
      replace_na(nonsynonymous_mutations, 0),
    
    synonymous_variants =
      replace_na(synonymous_variants, 0),
    
    nonsynonymous_variants =
      replace_na(nonsynonymous_variants, 0)
  )

################################################################################
# GENERATE WINDOW DATASETS
################################################################################

selection_w4s1.q20 <- selection_table.q20 %>%
  group_by(sample) %>%
  group_modify(
    ~ calculate_selection_windows(
      .x,
      window_size = 4,
      step_size = 1
    )
  ) %>%
  ungroup()

selection_w4s2.q20 <- selection_table.q20 %>%
  group_by(sample) %>%
  group_modify(
    ~ calculate_selection_windows(
      .x,
      window_size = 4,
      step_size = 2
    )
  ) %>%
  ungroup()

################################################################################
# GLOBAL TESTS (FIGURE 7A)
################################################################################

cat("\n")
cat("############################################################\n")
cat("FIGURE 7A\n")
cat("############################################################\n")

kw.q20 <- selection_w4s1.q20 %>%
  mutate(
    dNdS = dN - dS
  ) %>%
  kruskal.test(
    dNdS ~ protein,
    data = .
  )

print(kw.q20)

pw.q20 <- selection_w4s1.q20 %>%
  mutate(
    dNdS = dN - dS
  ) %>%
  with(
    pairwise.wilcox.test(
      x = dNdS,
      g = protein,
      p.adjust.method = "BH"
    )
  )

print(pw.q20)

################################################################################
# SAMPLE SUMMARIES
################################################################################

selection_sample_summary.q20 <- selection_w4s1.q20 %>%
  left_join(
    samples_info %>%
      dplyr::select(
        sample,
        hpi,
        replica
      ),
    by = "sample"
  ) %>%
  mutate(
    dNdS = dN - dS
  ) %>%
  group_by(
    sample,
    hpi,
    protein
  ) %>%
  summarise(
    mean_dNdS =
      mean(dNdS,
           na.rm = TRUE),
    
    median_dNdS =
      median(dNdS,
             na.rm = TRUE),
    
    .groups = "drop"
  )

################################################################################
# LINEAR MIXED MODEL (FIGURE 7B)
################################################################################

selection_w4s1_time.q20 <- selection_w4s1.q20 %>%
  left_join(
    samples_info %>%
      dplyr::select(
        sample,
        hpi,
        replica
      ),
    by = "sample"
  ) %>%
  mutate(
    dNdS = dN - dS
  )

m0.q20 <- lmer(
    dNdS ~ protein +
      hpi +
      (1 | sample),
    data = selection_w4s1_time.q20
  )

m1.q20 <- lmer(
    dNdS ~ protein *
      hpi +
      (1 | sample),
    data = selection_w4s1_time.q20
  )

cat("\n")
cat("############################################################\n")
cat("MODEL COMPARISON\n")
cat("############################################################\n")

print(
  anova(
    m0.q20,
    m1.q20
  )
)

cat("\n")
cat("############################################################\n")
cat("TYPE III ANOVA\n")
cat("############################################################\n")

print(
  car::Anova(
    m0.q20,
    type = 3
  )
)

cat("\n")
cat("############################################################\n")
cat("MODEL SUMMARY\n")
cat("############################################################\n")

print(
  summary(m0.q20)
)

cat("\n")
cat("############################################################\n")
cat("POST-HOC COMPARISONS\n")
cat("############################################################\n")

print(
  emmeans(
    m0.q20,
    ~ protein
  ) %>%
    pairs(
      adjust = "fdr"
    )
)

################################################################################
# PREDICTIONS FOR FIGURE 7B
################################################################################

newdat.q20 <- expand.grid(
    protein =
      unique(
        selection_w4s1_time.q20$protein
      ),
    hpi =
      sort(
        unique(
          selection_w4s1_time.q20$hpi
        )
      )
  )

newdat.q20$pred <- predict(
    m0.q20,
    newdata = newdat.q20,
    re.form = NA
  )

X.q20 <- model.matrix(
    ~ protein + hpi,
    data = newdat.q20
  )

V.q20 <- vcov(m0.q20)

newdat.q20$se <-sqrt(diag(X.q20 %*% V.q20 %*% t(X.q20)))

newdat.q20 <- newdat.q20 %>%
  mutate(
    lwr = pred - 1.96 * se,
    upr = pred + 1.96 * se
  )

################################################################################
# HOTSPOT ANALYSIS (FIGURE 7C)
################################################################################

selection_hotspots_w4s1.q20 <- selection_w4s1.q20 %>%
  mutate(
    dNdS = dN - dS
  ) %>%
  group_by(
    protein,
    window_start
  ) %>%
  summarise(
    
    n =
      sum(!is.na(dNdS)),
    
    median_dNdS =
      median(dNdS,
             na.rm = TRUE),
    
    mean_dNdS =
      mean(dNdS,
           na.rm = TRUE),
    
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
  mutate(
    
    q_positive =
      p.adjust(
        p_positive,
        method = "BH"
      ),
    
    q_negative =
      p.adjust(
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

################################################################################
# HOTSPOT SUMMARY
################################################################################

cat("\n")
cat("############################################################\n")
cat("HOTSPOT SUMMARY\n")
cat("############################################################\n")

print(
  selection_hotspots_w4s1.q20 %>%
    dplyr::count(hotspot)
)

cat("\n")
cat("############################################################\n")
cat("SIGNIFICANT WINDOWS\n")
cat("############################################################\n")

print(
  selection_hotspots_w4s1.q20 %>%
    filter(
      hotspot != "none"
    ) %>%
    arrange(
      protein,
      window_start
    )
)

################################################################################
# COMPARISON OF RAW COUNTS
################################################################################

cat("\n")
cat("############################################################\n")
cat("SYNONYMOUS VS NONSYNONYMOUS MUTATIONS\n")
cat("############################################################\n")

print(selection_table.q20 %>%
    summarise(
      synonymous =
        sum(
          synonymous_mutations
        ),
      nonsynonymous =
        sum(
          nonsynonymous_mutations
        )
    )
)


######

selec.A.q20 <-  ggplot(
  selection_sample_summary.q20 %>%
    mutate(
      protein = factor(
        protein,
        levels = c("RdRP", "CP", "delta")
      )
    ),
  aes(protein, mean_dNdS)
) +
  
  geom_hline(yintercept = 0) +
  
  geom_quasirandom(
    aes(
      size = hpi,
      fill = protein
    ),
    shape = 21,
    width = .15,
    alpha = .9
  ) +
  
  stat_summary(
    aes(color = protein),
    fun = median,
    shape = 21,
    fill = "white",
    size = 1
  ) +
  
  labs(
    x = NULL,
    y = expression("Mean " * (d[N] - d[S]))
  ) +
  
  scale_color_manual(values = color.proteins) +
  scale_fill_manual(values = color.proteins) +
  
  scale_size_continuous(
    range = c(1, 3)
  ) +
  
  ylim(c(-0.00765, 0.0025)) +
  
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = "transparent"),
    text = element_text(color = "black"),
    axis.text = element_text(
      size = 12,
      color = "black"
    ),
    strip.text = element_text(size = 14)
  )

selec.B.q20 <- selection_sample_summary.q20 %>%
  
  mutate(
    protein = factor(
      protein,
      levels = c(
        "RdRP",
        "CP",
        "delta"
      )
    )
  ) %>%
  
  arrange(hpi) %>%
  
  ggplot(
    aes(
      hpi,
      mean_dNdS,
      color = protein
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    color = "grey50"
  ) +
  
  geom_jitter(
    aes(
      fill = protein,
      size = hpi
    ),
    color = "black",
    shape = 21,
    width = 0.4,
    height = 0
  ) +
  
  geom_ribbon(
    data = newdat.q20,
    aes(
      x = hpi,
      ymin = lwr,
      ymax = upr,
      fill = protein
    ),
    alpha = .2,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = newdat.q20,
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
  
  scale_size_continuous(
    range = c(1, 3)
  ) +
  
  ylim(c(-0.00765, 0.0025)) +
  
  labs(
    x = "Time (hpi)",
    y = expression("Mean " * (d[N] - d[S]))
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

selection_w4s1_4plot.q20 <- selection_w4s1.q20 %>%
  
  mutate(
    dNdS = dN - dS
  ) %>%
  
  left_join(
    
    selection_hotspots_w4s1.q20 %>%
      
      mutate(
        significance = case_when(
          
          q_negative < 0.05 ~ "q_negative",
          
          q_positive < 0.05 ~ "q_positive",
          
          p_negative < 0.05 ~ "p_negative",
          
          p_positive < 0.05 ~ "p_positive",
          
          TRUE ~ "none"
        )
      ) %>%
      
      dplyr::select(
        protein,
        window_start,
        significance,
        p_positive,
        p_negative,
        q_positive,
        q_negative
      ),
    
    by = c(
      "protein",
      "window_start"
    )
  )

selection_w4s1_4hotspot_labels.q20 <- selection_w4s1_4plot.q20 %>%
  
  filter(
    p_positive < 0.05 |
      p_negative < 0.05
  ) %>%
  
  mutate(
    
    window_mid =
      window_start +
      (window_end - window_start) / 2,
    
    label =
      paste0(
        window_start,
        "-",
        window_end
      ),
    
    label_type = case_when(
      
      q_positive < 0.05 ~ "q_positive",
      
      q_negative < 0.05 ~ "q_negative",
      
      p_positive < 0.05 ~ "p_positive",
      
      p_negative < 0.05 ~ "p_negative"
    )
  ) %>%
  
  distinct(
    protein,
    window_start,
    .keep_all = TRUE
  )


################################################################################
# FIGURE 7C (Q20)
################################################################################

library(ggrepel)
library(patchwork)

################################################################################
# ADD WINDOW SIGNIFICANCE
################################################################################

selection_w4s1_4plot.q20 <-
  
  selection_w4s1.q20 %>%
  
  mutate(
    dNdS = dN - dS
  ) %>%
  
  left_join(
    
    selection_hotspots_w4s1.q20 %>%
      
      mutate(
        
        significance = case_when(
          
          q_negative < 0.05 ~ "q_negative",
          
          q_positive < 0.05 ~ "q_positive",
          
          p_negative < 0.05 ~ "p_negative",
          
          p_positive < 0.05 ~ "p_positive",
          
          TRUE ~ "none"
        )
      ) %>%
      
      select(
        protein,
        window_start,
        p_positive,
        p_negative,
        q_positive,
        q_negative,
        significance
      ),
    
    by = c(
      "protein",
      "window_start"
    )
  )

################################################################################
# LABELS
################################################################################

selection_w4s1_4hotspot_labels.q20 <-
  
  selection_w4s1_4plot.q20 %>%
  
  filter(
    p_positive < 0.05 |
      p_negative < 0.05
  ) %>%
  
  mutate(
    
    window_mid =
      window_start +
      (window_end - window_start) / 2,
    
    label =
      paste0(
        window_start,
        "-",
        window_end
      ),
    
    label_type = case_when(
      
      q_positive < 0.05 ~ "q_positive",
      
      q_negative < 0.05 ~ "q_negative",
      
      p_positive < 0.05 ~ "p_positive",
      
      p_negative < 0.05 ~ "p_negative"
    )
  ) %>%
  
  distinct(
    protein,
    window_start,
    .keep_all = TRUE
  )

################################################################################
# PANEL 1 (RdRP)
################################################################################

fig7C_top.q20 <-
  
  {
    selection_w4s1_4plot.q20 %>%
      
      left_join(
        samples_info %>%
          dplyr::select(
            sample,
            hpi,
            replica
          ),
        by = "sample"
      ) %>%
      
      rowwise() %>%
      
      mutate(
        
        window_mid =
          window_start +
          (window_end - window_start) / 2,
        
        `pN/pS` = pN / pS,
        
        `pN - pS` = pN - pS,
        
        `dN - dS` = dN - dS,
        
        `NS-S` =
          nonsynonymous_mutations -
          synonymous_mutations,
        
        protein = factor(
          protein,
          levels = c(
            "RdRP",
            "CP",
            "delta"
          )
        )
      ) %>%
      
      arrange(hpi) %>%
      
      filter(
        abs(`dN - dS`) > 0,
        protein == "RdRP"
      ) %>%
      
      ggplot(
        aes(
          x = window_mid,
          y = `dN - dS`,
          colour = protein
        )
      ) +
      
      geom_rect(
        
        data = palmprint.orv %>%
          
          mutate(
            protein = factor(
              "RdRP",
              levels = c(
                "RdRP",
                "CP",
                "delta"
              )
            )
          ),
        
        aes(
          xmin = aa.start,
          xmax = aa.end,
          ymin = -Inf,
          ymax = Inf
        ),
        
        inherit.aes = FALSE,
        
        alpha = .10,
        
        linewidth = .05,
        
        colour = "black"
      ) +
      
      geom_hline(
        yintercept = 0,
        linewidth = .2
      ) +
      
      geom_point(
        aes(
          colour = significance,
          alpha = significance,
          size = significance
        )
      ) +
      
      geom_label_repel(
        
        data =
          selection_w4s1_4hotspot_labels.q20 %>%
          filter(
            protein == "RdRP",
            label_type %in%
              c(
                "p_positive",
                "p_negative"
              )
          ),
        
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
      
      geom_label_repel(
        
        data =
          selection_w4s1_4hotspot_labels.q20 %>%
          filter(
            protein == "RdRP",
            label_type %in%
              c(
                "q_positive",
                "q_negative"
              )
          ),
        
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
      
      facet_grid(
        ~ protein,
        scales = "free_x",
        space = "free_x"
      ) +
      
      scale_x_continuous(
        expand = c(0,0)
      ) +
      
      coord_cartesian(
        clip = "off"
      ) +
      
      scale_colour_manual(
        values = c(
          none = "black",
          p_negative = "#ff1493",
          q_negative = "#ff1493",
          p_positive = "#6b8e23",
          q_positive = "#6b8e23"
        )
      ) +
      
      scale_alpha_manual(
        values = c(
          none = .5,
          p_negative = .5,
          q_negative = 1,
          p_positive = .5,
          q_positive = 1
        )
      ) +
      
      scale_size_manual(
        values = c(
          none = .5,
          p_negative = .5,
          q_negative = .8,
          p_positive = .5,
          q_positive = .8
        )
      ) +
      
      scale_fill_manual(
        values = c(
          q_negative = "#ff1493",
          q_positive = "#6b8e23"
        )
      ) +
      
      labs(x = "") +
      
      theme_classic() +
      
      theme(
        legend.position = "none",
        panel.border = element_rect(
          fill = "transparent"
        ),
        text = element_text(
          colour = "black"
        ),
        axis.text = element_text(
          size = 12,
          colour = "black"
        ),
        strip.text = element_text(
          size = 14
        )
      )
  }

################################################################################
# PANEL 2 (CP + delta)
################################################################################

fig7C_bottom.q20 <-
  
  {
    selection_w4s1_4plot.q20 %>%
      
      left_join(
        samples_info %>%
          dplyr::select(
            sample,
            hpi,
            replica
          ),
        by = "sample"
      ) %>%
      
      rowwise() %>%
      
      mutate(
        
        window_mid =
          window_start +
          (window_end - window_start) / 2,
        
        `pN/pS` = pN / pS,
        
        `pN - pS` = pN - pS,
        
        `dN - dS` = dN - dS,
        
        `NS-S` =
          nonsynonymous_mutations -
          synonymous_mutations,
        
        protein = factor(
          protein,
          levels = c(
            "RdRP",
            "CP",
            "delta"
          )
        )
      ) %>%
      
      arrange(hpi) %>%
      
      filter(
        abs(`pN - pS`) > 0,
        protein != "RdRP"
      ) %>%
      
      ggplot(
        aes(
          x = window_mid,
          y = `dN - dS`,
          colour = protein
        )
      ) +
      
      geom_hline(
        yintercept = 0,
        linewidth = .2
      ) +
      
      geom_point(
        aes(
          colour = significance,
          alpha = significance,
          size = significance
        )
      ) +
      
      geom_label_repel(
        
        data =
          selection_w4s1_4hotspot_labels.q20 %>%
          filter(
            protein != "RdRP",
            label_type %in%
              c(
                "p_positive",
                "p_negative"
              )
          ),
        
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
      
      geom_label_repel(
        
        data =
          selection_w4s1_4hotspot_labels.q20 %>%
          filter(
            protein != "RdRP",
            label_type %in%
              c(
                "q_positive",
                "q_negative"
              )
          ),
        
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
      
      facet_grid(
        ~ protein,
        scales = "free_x",
        space = "free_x"
      ) +
      
      scale_x_continuous(
        expand = c(0,0)
      ) +
      
      coord_cartesian(
        clip = "off"
      ) +
      
      scale_colour_manual(
        values = c(
          none = "black",
          p_negative = "#ff1493",
          q_negative = "#ff1493",
          p_positive = "#6b8e23",
          q_positive = "#6b8e23"
        )
      ) +
      
      scale_alpha_manual(
        values = c(
          none = .5,
          p_negative = .5,
          q_negative = 1,
          p_positive = .5,
          q_positive = 1
        )
      ) +
      
      scale_size_manual(
        values = c(
          none = .5,
          p_negative = .5,
          q_negative = .8,
          p_positive = .5,
          q_positive = .8
        )
      ) +
      
      scale_fill_manual(
        values = c(
          q_negative = "#ff1493",
          q_positive = "#6b8e23"
        )
      ) +
      
      labs(
        x = "Codon"
      ) +
      
      theme_classic() +
      
      theme(
        legend.position = "none",
        panel.border = element_rect(
          fill = "transparent"
        ),
        text = element_text(
          colour = "black"
        ),
        axis.text = element_text(
          size = 12,
          colour = "black"
        ),
        strip.text = element_text(
          size = 14
        )
      )
  }



fig7C.q20 <-
  
  fig7C_top.q20 /
  fig7C_bottom.q20

################################################################################
# FINAL FIGURE
################################################################################

(
  selec.A.q20 +
    selec.B.q20 +
    plot_layout(widths = c(1.5, 2))
) /
  fig7C.q20 +
  plot_layout(
    heights = c(0.75, 1.5)
  )



################################################################################

all.snv.byStrand.bed.q20 %>%
  filter(
    strand == "+",
    id_snv_chrm %in% c(
      "RNA2:631U>C",
      "RNA2:718C>U"
    )
  ) %>%
  group_by(id_snv_chrm) %>%
  summarise(
    n_samples_total = n_distinct(sample),
    n_samples_AF001 = n_distinct(sample[AF > 0.01]),
    .groups = "drop"
  )

