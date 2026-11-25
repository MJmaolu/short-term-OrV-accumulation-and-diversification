#===============================================================================
# MJ OLMO-UCEDA
# 2026/11/09
#===============================================================================
# KDE ENRICHMENT ANALYSIS v2
#
# Improvements over original version:
#   1. Weighted KDE maintained (weights = SNV recurrence among samples)
#   2. Permutation-based null model maintained
#   3. Genome-wide significance threshold added
#      (controls for multiple testing across positions)
#   4. Local envelope retained for visualization
#   5. Global significant hotspots identified
#   6. No existing objects overwritten (_v2.q20 suffix)
#===============================================================================

#===============================================================================
# 1. INPUT DATA
#===============================================================================

# Expected columns in snv_events.q20:
#
# POS
# sample_count
# rna
# category
#
# category = plus / minus / both

#===============================================================================
# 2. WEIGHTED KDE
#===============================================================================

weighted_kde_v2 <- function(
    pos,
    weights,
    genome_length,
    bw = 100,
    n = 512
){
  
  density(
    pos,
    weights = weights / sum(weights),
    bw = bw,
    from = 1,
    to = genome_length,
    n = n
  )
  
}

#===============================================================================
# 3. PERMUTATION NULL MODEL
#===============================================================================

permute_kde_v2 <- function(
    pos,
    weights,
    genome_length,
    bw = 100,
    nperm = 5000,
    ngrid = 512
){
  
  grid <- seq(
    1,
    genome_length,
    length.out = ngrid
  )
  
  perm_dens <- matrix(
    NA_real_,
    nrow = nperm,
    ncol = ngrid
  )
  
  for(i in seq_len(nperm)){
    
    pos_perm <- sample(
      seq_len(genome_length),
      size = length(pos),
      replace = FALSE
    )
    
    d <- density(
      pos_perm,
      weights = weights / sum(weights),
      bw = bw,
      from = 1,
      to = genome_length,
      n = ngrid
    )
    
    perm_dens[i, ] <- d$y
    
  }
  
  #---------------------------------------------------------------------------
  # Local envelope
  #---------------------------------------------------------------------------
  
  low_local <- apply(
    perm_dens,
    2,
    quantile,
    probs = 0.025
  )
  
  high_local <- apply(
    perm_dens,
    2,
    quantile,
    probs = 0.975
  )
  
  #---------------------------------------------------------------------------
  # Global permutation threshold
  #
  # Uses maximum KDE value from every permutation.
  # Controls family-wise error rate across the genome.
  #---------------------------------------------------------------------------
  
  max_perm <- apply(
    perm_dens,
    1,
    max
  )
  
  global_threshold <- quantile(
    max_perm,
    probs = 0.95
  )
  
  list(
    
    grid = grid,
    
    mean_expected = colMeans(perm_dens),
    
    low_local = low_local,
    
    high_local = high_local,
    
    global_threshold = global_threshold,
    
    perm_dens = perm_dens,
    
    max_perm = max_perm
    
  )
  
}

#===============================================================================
# 4. RUN COMPLETE KDE ENRICHMENT
#===============================================================================

run_kde_enrichment_v2 <- function(
    df,
    genome_length,
    bw = 100,
    nperm = 5000
){
  
  obs <- weighted_kde_v2(
    pos = df$POS,
    weights = df$sample_count,
    genome_length = genome_length,
    bw = bw
  )
  
  null <- permute_kde_v2(
    pos = df$POS,
    weights = df$sample_count,
    genome_length = genome_length,
    bw = bw,
    nperm = nperm,
    ngrid = length(obs$x)
  )
  
  data.frame(
    
    x = obs$x,
    
    observed = obs$y,
    
    expected = null$mean_expected,
    
    low_local = null$low_local,
    
    high_local = null$high_local,
    
    global_threshold = null$global_threshold
    
  )
  
}

#===============================================================================
# 5. RUN ANALYSIS
#===============================================================================

res_rna1_plus_v2.q20 <- run_kde_enrichment_v2(
  snv_events.q20 %>%
    filter(
      rna == "RNA 1",
      category == "plus"
    ),
  genome_length = 3421
)

res_rna1_minus_v2.q20 <- run_kde_enrichment_v2(
  snv_events.q20 %>%
    filter(
      rna == "RNA 1",
      category == "minus"
    ),
  genome_length = 3421
)

res_rna2_plus_v2.q20 <- run_kde_enrichment_v2(
  snv_events.q20 %>%
    filter(
      rna == "RNA 2",
      category == "plus"
    ),
  genome_length = 2574
)

res_rna2_minus_v2.q20 <- run_kde_enrichment_v2(
  snv_events.q20 %>%
    filter(
      rna == "RNA 2",
      category == "minus"
    ),
  genome_length = 2574
)

#===============================================================================
# 6. HOTSPOT CALLING
#
# Genome-wide significant hotspots
#===============================================================================

find_hotspots_v2 <- function(res){
  
  res$enriched_global <-
    res$observed > res$global_threshold
  
  rr <- rle(
    res$enriched_global
  )
  
  ends <- cumsum(rr$lengths)
  
  starts <- c(
    1,
    head(ends, -1) + 1
  )
  
  data.frame(
    enriched = rr$values,
    start = res$x[starts],
    end = res$x[ends]
  ) %>%
    filter(enriched)
  
}

hotspots_rna1_plus_v2.q20 <-
  find_hotspots_v2(
    res_rna1_plus_v2.q20
  )

hotspots_rna1_minus_v2.q20 <-
  find_hotspots_v2(
    res_rna1_minus_v2.q20
  )

hotspots_rna2_plus_v2.q20 <-
  find_hotspots_v2(
    res_rna2_plus_v2.q20
  )

hotspots_rna2_minus_v2.q20 <-
  find_hotspots_v2(
    res_rna2_minus_v2.q20
  )

#===============================================================================
# 7. COMBINED DATASET
#===============================================================================

kde_plot_v2.q20 <- bind_rows(
  
  res_rna1_plus_v2.q20 %>%
    mutate(
      rna = "RNA 1",
      strand = "+"
    ),
  
  res_rna1_minus_v2.q20 %>%
    mutate(
      rna = "RNA 1",
      strand = "-"
    ),
  
  res_rna2_plus_v2.q20 %>%
    mutate(
      rna = "RNA 2",
      strand = "+"
    ),
  
  res_rna2_minus_v2.q20 %>%
    mutate(
      rna = "RNA 2",
      strand = "-"
    )
  
) %>%
  
  mutate(
    
    enriched_local =
      observed > high_local,
    
    enriched_global =
      observed > global_threshold
    
  )

#===============================================================================
# 8. OPTIONAL QC PLOT
#
# Grey ribbon = local 95% envelope
# Red dashed = genome-wide threshold
#===============================================================================

ggplot(
  res_rna1_plus_v2.q20,
  aes(x = x)
) +
  
  geom_ribbon(
    aes(
      ymin = low_local,
      ymax = high_local
    ),
    fill = "grey80",
    alpha = 0.4
  ) +
  
  geom_hline(
    yintercept =
      unique(
        res_rna1_plus_v2.q20$global_threshold
      ),
    colour = "red",
    linetype = "dashed"
  ) +
  
  geom_line(
    aes(y = observed),
    linewidth = 0.5
  ) +
  
  theme_classic()

#===============================================================================
# 9. RESCALE KDE FOR OVERLAY WITH FIGURE 5
#===============================================================================

max_snv.q20 <- 13

max_kde_v2.q20 <- max(
  kde_plot_v2.q20$observed,
  na.rm = TRUE
)

kde_scale_v2.q20 <- max_snv.q20 / max_kde_v2.q20

kde_plot2_v2.q20 <- kde_plot_v2.q20 %>%
  
  mutate(
    
    observed_s =
      observed * kde_scale_v2.q20,
    
    expected_s =
      expected * kde_scale_v2.q20,
    
    high_local_s =
      high_local * kde_scale_v2.q20,
    
    low_local_s =
      low_local * kde_scale_v2.q20,
    
    global_threshold_s =
      global_threshold * kde_scale_v2.q20
    
  )

#===============================================================================
# 10. KDE LAYERS FOR FIGURE 5
#
# Replace your previous KDE layer by this:
#===============================================================================

geom_line(
  data = kde_plot2_v2.q20,
  aes(
    x = x,
    y = ifelse(
      strand == "+",
      observed_s,
      -observed_s
    ),
    color = strand,
    group = interaction(rna, strand)
  ),
  linewidth = 0.4,
  inherit.aes = FALSE
)

geom_line(
  data = kde_plot2_v2.q20,
  aes(
    x = x,
    y = ifelse(
      strand == "+",
      global_threshold_s,
      -global_threshold_s
    ),
    color = strand,
    group = interaction(rna, strand)
  ),
  linewidth = 0.25,
  linetype = "dashed",
  alpha = 0.6,
  inherit.aes = FALSE
)

#===============================================================================
# 11. SUMMARY TABLES
#===============================================================================

hotspot_summary_v2.q20 <- bind_rows(
  
  hotspots_rna1_plus_v2.q20 %>%
    mutate(
      rna = "RNA 1",
      strand = "+"
    ),
  
  hotspots_rna1_minus_v2.q20 %>%
    mutate(
      rna = "RNA 1",
      strand = "-"
    ),
  
  hotspots_rna2_plus_v2.q20 %>%
    mutate(
      rna = "RNA 2",
      strand = "+"
    ),
  
  hotspots_rna2_minus_v2.q20 %>%
    mutate(
      rna = "RNA 2",
      strand = "-"
    )
  
)

hotspot_summary_v2.q20


#==============================================================================
# FINAL FIGURE
#   - SNV recurrence
#   - Global KDE enrichment analysis
#   - Genome-wide permutation threshold
#   - Significant hotspot regions
#   - Protein annotation
#   - Palmprint motifs
#==============================================================================

hotspot_rects_v2.q20 <- hotspot_summary_v2.q20 %>%
  mutate(
    ymin = -Inf,
    ymax = Inf
  )


ggplot() +
  
  #--------------------------------------------------------------------------
# zero line
#--------------------------------------------------------------------------

geom_hline(
  yintercept = 0
) +
  
  #--------------------------------------------------------------------------
# significant hotspot regions (global test)
#--------------------------------------------------------------------------

# geom_rect(
#   data = hotspot_rects_v2.q20,
#   aes(
#     xmin = start,
#     xmax = end,
#     ymin = ymin,
#     ymax = ymax
#   ),
#   fill = "red",
#   alpha = 0.08,
#   inherit.aes = FALSE
# ) +
  
  #--------------------------------------------------------------------------
# proteins
#--------------------------------------------------------------------------

geom_rect(
  data = orv.proteins %>% filter(predicted == "No") ,
  aes(
    xmin = start,
    xmax = end,
    ymin = -Inf,
    ymax = Inf,
    fill = protein
  ),
  alpha = 0.10,
  inherit.aes = FALSE
) +
  
  scale_fill_manual(values = color.proteins) +
  
  ggnewscale::new_scale_fill() +
  
  #--------------------------------------------------------------------------
# segment borders
#--------------------------------------------------------------------------

geom_rect(
  data = data.frame(
    rna = c("RNA 1", "RNA 2"),
    start = c(1, 1),
    end = c(3421, 2574)
  ),
  aes(
    xmin = start,
    xmax = end,
    ymin = -Inf,
    ymax = Inf,
    color = as.factor(rna)
  ),
  fill = "transparent",
  linewidth = 0.1,
  inherit.aes = FALSE
) +
  
  scale_color_manual(values = color.rna) +
  
  ggnewscale::new_scale_color() +
  
  #--------------------------------------------------------------------------
# KDE observed
#--------------------------------------------------------------------------

geom_line(
  data = kde_plot2_v2.q20,
  aes(
    x = x,
    y = ifelse(
      strand == "+",
      observed_s,
      -observed_s
    ),
    color = strand,
    group = interaction(rna, strand)
  ),
  linewidth = 0.45,
  inherit.aes = FALSE
) +
  
  scale_color_manual(values = color.strand) +
  
  ggnewscale::new_scale_color() +
  
  
  #--------------------------------------------------------------------------
# GLOBAL permutation threshold
#--------------------------------------------------------------------------

geom_line(
  data = kde_plot2_v2.q20,
  aes(
    x = x,
    y = ifelse(
      strand == "+",
      global_threshold_s,
      -global_threshold_s
    ),
    color = strand,
    group = interaction(rna, strand)
  ),
  linetype = "dashed",
  linewidth = 0.35,
  alpha = 0.8,
  inherit.aes = FALSE
) +
  
  scale_color_manual(values = color.strand) +
  
  ggnewscale::new_scale_color() +
  
  #--------------------------------------------------------------------------
# significant KDE regions highlighted
#--------------------------------------------------------------------------

geom_line(
  data = kde_plot2_v2.q20 %>%
    filter(enriched_global),
  aes(
    x = x,
    y = ifelse(
      strand == "+",
      observed_s,
      -observed_s
    ),
    group = interaction(rna, strand),
    color = strand
  ),
  linewidth = 0.9,
  inherit.aes = FALSE
) +
  
  scale_color_manual(
    values = color.strand
  ) +

  ggnewscale::new_scale_color() +
  # 
  #--------------------------------------------------------------------------
# palmprint motifs
#--------------------------------------------------------------------------

geom_rect(
  data = palmprint.orv,
  aes(
    xmin = nt.start,
    xmax = nt.end,
    ymin = -Inf,
    ymax = Inf
  ),
  color = "black",
  alpha = 0.1,
  linewidth = 0.05,
  inherit.aes = FALSE
) +
  
  #==========================================================================
# SNVs ONLY IN PLUS
#==========================================================================

geom_segment(
  data = all.snv.byStrand.bed.q20 %>%
    filter(
      AF > af,
      id_snv_chrm %in% c(
        unique2Plus.rna1.q20,
        unique2Plus.rna2.q20
      )
    ) %>%
    group_by(id_snv_chrm) %>%
    mutate(
      sample_count = n_distinct(sample)
    ) %>%
    distinct(id_snv, .keep_all = TRUE),
  
  aes(
    x = POS,
    xend = POS,
    y = sample_count,
    yend = 0,
    color = strand
  ),
  
  linewidth = 0.1
) +
  
  geom_point(
    data = all.snv.byStrand.bed.q20 %>%
      filter(
        AF > af,
        id_snv_chrm %in% c(
          unique2Plus.rna1.q20,
          unique2Plus.rna2.q20
        )
      ) %>%
      group_by(id_snv_chrm) %>%
      mutate(
        sample_count = n_distinct(sample)
      ) %>%
      distinct(id_snv, .keep_all = TRUE),
    
    aes(
      x = POS,
      y = sample_count,
      shape = aa_change,
      color = strand
    ),
    
    size = 2
  ) +
  
  geom_text(
    data = all.snv.byStrand.bed.q20 %>%
      filter(
        AF > af,
        id_snv_chrm %in% c(
          unique2Plus.rna1.q20,
          unique2Plus.rna2.q20
        )
      ) %>%
      group_by(id_snv_chrm) %>%
      mutate(
        sample_count = n_distinct(sample)
      ) %>%
      distinct(id_snv, .keep_all = TRUE) %>%
      filter(sample_count > 1),
    
    aes(
      x = POS,
      y = 0.6 * sample_count,
      label = ifelse(
        aa_change == "non-synonymous",
        paste(SNV, "(", sample_count, ")", sep = ""),
        paste(id_snv, "(", sample_count, ")", sep = "")
      ),
      color = strand
    ),
    
    size = 3,
    hjust = 0,
    angle = 45,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(values = color.strand) +
  ggnewscale::new_scale_color() +
  
  #==========================================================================
# SNVs ONLY IN MINUS
#==========================================================================

geom_segment(
  data = all.snv.byStrand.bed.q20 %>%
    filter(
      AF > af,
      id_snv_chrm %in% c(
        unique2minus.rna1.q20,
        unique2minus.rna2.q20
      )
    ) %>%
    group_by(id_snv_chrm) %>%
    mutate(
      sample_count = n_distinct(sample)
    ) %>%
    distinct(id_snv, .keep_all = TRUE),
  
  aes(
    x = POS,
    xend = POS,
    y = -sample_count,
    yend = 0,
    color = strand
  ),
  
  linewidth = 0.1
) +
  
  geom_point(
    data = all.snv.byStrand.bed.q20 %>%
      filter(
        AF > af,
        id_snv_chrm %in% c(
          unique2minus.rna1.q20,
          unique2minus.rna2.q20
        )
      ) %>%
      group_by(id_snv_chrm) %>%
      mutate(
        sample_count = n_distinct(sample)
      ) %>%
      distinct(id_snv, .keep_all = TRUE),
    
    aes(
      x = POS,
      y = -sample_count,
      shape = aa_change,
      color = strand
    ),
    
    size = 2
  ) +
  
  geom_text(
    data = all.snv.byStrand.bed.q20 %>%
      filter(
        AF > af,
        id_snv_chrm %in% c(
          unique2minus.rna1.q20,
          unique2minus.rna2.q20
        )
      ) %>%
      group_by(id_snv_chrm) %>%
      mutate(
        sample_count = n_distinct(sample)
      ) %>%
      distinct(id_snv, .keep_all = TRUE) %>%
      filter(sample_count > 1),
    
    aes(
      x = POS,
      y = -0.8 * sample_count,
      label = paste(
        id_snv,
        "(",
        sample_count,
        ")",
        sep = ""
      ),
      color = strand
    ),
    
    size = 3,
    hjust = 1,
    angle = 45,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(values = color.strand) +
  ggnewscale::new_scale_color() +
  
  #==========================================================================
# SNVs FOUND IN BOTH STRANDS
#==========================================================================

geom_segment(
  data = all.snv.byStrand.bed.q20 %>%
    filter(
      AF > af,
      id_snv_chrm %in% c(
        realBoth.rna1.q20,
        realBoth.rna2.q20
      )
    ),
  
  aes(
    x = POS,
    xend = POS,
    y = -0.2,
    yend = 0.2
  ),
  
  color = "#803058",
  linewidth = 0.2,
  alpha = 0.7
) +
  
  #==========================================================================
# FACETS
#==========================================================================

facet_grid(
  rna ~ .,
  scales = "free_x",
  space = "free_x"
) +
  
  scale_shape_manual(
    values = c(
      "synonymous" = 20,
      "non-synonymous" = 18,
      "NA" = 1
    )
  ) +
  
  scale_y_continuous(
    name = "SNV recurrence (samples)",
    
    sec.axis = sec_axis(
      ~ . / kde_scale_v2.q20,
      name = "Weighted SNV density\n(global permutation test)"
    )
  ) +
  
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  
  labs(
    x = "Genomic position (nt)",
    y = NULL
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "bottom",
    
    axis.line.x = element_line(colour = "black"),
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    
    panel.background = element_blank(),
    
    strip.text = element_blank(),
    
    panel.border = element_rect(
      colour = "black",
      fill = "transparent"
    ),
    
    panel.spacing.y = unit(
      1.2,
      "lines"
    )
  )
