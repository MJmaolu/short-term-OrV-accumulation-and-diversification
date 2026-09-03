# Fig 6A: AF CORRELATION
color.proteins 
all.snv.byStrand.bed_wider %>%
  filter(`AF_+` > .01 & `AF_-` > 0.01) %>%
  ggplot(.,
         aes(`AF_-`,
             `AF_+`,
             fill = region)) +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dotted",
              linewidth = .2) +
  geom_smooth(method = "lm",
              linewidth = .6,
              aes(color = region)) +
  
  geom_point(shape = 21,
             size = 2,
             alpha = .8) +
  
  geom_text_repel(
    data = ~ .x %>%
      filter(`AF_+` > 0.2 | `AF_-` > 0.2) %>%
      mutate(
        label = if_else(
          aa_change == "synonymous",
          id_snv,
          SNV
        )
      ),
    aes(label = label,
        color = region),
    size = 3,
    max.overlaps = Inf
  ) +
 
  scale_fill_manual(values = c(color.proteins,
                              "5'UTR" = "gray30",
                              "3'UTR" = "white")) +
  scale_color_manual(values = color.proteins) +
  
  facet_grid(rna ~ aa_change) +
  labs(x = "AF(-)",
       y = "AF(+)") +
  theme_classic() +
  coord_equal() +
  lims(x = c(0,1),
       y = c(0,1)) +
  theme(legend.position = "bottom",
        axis.line.x = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 10),
        panel.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        # panel.spacing.y = unit(1.2, 
        #                        "lines")
  )

## ADD n = number of SNVs and N = number of distinct SNVs
all.snv.byStrand.bed_wider %>%
  filter(`AF_+` > .01 & `AF_-` > .01) %>%
  # # Filter used in the previous version
  # filter(`AF_+` >= .01 | `AF_-` >= 0.01, 
  #        present == "both") %>%
  
  group_by(rna, region, aa_change) %>%
  summarise(N = n(),
            nSNVs = n_distinct(id_snv_chrm))

### RECURRENCE PATTERN
all.snv.byStrand.bed_wider %>%
  filter(aa_change %in% c("synonymous", "non-synonymous")) %>%
  group_by(region, aa_change) %>%
  summarise(
    #recurrent = prop(present == "recurrent", na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# FIG 6B)
# DISTRIBUTION OF THE RECURRENT SNVS (ONLY BOTH)
# with labels and space adapted to the distribution
all.snv.byStrand.bed_wider %>%
  filter(
    `AF_+` > 0.01,
    `AF_-` > 0.01,
    region %in% c("RdRP", "CP", "delta P")
  ) %>%
  
  # 1. recurrencia por SNV
  group_by(region, aa_change, id_snv_chrm) %>%
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
  
  facet_grid(~region, scales = "free_x", space = "free_x") +
  
  scale_y_continuous(
    breaks = seq(0, 50, by = 5),
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

# DISTRIBUTION OF THE RECURRENT SNVS, but inherits the sample_count column that
# doesn't filter by present == "both". So it is the SNVs that are in some samples
# in both but also counted when they appear only in one of the strands in some samples
# with labels and space adapted to the distribution
all.snv.byStrand.bed_wider %>%
  filter(
    `AF_+` > 0.01,
    `AF_-` > 0.01,
    region %in% c("RdRP", "CP", "delta P")
  ) %>%
  group_by(region, aa_change, sample_count) %>%
  summarise(
    n = n_distinct(id_snv_chrm),
    .groups = "drop"
  ) %>%
  group_by(region, aa_change) %>%
  complete(
    sample_count = full_seq(sample_count, 1),
    fill = list(n = 0)
  ) %>%
  ungroup() %>%
  mutate(
    region = factor(
      region,
      levels = c("5'UTR", "RdRP", "CP", "delta P", "3'UTR")
    ),
    sample_count = factor(sample_count)
    ) %>% 
  
  ggplot(
    aes(
      group = region,
      fill = region,
      linewidth = aa_change,
      linetype = aa_change,
      color = aa_change
    )
  ) +
  
  geom_col(
    aes(
      x = sample_count,
      y = n
    ),
    position = "dodge2"
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
    data = \(x) filter(x, n > 0),
    aes(
      x = sample_count,
      y = n,
      label = n,
      fill = aa_change
    ),
    position = position_stack(vjust = 1.05),
    color = "black",
    size = 3,
    label.size = 0.25,
    show.legend = FALSE,
    inherit.aes = FALSE
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
  
  facet_grid(~ region,
             scales = "free_x",
             space = "free_x") +
  
  # coord_cartesian(
  #   ylim = c(0, 38)
  # ) +
  # 
  # scale_x_continuous(
  #   breaks = seq(1, 30, by = 2)
  # ) +
  # 
  scale_y_continuous(
    breaks = seq(0, 50, by = 5),
    expand = c(0, 0)
  ) +
  labs(x = "Samples recurrence",
       y = "Number of SNVs \n(both strands AF > 0.01)") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.line.x = element_line(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    panel.background = element_blank(),
    strip.text = element_text(size = 12),
    panel.border = element_rect(
      colour = "black",
      fill = "transparent"
    )
  )
  
# FIG 6C)
### BOTH AND RECURRENT (n_samples > 3)

all.snv.af0.01.AFmat[rev(rownames(all.snv.af0.01.binMat.df[all.snv.af0.01.binMat.df$N > 3,])),
                     samples_info$sample][(all.snv.af0.01.AFmat[rev(rownames(all.snv.af0.01.binMat.df[all.snv.af0.01.binMat.df$N > 3,])),
                     samples_info$sample] %>% rownames()) %in% (all.snv.byStrand.bed_wider %>%
  filter(#present == "both",
         `AF_+` > .01 & `AF_-` > 0.01
         ) %>% pull(id_snv_chrm)),]

# FIG 6C
# WITH values FROM (+): But only SNVs that appear in both strands at least in 4 samples
mat.both.af0.01.pAF.new <- both.snv.af0.01.pAFMat[
  rev(rownames(both.snv.af0.01.binMat[
    both.snv.af0.01.binMat$N > 3,
  ])),
  samples_info$sample
]

both.snv.up3samples_annot <- mat.both.af0.01.pAF.new %>%
  as.data.frame() %>%
  rownames_to_column(var = "id_snv_chrm") %>%
  left_join(
    all.snv %>%
      dplyr::select(id_snv_chrm, rna, region, aa_change) %>%
      distinct(),
    by = "id_snv_chrm"
  ) %>%
  dplyr::select(id_snv_chrm, rna, region, aa_change)

# # SHOWING ONLY THE AF+ OF THE SAMPLES WHERE THE SNV APPEAR IN BOTH STRANDS
# reshape2::melt(as.matrix(mat.both.af0.01.pAF.new)) %>%
#   rowwise() %>%
#   mutate(
#     timepoint = strsplit(as.character(Var2), "_")[[1]][2],
#     hpi = as.numeric(str_replace(timepoint, "h", "")),
#     replica = strsplit(as.character(Var2), "_")[[1]][3],
#     sample = paste0(timepoint, " r", replica),
#     sample = factor(
#       sample,
#       levels = c(
#         "0h r1", "0h r2", "0h r3",
#         "2h r1", "2h r2", "2h r3",
#         "4h r1", "4h r2", "4h r3",
#         "6h r1", "6h r2", "6h r3",
#         "8h r1", "8h r2", "8h r3",
#         "12h r1", "12h r2", "12h r3",
#         "18h r1", "18h r2", "18h r3",
#         "22h r1", "22h r2", "22h r3",
#         "26h r1", "26h r2", "26h r3",
#         "32h r1", "32h r2", "32h r3",
#         "38h r1", "38h r2", "38h r3",
#         "44h r1", "44h r2", "44h r3"
#       )
#     )
#   ) %>%
#   arrange(hpi) %>%
#   merge(
#     all.snv[, c(
#       "rna", "id_snv", "id_snv_chrm", "cds_region",
#       "aa_change", "aa_ref", "aa_alt", "SNV", "region"
#     )],
#     by.x = "Var1",
#     by.y = "id_snv_chrm"
#   ) %>%
#   distinct() %>%
#   
#   ggplot(aes(
#     x = sample,
#     y = Var1,
#     fill = value
#   )) +
#   
#   geom_tile(
#     color = "black",
#     linewidth = 0.1
#   ) +
#   
#   geom_text(
#     data = . %>% filter(value >= 0.1),
#     aes(
#       x = sample,
#       y = Var1,
#       label = sprintf("%.2f", value)
#     ),
#     size = 2.5,
#     color = "white"
#   ) +
#   
#   scale_fill_gradient(
#     low = "white",
#     high = "black"
#   ) +
#   
#   ggnewscale::new_scale_fill() +
#   
#   geom_tile(
#     data = both.snv.up3samples_annot,
#     aes(
#       x = "rna",
#       y = id_snv_chrm,
#       fill = rna
#     ),
#     color = "black",
#     linewidth = 0.1,
#     inherit.aes = FALSE
#   ) +
#   
#   geom_tile(
#     data = both.snv.up3samples_annot,
#     aes(
#       x = "region",
#       y = id_snv_chrm,
#       fill = region
#     ),
#     color = "black",
#     linewidth = 0.1,
#     inherit.aes = FALSE
#   ) +
#   
#   scale_fill_manual(
#     values = c(color.proteins, color.rna),
#     name = "Region"
#   ) +
#   
#   theme_minimal() +
#   
#   labs(
#     x = "Samples",
#     y = "",
#     fill = "AF"
#   ) +
#   
#   coord_equal() +
#   
#   theme(
#     legend.position = "none",
#     axis.text.x = element_text(
#       size = 8,
#       angle = 90,
#       hjust = 1,
#       color = "black"
#     ),
#     axis.text.y = element_text(
#       face = ifelse(
#         both.snv.up3samples_annot[match(both.snv.up3samples_annot$id_snv_chrm, rownames(mat.both.af0.01.pAF.new)), "aa_change"] == "non-synonymous",
#         "bold",
#         "italic"
#       ),
#       color = "black"
#     )
#   )

# SHOWING ALSO THE AF+ OF THE SAMPLES WHERE THE SNV APPEAR ONLY IN THE + STRANDS
reshape2::melt(as.matrix(both.snv.af0.01.pAFMat_inAllPlus)) %>%
  rowwise() %>%
  mutate(
    timepoint = strsplit(as.character(Var2), "_")[[1]][2],
    hpi = as.numeric(str_replace(timepoint, "h", "")),
    replica = strsplit(as.character(Var2), "_")[[1]][3],
    sample = paste0(timepoint, " r", replica),
    sample = factor(
      sample,
      levels = c(
        "0h r1", "0h r2", "0h r3",
        "2h r1", "2h r2", "2h r3",
        "4h r1", "4h r2", "4h r3",
        "6h r1", "6h r2", "6h r3",
        "8h r1", "8h r2", "8h r3",
        "12h r1", "12h r2", "12h r3",
        "18h r1", "18h r2", "18h r3",
        "22h r1", "22h r2", "22h r3",
        "26h r1", "26h r2", "26h r3",
        "32h r1", "32h r2", "32h r3",
        "38h r1", "38h r2", "38h r3",
        "44h r1", "44h r2", "44h r3"
      )
    )
  ) %>%
  arrange(hpi) %>%
  merge(
    all.snv[, c(
      "rna", "id_snv", "id_snv_chrm", "cds_region",
      "aa_change", "aa_ref", "aa_alt", "SNV", "region"
    )],
    by.x = "Var1",
    by.y = "id_snv_chrm"
  ) %>%
  distinct() %>%
  
  ggplot(aes(
    x = sample,
    y = Var1,
    fill = value
  )) +
  
  geom_tile(
    color = "black",
    linewidth = 0.1
  ) +
  
  geom_text(
    data = . %>% filter(value >= 0.1),
    aes(
      x = sample,
      y = Var1,
      label = sprintf("%.2f", value)
    ),
    size = 2.5,
    color = "white"
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "black"
  ) +
  
  ggnewscale::new_scale_fill() +
  
  geom_tile(
    data = both.snv.up3samples_annot,
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
    data = both.snv.up3samples_annot,
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
    values = c(color.proteins, color.rna),
    name = "Region"
  ) +
  
  theme_minimal() +
  
  labs(
    x = "Samples",
    y = "",
    fill = "AF"
  ) +
  
  coord_equal() +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      size = 8,
      angle = 90,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      face = ifelse(
        both.snv.up3samples_annot[match(both.snv.up3samples_annot$id_snv_chrm, rownames(both.snv.af0.01.pAFMat_inAllPlus)), "aa_change"] == "non-synonymous",
        "bold",
        "italic"
      ),
      color = "black"
    )
  )

# Number of samples where appear the selected SNVs
################################################################################
## BOTH
all.snv.byStrand.bed_wider %>%
  filter(`AF_+` > 0.01 & `AF_-` > 0.01) %>%
  dplyr::select(c(id_snv_chrm, sample, `AF_+`, `AF_-`, region, aa_change)) %>%
  group_by(id_snv_chrm, region) %>%
  summarise(recurrence = n_distinct(sample)) %>% 
  filter(recurrence > 3) %>% 
  arrange(desc(recurrence))

snv.both.up0.01 <- all.snv.byStrand.bed_wider %>%
  filter(`AF_+` > 0.01 & `AF_-` > 0.01) %>%
  dplyr::select(c(id_snv_chrm, sample, `AF_+`, `AF_-`, region, aa_change)) %>%
  group_by(id_snv_chrm, region) %>%
  summarise(recurrence = n_distinct(sample)) %>% 
  filter(recurrence > 3) %>% 
  arrange(desc(recurrence)) %>% pull(id_snv_chrm)

## (+) (COULD IN - BUT WITH AF_- < 0.01)
all.snv.byStrand.bed_wider %>%
  filter(`AF_+` > 0.01 & `AF_-` < 0.01) %>%
  dplyr::select(c(id_snv_chrm, sample, `AF_+`, `AF_-`, region, aa_change)) %>%
  group_by(id_snv_chrm, region) %>%
  summarise(recurrence = n_distinct(sample)) %>%
  filter(id_snv_chrm %in% snv.both.up0.01)

## (-) (COULD IN + BUT WITH AF_+ < 0.01)
all.snv.byStrand.bed_wider %>%
  filter(`AF_+` < 0.01 & `AF_-` > 0.01) %>%
  dplyr::select(c(id_snv_chrm, sample, `AF_+`, `AF_-`, region, aa_change)) %>%
  group_by(id_snv_chrm, region) %>%
  summarise(recurrence = n_distinct(sample)) %>%
  filter(id_snv_chrm %in% snv.both.up0.01)




################################################################################
## BOTH
all.snv.byStrand.bed_wider %>%
  filter(`AF_+` > 0.01 & `AF_-` > 0.01,
         cds_region == "coding") %>%
  dplyr::select(c(id_snv_chrm, sample, `AF_+`, `AF_-`, region, aa_change)) %>%
  group_by(region, aa_change) %>%
  summarise(n = n()) 
