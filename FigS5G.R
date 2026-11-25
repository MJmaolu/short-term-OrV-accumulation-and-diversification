################################################################################
# FIGURE S5G (Q20)
# 
# MJ OLMO-UCEDA
# 2026/11/23
################################################################################

library(tidyverse)
library(ggrepel)

################################################################################
# HOTSPOT REGIONS IDENTIFIED IN FIGURE 7
################################################################################

# hotspot_regions.q20 <- tibble(
#   protein = c(
#     "RdRP",
#     "CP",
#     "CP"
#   ),
#   codon_start = c(
#     622,
#     176,
#     204
#   ),
#   codon_end = c(
#     628,
#     182,
#     210
#   )
# )
hotspot_regions.q20
################################################################################
# SNVs USED FOR dN-dS ANALYSIS
################################################################################

snv.selection.q20 <- all.snv.byStrand.bed.q20 %>%
  
  filter(
    strand == "+",
    AF > 0.01,
    DP > 20,
    cds_region == "coding"
  ) %>%
  
  left_join(
    genomic_position_annotation,
    by = c("rna","POS")
  )

################################################################################
# IDENTIFY VARIANTS INSIDE HOTSPOTS
################################################################################

snv.selection.q20 <- snv.selection.q20 %>%
  rowwise() %>%
  mutate(
    
    hotspot = any(
      protein == hotspot_regions.q20$protein &
        codon_position >= hotspot_regions.q20$codon_start &
        codon_position <= hotspot_regions.q20$codon_end
    )
  ) %>%
  ungroup()

################################################################################
# REPRESENTATION USED IN ORIGINAL FIGURE
################################################################################

snv.selection.q20 <- snv.selection.q20 %>%
  mutate(
    
    AF.plot = case_when(
      aa_change == "synonymous" ~ -AF,
      aa_change == "non-synonymous" ~ AF
    )
    
  )

################################################################################
# LABELS (OPTIONAL)
################################################################################

labels.q20 <- snv.selection.q20 %>%
  group_by(id_snv_chrm) %>%
  slice_max(AF, n = 1) %>%
  ungroup() %>%
  filter(
    AF > 0.05 | hotspot
  )

################################################################################
# FIGURE S5G
################################################################################

figS5G_top.q20 <- snv.selection.q20 %>%
  
  filter(protein == "RdRP") %>%
  
  ggplot(
    aes(
      POS,
      AF.plot
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = .2
  ) +
  
  geom_point(
    aes(shape = aa_change,
        colour = hotspot
    ),
    size = 1.3
  ) +
  facet_grid(
    ~ protein,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_colour_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "#ff1493"
    )
  ) +
  
  geom_text(data = snv.selection.q20 %>%
              filter(protein == "RdRP",
                     abs(AF.plot) > 0.1),
            aes(x = POS,
                y = 1.1 * AF.plot,
                label = ifelse(
                  aa_change == "non-synonymous",
                  SNV,
                  id_snv)
                ),
            inherit.aes = FALSE) +
  
  scale_shape_manual(values = c("synonymous" = 20,
                                "non-synonymous" = 18,
                                "NA" = 1)) +
  labs(
    x = "",
    y = "AF"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = "transparent"),
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black")
  )

figS5G_bottom.q20 <- snv.selection.q20 %>%

  filter(
    protein %in% c(
      "CP",
      "delta"
    )
  ) %>%

  ggplot(
    aes(
      POS,
      AF.plot
    )
  ) +

  geom_hline(
    yintercept = 0,
    linewidth = .2
  ) +

  geom_point(
    aes(shape = aa_change,
        colour = hotspot
    ),
    size = 1.3
  ) +
  geom_text(data = snv.selection.q20 %>%
              filter(protein != "RdRP",
                     abs(AF.plot) > 0.1),
            aes(x = POS,
                y = 1.1 * AF.plot,
                label = ifelse(
                  aa_change == "non-synonymous",
                  SNV,
                  id_snv)
            ),
            inherit.aes = FALSE) +
  
  facet_grid(
    ~ protein,
    scales = "free_x",
    space = "free_x"
  ) +

  scale_colour_manual(
    values = c(
      "FALSE" = "black",
      "TRUE" = "#ff1493"
    )
  ) +
  
  scale_shape_manual(values = c("synonymous" = 20,
                                "non-synonymous" = 18,
                                "NA" = 1)) +

  labs(
    x = "Nucleotide position",
    y = "AF"
  ) +

  theme_classic() +

  theme(
    legend.position = "none",
    panel.border = element_rect(fill = "transparent"),
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black")
  )


figS5G.q20 <- figS5G_top.q20 /
  figS5G_bottom.q20 +
  plot_layout(
    heights = c(1, 1)
  )

figS5G.q20
