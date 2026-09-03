
# TO ADD IN FIG 5
# FIXED SCALING FOR FIG 5 + KDE OVERLAY

# 1. known from your plot
max_snv <- 13

# 2. KDE scale reference
max_kde <- max(kde_plot$observed, na.rm = TRUE)

# 3. scaling factor
kde_scale <- max_snv / max_kde

# 4. rescale KDE
kde_plot2 <- kde_plot %>%
  mutate(
    observed_s = observed * kde_scale,
    high_s = high * kde_scale,
    low_s  = low * kde_scale
  )

max(kde_plot2$observed_s)


# EXTRACT KERNEL PEAKS
get_hotspot_stats <- function(res, hotspots){
  
  do.call(
    rbind,
    lapply(seq_len(nrow(hotspots)), function(i){
      
      tmp <- res %>%
        filter(
          x >= hotspots$start[i],
          x <= hotspots$end[i]
        )
      
      data.frame(
        start = hotspots$start[i],
        end = hotspots$end[i],
        peak_position = tmp$x[which.max(tmp$observed)],
        peak_density = max(tmp$observed)
      )
      
    })
  )
}

hotspots_summary <- bind_rows(
  
  get_hotspot_stats(res_rna1_plus, hotspots_rna1_plus) %>%
    mutate(rna="RNA 1", strand="+"),
  
  get_hotspot_stats(res_rna1_minus, hotspots_rna1_minus) %>%
    mutate(rna="RNA 1", strand="-"),
  
  get_hotspot_stats(res_rna2_plus, hotspots_rna2_plus) %>%
    mutate(rna="RNA 2", strand="+"),
  
  get_hotspot_stats(res_rna2_minus, hotspots_rna2_minus) %>%
    mutate(rna="RNA 2", strand="-")
  
) %>%
  dplyr::select(
    rna,
    strand,
    start,
    end,
    peak_position,
    peak_density
  )

#===============================================================================
## WITH THE KDE INFO
#===============================================================================

ggplot() +
  geom_hline(yintercept = 0) +
  # Proteins LIMITS
  geom_rect(data = orv.proteins,
            aes(xmin = start,
                xmax = end,
                ymin = -Inf,
                ymax = Inf,
                fill = protein),
            alpha = 0.1,
            
            inherit.aes = F) +
  scale_fill_futurama() +
  ggnewscale::new_scale_fill() +
  
  # segments limits
  geom_rect(data = data.frame(rna = c("RNA 1", "RNA 2"),
                              start = c(1,1),
                              end = c(3421,2574)),
            aes(xmin = start,
                xmax = end,
                ymin = -Inf,
                ymax = Inf,
                color = as.factor(rna)),
            fill = "transparent",
            alpha = 0.1,
            linewidth = .1,
            inherit.aes = F) +

# KDE ENRICHMENT LAYER (SCALED TO SNV AXIS)
geom_line(
  data = kde_plot2,
  aes(
    x = x,
    y = ifelse(strand == "+", observed_s, -observed_s),
    color = strand,
    group = interaction(rna, strand)
  ),
  linewidth = 0.4,
  inherit.aes = FALSE
) +
  
  geom_line(
    data = kde_plot2,
    aes(
      x = x,
      y = ifelse(strand == "+", high_s, -high_s),
      color = strand,
      group = interaction(rna, strand)
    ),
    linetype = "dashed",
    linewidth = 0.2,
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = color.strand) +
  ggnewscale::new_scale_color() +
  # palmprint motifs
  geom_rect(data = palmprint.orv,
            aes(xmin = nt.start,
                xmax = nt.end,
                ymin = -Inf,
                ymax = Inf),
            color = "black",
            alpha = 0.1,
            linewidth = .05,
            inherit.aes = F) +
  ## Add info SNV (only in )
  ##########################
# only in Plus
geom_segment(data = all.snv.byStrand.bed %>%
               filter(AF > af,
                      id_snv_chrm %in% c(unique2Plus.rna1, unique2Plus.rna2)) %>%
               group_by(id_snv_chrm) %>%
               mutate(sample_count = n_distinct(sample)) %>%  # Contar número de ocurrencias de cada SNV
               distinct(id_snv, .keep_all = T),
             aes(x = POS,
                 xend = POS,
                 y = sample_count,
                 yend = 0,
                 #size = AF,
                 #alpha = AF,
                 #alpha = strand,
                 color = strand,
             ),
             linewidth = .1,
             #shape = 20
) +
  # only in PLUS
  # POINT
  geom_point(data = all.snv.byStrand.bed %>%
               filter(AF > af,
                      id_snv_chrm %in% c(unique2Plus.rna1, unique2Plus.rna2)) %>%
               group_by(id_snv_chrm) %>%
               mutate(sample_count = n_distinct(sample)) %>%  # Contar número de ocurrencias de cada SNV
               distinct(id_snv, .keep_all = T),
             aes(x = POS,
                 y = sample_count,
                 #size = AF,
                 #alpha = AF,
                 #alpha = strand,
                 shape = aa_change,
                 color = strand,
             ),
             size = 2,
             #shape = 20
  ) +
  
  # Etiquetas para los SNVs (únicamente una vez por SNV con su número de ocurrencias)
  ## (+)
  geom_text(data = all.snv.byStrand.bed %>%
              filter(AF > af,
                     id_snv_chrm %in% c(unique2Plus.rna1, unique2Plus.rna2)) %>%
              group_by(id_snv_chrm) %>%
              mutate(sample_count = n_distinct(sample)) %>%  # Contar número de ocurrencias de cada SNV
              distinct(id_snv, .keep_all = T) %>%
              filter(sample_count > 1),
            # Mantener solo una fila por SNV
            aes(x = POS,
                y = .6*sample_count,  # Ajusta la posición de la etiqueta en el eje Y
                label = ifelse(aa_change == "non-synonymous",
                               paste(SNV, "(", sample_count, ")", sep = ""),
                               paste(id_snv,  "(", sample_count, ")", sep = ""))),  # Etiqueta: SNV (número de ocurrencias)
            size = 3,  # Ajusta el tamaño del texto
            hjust = 0,
            angle = 45,  # Ajusta el ángulo si es necesario
            color = "black",  # Color de la etiqueta
            inherit.aes = F
  ) +
  ## MINUS
  # only in minus
  geom_segment(data = all.snv.byStrand.bed %>%
                 filter(AF > af,
                        id_snv_chrm %in% c(unique2minus.rna1, unique2minus.rna2)) %>%
                 group_by(id_snv_chrm) %>%
                 mutate(sample_count = n_distinct(sample)) %>%  # Contar número de ocurrencias de cada SNV
                 distinct(id_snv, .keep_all = T),
               aes(x = POS,
                   xend = POS,
                   y = - sample_count,
                   yend = 0,
                   #alpha = AF,
                   #alpha = strand,
                   color = strand,
               ),
               linewidth = .1,
               #shape = 20
  ) +
  # POINT
  geom_point(data = all.snv.byStrand.bed %>%
               filter(AF > af,
                      id_snv_chrm %in% c(unique2minus.rna1, unique2minus.rna2)) %>%
               group_by(id_snv_chrm) %>%
               mutate(sample_count = n_distinct(sample)) %>%  # Contar número de ocurrencias de cada SNV
               distinct(id_snv, .keep_all = T),
             aes(x = POS,
                 y = -sample_count,
                 #size = AF,
                 #alpha = AF,
                 #alpha = strand,
                 shape = aa_change,
                 color = strand,
             ),
             size = 2,
             #shape = 20
  ) +
  # Labeling the SNVs in the MINUS category
  # Etiquetas para los SNVs (únicamente una vez por SNV con su número de ocurrencias)
  geom_text(data = all.snv.byStrand.bed %>%
              filter(AF > af,
                     id_snv_chrm %in% c(unique2minus.rna1, unique2minus.rna2)) %>%
              group_by(id_snv_chrm) %>%
              mutate(sample_count = n_distinct(sample)) %>%  # Contar número de ocurrencias de cada SNV
              distinct(id_snv, .keep_all = T) %>%
              filter(sample_count > 1),
            # Mantener solo una fila por SNV
            aes(x = POS,
                y = -.8 * sample_count,  # Ajusta la posición de la etiqueta en el eje Y
                # label = ifelse(aa_change == "non-synonymous",
                #                paste(SNV, "(", sample_count, ")", sep = ""),
                #                paste(id_snv,  "(", sample_count, ")", sep = "")),
                label = paste(id_snv,  "(", sample_count, ")", sep = ""),
                color = strand,  # Color de la etiqueta
            ),  # Etiqueta: SNV (número de ocurrencias)
            size = 3,  # Ajusta el tamaño del texto
            hjust = 1,
            angle = 45,  # Ajusta el ángulo si es necesario
            inherit.aes = F
  ) +
  
  # # IN BOTH STRANDs BUT NOT IN THE SAME SAMPLE
  # geom_segment(data = all.snv.byStrand.bed %>% 
  #                filter(AF > af,
  #                       #replica == 1,
  #                       id_snv_chrm %in% c(inBothButDiffered.rna1, inBothButDiffered.rna2)),
  #              aes(x = POS,
  #                  xend = POS,
  #                  y = -.2,
  #                  yend = .2,
  #                  alpha = 0.7,
#                  #alpha = strand
#              ),
#              color = "white",
#              linewidth = .1,
#              #shape = 20
# ) +
# IN BOTH
geom_segment(data = all.snv.byStrand.bed %>% 
               filter(AF > af,
                      #replica == 1,
                      id_snv_chrm %in% c(realBoth.rna1, realBoth.rna2)),
             aes(x = POS,
                 xend = POS,
                 y = -.2,
                 yend = .2,
                 alpha = 0.7,
                 #alpha = strand
             ),
             color = "#803058",
             linewidth = .2,
             #shape = 20
) +
  #scale_color_viridis_c(direction = -1) +
  #coord_cartesian(ylim = c(-0.2, 1)) +
  facet_grid(rna ~ .,
             scales = "free_x",
             space = "free_x"
  ) +
  #scale_color_gradientn(colors = c("white", "#A9D6E5","yellow", "#FFB800", "#FF6A00", "#FF2A00", "black")) +
  scale_color_manual(values = color.strand) +
  # scale_alpha_manual(values = c("+" = .6,
  #                               "-" = .4)) +
  scale_shape_manual(values = c("synonymous" = 20,
                                "non-synonymous" = 18,
                                "NA" = 1)) +
  #scale_y_continuous(expand = c(.1,.1)) +
  scale_y_continuous(
    name = "SNV recurrence (samples)",
    sec.axis = sec_axis(~ . / kde_scale,
                        name = "Weighted SNV density")
  ) + 
  xlab("Genomic position (nt)") +
  ylab("") +
  #xlim(c(3150,3300)) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.line.x = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 10),
        panel.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.spacing.y = unit(1.2, 
                               "lines")  # Adjust the spacing between facets vertically
  ) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.line.x = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 10),
        panel.background = element_blank(),
        #strip.text = element_text(size = 12),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"),
        # panel.spacing.y = unit(1.2, 
        #                        "lines")
  )
