# Fig S1
# SHANNON ENTROPY PER SITE
################################################################################
df_complete.normH %>%
  filter(norm.H > 0) %>%
  ggplot(.) +
  geom_hline(yintercept = 0,
             linewidth = .2) +
  # Proteins
  geom_rect(data = orv.proteins %>%
              filter(protein != "putative"),
            aes(xmin = start,
                xmax = end,
                ymin = -Inf,
                ymax = Inf,
                fill = protein),
            alpha = 0.1,
            inherit.aes = F) +
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
  
  # Shannon +
  geom_point(data = . %>%
               filter(strand == "plus",
                      reads_pp > 10),
             aes(x = pos, 
                 y = norm.H,
                 color = norm.H,
                 # alpha = ifelse(norm.H_plus < 0.1,
                 #                0.7,
                 #                1)
  ),
  size = 0.5) +
  # shannon in -
  geom_point(data = . %>%
               filter(strand == "minus",
                      reads_pp > 10),
             aes(x = pos, 
                 y = - norm.H,
                 color = norm.H,
                 # alpha = ifelse(norm.H_minus < 0.1,
                 #                0.7,
                 #                1)
  ),
  size = 0.5) +
  geom_segment(data = . %>%
                 filter(strand == "plus", 
                        reads_pp > 10),
               aes(x = pos, 
                   xend = pos, 
                   y = 0, 
                   yend = norm.H, 
                   color = norm.H),
               linewidth = 0.1) +
  
  geom_segment(data = . %>%
                 filter(strand == "minus", 
                        reads_pp > 10),
               aes(x = pos, 
                   xend = pos, 
                   y = 0, 
                   yend = -norm.H, 
                   color = norm.H),
               linewidth = 0.1) +
  # #Add info SNV
  # geom_hline(yintercept = - 0.1,
  #            linewidth = 0.3) +
  scale_color_viridis_c(option = "magma",
                        direction = -1) +
  scale_fill_manual(values = color.proteins) +
  facet_grid(hpi ~ rna,
             scales = "free_x",
             space = "free_x"
             ) +
  # ggrepel::geom_label_repel(label = ifelse(abs(Sn_plus - Sn_minus) > 0.4),
  #                           pos,
  #                           NA) +
  xlab("Genomic position (nt)") +
  ylab("Normalized Shannon Entropy") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.x = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 10),
        panel.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent")
  )
