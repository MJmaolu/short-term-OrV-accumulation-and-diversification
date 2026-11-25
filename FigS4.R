# FigS4
# DATA Q20
#
# 2026/11/20
# MJ OLMO-UCEDA
################################################################################

snv_plus_df_rna1.q20 <- snv_events.q20 %>%
  filter(category == "plus",
         rna == "RNA 1") %>%
  mutate(mut = sub(".*([A-Z]>[A-Z]).*", "\\1", id_snv_chrm)) %>%
  group_by(mut) %>%
  summarise(
    n = sum(sample_count),
    n_events = n(),
    .groups = "drop"
  ) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(n))

snv_plus_df_rna1.q20 <- snv_plus_df_rna1.q20 %>%
  mutate(
    from = substr(mut, 1, 1),
    to   = substr(mut, 3, 3),
  ) %>%
  dplyr::rename("mut(+)" = "mut")

snv_plus_df_rna2.q20 <- snv_events.q20 %>%
  filter(category == "plus",
         rna == "RNA 2") %>%
  mutate(mut = sub(".*([A-Z]>[A-Z]).*", "\\1", id_snv_chrm)) %>%
  group_by(mut) %>%
  summarise(
    n = sum(sample_count),
    n_events = n(),
    .groups = "drop"
  ) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(n))

snv_plus_df_rna2.q20 <- snv_plus_df_rna2.q20 %>%
  mutate(
    from = substr(mut, 1, 1),
    to   = substr(mut, 3, 3),
  ) %>%
  dplyr::rename("mut(+)" = "mut")

snv_plus_df.q20 <- rbind(snv_plus_df_rna1.q20 %>%
                       mutate(rna = "RNA 1",
                              presence = "+"),
                     snv_plus_df_rna2.q20 %>%
                       mutate(rna = "RNA 2",
                              presence = "+"))


snv_minus_df_rna1.q20 <- snv_events.q20 %>%
  filter(category == "minus",
         rna == "RNA 1") %>%
  mutate(mut = sub(".*([A-Z]>[A-Z]).*", "\\1", id_snv_chrm)) %>%
  group_by(mut) %>%
  summarise(
    n = sum(sample_count),
    n_events = n(),
    .groups = "drop"
  ) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(n))

snv_minus_df_rna1.q20 <- snv_minus_df_rna1.q20 %>%
  mutate(
    from = substr(mut, 1, 1),
    to   = substr(mut, 3, 3),
    
    mut_minus = paste0(
      comp(from),
      ">",
      comp(to)
    )
  ) %>%
  dplyr::rename("mut(+)" = "mut",
         "mut(-)" = "mut_minus")

snv_minus_df_rna2.q20 <- snv_events.q20 %>%
  filter(category == "minus",
         rna == "RNA 2") %>%
  mutate(mut = sub(".*([A-Z]>[A-Z]).*", "\\1", id_snv_chrm)) %>%
  group_by(mut) %>%
  summarise(
    n = sum(sample_count),
    n_events = n(),
    .groups = "drop"
  ) %>%
  mutate(freq = n / sum(n)) %>%
  arrange(desc(n))

snv_minus_df_rna2.q20 <- snv_minus_df_rna2.q20 %>%
  mutate(
    from = substr(mut, 1, 1),
    to   = substr(mut, 3, 3),
    
    mut_minus = paste0(
      comp(from),
      ">",
      comp(to)
    )
  ) %>%
  dplyr::rename("mut(+)" = "mut",
         "mut(-)" = "mut_minus")


snv_minus_df.q20 <- rbind(snv_minus_df_rna1.q20 %>%
                        mutate(rna = "RNA 1",
                               presence = "-"),
                      snv_minus_df_rna2.q20 %>%
                        mutate(rna = "RNA 2",
                               presence = "-"))


snv_changes_df.q20 <- rbind(snv_minus_df.q20 %>%
                          dplyr::select(-`mut(+)`) %>%
                          dplyr::rename(`mut` = `mut(-)`),
                        snv_plus_df.q20 %>%
                          dplyr::rename(`mut` = `mut(+)`))


# CALCULATE TRANSITION & TRANSVERSIONS

transition <- c("A>G","G>A","C>U","U>C")

snv_changes_trans.q20 <- snv_changes_df.q20 %>%
  mutate(type = ifelse(mut %in% transition,
                       "transition",
                       "transversion")) %>%
  group_by(presence, rna, type) %>%
  summarise(n = sum(freq), .groups = "drop") %>%
  group_by(presence, rna) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()




mut_order.q20 <- snv_changes_df.q20 %>%
  filter(rna == "RNA 1",
         presence == "-") %>%
  arrange(desc(freq)) %>%
  pull(mut)


# FIG S4A)
ggplot(snv_changes_trans.q20, 
       aes(x = presence, 
           y = prop, 
           fill = type)) +
  geom_col(color = "black", 
           alpha = 0.9,
           width = 1) +
  facet_wrap(~ rna) +
  scale_fill_manual(values = c("transition" = "#B69ECD", # BLUE
                               "transversion" = "#5C7BB4")) + # PURPLE
  theme_classic() +
  labs(y = "Proportion",
       x = "Substitution type") +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black"),
    axis.line = element_blank(),
    panel.border = element_rect(colour = "black",
                                fill = "transparent")
  ) +

  
# FIG S4B) TILE WITH ALL SUSTITUTIONS
ggplot(snv_changes_df.q20 %>%
         mutate(mut = factor(mut, 
                             levels = rev(mut_order.q20))),
       aes(rna, 
           mut, 
           fill = freq)) +
  geom_tile(color = "black") +
  geom_text(aes(label = sprintf("%.2f", freq),
                colour = freq > max(freq)/2),
            size = 3,
            show.legend = FALSE) +
  facet_grid(~ presence) +
  scale_fill_gradient(low = "white", 
                      high = "black") +
  scale_colour_manual(values = c("black", 
                                 "white")) +
  labs(x = "RNA",
       y = "Substitution") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

snv_changes_df.q20$mut %>% table()
(snv_changes_df.q20 %>%
  mutate(mut = factor(mut, 
                      levels = rev(mut_order))))$mut %>% table()
