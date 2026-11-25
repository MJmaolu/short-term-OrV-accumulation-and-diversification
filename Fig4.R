# Fig4.R
# SHORT-TERM ACCUMULATION AND DIVERSITY DYNAMICS
# MS + supplentary FIGURES and STATS
# MJ OLMO-UCEDA
################################################################################
# FIG 4A
################################################################################
################################################################################
# SNV RICHNESS AND ABUNDANCE (Q20)
################################################################################

snv.summary.q20 <- all.snv.byStrand.bed.q20 %>%
  filter(AF > 0.01) %>%
  group_by(
    sample,
    rna,
    strand,
    hpi,
    replica
  ) %>%
  summarise(
    
    richness = n(),
    
    abundance = sum(
      AF,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

snv.summary.time.q20 <- snv.summary.q20 %>%
  group_by(
    rna,
    strand,
    hpi
  ) %>%
  summarise(
    
    mean_richness =
      mean(richness, na.rm = TRUE),
    
    sd_richness =
      sd(richness, na.rm = TRUE),
    
    mean_abundance =
      mean(abundance, na.rm = TRUE),
    
    sd_abundance =
      sd(abundance, na.rm = TRUE),
    
    .groups = "drop"
  )


snv.summary.q20 %>%
  group_by(
    rna,
    strand
  ) %>%
  summarise(
    mean_richness = mean(richness),
    mean_abundance = mean(abundance)
  )

# Fig4.R

# SHORT-TERM ACCUMULATION AND DIVERSITY DYNAMICS
# MS FIGURES
# MJ OLMO-UCEDA
################################################################################
# FIG 4A
################################################################################
af = .01
p.snv.richness.q20 <- all.snv.byStrand.bed.q20 %>%
  filter(AF > af) %>%
  group_by(sample, rna, strand, hpi) %>%
  summarise(`Richness (AF > 0.01)` = n_distinct(id_snv)
  ) %>%
  mutate(viralSpecie = case_when((rna == "RNA 1" & strand == "+") ~  "RNA1 +",
                                 (rna == "RNA 1" & strand == "-") ~  "RNA1 -",
                                 (rna == "RNA 2" & strand == "+") ~  "RNA2 +",
                                 (rna == "RNA 2" & strand == "-") ~  "RNA2 -")) %>%
  ggplot(.,
         aes(hpi,
             `Richness (AF > 0.01)`,
             group = interaction(rna, strand),
             fill = viralSpecie,
             color = rna)) +
  
  # SD ribbon
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      
      data.frame(
        y = m,
        ymin = m - s,
        ymax = m + s
      )
    },
    geom = "ribbon",
    alpha = .15,
    color = NA
  ) +
  
  # points
  geom_point(aes(shape = strand,
                 fill = case_when((rna == "RNA 1" & strand == "+") ~  "#4E3BA7",
                                  (rna == "RNA 1" & strand == "-") ~  "white",
                                  (rna == "RNA 2" & strand == "+") ~  "#90B083",
                                  (rna == "RNA 2" & strand == "-") ~  "white",
                 )),
             position = position_jitter(height = 0,
                                        width = .2),
             size = 2,
             shape = 21,
             alpha = .8) +
  
  # # mean points
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 21,
  #              size = 2,
  #              alpha = .9) +
  # 
  # mean lines
  stat_summary(fun = mean,
               geom = "line",
               aes(linetype = strand),
               linewidth = .4) +
  
  scale_color_manual(values = c("RNA 1" = "#4E3BA7",
                                "RNA 2" = "#90B083")) +
  
  scale_fill_manual(values = c("RNA1 +" = "#4E3BA7",
                               "RNA2 +" = "#90B083",
                               "RNA1 -" = "#4E3BA7",
                               "RNA2 -" = "#90B083",
                               "#4E3BA7" = "#4E3BA7",
                               "#90B083" = "#90B083",
                               "white" = "white"
  )) +
  
  scale_linetype_manual(values = c("+" = "solid",
                                   "-" = "dashed")) +
  
  # facet_wrap(~strand,
  #            scales = "free") +
  
  ylab("Richness (number SNVs AF > 0.01)") +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

# FIG 4A normalized by segment length
################################################################################
af = .01
all.snv.byStrand.bed.q20 %>%
  filter(AF > af,
         DP > 20) %>%
  group_by(sample, rna, strand, hpi) %>%
  summarise(`Richness (AF > 0.01)` = n(),
  ) %>%
  mutate(viralSpecie = case_when((rna == "RNA 1" & strand == "+") ~  "RNA1 +",
                                 (rna == "RNA 1" & strand == "-") ~  "RNA1 -",
                                 (rna == "RNA 2" & strand == "+") ~  "RNA2 +",
                                 (rna == "RNA 2" & strand == "-") ~  "RNA2 -"),
         `Richness normalized` = (`Richness (AF > 0.01)` / (ifelse(rna == "RNA 1",
                                                                   3421,
                                                                   2574)))) %>%
  ggplot(.,
         aes(hpi,
             `Richness normalized`,
             group = interaction(rna, strand),
             fill = viralSpecie,
             color = rna)) +
  
  # SD ribbon
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      
      data.frame(
        y = m,
        ymin = m - s,
        ymax = m + s
      )
    },
    geom = "ribbon",
    alpha = .25,
    color = NA
  ) +
  
  # points
  geom_point(aes(shape = strand,
                 fill = case_when((rna == "RNA 1" & strand == "+") ~  "#4E3BA7",
                                  (rna == "RNA 1" & strand == "-") ~  "white",
                                  (rna == "RNA 2" & strand == "+") ~  "#90B083",
                                  (rna == "RNA 2" & strand == "-") ~  "white",
                 )),
             position = position_jitter(height = 0,
                                        width = .2),
             size = 2,
             shape = 21,
             alpha = .8) +
  
  # # mean points
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 21,
  #              size = 2,
  #              alpha = .9) +
  # 
  # mean lines
  stat_summary(fun = mean,
               geom = "line",
               aes(linetype = strand),
               linewidth = .4) +
  
  scale_color_manual(values = c("RNA 1" = "#4E3BA7",
                                "RNA 2" = "#90B083")) +
  
  scale_fill_manual(values = c("RNA1 +" = "#4E3BA7",
                               "RNA2 +" = "#90B083",
                               "RNA1 -" = "#4E3BA7",
                               "RNA2 -" = "#90B083",
                               "#4E3BA7" = "#4E3BA7",
                               "#90B083" = "#90B083",
                               "white" = "white"
  )) +
  
  scale_linetype_manual(values = c("+" = "solid",
                                   "-" = "dashed")) +
  
  # facet_wrap(~strand,
  #            scales = "free") +
  
  ylab("Richness (number SNVs AF > 0.01") +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))


################################################################################
# FIG 4B (ABUNDANCE AS SUM OF SNVs'AF)
################################################################################
af = .01
p.snv.abundance.q20 <- all.snv.byStrand.bed.q20 %>%
  filter(AF > af) %>%
  group_by(sample, rna, strand, hpi) %>%
  summarise(`Abundance (AF > 0.01)` = sum(AF),
  ) %>%
  mutate(viralSpecie = case_when((rna == "RNA 1" & strand == "+") ~  "RNA1 +",
                                 (rna == "RNA 1" & strand == "-") ~  "RNA1 -",
                                 (rna == "RNA 2" & strand == "+") ~  "RNA2 +",
                                 (rna == "RNA 2" & strand == "-") ~  "RNA2 -")) %>%
  ggplot(.,
         aes(hpi,
             `Abundance (AF > 0.01)`,
             group = interaction(rna, strand),
             fill = viralSpecie,
             color = rna)) +
  
  # SD ribbon
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      
      data.frame(
        y = m,
        ymin = m - s,
        ymax = m + s
      )
    },
    geom = "ribbon",
    alpha = .15,
    color = NA
  ) +
  
  # points
  geom_point(aes(shape = strand,
                 fill = case_when((rna == "RNA 1" & strand == "+") ~  "#4E3BA7",
                                  (rna == "RNA 1" & strand == "-") ~  "white",
                                  (rna == "RNA 2" & strand == "+") ~  "#90B083",
                                  (rna == "RNA 2" & strand == "-") ~  "white",
                 )),
             position = position_jitter(height = 0,
                                        width = .2),
             size = 2,
             shape = 21,
             alpha = .8) +
  
  # # mean points
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 21,
  #              size = 2,
  #              alpha = .9) +
  # 
  # mean lines
  stat_summary(fun = mean,
               geom = "line",
               aes(linetype = strand),
               linewidth = .4) +
  
  scale_color_manual(values = c("RNA 1" = "#4E3BA7",
                                "RNA 2" = "#90B083")) +
  
  scale_fill_manual(values = c("RNA1 +" = "#4E3BA7",
                               "RNA2 +" = "#90B083",
                               "RNA1 -" = "#4E3BA7",
                               "RNA2 -" = "#90B083",
                               "#4E3BA7" = "#4E3BA7",
                               "#90B083" = "#90B083",
                               "white" = "white"
  )) +
  
  scale_linetype_manual(values = c("+" = "solid",
                                   "-" = "dashed")) +
  
  # facet_wrap(~strand,
  #            scales = "free") +
  
  ylab(expression(paste("Abundance (", sum(AF > 0.01), ")"))) +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

p.snv.richness.q20 + p.snv.abundance.q20

################################################################################
# SUPPLEMENTARY
################################################################################
################################################################################
# FIG S2A
################################################################################
af = .01
p.snv.richness.global.q20 <- all.snv.byStrand.bed.q20 %>%
  group_by(sample, rna, hpi) %>%
  filter(AF > af) %>%
  summarise(`Richness (AF > 0.01)` = n_distinct(id_snv),
            `Abundance (AF > 0.01)` = sum(AF),
            DP = sum(DP)) %>%
  ggplot(.,
         aes(hpi,
             `Richness (AF > 0.01)`,
             group = rna,
             fill = rna,
             color = rna)) +
  
  # SD ribbon
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      
      data.frame(
        y = m,
        ymin = m - s,
        ymax = m + s
      )
    },
    geom = "ribbon",
    alpha = .25,
    color = NA
  ) +
  
  # points
  geom_point(aes(fill = rna,
  ),
  position = position_jitter(height = 0,
                             width = .2),
  size = 2,
  shape = 21,
  alpha = .8) +
  
  # # mean points
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 21,
  #              size = 2,
  #              alpha = .9) +
  # 
  # mean lines
  stat_summary(fun = mean,
               geom = "line",
               linewidth = .4) +
  
  scale_color_manual(values = c("RNA 1" = "#4E3BA7",
                                "RNA 2" = "#90B083")) +
  
  scale_fill_manual(values = c("RNA 1" = "#4E3BA7",
                               "RNA 2" = "#90B083")) +
  
  # facet_wrap(~strand,
  #            scales = "free") +
  
  ylab("Richness (number SNVs AF > 0.01") +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))


# NORMALIZING BY SEGMENT SIZE

all.snv.byStrand.bed.q20 %>%
  filter(AF > af) %>%
  group_by(sample, rna, hpi) %>%
  summarise(`Richness (AF > 0.01)` = n()) %>%
  mutate(`Richness normalized` = (`Richness (AF > 0.01)` / (ifelse(rna == "RNA 1",
                                                                   3421,
                                                                   2574)))
  ) %>%
  
  ggplot(.,
         aes(hpi,
             `Richness normalized` ,
             group = rna,
             fill = rna,
             color = rna)) +
  
  # SD ribbon
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      
      data.frame(
        y = m,
        ymin = m - s,
        ymax = m + s
      )
    },
    geom = "ribbon",
    alpha = .25,
    color = NA
  ) +
  
  # points
  geom_point(aes(fill = rna,
  ),
  position = position_jitter(height = 0,
                             width = .2),
  size = 2,
  shape = 21,
  alpha = .8) +
  
  # # mean points
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 21,
  #              size = 2,
  #              alpha = .9) +
  # 
  # mean lines
  stat_summary(fun = mean,
               geom = "line",
               linewidth = .4) +
  
  scale_color_manual(values = c("RNA 1" = "#4E3BA7",
                                "RNA 2" = "#90B083")) +
  
  scale_fill_manual(values = c("RNA 1" = "#4E3BA7",
                               "RNA 2" = "#90B083")) +
  
  # facet_wrap(~strand,
  #            scales = "free") +
  
  ylab("Normalized richness") +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

################################################################################
# FIG S2B (ABUNDANCE AS SUM OF SNVs'AF)
################################################################################
af = .01
all.snv.byStrand.bed.q20 %>%
  filter(AF > af) %>%
  group_by(sample, rna, hpi) %>%
  summarise(`Abundance (AF > 0.01)` = sum(AF)) %>%
  
  mutate(`Abundance normalized` = (`Abundance (AF > 0.01)` / (ifelse(rna == "RNA 1", 
                                                                     3421,
                                                                     2574)))
  ) %>%
  
  ggplot(.,
         aes(hpi,
             `Abundance normalized`,
             group = rna,
             fill = rna,
             color = rna)) +
  
  # SD ribbon
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      
      data.frame(
        y = m,
        ymin = m - s,
        ymax = m + s
      )
    },
    geom = "ribbon",
    alpha = .25,
    color = NA
  ) +
  
  # points
  geom_point(aes(fill = rna,
  ),
  position = position_jitter(height = 0,
                             width = .2),
  size = 2,
  shape = 21,
  alpha = .8) +
  
  # # mean points
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 21,
  #              size = 2,
  #              alpha = .9) +
  # 
  # mean lines
  stat_summary(fun = mean,
               geom = "line",
               linewidth = .4) +
  
  scale_color_manual(values = c("RNA 1" = "#4E3BA7",
                                "RNA 2" = "#90B083")) +
  
  scale_fill_manual(values = c("RNA 1" = "#4E3BA7",
                               "RNA 2" = "#90B083")) +
  
  # facet_wrap(~strand,
  #            scales = "free") +
  
  ylab(expression(paste("Abundance normalized (", sum(AF > 0.01), "/ RNA segment length)"))) +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

################################################################################
# FIG S2B (ABUNDANCE AS SUM OF SNVs'AF)
################################################################################
af = .01
p.snv.abundance.global.q20 <- all.snv.byStrand.bed.q20 %>%
  filter(AF > af) %>%
  group_by(sample, rna, hpi) %>%
  summarise(`Abundance (AF > 0.01)` = sum(AF),
            DP = sum(DP)
  ) %>%
  ggplot(.,
         aes(hpi,
             `Abundance (AF > 0.01)`,
             group = rna,
             fill = rna,
             color = rna)) +
  
  # SD ribbon
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      s <- sd(x, na.rm = TRUE)
      
      data.frame(
        y = m,
        ymin = m - s,
        ymax = m + s
      )
    },
    geom = "ribbon",
    alpha = .25,
    color = NA
  ) +
  
  # points
  geom_point(aes(fill = rna,
  ),
  position = position_jitter(height = 0,
                             width = .2),
  size = 2,
  shape = 21,
  alpha = .8) +
  
  # # mean points
  # stat_summary(fun = mean,
  #              geom = "point",
  #              shape = 21,
  #              size = 2,
  #              alpha = .9) +
  # 
  # mean lines
  stat_summary(fun = mean,
               geom = "line",
               linewidth = .4) +
  
  scale_color_manual(values = c("RNA 1" = "#4E3BA7",
                                "RNA 2" = "#90B083")) +
  
  scale_fill_manual(values = c("RNA 1" = "#4E3BA7",
                               "RNA 2" = "#90B083")) +
  
  scale_y_continuous(expand = c(0,0)) +
  
  # facet_wrap(~strand,
  #            scales = "free") +
  
  ylab(expression(paste("Abundance (", sum(AF > 0.01), ")"))) +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))


p.snv.richness.global.q20 + p.snv.abundance.global.q20

################################################################################
## STATS
################################################################################
### PER SEGMENT

summary.snvs.segment.q20 <- all.snv.byStrand.bed.q20 %>%
  filter(AF > af) %>%
  group_by(sample, rna, hpi) %>%
  summarise(`Richness (AF > 0.01)` = n_distinct(id_snv),
            `Abundance (AF > 0.01)` = sum(AF)
            
  ) %>%
  mutate(segment_length = ifelse(rna == "RNA 1",
                                 3421,
                                 2574))

summary.snvs.q20 <- all.snv.byStrand.bed.q20 %>%
  filter(AF > af) %>%
  group_by(sample, rna, strand, hpi) %>%
  summarise(`Richness (AF > 0.01)` = n_distinct(id_snv),
            `Abundance (AF > 0.01)` = sum(AF)
            
  ) %>%
  mutate(viralSpecie = case_when((rna == "RNA 1" & strand == "+") ~  "RNA1 +",
                                 (rna == "RNA 1" & strand == "-") ~  "RNA1 -",
                                 (rna == "RNA 2" & strand == "+") ~  "RNA2 +",
                                 (rna == "RNA 2" & strand == "-") ~  "RNA2 -"),
         `Richness normalized` = (`Richness (AF > 0.01)` / (ifelse(rna == "RNA 1",
                                                                   3421,
                                                                   2574))),
         `Abundance normalized` = (`Abundance (AF > 0.01)` / (ifelse(rna == "RNA 1",
                                                                     3421,
                                                                     2574))),
         segment_length = ifelse(rna == "RNA 1",
                                 3421,
                                 2574))

glmer(
  `Richness (AF > 0.01)` ~ hpi + rna +
    offset(log(segment_length)) +
    (1 | sample),
  family = poisson,
  data = summary.snvs.segment.q20
) %>%
  summary()

glmer(
  `Abundance (AF > 0.01)` ~ hpi + rna +
    offset(log(segment_length)) +
    (1 | sample),
  family = Gamma(link="log"),
  data = summary.snvs.segment.q20
) %>%
  summary()

# STRAND-RESOLVED ANALYSIS
performance::check_overdispersion(
  glmer(
    `Richness (AF > 0.01)` ~
      hpi + rna + strand +
      offset(log(segment_length)) +
      (1|sample),
    family = poisson,
    data = summary.snvs.q20
  )
)
# Poisson
glmer(
  `Richness (AF > 0.01)` ~ hpi + rna + strand +
    offset(log(segment_length)) +
    (1 | sample),
  family = poisson,
  data = summary.snvs.q20
) %>% summary()

glmer(
  `Abundance (AF > 0.01)` ~ hpi + rna + strand +
    offset(log(segment_length)) +
    (1 | sample),
  family = Gamma(link="log"),
  data = summary.snvs.q20
) %>% summary()

###############
# THIS CODE
##############
# PER SEGMENT

## richness
mod.rich.seg.q20 <- glmer(
  `Richness (AF > 0.01)` ~
    hpi +
    rna +
    offset(log(segment_length)) +
    (1 | sample),
  family = poisson,
  data = summary.snvs.segment.q20
)

summary(mod.rich.seg.q20)

performance::check_overdispersion(
  mod.rich.seg.q20
)

## abundance
mod.abund.seg.q20 <- glmmTMB::glmmTMB(
  `Abundance (AF > 0.01)` ~
    hpi +
    rna +
    offset(log(segment_length)) +
    (1 | sample),
  family = Gamma(link = "log"),
  data = summary.snvs.segment.q20
)

summary(mod.abund.seg.q20)

# STRAND-RESOLVED
## richness
mod.rich.strand.q20 <- glmer(
  `Richness (AF > 0.01)` ~
    hpi +
    rna +
    strand +
    offset(log(segment_length)) +
    (1 | sample),
  family = poisson,
  data = summary.snvs.q20
)

summary(mod.rich.strand.q20)

performance::check_overdispersion(
  mod.rich.strand.q20
)

# abundance
mod.abund.strand.q20 <- glmmTMB::glmmTMB(
  `Abundance (AF > 0.01)` ~
    hpi +
    rna +
    strand +
    offset(log(segment_length)) +
    (1 | sample),
  family = Gamma(link = "log"),
  data = summary.snvs.q20
)

summary(mod.abund.strand.q20)

#### results
summary(mod.rich.seg.q20)$coefficients
summary(mod.abund.seg.q20)

summary(mod.rich.strand.q20)$coefficients
summary(mod.abund.strand.q20)

performance::check_overdispersion(mod.rich.seg.q20)

performance::check_overdispersion(mod.rich.strand.q20)





# ABUNDANCE VS RICHNESS

snv.mean.q20 <- snv.summary.q20 %>%
  group_by(rna, strand, hpi) %>%
  summarise(
    mean_richness  = mean(richness),
    sd_richness    = sd(richness),
    mean_abundance = mean(abundance),
    sd_abundance   = sd(abundance),
    .groups = "drop"
  )

ggplot() +
  
  ## puntos individuales
  geom_point(
    data = snv.summary.q20,
    aes(
      richness,
      abundance,
      color = rna,
      shape = strand
    ),
    alpha = .7,
    size = 2
  ) +
  
  # ## SD horizontal (richness)
  # geom_segment(
  #   data = snv.mean.q20,
  #   aes(
  #     x = mean_richness - sd_richness,
  #     xend = mean_richness + sd_richness,
  #     y = mean_abundance,
  #     yend = mean_abundance,
  #     color = rna
  #   ),
  #   linewidth = .4
  # ) +
  # 
  # ## SD vertical (abundance)
  # geom_segment(
  #   data = snv.mean.q20,
  #   aes(
  #     x = mean_richness,
  #     xend = mean_richness,
  #     y = mean_abundance - sd_abundance,
  #     yend = mean_abundance + sd_abundance,
  #     color = rna
  #   ),
  #   linewidth = .4
  # ) +
  
  ## trayectorias temporales
  geom_path(
    data = snv.mean.q20,
    aes(
      mean_richness,
      mean_abundance,
      color = rna,
      linetype = strand,
      group = interaction(rna, strand)
    ),
    linewidth = .5
  ) +
  
  ## puntos medios
  geom_point(
    data = snv.mean.q20,
    aes(
      mean_richness,
      mean_abundance,
      color = rna,
      fill = strand
    ),
    shape = 21,
    size = 4,
    stroke = 1
  ) +
  
  ## tiempo
  geom_text(
    data = snv.mean.q20,
    aes(
      mean_richness,
      mean_abundance,
      label = hpi,
      color = rna
    ),
    size = 3.5,
    nudge_y = 0.04
  ) +
  
  scale_color_manual(values = c(
    "RNA 1" = "#4E3BA7",
    "RNA 2" = "#90B083"
  )) +
  
  scale_fill_manual(values = c(
    "+" = "black",
    "-" = "pink"
  )) +
  
  scale_linetype_manual(values = c(
    "+" = "solid",
    "-" = "dashed"
  )) +
  
  theme_classic() +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = "transparent"),
        axis.line = element_blank())
