# Fig4.R

# SHORT-TERM ACCUMULATION AND DIVERSITY DYNAMICS
# MS FIGURES
# MJ OLMO-UCEDA
################################################################################
# FIG 4A
################################################################################
af = .01
all.snv.byStrand.bed %>%
  filter(AF > af) %>%
  group_by(sample, rna, strand, hpi) %>%
  summarise(`Richness (AF > 0.01)` = n(),
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

# FIG 4A normalized by segment length
################################################################################
af = .01
all.snv.byStrand.bed %>%
  filter(AF > af) %>%
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
all.snv.byStrand.bed %>%
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
  
  ylab(expression(paste("Abundance (", sum(AF > 0.01), ")"))) +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

################################################################################
# SUPPLEMENTARY
################################################################################
################################################################################
# FIG S2A
################################################################################
af = .01
all.snv.byStrand.bed %>%
  filter(AF > af) %>%
  group_by(sample, rna, hpi) %>%
  summarise(`Richness (AF > 0.01)` = n_distinct(id_snv)
  ) %>%
  
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

all.snv.byStrand.bed %>%
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
all.snv.byStrand.bed %>%
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
  
  ylab("Richness (number SNVs AF > 0.01") +
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
all.snv.byStrand.bed %>%
  filter(AF > af) %>%
  group_by(sample, rna, hpi) %>%
  summarise(`Abundance (AF > 0.01)` = sum(AF),
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
## STATS
################################################################################
glmer(
  `Richness (AF > 0.01)` ~ hpi + rna +
    offset(log(segment_length)) +
    (1 | sample),
  family = poisson,
  data = summary.snvs.segment
) %>%
  summary()

glmer(
  `Abundance (AF > 0.01)` ~ hpi + rna +
    offset(log(segment_length)) +
    (1 | sample),
  family = Gamma(link="log"),
  data = summary.snvs.segment
) %>%
  summary()

# STRAND-RESOLVED ANALYSIS
# Poisson
glmer(
  `Richness (AF > 0.01)` ~ hpi + rna + strand +
    offset(log(segment_length)) +
    (1 | sample),
  family = poisson,
  data = summary.snvs
) %>% summary()

glmer(
  `Abundance (AF > 0.01)` ~ hpi + rna + strand +
    offset(log(segment_length)) +
    (1 | sample),
  family = Gamma(link="log"),
  data = summary.snvs
) %>% summary()
