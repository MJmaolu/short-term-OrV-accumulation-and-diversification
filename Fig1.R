# SHORT-TERM ACCUMULATION AND DIVERSITY DYNAMICS
# MS FIGURES
# MJ OLMO-UCEDA
################################################################################
# FIG 1A with SD as shadow 
################################################################################
VPKB.bed %>%
  pivot_longer(cols = c(`RNA1 +`, 
                        `RNA2 +`,
                        `RNA1 -`,
                        `RNA2 -`
  ),
  names_to = "RNA",
  values_to = "VPKB") %>%
  rowwise() %>%
  mutate(timepoint = strsplit(sample, "_")[[1]][2],
         hpi = as.numeric(str_replace(timepoint, "h", "")),
         replica = as.numeric(strsplit(strsplit(sample, "_")[[1]][3], "\\.")[[1]][1]),
         strand = ifelse(grepl("-", RNA),
                         "-",
                         "+"),
         chrom = ifelse(grepl("1", RNA),
                        "RNA 1",
                        "RNA 2")) %>%
  ungroup() %>%
  ggplot(.,
         aes(x = hpi,
             y = VPKB,
             group = RNA,
             color = chrom,
             fill = RNA
         )) +
  
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
                 fill = case_when((chrom == "RNA 1" & strand == "+") ~  "#4E3BA7",
                                  (chrom == "RNA 1" & strand == "-") ~  "white",
                                  (chrom == "RNA 2" & strand == "+") ~  "#90B083",
                                  (chrom == "RNA 2" & strand == "-") ~  "white",
                                  )),
             position = position_jitter(height = 0,
                                        width = .1),
             size = 1.5,
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
  
  ylab("VPMB") +
  xlab("Time (hpi)") +
  
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

# FIG 1B
# STATS SATURATION CURVE
################################################################################
# FITING `MICHAELIS-MENTEN` RELATION
################################################################################
# Fit the Michaelis-Menten model for RNA1 + vs RNA1 -
fit_rna1.MM <- nls(`RNA1 +` ~ Vmax * (`RNA1 -`) / (Km + `RNA1 -`),
                   data = VPKB.bed,
                   start = list(Vmax = max(VPKB.bed$`RNA1 +`, na.rm = TRUE), 
                                Km = median(VPKB.bed$`RNA1 -`, na.rm = TRUE)))

# Get the summary of the fitted model
summary(fit_rna1.MM)
# Predict values using the fitted model
VPKB.bed$`RNA1 +.pred` <- predict(fit_rna1.MM)

# IN THE TEXT WE TALK SATURATION MODEL ADN VMAX = A AND KM = B
# Fit the Michaelis-Menten model for RNA1 + vs RNA1 -
fit_rna2.MM <- nls(`RNA2 +` ~ Vmax * (`RNA2 -`) / (Km + `RNA2 -`),
                   data = VPKB.bed,
                   start = list(Vmax = max(VPKB.bed$`RNA2 +`, na.rm = TRUE), 
                                Km = median(VPKB.bed$`RNA2 -`, na.rm = TRUE)))

# Get the summary of the fitted model
summary(fit_rna2.MM)
VPKB.bed$`RNA2 +.pred` <- predict(fit_rna2.MM)

#=======================================
# COMPARISON WITH A LINEAL MODEL
# LINEAR MODEL RNA1
fit_rna1.lin <- lm(`RNA1 +` ~ `RNA1 -`, data = VPKB.bed)
summary(fit_rna1.lin)
# LINEAR MODEL RNA2
fit_rna2.lin <- lm(`RNA2 +` ~ `RNA2 -`, data = VPKB.bed)
summary(fit_rna2.lin)

# RNA1
AIC(fit_rna1.MM, fit_rna1.lin)

# RNA2
AIC(fit_rna2.MM, fit_rna2.lin)
#=======================================
# COMPARISON WITH EXPONENTIAL MODEL ~ GR
# GR proxy (log-linear exponential)
fit_rna1.exp <- lm(log(`RNA1 +`) ~ `RNA1 -`, data = VPKB.bed)
fit_rna2.exp <- lm(log(`RNA2 +`) ~ `RNA2 -`, data = VPKB.bed)

summary(fit_rna1.exp)
summary(fit_rna2.exp)


# RNA1
AIC(fit_rna1.MM, fit_rna1.exp)
# RNA2
AIC(fit_rna2.MM, fit_rna2.exp)


################################################################################
## Adding the IC as ribbon
################################################################################

VPKB.clean <- VPKB.bed %>%
  mutate(
    RNA1_plus  = `RNA1 +`,
    RNA1_minus = `RNA1 -`,
    RNA2_plus  = `RNA2 +`,
    RNA2_minus = `RNA2 -`
  )

fit_rna1.MM <- nls(
  RNA1_plus ~ Vmax * RNA1_minus / (Km + RNA1_minus),
  data = VPKB.clean,
  start = list(
    Vmax = max(VPKB.clean$RNA1_plus, na.rm = TRUE),
    Km = median(VPKB.clean$RNA1_minus, na.rm = TRUE)
  )
)

fit_rna2.MM <- nls(
  RNA2_plus ~ Vmax * RNA2_minus / (Km + RNA2_minus),
  data = VPKB.clean,
  start = list(
    Vmax = max(VPKB.clean$RNA2_plus, na.rm = TRUE),
    Km = median(VPKB.clean$RNA2_minus, na.rm = TRUE)
  )
)
new_rna1 <- data.frame(
  RNA1_minus = seq(
    min(VPKB.clean$RNA1_minus),
    max(VPKB.clean$RNA1_minus),
    length.out = 200
  )
)

pred_rna1 <- predFit(
  fit_rna1.MM,
  newdata = new_rna1,
  interval = "confidence",
  level = 0.95
)

pred_rna1_df <- cbind(
  new_rna1,
  fit   = pred_rna1[,1],
  lower = pred_rna1[,2],
  upper = pred_rna1[,3]
)

new_rna2 <- data.frame(
  RNA2_minus = seq(
    min(VPKB.clean$RNA2_minus),
    max(VPKB.clean$RNA2_minus),
    length.out = 200
  )
)

pred_rna2 <- predFit(
  fit_rna2.MM,
  newdata = new_rna2,
  interval = "confidence",
  level = 0.95
)

pred_rna2_df <- cbind(
  new_rna2,
  fit   = pred_rna2[,1],
  lower = pred_rna2[,2],
  upper = pred_rna2[,3]
)

# New Fig2B

ggplot() +
  geom_ribbon(data = pred_rna1_df,
              aes(x = RNA1_minus,
                  ymin = lower,
                  ymax = upper),
              fill = "#4E3BA7",
              alpha = 0.15) +
  
  geom_ribbon(data = pred_rna2_df,
              aes(x = RNA2_minus,
                  ymin = lower,
                  ymax = upper),
              fill = "#90B083",
              alpha = 0.15) +
  
  geom_line(data = pred_rna1_df,
            aes(x = RNA1_minus,
                y = fit),
            color = "#4E3BA7",
            linewidth = .4) +
  
  geom_line(data = pred_rna2_df,
            aes(x = RNA2_minus,
                y = fit),
            color = "#90B083",
            linewidth = .4) +
  
  
  geom_point(data = VPKB.bed,
             aes(x = `RNA1 -`,
                 y = `RNA1 +`),
             fill = "#4E3BA7",
             shape = 21,
             size = 2.3,
             alpha = .8) +
  
  geom_point(data = VPKB.clean,
             aes(x = `RNA2 -`,
                 y = `RNA2 +`),
             fill = "#90B083",
             shape = 21,
             size = 2.3,
             alpha = .8) +
  
  
  theme_classic() +
  ylab("VPMB (+)-strand") +
  xlab("VPMB (-)-strand") +
  scale_y_continuous(expand = c(.01,0)) +
  scale_x_continuous(expand = c(.01,0)) +
  coord_cartesian(clip = "off") +
  #scale_y_log10() +
  #scale_x_log10() +
  theme_classic() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = "transparent"))

cor.test(VPKB.bed$`RNA1 +`, VPKB.bed$`RNA1 -`, method = "spearman")
cor.test(VPKB.bed$`RNA2 +`, VPKB.bed$`RNA2 -`, method = "spearman")
