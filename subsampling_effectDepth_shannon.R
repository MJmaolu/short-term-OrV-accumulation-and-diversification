################################################################################
# DEPTH-MATCHED SUBSAMPLING
################################################################################
# Sample: orv_12h_1
#
# Q20 dataset
#
# MJ OLMO-UCEDA
# 2026/11/16
################################################################################
library(QSutils)
set.seed(123)
################################################################################
# JOIN STRANDS
################################################################################
subsample.df <- orv12h1_plus_q20 %>%
  dplyr::select(
    chrom,
    pos,
    A_pp,
    C_pp,
    G_pp,
    T_pp,
    reads_pp,
    norm.H,
    `D_j/site_pairs`
  ) %>%
  dplyr::rename(
    A_plus = A_pp,
    C_plus = C_pp,
    G_plus = G_pp,
    T_plus = T_pp,
    depth_plus = reads_pp,
    H_plus = norm.H,
    APD_plus = `D_j/site_pairs`
  ) %>%
  left_join(orv12h1_minus_q20 %>%
       dplyr::select(
        chrom,
        pos,
        A_pp,
        C_pp,
        G_pp,
        T_pp,
        reads_pp,
        norm.H,
        `D_j/site_pairs`
      ) %>%
      dplyr::rename(
        A_minus = A_pp,
        C_minus = C_pp,
        G_minus = G_pp,
        T_minus = T_pp,
        depth_minus = reads_pp,
        H_minus = norm.H,
        APD_minus = `D_j/site_pairs`
      ),
    by = c("chrom","pos")
  ) %>%
  mutate(
    min_depth = pmin(depth_plus, depth_minus)) %>%
  filter(
    min_depth > 10
  )

################################################################################
# FUNCTIONS
# APD is calculated as indicated in Shao et al, 2014
################################################################################

calc_metrics <- function(
    A,
    C,
    G,
    T
){
  counts <- c(
    A,
    C,
    G,
    T
  )
  total <- sum(counts)
  if(total <= 1){
    return(
      c(
        H = NA,
        APD = NA
      )
    )
  }
  H <- QSutils::NormShannon(
    counts
  )
  D_j <-
    A*(C+G+T) +
    C*(G+T) +
    G*T
  APD <-
    D_j /
    (total*(total-1)/2)
  c( H = H,
    APD = APD
  )
}

################################################################################
# BOOTSTRAP SUBSAMPLING
################################################################################

n.boot <- 100

boot.results <- vector(
  "list",
  n.boot
)

for(i in seq_len(n.boot)){
  tmp <- subsample.df
  plus.metrics <-
    map2_dfr(
      seq_len(nrow(tmp)),
      tmp$min_depth,
      
      function(j,N){
        counts <- c(
          tmp$A_plus[j],
          tmp$C_plus[j],
          tmp$G_plus[j],
          tmp$T_plus[j]
        )
        
        total <- sum(counts)
        probs <- counts/total
        sampled <-
          as.vector(
            rmultinom(
              1,
              size = N,
              prob = probs
            )
          )
        metrics <-
          calc_metrics(
            sampled[1],
            sampled[2],
            sampled[3],
            sampled[4]
          )
        
        tibble(
          H = metrics["H"],
          APD = metrics["APD"]
        )
      }
    )
  
  minus.metrics <-
    map2_dfr(
      seq_len(nrow(tmp)),
      tmp$min_depth,
      
      function(j,N){
        counts <- c(
          tmp$A_minus[j],
          tmp$C_minus[j],
          tmp$G_minus[j],
          tmp$T_minus[j]
        )
        
        total <- sum(counts)
        probs <- counts/total
        sampled <-
          as.vector(
            rmultinom(
              1,
              size = N,
              prob = probs
            )
          )
        metrics <-
          calc_metrics(
            sampled[1],
            sampled[2],
            sampled[3],
            sampled[4]
          )
        
        tibble(
          H = metrics["H"],
          APD = metrics["APD"]
        )
      }
    )
  
  boot.results[[i]] <-
    tibble(
      iteration = i,
      mean_H_plus =
        mean(
          plus.metrics$H,
          na.rm = TRUE
        ),
      mean_H_minus =
        mean(
          minus.metrics$H,
          na.rm = TRUE
        ),
      mean_APD_plus =
        mean(
          plus.metrics$APD,
          na.rm = TRUE
        ),
      mean_APD_minus =
        mean(
          minus.metrics$APD,
          na.rm = TRUE
        )
    )
}

boot.results <- bind_rows(
    boot.results
  )

################################################################################
# SUMMARY
################################################################################

boot.results %>%
  summarise(
    H_plus = mean(mean_H_plus),
    H_minus = mean(mean_H_minus),
    APD_plus =
      mean(
        mean_APD_plus
      ),
    APD_minus =
      mean(
        mean_APD_minus
      )
  )

################################################################################
# DELTA DISTRIBUTIONS
################################################################################

boot.results <- boot.results %>%
  
  mutate(
    
    delta_H =
      mean_H_minus -
      mean_H_plus,
    
    delta_APD =
      mean_APD_minus -
      mean_APD_plus
    
  )

################################################################################
# PLOT DELTA DISTRIBUTIONS (100 BOOTSTRAP)
################################################################################

ggplot(
  boot.results,
  aes(delta_H)
) +
  geom_histogram(
    bins = 70,
    fill = "grey70",
    colour = "black"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "red"
  ) +
  
  theme_classic() +
  
  labs(
    x = expression(Delta*Shannon~"(" * "\u2212" * ")" - "(+)"),
    y = "100 bootstrap replicates"
  ) +
  scale_y_continuous(expand = c(0,0)) +
  theme(panel.border = element_rect(fill = "transparent",
                                    color = "black"),
        axis.line = element_blank()) +

################################################################################
# CORRELATION SHANNON-APD
################################################################################

ggplot(
  boot.results,
  aes(
    delta_H,
    delta_APD
  )
) +
  geom_point() +
  geom_smooth(
    method = "lm"
  ) +
  labs(
    x = expression(Delta*Shannon~"(" * "\u2212" * ")" - "(+)"),
    y = expression(Delta*APD~"(" * "\u2212" * ")" - "(+)")
  ) +
  theme_classic() +
  theme(panel.border = element_rect(fill = "transparent",
                                    color = "black"),
        axis.line = element_blank())


cor.test(boot.results$delta_H,
         boot.results$delta_APD,
         method = "spearman")
