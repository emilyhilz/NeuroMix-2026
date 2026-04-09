library("doBy")
library("car")
library("DescTools")
library("tidyverse")
library("ggplot2")
library("reshape2")
library("readr")
library("dplyr")


################################################################################
#                      Trait Data Prediction by GO - HFFP
################################################################################

# Emily N. Hilz, Phd. 
# The University of Texas at Austin
# ehilz@utexas.edu
# please reference _________________________________ when using this code

################################################################################
# Example code for one sex/brain region and outcome variable (high fat food):
################################################################################

PCs<- read.csv("MaleNAC_PCA_results_IndvGO.csv")
HFFP <- read.csv("NMX_HFFP.csv") %>%
  group_by(ID) %>%  
  summarise(
    HF_g = mean(HF_g, na.rm = TRUE)
  ) %>%
  ungroup()

data <- left_join(HFFP,PCs) %>%
  mutate(HF_g=scale(HF_g))%>%
  group_by(Treatment, name) %>%  
  mutate(
    PC1_z = as.numeric(scale(PC1))  # z-score within group, unused here
  ) %>%
  ungroup() %>%
  na.omit() #%>% filter(between(PC1_z, -3, 3),between(HF_g, -3, 3))

custom_order <- c("Veh", "NMX")

data$Treatment <- factor(data$Treatment, levels = custom_order)

################################################################################
#                   HF_g permutation across samples
################################################################################
# ------------------------------------------------------------
# Permutation p-values for a single GO term, using global shuffle of HF_g
# Model: HF_g ~ PC1 * Treatment
#
# Returns permutation p-values for:
#   - interaction (as a whole): nested-model F (full vs no-interaction)
#   - overall model: model F-statistic
#   - PC1 main effect: |t|
#   - Treatment main effect: |t| (2-level Treatment assumed)
#
# Selection rule:
#   keep if (p_int < alpha) AND (p_model < alpha) AND (p_PC1 < alpha OR p_Treatment < alpha)
#
# Uses +1 correction so p-values never equal 0.
# ------------------------------------------------------------
perm_pvals_and_keep <- function(df, n_perm = 1000, alpha = 0.05) {
  
  df <- df %>% mutate(Treatment = factor(Treatment))
  
  # observed fits
  mod_full <- lm(HF_g ~ PC1 * Treatment, data = df)
  sm_full  <- summary(mod_full)
  
  # adjusted R^2 (observed)
  adjR2 <- sm_full$adj.r.squared
  
  # overall model F
  F_model_obs <- unname(sm_full$fstatistic[1])
  
  # interaction as a whole: compare full vs reduced (no interaction)
  mod_red <- lm(HF_g ~ PC1 + Treatment, data = df)
  F_int_obs <- unname(anova(mod_red, mod_full)$F[2])
  
  # main-effect t stats (2-level Treatment assumed)
  rn <- rownames(sm_full$coefficients)
  if (!("PC1" %in% rn)) return(NULL)
  
  trt_term <- grep("^Treatment", rn, value = TRUE)[1]
  if (is.na(trt_term)) return(NULL)
  
  t_PC1_obs <- unname(sm_full$coefficients["PC1", "t value"])
  t_Trt_obs <- unname(sm_full$coefficients[trt_term, "t value"])
  
  # permutation containers
  F_model_perm <- numeric(n_perm)
  F_int_perm   <- numeric(n_perm)
  t_PC1_perm   <- numeric(n_perm)
  t_Trt_perm   <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    dfp <- df
    dfp$HF_g <- sample(dfp$HF_g)
    
    mf <- lm(HF_g ~ PC1 * Treatment, data = dfp)
    sf <- summary(mf)
    
    mr <- lm(HF_g ~ PC1 + Treatment, data = dfp)
    
    F_model_perm[b] <- unname(sf$fstatistic[1])
    F_int_perm[b]   <- unname(anova(mr, mf)$F[2])
    
    # coefficient t's (guard in case of odd singularities)
    co <- sf$coefficients
    if ("PC1" %in% rownames(co) && trt_term %in% rownames(co)) {
      t_PC1_perm[b] <- unname(co["PC1", "t value"])
      t_Trt_perm[b] <- unname(co[trt_term, "t value"])
    } else {
      t_PC1_perm[b] <- NA_real_
      t_Trt_perm[b] <- NA_real_
    }
  }
  
  # +1 corrected permutation p-values
  p_model <- (1 + sum(F_model_perm >= F_model_obs, na.rm = TRUE)) / (n_perm + 1)
  p_int   <- (1 + sum(F_int_perm   >= F_int_obs,   na.rm = TRUE)) / (n_perm + 1)
  p_PC1   <- (1 + sum(abs(t_PC1_perm) >= abs(t_PC1_obs), na.rm = TRUE)) / (n_perm + 1)
  p_Trt   <- (1 + sum(abs(t_Trt_perm) >= abs(t_Trt_obs), na.rm = TRUE)) / (n_perm + 1)
  
  keep <- (p_int < alpha) && (p_model < alpha) && (p_PC1 < alpha || p_Trt < alpha)
  
  tibble(
    n = nrow(df),
    adjR2 = adjR2,
    F_model_obs = F_model_obs, p_model = p_model,
    F_int_obs   = F_int_obs,   p_int   = p_int,
    t_PC1_obs   = t_PC1_obs,   p_PC1   = p_PC1,
    t_Trt_obs   = t_Trt_obs,   p_Treatment = p_Trt,
    keep = keep
  )
}


# ------------------------------------------------------------
# Run across ALL GO terms
# ------------------------------------------------------------
set.seed(123)
n_perm <- 1000
alpha  <- 0.05

all_clusters <- unique(data$name)

perm_results_all <- lapply(all_clusters, function(clust) {
  
  df <- data %>% filter(name == clust) %>% na.omit()
  
  # basic requirements
  if (nrow(df) < 5) return(NULL)
  if (n_distinct(df$Treatment) < 2) return(NULL)
  if (any(is.na(df$PC1)) || any(is.na(df$HF_g))) return(NULL)
  
  out <- perm_pvals_and_keep(df, n_perm = n_perm, alpha = alpha)
  if (is.null(out)) return(NULL)
  
  mutate(out, name = clust, .before = 1)
}) %>%
  bind_rows()

# GO terms that pass permutation
kept_terms <- perm_results_all %>%
  filter(keep) %>%
  arrange(p_int, p_model)
kept_terms %>% print(n = 100)
kept_terms$name


# Perform same analysis with z-score outlier filtering
data <- left_join(HFFP,PCs) %>%
  mutate(HF_g=scale(HF_g))%>%
  group_by(Treatment, name) %>%  
  mutate(
    PC1_z = as.numeric(scale(PC1))  
  ) %>%
  ungroup() %>%
  na.omit() %>% filter(between(PC1_z, -3, 3),between(HF_g, -3, 3))

perm_results_all <- lapply(all_clusters, function(clust) {
  
  df <- data %>% filter(name == clust) %>% na.omit()
  
  # basic requirements
  if (nrow(df) < 5) return(NULL)
  if (n_distinct(df$Treatment) < 2) return(NULL)
  if (any(is.na(df$PC1)) || any(is.na(df$HF_g))) return(NULL)
  
  out <- perm_pvals_and_keep(df, n_perm = n_perm, alpha = alpha)
  if (is.null(out)) return(NULL)
  
  mutate(out, name = clust, .before = 1)
}) %>%
  bind_rows()

# GO terms that pass permutation after outlier removal
kept_terms2 <- perm_results_all %>%
  filter(keep) %>%
  arrange(p_int, p_model)
kept_terms2 %>% print(n = 100)
kept_terms2$name

# Create a dataframe that contains only models that meet significance requirements
# both with and without outlier removal.
df3 <- inner_join(kept_terms1, kept_terms2, by = c("name"))