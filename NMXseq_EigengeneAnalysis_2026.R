
################################################################################
#                         Gene Ontology Eigengene Scores
################################################################################

# Emily N. Hilz, Phd. 
# The University of Texas at Austin
# ehilz@utexas.edu
# please reference _________________________________ when using this code

################################################################################
# PCA by sample - loop using linear mixed model - individual GO terms
################################################################################
# The following code provides example flow for one sex and brain region. 

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(ggplot2)

# Step 1: Load normalized gene expression and GO data files
IDs <- read.csv("NMXseq_IDs.csv")

norm <- read.csv("NMXseq_NormalizedGenes.csv") 

geneGOs <- read.csv("GeneGOs.csv") # Derived from Msibdbr package version 7.1
GOresults <- read.csv("GOall_Fem_OFC.csv") %>% filter(fdr < 0.05)

# Step 2: Match genes to GO terms and clusters
matched_genes <- GOresults %>%
  inner_join(geneGOs, by = c("id" = "GO")) %>%
  select(name, gene)  # Keep only relevant columns

matched_summary <- matched_genes %>%
  group_by(name) %>%
  summarise(Matching_Genes = list(unique(gene)),  # Store genes as a list
            Gene_Count = n())  # Count matching genes per summary

# Step 3: Initialize lists for results
lmm_results_list <- list()
pca_plots <- list()
pc_results_list <- list()  # Stores PC1 & PC2 results for all 

# Step 4: Loop through each unique name (i.e., GO term)
for (summary_category in unique(matched_genes$name)) {
  
  # Filter genes for the current summary category
  summary_genes <- matched_genes %>%
    filter(name == summary_category) %>%
    distinct(gene) %>%
    pull(gene)  # Extract gene names as a vector
  
  # Filter expression data for relevant genes; mind the filter
  filtered_norm_expr <- norm %>%
    filter(gene %in% summary_genes, Sex == "Female", Brain_Region == "OFC")  
  
  # Aggregate duplicate genes per Sample.ID
  filtered_genes <- filtered_norm_expr %>%
    group_by(gene, Sample.ID) %>%
    summarize(Expression = mean(Expression, na.rm = TRUE), .groups = "drop")
  
  # Transform to wide format
  filtered_wide <- filtered_genes %>%
    pivot_wider(names_from = gene, values_from = Expression) %>%
    column_to_rownames(var = "Sample.ID")  
  
  # Skip if not enough genes for PCA
  if (ncol(filtered_wide) < 2) next  
  
  # Perform PCA
  pca_results <- prcomp(filtered_wide, scale. = F, center = TRUE)
  
  # Convert PCA results to dataframe
  pca_df <- as.data.frame(pca_results$x) %>%
    rownames_to_column(var = "Sample.ID") %>%
    left_join(IDs, by = "Sample.ID")  # Add Treatment info
  
  # Skip if fewer than 2 PCs
  if (ncol(pca_df) < 3) next  
  
  # Fit Linear Model for PC1
  lm_PC1 <- lm(PC1 ~ Treatment, data = pca_df)
  
  # Fit Linear Model for PC2
  lm_PC2 <- lm(PC2 ~ Treatment, data = pca_df)
  
  # Store LMM results
  lmm_results_list[[summary_category]] <- list(
    PC1_LM = summary(lm_PC1),
    PC2_LM = summary(lm_PC2)
  )
  
  # Store PCA results with name
  pc_results_list[[summary_category]] <- pca_df %>%
    select(ID, Brain_Region, Sex, Treatment, PC1, PC2) %>%
    mutate(name = summary_category,  # Add cluster name
           ID = paste0("MDR ", ID))  # Modify ID
  
  # Generate PCA plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample.ID, color = Treatment)) +
    geom_point(size = 4, alpha = 0.8, show.legend = TRUE) +
    geom_text(vjust = 1.5, size = 3, check_overlap = TRUE, color = "black") +
    labs(title = paste("PCA of", summary_category, "cluster"),
         x = "PC1",
         y = "PC2",
         color = "GO Term") +
    theme_classic()
  
  # Store the plot in a list
  pca_plots[[summary_category]] <- p
}

# Optional: Print LMM results to look at treatment effects on PC scores per term
lmm_results_list

# Optional: Check for outliers by viewing PCA plots, perform exclusion as needed
pca_plots


# Step 5: Combine all PC1 & PC2 results into one dataframe
PCS <- bind_rows(pc_results_list)

# Step 6: Save PCA results as a CSV file
write.csv(PCS, "FemOFC_PCA_results_IndvGO.csv", row.names = FALSE)


#Optional: Filter and view significant models and plots
extract_significant_models <- function(lmm_results, threshold = 0.05) {
  significant_results <- list()
  
  for (cluster in names(lmm_results)) {
    model_PC1 <- lmm_results[[cluster]]$PC1_LM
    model_PC2 <- lmm_results[[cluster]]$PC2_LM
    
    # Extract p-values from coefficients
    pvals_PC1 <- coef(model_PC1)[, "Pr(>|t|)"]
    pvals_PC2 <- coef(model_PC2)[, "Pr(>|t|)"]
    
    # Check if any term in PC1 or PC2 is below threshold
    if (any(pvals_PC1 < threshold, na.rm = TRUE) || any(pvals_PC2 < threshold, na.rm = TRUE)) {
      significant_results[[cluster]] <- lmm_results[[cluster]]
    }
  }
  
  return(significant_results)
}
significant_lmm_results <- extract_significant_models(lmm_results_list, threshold = 0.055)
print(significant_lmm_results)


# Create an empty list to store plots
significant_lm_plots <- list()

# Extract only significant clusters
significant_clusters <- names(significant_lmm_results)

# Loop through each significant cluster and generate regression plots
for (cluster in significant_clusters) {
  
  # Retrieve PCA scores for this cluster
  cluster_data <- data %>% filter(name == cluster)
  
  # Skip if not enough data
  if (nrow(cluster_data) < 5) next
  
  # Generate regression plot for PC1
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample.ID, color = Treatment)) +
    geom_point(size = 4, alpha = 0.8, show.legend = TRUE) +
    geom_text(vjust = 1.5, size = 3, check_overlap = TRUE, color = "black") +
    labs(title = paste("PCA of", summary_category, "cluster"),
         x = "PC1",
         y = "PC2",
         color = "GO Term") +
    theme_classic() +
    scale_color_manual(values = c("NMX" = "orange", 
                                  "Veh" = "cornflowerblue")) 
  
  # Store only the plots from significant models
  significant_lm_plots[[paste(cluster, "PC1")]] <- p1
}

# View only significant model plots
significant_lm_plots

################################################################################
#                   PCA by GO term for module membership
################################################################################

library(dplyr)
library(tidyr)
library(tibble)


# Step 1: Load data files
IDs <- read.csv("NMXseq_IDs.csv")
norm <- read.csv("NMXseq_NormalizedGenes.csv")
geneGOs <- read.csv("GeneGOs.csv")
GOresults <- read.csv("GOall_Fem_NAC.csv")  %>% na.omit()

# Step 2: Match genes to GO terms and clusters
matched_genes <- GOresults %>%
  inner_join(geneGOs, by = c("id" = "GO")) %>%
  select(name, gene)  # Keep only relevant columns

matched_summary <- matched_genes %>%
  group_by(name) %>%
  summarise(Matching_Genes = list(unique(gene)),  # Store genes as a list
            Gene_Count = n())  # Count matching genes per summary

run_go_pca <- function(go_name,
                       sex = "Female",
                       brain_region = "NAC",
                       pc_keep = 1:2,
                       scale = FALSE,
                       center = TRUE,
                       membership_pc = 1,          # <-- which PC defines "eigengene"
                       membership_method = "pearson",
                       membership_use = "pairwise.complete.obs") {
  
  summary_genes <- matched_genes %>%
    filter(name == go_name) %>%
    distinct(gene) %>%
    pull(gene)
  
  filtered_norm_expr <- norm %>%
    filter(gene %in% summary_genes,
           Sex == sex,
           Brain_Region == brain_region)
  
  filtered_genes <- filtered_norm_expr %>%
    group_by(gene, Sample.ID) %>%
    summarize(Expression = mean(Expression, na.rm = TRUE), .groups = "drop")
  
  filtered_wide <- filtered_genes %>%
    pivot_wider(names_from = gene, values_from = Expression) %>%
    column_to_rownames("Sample.ID")
  
  if (ncol(filtered_wide) < 2) return(NULL)
  if (nrow(filtered_wide) < 3) return(NULL)
  
  pca <- prcomp(filtered_wide, scale. = scale, center = center)
  
  # ---- Loadings ----
  loadings <- as.data.frame(pca$rotation) %>%
    rownames_to_column("gene")
  
  pcs_exist <- paste0("PC", pc_keep)
  pcs_exist <- pcs_exist[pcs_exist %in% colnames(loadings)]
  loadings <- loadings %>%
    select(gene, all_of(pcs_exist)) %>%
    mutate(go_name = go_name,
           sex = sex,
           brain_region = brain_region)
  
  # ---- Variance explained ----
  var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  var_explained <- tibble(
    go_name = go_name,
    sex = sex,
    brain_region = brain_region,
    PC = paste0("PC", seq_along(var_explained)),
    pct = var_explained
  )
  
  # ---- Module membership (kME-like): cor(gene expression, eigengene) ----
  eigengene <- pca$x[, paste0("PC", membership_pc)]
  
  # Make sure columns are genes
  expr_mat <- as.matrix(filtered_wide)  # rows=samples, cols=genes
  
  kME <- apply(expr_mat, 2, function(gene_vec) {
    cor(gene_vec, eigengene,
        method = membership_method,
        use = membership_use)
  })
  
  membership <- tibble(
    gene = names(kME),
    kME = unname(kME),
    abs_kME = abs(kME),
    go_name = go_name,
    sex = sex,
    brain_region = brain_region
  ) %>%
    arrange(desc(abs_kME))
  
  list(
    pca = pca,
    loadings = loadings,
    var_explained = var_explained,
    membership = membership,
    eigengene = tibble(Sample.ID = rownames(filtered_wide),
                       eigengene = as.numeric(eigengene),
                       go_name = go_name,
                       sex = sex,
                       brain_region = brain_region)
  )
}

