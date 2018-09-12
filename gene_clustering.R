# This Rscript uses gene expression data of different timepoints
# throughout development of the mouse retina. The data can be found here:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101986

library(tidyverse)
library(ggplot2)
#library(edgeR)
#library(xlsx)
library(apcluster)

read_in_file <- function(file_name) {
  read.table(paste0("~/Documents/GIT/gene_analysis/GSE101986_", file_name, ".txt"), 
             sep ="\t",
             quote = "",
             stringsAsFactors = F,
             header = T,
             row.names = 1
  )
}

gene_normalized <- read_in_file("Gene_Normalized-CPM")

# The data are available in different stages of processing.
# If desired the can be read in using the functions below

# transcript_counts <- read_in_file("Transcript_Counts")
# differential_expression <- read_in_file("DifferentialExpression")
# gene_counts <- read_in_file("Gene_Counts") 


# The dataframe contains unused columns and multiple measures for the same timepoint.
# Create a new dataframe using only the required information and combine similar timepoints
# Can function-ize this, refactor later...

gene_normalized_clean <- data.frame(gene = gene_normalized$external_gene_name,
                                    ensemble_id = rownames(gene_normalized),
                                    E11 = (gene_normalized$E11.1 + gene_normalized$E11.2)/2,
                                    E12 = (gene_normalized$E12.1 + gene_normalized$E12.2)/2,
                                    E14 = (gene_normalized$E14.1 + gene_normalized$E14.2)/2,
                                    E16 = (gene_normalized$E16.1 + gene_normalized$E16.2)/2,
                                    P0 = (gene_normalized$P0.1 + gene_normalized$P0.2)/2,
                                    P2 = (gene_normalized$P2.1 + gene_normalized$P2.2)/2,
                                    P4 = (gene_normalized$P4.1 + gene_normalized$P4.2)/2,
                                    P6 = (gene_normalized$P6.1 + gene_normalized$P6.2)/2,
                                    P10 = (gene_normalized$P10.1 + gene_normalized$P10.2)/2,
                                    P14 = (gene_normalized$P14.1 + gene_normalized$P14.2)/2,
                                    P21 = (gene_normalized$P21.1 + gene_normalized$P21.2)/2,
                                    P28 = (gene_normalized$P28.1 + gene_normalized$P28.2)/2)

####################

# The dataframe contains data for all detected genes (~50,000) at every timepoint (12)
# This is a lot of information to plow through. Based on previous literature, groups of
# genes exist. Check if can apply cut-offs to identify genes that follow patterns in expression.

# This graph will show genes that decrease 10-fold from E11 to E14
# and then don't change from E14 to P2.
# It also cuts out genes that are expressed below 100 CPM at E11
gene_normalized_clean %>% 
  filter(E12 <= E11*0.15 & abs(P2 - E12) < abs(E12 - E11)/2 & E11 > 25) %>% 
  gather(key = "tp", value = "expression", E11:P2) %>% 
  ggplot(aes(x = factor(tp, levels = c(colnames(gene_normalized_clean[,2:7]))), y = expression, group = gene, color = gene)) + 
  geom_line(size = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


# This graph will show genes that decrease 6.66-fold from E11 to E14
# and then don't change from E12 to P2.
# It also cuts out genes that are expressed below 25 CPM at E11
# The intention of this graph was to capture Mitf, a gene of interest

gene_normalized_clean %>% 
  filter(E12 <= E11*0.25 & abs(P2 - E12) < abs(E12 - E11)/2 & E11 > 25) %>% 
  gather(key = "tp", value = "expression", E11:P2) %>%
  ggplot(aes(x = factor(tp, levels = c(colnames(gene_normalized_clean[,2:7]))), y = expression, group = gene, color = gene)) +# %in% c("Mitf"))) + 
  geom_line(size = 1) +
  theme_classic() +
  labs(colour = "Gene", x = "Time Point", y = "Expression")

# Identify genes that go up from E11 to E14 and stay up through P2

gene_normalized_clean %>% 
  filter((E14 - E11) > abs(P2 - E14) & P2 > 400) %>% 
  select(gene:P2) %>% 
  gather(key = "tp", value = "expression", E11:P2) %>%
  ggplot(aes(x = factor(tp, levels = c(colnames(gene_normalized_clean[,2:7]))), y = expression, group = gene)) +#, color = gene)) +# %in% c("Mitf"))) + 
  geom_line(size = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Somewhat the opposite of previous graph. E14 is lower than E11 with at least 50 counts at E11

gene_normalized_clean %>% 
  filter(E14 <= E11*0.5& E11 > 50) %>% 
  select(gene:P2) %>%
  gather(key = "tp", value = "expression", E11:P2) %>%
  ggplot(aes(x = factor(tp, levels = c(colnames(gene_normalized_clean[,2:7]))), y = expression, group = gene)) +#, color = gene)) + 
  geom_line(size = 1) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

####################

# Looking at straight counts of genes is difficult to compare across timepoints since 
# gene can increase in relative frequency while the tissue also grows.
# Calculate log2 fold change to observe relative changes between sequential timepoints

calculate_lg2_fc <- function(first_tp, second_tp) {
  col_length <- length(first_tp)
  fc_col = numeric(col_length)
  for (row in seq(1, col_length, by = 1)) {
    if (first_tp[row] == 0) {
      fc  <- second_tp[row]
    }
    else if (second_tp[row] == 0) {
      fc  <- -first_tp[row]
    }
    else {
      fc  <- (log2(second_tp[row]) - log2(first_tp[row]))
    }
    fc_col[row] = fc
  }
  return(fc_col)
}


gene_lg2_fc <- data.frame(E11_E12 = calculate_lg2_fc(gene_normalized_clean$E11, gene_normalized_clean$E12),
                          E12_E14 = calculate_lg2_fc(gene_normalized_clean$E12, gene_normalized_clean$E14),
                          E14_E16 = calculate_lg2_fc(gene_normalized_clean$E14, gene_normalized_clean$E16),
                          E16_P0 = calculate_lg2_fc(gene_normalized_clean$E16, gene_normalized_clean$P0),
                          P0_P2 = calculate_lg2_fc(gene_normalized_clean$P0, gene_normalized_clean$P2),
                          P2_P4 = calculate_lg2_fc(gene_normalized_clean$P2, gene_normalized_clean$P4),
                          P4_P6 = calculate_lg2_fc(gene_normalized_clean$P4, gene_normalized_clean$P6),
                          P6_P10 = calculate_lg2_fc(gene_normalized_clean$P6, gene_normalized_clean$P10),
                          P10_P14 = calculate_lg2_fc(gene_normalized_clean$P10, gene_normalized_clean$P14),
                          P14_P21 = calculate_lg2_fc(gene_normalized_clean$P14, gene_normalized_clean$P21),
                          P21_P28 = calculate_lg2_fc(gene_normalized_clean$P21, gene_normalized_clean$P28))

####################

# Previous graphs suggest that there are sets of genes that are similar to each other
# and distinct from other genes.
# Use clustering to group similar genes.
# Agglomerative clustering has been used previously to group genes so start with that.

# Create similarity matrix (Euclidian Distance) for input to clustering algorithm.
# This function is slow to run. Not surprising for calculating similarity matrix
# but may be able to refactor (sub for-loops with apply/map)

l2 = length(gene_lg2_fc)
# Calculating similarity matrix for entire dataset will take days, run on subset.
rows = c(1:3000)
l1 = length(rows)
diff_mat_test = list(l1)
start_time <- Sys.time()
for (row in rows) {
  row_comp = 0
  diffs = numeric(l1)
  while (row+row_comp <= l1) {
    column = 1
    sum = 0
    while (column <= l2) {
      diff = gene_lg2_fc[row, column] - gene_lg2_fc[row+row_comp, column]
      sum = sum + diff**2
      column = column + 1
    } 
    diffs[row+row_comp] = sum
    row_comp = row_comp + 1
  }
  diff_mat_test[[row]] <- diffs
}
end_time <- Sys.time()
end_time-start_time # Calculate time it takes to run code chunk if interested
diff_mat_test <- do.call(rbind, diff_mat_test) # Combine list of similarity measures into matrix
rownames(diff_mat_test) <- gene_normalized$external_gene_name[1:3000] # gene names to rownames
colnames(diff_mat_test) <- gene_normalized$external_gene_name[1:3000] # gene names to colnames

# The way similarity matrix was calculated somewhat reverses conventional similarity measures.
# Add new calculation step to get it closer to true similarity measure

proportional_similarity <- function (row_name) {
  diff_mat_test[row_name, ]/rowSums(diff_mat_test)*-1
}

diff_mat_test_3 <- sapply(rownames(diff_mat_test), proportional_similarity)

# Similarity matrix is just one triangle of whole square matrix, copy values to empty triange

diff_mat_test_3[upper.tri(diff_mat_test_3)] <- t(diff_mat_test_3)[upper.tri(diff_mat_test_3)]

# Can save similarity matrix if desired.
#Recommended if going to be using the matrix for many downstream purposes

# save(diff_mat_test_3, file = "~/Documents/GIT/gene_analysis/diff_mat_test_3")

# apcluster is a package that performs agglomerative clustering on similarity matrix.
# Main function is 'apcluster'. Try it out of the box, accepting any default arg values.

ap_test <- apcluster(diff_mat_test_3)

# apcluster returns cluster numbers with rownames that belong to that cluster

ap_test

# Save rownames of reasonable sized clusters to own variables

cluster_72 <- names(ap_test@clusters[[72]]) # ok
cluster_57 <- names(ap_test@clusters[[57]]) # ok
cluster_55 <- names(ap_test@clusters[[55]]) # good
cluster_49 <- names(ap_test@clusters[[49]]) # good
cluster_43 <- names(ap_test@clusters[[43]]) # good
cluster_42 <- names(ap_test@clusters[[42]]) # good
cluster_41 <- names(ap_test@clusters[[41]]) # very good
cluster_40 <- names(ap_test@clusters[[40]]) # ok/good
cluster_37 <- names(ap_test@clusters[[37]]) # good
cluster_33 <- names(ap_test@clusters[[33]]) # ok
cluster_25 <- names(ap_test@clusters[[25]]) # ok
cluster_23 <- names(ap_test@clusters[[23]]) # good
cluster_19 <- names(ap_test@clusters[[19]]) # good
cluster_13 <- names(ap_test@clusters[[13]]) # good
cluster_7 <- names(ap_test@clusters[[7]]) # ok
cluster_4 <- names(ap_test@clusters[[4]]) # ok

# Create subset df that contains only the genes included in similarity matrix

gene_lg2_fc_3000 <- gene_lg2_fc[1:3000, ]

gene_lg2_fc_3000$gene <- gene_normalized$external_gene_name[1:3000]

# Plot one cluster to see how good the algorithm was. Remember to use log-fold change matrix
# Straight counts can work but grouping may vary in how tight it is.

gene_lg2_fc_3000 %>% 
  filter(gene %in% cluster_2_4) %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  ggplot(aes(x = factor(tp, levels = c(colnames(gene_lg2_fc_3000[,1:11]))), y = fold_change, group = gene)) + 
  geom_line() +
  theme_classic() +
  theme(legend.position="none") +
  labs(colour = "Gene", x = "Time Point", y = "Fold Change", title = "4")

# Create dfs that contain genes for good clusters

gene_lg2_fc_3000_72 <- gene_lg2_fc_3000 %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  filter(gene %in% cluster_2_72)

gene_lg2_fc_3000_41 <- gene_lg2_fc_3000 %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  filter(gene %in% cluster_2_41)

gene_lg2_fc_3000_33 <- gene_lg2_fc_3000 %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  filter(gene %in% cluster_2_33)

gene_lg2_fc_3000_23 <- gene_lg2_fc_3000 %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  filter(gene %in% cluster_2_23)

gene_lg2_fc_3000_19 <- gene_lg2_fc_3000 %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  filter(gene %in% cluster_2_19)

gene_lg2_fc_3000_7 <- gene_lg2_fc_3000 %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  filter(gene %in% cluster_2_7)

gene_lg2_fc_3000_4 <- gene_lg2_fc_3000 %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  filter(gene %in% cluster_2_4)

# Create df containing all remaining genes
# Not all genes are included in downstream graph so only remove genes in the clusters considered

gene_lg2_fc_3000_remainder <- gene_lg2_fc_3000 %>% 
  gather(key = "tp", value = "fold_change", 1:11) %>%
  filter(#!(gene %in% cluster_2_72)&
           !(gene %in% cluster_2_41) &
           #!(gene %in% cluster_2_33) &
           #!(gene %in% cluster_2_23) &
           !(gene %in% cluster_2_19) &
           #!(gene %in% cluster_2_7) &
           #!(gene %in% cluster_2_4) &
           (fold_change >= -5) &
           (fold_change <= 5))

gene_lg2_fc_3000_remainder %>% 
  ggplot(aes(x = factor(tp, levels = c(colnames(gene_lg2_fc_1000[,1:11]))), y = fold_change, group = gene)) + 
  geom_line(color = "gray", alpha = 0.25) +
  # This is a third cluster is want to observe. If included be sure to create new 'remainder' df
  # to also include genes from cluster 7
  # geom_line(data = gene_lg2_fc_3000_7,
  #           aes(x = factor(tp, levels = c(colnames(gene_lg2_fc_3000[,1:11]))),
  #               y = fold_change,
  #               group = gene),
  #           color = "#009E73",
  #           alpha = 0.5) +
  geom_line(data = gene_lg2_fc_3000_19,
            aes(x = factor(tp, levels = c(colnames(gene_lg2_fc_3000[,1:11]))),
                y = fold_change,
                group = gene),
            color = "#D55E00",
            alpha = 0.4) +
  geom_line(data = gene_lg2_fc_3000_41,
            aes(x = factor(tp, levels = c(colnames(gene_lg2_fc_3000[,1:11]))),
                y = fold_change,
                group = gene),
            color = "#0072B2",
            alpha = 0.4) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(colour = "grey38"),
        axis.title.y = element_text(colour = "grey38")) +
  labs(x = "Time", y = "Fold Change")

# Clustering seems to have performed relatively well.
# Of the two clusters considered both have distinct expression patterns
# and each cluster is tightly grouped for the most part.
# This can be used for downstream analysis to identify genes that fit a particular pattern
# or that are similar to another gene(s)
