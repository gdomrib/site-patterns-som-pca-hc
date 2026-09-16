"""

CODE FOR THE PAPER: Uncovering Spatial Patterns in Archaeological Data through the Integration of Self-Organizing Maps, Principal Component Analysis and Hierarchical Clustering

Developed by Guillem Domingo-Ribas

This code is divided in three parts, corresponding to the approaches explored in the paper: principal component analysis, hierarchical clustering and self-organising maps

Development of this code is by the author, and it has only been refined with the assistance of AI tools (Claude).

"""


### PRINCIPAL COMPONENT ANALYSIS ### 

library(readr)
library(ggplot2)
library(ggfortify)
library(FactoMineR)
library(factoextra)

setwd("your_directory")

# 1. Load your data
aa_samp <- read.csv("yourSamplingFile.csv", encoding = "UTF-8")

# 2. Select and rename variables - Change according to your variable names
vars_keep <- c("HubDist", "dem05_med", "current_med", "aspect05_med",
               "slope05_med", "distRiver_med", "visInOut_med", "visInIng_med",
               "prom10k_med", "rough_med", "globRad_med")
  
var_labels <- c(
  "HubDist"       = "Distance nearest site",
  "dem05_med"     = "Altitude",
  "current_med"   = "Connectivity",
  "aspect05_med"  = "Aspect",
  "slope05_med"   = "Slope",
  "distRiver_med" = "Distance nearest river",
  "visInOut_med"  = "Visibility index (out)",
  "visInIng_med"  = "Visibility index (in)",
  "prom10k_med"   = "Prominence",
  "rough_med"     = "Roughness",
  "globRad_med"   = "Solar radiance"
)

# Subset to selected variables and rename columns
pca_data <- aa_samp[, vars_keep]
colnames(pca_data) <- var_labels[vars_keep]

# 3. Run PCA
pca_result <- prcomp(pca_data, scale. = TRUE)

# 4. Shared theme for publication-ready figures 
pub_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 13, hjust = 0.5),
    axis.title   = element_text(size = 11),
    axis.text    = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    plot.margin  = margin(10, 15, 10, 10)
  )

# 5. Scree plot - FIGURE 3
p1 <- fviz_eig(pca_result,
               addlabels  = TRUE,
               ylim       = c(0, 50),
               barfill    = "steelblue",
               barcolor   = "steelblue",
               linecolor  = "grey30") +
  labs(title = "Scree Plot",
       x     = "Principal Component",
       y     = "Variance Explained (%)") +
  pub_theme

ggsave("figs/pca/fig3_pca_scree.png", p1, width = 6,  height = 4,   dpi = 300)
print(p1)
dev.off()

# 6. Variable loadings plot (variable circle) 
p2 <- fviz_pca_var(pca_result,
                   repel       = TRUE,
                   col.var     = "contrib",          # colour by contribution
                   gradient.cols = c("#4393c3", "#f7f7f7", "#d6604d"),
                   legend.title = "Contribution (%)") +
  labs(title = "PCA – Variable Loadings",
       x     = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y     = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  pub_theme

print(p2)

# 7. Biplot (sites + variables) - FIGURE 4
p3 <- fviz_pca_biplot(pca_result,
                      repel        = TRUE,
                      col.var      = "#d6604d",
                      col.ind      = "grey40",
                      alpha.ind    = 0.6,
                      label        = "var",           # label variables only
                      pointsize    = 1.5,
                      arrowsize    = 0.6) +
  labs(title = "PCA Biplot – Sites and Variable Loadings",
       x     = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y     = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  pub_theme

ggsave("figs/pca/fig4_pca_biplot.png", p3, width = 7,  height = 6,   dpi = 300)
print(p3 + theme(panel.grid = element_blank()) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60")
)
dev.off()


# 8. Variable contributions to PC1 and PC2 (side-by-side) 
p4a <- fviz_contrib(pca_result, choice = "var", axes = 1) +
  labs(title = "Contributions to PC1", x = NULL, y = "Contribution (%)") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

p4b <- fviz_contrib(pca_result, choice = "var", axes = 2) +
  labs(title = "Contributions to PC2", x = NULL, y = "Contribution (%)") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

# Combine side by side (requires patchwork)
# install.packages("patchwork")  # run once if needed
library(patchwork)
p4 <- p4a + p4b +
  plot_annotation(title = "Variable Contributions by Principal Component",
                  theme = theme(plot.title = element_text(face = "bold",
                                                          hjust = 0.5, size = 13)))
print(p4)

# 9. Individuals plot coloured by cos2 (quality of representation)
p5 <- fviz_pca_ind(pca_result,
                   geom.ind     = "point",
                   col.ind      = "cos2",
                   gradient.cols = c("#4393c3", "#fee090", "#d73027"),
                   alpha.ind    = 0.7,
                   pointsize    = 2,
                   repel        = FALSE,
                   legend.title = "cos²") +
  labs(title = "PCA – Site Distribution (PC1–PC2)",
       x     = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
       y     = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")) +
  pub_theme

print(p5)


#############################################################################################
#############################################################################################


### HIERARCHICAL CLUSTERING ### 

library(dplyr)
library(dendextend)
library(ggdendro)
library(patchwork)

# 1. Prepare data

aa_data_wName <- aa_samp[, c("id", "name", vars_keep)]

# 2. Scale and cluster
data_scaled  <- scale(aa_data_wName[, vars_keep])
rownames(data_scaled) <- aa_data_wName$name   # attach names to rows

dist_matrix  <- dist(data_scaled, method = "euclidean")
hc           <- hclust(dist_matrix, method = "ward.D2")


# 3. Choose number of clusters - FIGURE 5
wss <- sapply(1:10, function(k) {
  cl  <- cutree(hc, k = k)
  sum(sapply(unique(cl), function(c) {
    members <- data_scaled[cl == c, , drop = FALSE]
    if (nrow(members) < 2) return(0)
    sum(dist(members)^2) / nrow(members)
  }))
})

elbow_df <- data.frame(k = 1:10, wss = wss)

p_elbow <- ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_line(colour = "grey40", linewidth = 0.8) +
  geom_point(size = 3, colour = "steelblue") +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Elbow Plot – Optimal Number of Clusters",
       x = "Number of Clusters (k)", y = "Within-cluster Sum of Squares") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("figs/dend/fig5_hc_elbow.png", p_elbow, width = 6, height = 4, dpi = 300)
print(p_elbow)
dev.off()


# 4. Define cluster palette (k = 4) - FIGURE 6
k          <- 4
clust_cols <- c("#4393c3", "#d6604d", "#4dac26", "#8073ac")  
# blue, red, green, purple — colourblind-safer, prints well in greyscale

# ── 5. Publication dendrogram (base R + dendextend) ───────────────────────────
dend <- as.dendrogram(hc)

dend <- dend %>%
  set("branches_k_color", k = k, value = clust_cols) %>%
  set("branches_lwd", 2.5) %>%
  set("labels_cex", 0.55) %>%
  set("labels_col", "grey20")

# Save as high-res PNG
png("figs/dend/fig6_hc_dendrogram.png", width = 3000, height = 1800, res = 300)

# Larger bottom margin to prevent label clipping, extra space for cluster labels
par(mar = c(8, 4.5, 3, 1))

plot(dend,
     main     = "Hierarchical Clustering of Early Medieval Sites\n(Ward's method, Euclidean distance)",
     ylab     = "Height",
     xlab     = "",
     cex.main = 1.0,
     cex.axis = 0.8,
     cex.lab  = 0.9,
     ylim     = c(-1, max(hc$height) + 1))  # extend y-axis downward to make room

# ── Cluster labels only, no rectangles ───────────────────────────────────────
leaf_order     <- order.dendrogram(dend)
cluster_assign <- cutree(hc, k = k)
cluster_by_pos <- cluster_assign[leaf_order]

n       <- length(leaf_order)
label_x <- tapply(seq_len(n), cluster_by_pos, median)

clusters_lr <- as.integer(names(sort(label_x)))
x_pos_lr    <- sort(
  label_x)

mtext(paste("Cluster", clusters_lr),
      side = 1,
      line = 7,
      at   = as.numeric(x_pos_lr),
      col  = clust_cols,
      cex  = 0.9,
      font = 2)



dev.off()


# 6. Cluster profile heatmap
# Shows what characterises each cluster — very useful for interpretation
cluster_id <- cutree(hc, k = k)

profile_df <- as.data.frame(data_scaled) %>%
  mutate(cluster = factor(cluster_id, labels = paste("Cluster", 1:k))) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean)) %>%
  tidyr::pivot_longer(-cluster, names_to = "variable", values_to = "mean_z")

# Apply readable variable labels
profile_df$variable <- var_labels[profile_df$variable]

p_heat <- ggplot(profile_df, aes(x = variable, y = cluster, fill = mean_z)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = round(mean_z, 2)),
            size = 3, colour = "grey10") +
  scale_fill_gradient2(low = "#4393c3", mid = "white", high = "#d6604d",
                       midpoint = 0,
                       name = "Mean\nz-score") +
  labs(title = "Cluster Profiles – Mean Standardised Values",
       x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 12),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  )

print(p_heat)


# 7. Cluster size summary (quick check)
cat("\nCluster sizes:\n")
print(table(cluster_id))


#############################################################################################
#############################################################################################


### SELF-ORGANISING MAPS ###

# ══════════════════════════════════════════════════════════════════════════════

library(kohonen)
library(dplyr)
library(RColorBrewer)
library(tidyr)


# 1. Prepare data

aa_data <- aa_samp[, c("id", "name", vars_keep)]

# Scale only the numeric variables
aa_data_matrix <- as.matrix(scale(aa_data[, vars_keep]))
rownames(aa_data_matrix) <- aa_data$name

# 2. Shared colour palette
# Blue–white–red diverging: colourblind-safe, prints in greyscale
coolBlueHotRed <- function(n, alpha = 1) {
  colorRampPalette(
    c("#313695", "#4575b4", "#74add1", "#abd9e9",
      "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027")
  )(n)
}

# Cluster colours (consistent with PCA and dendrogram)
clust_cols <- c("#4393c3", "#d6604d", "#4dac26", "#8073ac", "#e08214", "#01665e")

pub_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", size = 11, hjust = 0.5),
    axis.title   = element_text(size = 10),
    axis.text    = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.margin  = margin(8, 10, 8, 8)
  )

# 3. Train SOM
set.seed(2026)

samp_grid <- somgrid(xdim = 4, ydim = 5,
                     topo = "hexagonal",
                     neighbourhood.fct = "gaussian")

samp_SOM_model <- som(aa_data_matrix,
                      grid     = samp_grid,
                      rlen     = 400,
                      alpha    = c(0.08, 0.01),
                      keep.data = TRUE)


# 4. Training progress - FIGURE 7
png("figs/som/fig7_som_training.png", width = 1800, height = 1200, res = 300)
par(mar = c(4, 4, 3, 1), cex.main = 1, cex.lab = 0.9, cex.axis = 0.8)
plot(samp_SOM_model, type = "changes",
     main = "SOM Training Progress")
dev.off()


# 5. Node counts - FIGURE 8
png("figs/som/fig8_som_counts.png", width = 1800, height = 2000, res = 300)
par(mar = c(1, 1, 3, 4), cex.main = 1)
plot(samp_SOM_model,
     type         = "counts",
     palette.name = coolBlueHotRed,
     main         = "SOM Node Counts – Sites per Node")
#add.cluster.boundaries(samp_SOM_model, som_cluster, lwd = 2.5, col = "firebrick")
dev.off()


# 6. Determine optimal number of SOM clusters (elbow) - FIGURE 9
# Use codebook vectors (there is one per node) as the clustering input
codes <- getCodes(samp_SOM_model)

wss <- sapply(1:10, function(k) {
  km <- kmeans(codes, centers = k, nstart = 25)
  km$tot.withinss
})

elbow_df <- data.frame(k = 1:10, wss = wss)

p_elbow <- ggplot(elbow_df, aes(x = k, y = wss)) +
  geom_line(colour = "grey40", linewidth = 0.8) +
  geom_point(size = 3, colour = "steelblue") +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Elbow Plot – SOM Node Clusters",
       x = "Number of Clusters (k)",
       y = "Within-cluster Sum of Squares") +
  pub_theme

ggsave("figs/som/fig9_som_elbow.png", p_elbow, width = 5, height = 3.5, dpi = 300)
print(p_elbow)

# 7. Assign clusters to SOM nodes (k = 4, consistent with dendrogram)
k_som       <- 4
som_hc      <- hclust(dist(codes), method = "ward.D2")
som_cluster <- cutree(som_hc, k = k_som)

# Cluster membership for each SITE (via its winning node)
site_cluster <- som_cluster[samp_SOM_model$unit.classif]
names(site_cluster) <- aa_data$name

# 8. Combined figure: U-Matrix + Cluster Map - FIGURE 10
png("figs/som/fig10_som_umatrix_clusters_combined.png", width = 3800, height = 2200, res = 300)

layout(matrix(c(1, 2), nrow = 1, ncol = 2), widths = c(1, 1))

# ── Left panel: U-Matrix ──────────────────────────────────────────────────────
par(mar = c(1, 1, 3, 4))
plot(samp_SOM_model,
     type         = "dist.neighbours",
     palette.name = grey.colors,
     main         = "A)  U-Matrix – Node Distances")
#add.cluster.boundaries(samp_SOM_model, som_cluster, lwd = 3, col = "firebrick")

# ── Right panel: Cluster Map ──────────────────────────────────────────────────
par(mar = c(1, 1, 1.8, 4), cex.main = 1)
plot(samp_SOM_model,
     type   = "mapping",
     col    = clust_cols[site_cluster],
     pch    = 19,
     #cex    = 0.9,
     bgcol  = clust_cols[som_cluster],
     main   = "B)  SOM Cluster Map")
#add.cluster.boundaries(samp_SOM_model, som_cluster, lwd = 3, col = "white")
legend("bottom",
       legend = paste("Cluster", 1:k_som),
       fill   = clust_cols[1:k_som],
       border = NA,
       bty    = "n",
       cex    = 0.85)

dev.off()



# 9. Variable heatmaps (all 11 variables, unscaled) - FIGURE 11
# Loop replaces the repetitive block-copy approach in the original code
png("figs/som/fig11_som_heatmaps_unscaled.png", width = 3600, height = 2800, res = 300)
par(mfrow = c(3, 4),
    mar   = c(1, 1, 2.5, 3),
    oma   = c(0, 0, 3, 0))

for (i in seq_along(vars_keep)) {
  var_col   <- vars_keep[i]
  orig_vals <- aa_samp[[var_col]]           # unscaled original column
  
  var_unscaled <- aggregate(orig_vals,
                            by  = list(samp_SOM_model$unit.classif),
                            FUN = mean,
                            simplify = TRUE)[, 2]
  
  plot(samp_SOM_model,
       type         = "property",
       property     = var_unscaled,
       main         = var_labels[var_col],
       palette.name = coolBlueHotRed,
       cex.main     = 0.9)
  
  add.cluster.boundaries(samp_SOM_model, som_cluster, lwd = 2, col = "black")
}

mtext("SOM Variable Heatmaps – Unscaled Mean Values per Node",
      outer = TRUE, cex = 1, font = 2, line = 1)
dev.off()



# 10. Variable heatmaps (scaled / codebook values)
par(mfrow = c(3, 4),
    mar   = c(1, 1, 2.5, 3),
    oma   = c(0, 0, 3, 0))

for (i in seq_along(vars_keep)) {
  plot(samp_SOM_model,
       type         = "property",
       property     = getCodes(samp_SOM_model)[, i],
       main         = var_labels[vars_keep[i]],
       palette.name = coolBlueHotRed,
       cex.main     = 0.9)
  
  add.cluster.boundaries(samp_SOM_model, som_cluster, lwd = 2, col = "black")
}

mtext("SOM Variable Heatmaps – Scaled Codebook Values per Node",
      outer = TRUE, cex = 1, font = 2, line = 1)
dev.off()


# 11. Codebook fan / radar plot - FIGURE 12
png("figs/som/fig12_som_codes.png", width = 2400, height = 2800, res = 300)
par(mar = c(0.5, 0.5, 2, 0.5), cex.main = 1)

# Temporarily rename codebook columns for the legend
colnames(samp_SOM_model$codes[[1]]) <- var_labels[vars_keep]

plot(samp_SOM_model,
     type = "codes",
     main = "SOM Codebook Vectors")
add.cluster.boundaries(samp_SOM_model, som_cluster, lwd = 2.5, col = "firebrick")

# Restore original column names so nothing else in the script is affected
colnames(samp_SOM_model$codes[[1]]) <- vars_keep
dev.off()



# 12. Cluster profile heatmap (with ggplot2)
# Shows what characterises each SOM cluster — key for archaeological interpretation
codes_df <- as.data.frame(codes)
colnames(codes_df) <- var_labels[vars_keep]
codes_df$cluster   <- factor(som_cluster, labels = paste("Cluster", 1:k_som))

profile_df <- codes_df %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(-cluster, names_to = "variable", values_to = "mean_z")

p_profile <- ggplot(profile_df, aes(x = variable, y = cluster, fill = mean_z)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = round(mean_z, 2)), size = 2.8, colour = "grey15") +
  scale_fill_gradient2(low      = "#4393c3",
                       mid      = "white",
                       high     = "#d6604d",
                       midpoint = 0,
                       name     = "Mean\nz-score") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title = "SOM Cluster Profiles – Mean Standardised Codebook Values",
       x = NULL, y = NULL) +
  pub_theme +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        panel.border = element_rect(colour = "grey80", fill = NA))

print(p_profile)



# 13. SOM vs HC cluster agreement - FIGURE 13
# Demonstrates added value of SOM: compare cluster assignments with the dendrogram
hc_clusters  <- cutree(hclust(dist(aa_data_matrix), method = "ward.D2"), k = k_som)
som_clusters <- site_cluster

comparison_df <- data.frame(
  site       = aa_data$name,
  HC_cluster = factor(hc_clusters,  labels = paste("HC Cluster",  1:k_som)),
  SOM_cluster= factor(som_clusters, labels = paste("SOM Cluster", 1:k_som))
)

# Confusion-style tile plot
agree_df <- comparison_df %>%
  count(HC_cluster, SOM_cluster)

p_agree <- ggplot(agree_df, aes(x = SOM_cluster, y = HC_cluster, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = n), size = 4, fontface = "bold") +
  scale_fill_gradient(low = "#f7f7f7", high = "#4393c3", name = "Sites (n)") +
  labs(title = "Cluster Agreement – Hierarchical Clustering vs SOM",
       subtitle = "Cell values = number of sites assigned to each cluster combination",
       x = "SOM Cluster", y = "HC Cluster") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("figs/som/fig13_som_hc_agreement.png", p_agree, width = 6, height = 5, dpi = 300)
print(p_agree)



# 14. PCA space coloured by SOM cluster - FIGURE 14
# The key "bridge" figure: shows SOM clusters in PCA space,
# demonstrating that SOM captures structure the PCA biplot hints at
pca_result <- prcomp(aa_data[, vars_keep], scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$site        <- aa_data$name
pca_scores$SOM_cluster <- factor(som_clusters, labels = paste("Cluster", 1:k_som))
pca_scores$HC_cluster  <- factor(hc_clusters,  labels = paste("Cluster", 1:k_som))

pct_var <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)

p_pca_som <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = SOM_cluster)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(group = SOM_cluster), level = 0.75,
               linewidth = 0.6, linetype = "dashed") +
  scale_colour_manual(values = clust_cols[1:k_som], name = NULL) +
  labs(title = "PCA Space – Sites Coloured by SOM Cluster",
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)")) +
  pub_theme +
  theme(legend.position = "right")

p_pca_hc <- ggplot(pca_scores, aes(x = PC1, y = PC2, colour = HC_cluster)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(aes(group = HC_cluster), level = 0.75,
               linewidth = 0.6, linetype = "dashed") +
  scale_colour_manual(values = clust_cols[1:k_som], name = NULL) +
  labs(title = "PCA Space – Sites Coloured by HC Cluster",
       x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)")) +
  pub_theme +
  theme(legend.position = "right")

p_pca_compare <- p_pca_som + p_pca_hc +
  plot_annotation(
    title = "Cluster Assignments in PCA Space: SOM vs Hierarchical Clustering",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
  )

ggsave("figs/som/fig14_som_pca_comparison.png", p_pca_compare, width = 12, height = 5, dpi = 300)
print(p_pca_compare)

# 16. Print cluster membership summary
cat("\n── SOM cluster sizes ──\n")
print(table(site_cluster))

cat("\n── Sites per SOM cluster ──\n")
for (cl in sort(unique(site_cluster))) {
  cat(sprintf("\nCluster %d:\n", cl))
  cat(paste(" ", aa_data$name[site_cluster == cl], collapse = "\n"), "\n")
}



# EXPORT FOR SPATIAL EXPLORATION OUTSIDE R — sites with PCA / HC / SOM results attached
# Employed for FIGURE 15 (a, b, c)
library(sf)

# ── Recover coordinates from the raw sampling table ───────────────────────────
# aa_data dropped x/y, so pull them from aa_samp. y imported as character in
# some rows, so coerce and strip any stray non-numeric characters.
x_num <- as.numeric(aa_samp$x)
y_num <- as.numeric(gsub("[^0-9.\\-]", "", trimws(aa_samp$y)))

# ── Assemble one attribute table ──────────────────────────────────────────────
# All results are in aa_samp / aa_data row order, so binding by position is safe.
export_df <- data.frame(
  id          = aa_samp$id,
  name        = aa_samp$name,
  x           = x_num,
  y           = y_num,
  PC1         = pca_result$x[, 1],   # continuous gradient (dominant axis)
  PC2         = pca_result$x[, 2],
  HC_clust    = hc_clusters,         # integer 1–4  (hierarchical clustering)
  SOM_clust   = som_clusters,        # integer 1–4  (SOM)
  SOM_node    = samp_SOM_model$unit.classif
)

# ── Guard: no missing coordinates allowed ─────────────────────────────────────
bad <- which(is.na(export_df$x) | is.na(export_df$y))
if (length(bad) > 0) {
  cat("Rows with missing coordinates (fix before export):\n")
  print(export_df[bad, c("id", "name", "x", "y")])
  stop("Missing coordinates — see rows above.")
}

# ── Build sf object and write ─────────────────────────────────────────────────
# CRS 25830 = ETRS89 / UTM 30N (Burgos). Change if your project differs.
sites_sf <- st_as_sf(export_df, coords = c("x", "y"), crs = 25830, remove = FALSE)

st_write(sites_sf, "figs/5_comparative_map/sites_classified.gpkg",
         layer = "sites", delete_dsn = TRUE)

cat("\nExported", nrow(sites_sf), "sites to figs/5_comparative_map/sites_classified.gpkg\n")
