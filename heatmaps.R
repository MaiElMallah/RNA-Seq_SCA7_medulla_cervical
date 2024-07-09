# Heatmap 1



# for annotation colors
ann_colors = list(
  Age = c("9w" = "#7570B3", "5w" = "#E7298A", "p" = "#66A61E"),
  Tissue = c("csc" = "#1B9E77", "med" = "#D95F02"),
  Genotype = c("WT" = "yellow", "SCA7" = "green")
)

select <- order(rowMeans(counts(dds, normalized=TRUE)),
                decreasing=TRUE)[1:40] #picks genes based on counts
df <- as.data.frame(colData(dds)[c("Genotype", "Tissue", "Age")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames = TRUE,
         show_colnames = TRUE, cluster_cols=FALSE, annotation_col = df,
         color = colorRampPalette(c("red", "white", "blue"))(100),
         annotation_colors = ann_colors,
         main = "Heatmap of count matrix")


# Heatmap 2


# for annotation colors
ann_colors = list(
  Age = c("9w" = "#7570B3", "5w" = "#E7298A", "p" = "#66A61E"),
  Tissue = c("csc" = "#1B9E77", "med" = "#D95F02"),
  Genotype = c("WT" = "yellow", "SCA7" = "green")
)

select <- order(res$padj)[1:40] #picks genes based on adjusted p-value
df <- as.data.frame(colData(dds)[c("Genotype", "Tissue", "Age")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames = TRUE,
         show_colnames = TRUE, cluster_cols=FALSE, annotation_col = df,
         color = colorRampPalette(c("red", "white", "blue"))(1000),
         annotation_colors = ann_colors,
         main = "Heatmap of normalized counts")



# Heatmap 3



#for annotation colors
ann_colors = list(
  Age = c("9w" = "#7570B3", "5w" = "#E7298A", "p" = "#66A61E"),
  Tissue = c("csc" = "#1B9E77", "med" = "#D95F02"),
  Genotype = c("WT" = "yellow", "SCA7" = "green")
)

df <- as.data.frame(colData(dds)[c("Genotype", "Tissue", "Age")])

select <- order(res$padj)[1:40] #picks genes based on adjusted padj
mat <- assay(vsd)[select,]
zscore <- (mat - rowMeans(mat)) / rowSds(mat)
pheatmap(zscore, cluster_rows=FALSE, show_rownames = TRUE,
         show_colnames = TRUE, cluster_cols=FALSE, annotation_col = df,
         color = colorRampPalette(c("red", "white", "blue"))(100),
         annotation_colors = ann_colors,
         main = "Heatmap of z-score")



# Heatmap 4

#for annotation colors
ann_colors = list(
  Age = c("9w" = "#7570B3", "5w" = "#E7298A", "p" = "#66A61E"),
  Tissue = c("csc" = "#1B9E77", "med" = "#D95F02"),
  Genotype = c("WT" = "yellow", "SCA7" = "green")
)

df <- as.data.frame(colData(med_dds)[c("Genotype", "Age")])

select <- order(med_res$padj)[1:40] #picks genes based on adjusted padj
mat <- assay(med_vsd)[select,]
zscore <- (mat - rowMeans(mat)) / rowSds(mat)
pheatmap(zscore, cluster_rows=FALSE, show_rownames = TRUE,
         show_colnames = TRUE, cluster_cols=FALSE, annotation_col = df,
         color = colorRampPalette(c("red", "white", "blue"))(100),
         annotation_colors = ann_colors,
         main = "Heatmap of z-score (med)")


# Heatmap 5



# for annotation colors
ann_colors = list(
  Age = c("9w" = "#7570B3", "5w" = "#E7298A", "p" = "#66A61E"),
  Tissue = c("csc" = "#1B9E77", "med" = "#D95F02"),
  Genotype = c("WT" = "yellow", "SCA7" = "green")
)

df <- as.data.frame(colData(csc_dds)[c("Genotype", "Age")])

select <- order(csc_res$padj)[1:40] #picks genes based on adjusted padj
mat <- assay(csc_vsd)[select,]
zscore <- (mat - rowMeans(mat)) / rowSds(mat)
pheatmap(zscore, cluster_rows=FALSE, show_rownames = TRUE,
         show_colnames = TRUE, cluster_cols=FALSE, annotation_col = df,
         color = colorRampPalette(c("red", "white", "blue"))(100),
         annotation_colors = ann_colors,
         main = "Heatmap of z-score (csc)")

# Heatmap 6


#for annotation colors
ann_colors = list(
  Age = c("9w" = "#7570B3", "5w" = "#E7298A", "p" = "#66A61E"),
  Tissue = c("csc" = "#1B9E77", "med" = "#D95F02"),
  Genotype = c("WT" = "yellow", "SCA7" = "green")
)

df <- as.data.frame(colData(med_dds)[c("Genotype", "Age")])

select <- order(med_res$padj)[1:100] #picks genes based on adjusted padj
mat <- assay(med_vsd)[select,]
zscore <- (mat - rowMeans(mat)) / rowSds(mat)
pheatmap(zscore, cluster_rows=FALSE, show_rownames = TRUE,
         show_colnames = TRUE, cluster_cols=FALSE, annotation_col = df,
         color = colorRampPalette(c("red", "white", "blue"))(100),
         annotation_colors = ann_colors,
         main = "Heatmap of z-score (med)")

# Heatmap 7


#for annotation colors
ann_colors = list(
  Age = c("9w" = "#7570B3", "5w" = "#E7298A", "p" = "#66A61E"),
  Tissue = c("csc" = "#1B9E77", "med" = "#D95F02"),
  Genotype = c("WT" = "yellow", "SCA7" = "green")
)

df <- as.data.frame(colData(csc_dds)[c("Genotype", "Age")])

select <- order(csc_res$padj)[1:100] #picks genes based on adjusted padj
mat <- assay(csc_vsd)[select,]
zscore <- (mat - rowMeans(mat)) / rowSds(mat)
pheatmap(zscore, cluster_rows=FALSE, show_rownames = TRUE,
         show_colnames = TRUE, cluster_cols=FALSE, annotation_col = df,
         color = colorRampPalette(c("red", "white", "blue"))(100),
         annotation_colors = ann_colors,
         main = "Heatmap of z-score (csc)")
