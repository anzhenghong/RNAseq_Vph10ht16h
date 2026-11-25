# ============================================================
# WGCNA pipeline
#   Network: signed, biweight midcorrelation (bicor)
#   kME: bicor
#   Module-trait: bicor + BH (ggplot heatmap, RdBu)
#   Hub: kME >= 0.8 (positive only, signed network)
#   Edge: intramodular TOM >= 0.1 (module_edges.txt AND Cytoscape export)
#   Sensitivity (OFAT): power even 18-30; deepSplit 1-4;
#                       minModuleSize 20/30/50;
#                       mergeCutHeight 0.10/0.15/0.20/0.25/0.30
#   Stability: bootstrap Jaccard overlap (100 iterations)
# ============================================================

rm(list = ls())   # uncomment only in a fresh R session
##BiocManager::install()  ###software installation
set.seed(12345)

library(WGCNA)
library(openxlsx)
library(tidyverse)
library(ggtext)

enableWGCNAThreads(nThreads = 4)

# ============================================================
# 0. Global parameters
# ============================================================
maxPOutliers    <- 0.05
main_dir        <- "."
if (!dir.exists(main_dir)) dir.create(main_dir, recursive = TRUE)

networkType     <- "signed"
TOMType         <- "signed"
corType         <- "bicor"
deepSplit       <- 1
minModuleSize   <- 30
mergeCutHeight  <- 0.25

kME_hub_cutoff  <- 0.8
edge_tom_cutoff <- 0.1   # used for module_edges.txt AND Cytoscape export

run_sensitivity <- TRUE
run_stability   <- TRUE
n_bootstrap     <- 100

# ============================================================
# 1. Load data
# ============================================================
setwd("D:/third/Vph10t16/stringtie/Recent_results/WGCNA")
dataExpr <- read.xlsx("./Vph10T16_All_DEG_TPM.xlsx", rowNames = TRUE)
dataExpr <- log2(dataExpr + 1)
datExpr0 <- as.data.frame(t(dataExpr))

stopifnot(!any(duplicated(colnames(datExpr0))))

gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
    datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}
nSamples <- nrow(datExpr0)
nGenes   <- ncol(datExpr0)
cat("Samples:", nSamples, " Genes:", nGenes, "\n")

write.table(data.frame(Sample = rownames(datExpr0)),
            file.path(main_dir, "sample_list.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ============================================================
# 2. Soft-threshold diagnostics (signed + bicor)
# ============================================================
powers <- c(1:10, seq(from = 12, to = 30, by = 2))
rsquare.cut <- 0.85

sft <- pickSoftThreshold(
    datExpr0,
    powerVector = powers,
    networkType = networkType,
    corFnc      = corType,
    corOptions  = list(maxPOutliers = maxPOutliers),
    RsquaredCut = rsquare.cut,
    removeFirst = FALSE,
    verbose     = 5
)

sft_table <- data.frame(
    Power      = sft$fitIndices[, 1],
    Raw_SFT_R2 = sft$fitIndices[, 2],
    Slope      = sft$fitIndices[, 3],
    Signed_R2  = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
    Mean_k     = sft$fitIndices[, 5],
    Median_k   = sft$fitIndices[, 6],
    Max_k      = sft$fitIndices[, 7]
)
write.table(sft_table,
            file.path(main_dir, "softThreshold_diagnostics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

softPower <- sft$powerEstimate
if (is.na(softPower)) {
    idx <- which.max(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2])
    softPower <- sft$fitIndices$Power[idx]
}
cat("softPower =", softPower, "\n")

sft_df <- sft$fitIndices %>%
    mutate(signed_R2 = -sign(slope) * SFT.R.sq, Power = Power)
beta_choose <- softPower
k_val <- sft_df %>% filter(Power == beta_choose) %>% pull(mean.k.)

pA <- ggplot(sft_df, aes(x = Power)) +
    geom_line(aes(y = signed_R2), color = "#92c5de", linewidth = 1.2) +
    geom_point(aes(y = signed_R2), color = "#1f78b4", size = 3) +
    geom_hline(yintercept = rsquare.cut, linetype = 2, color = "gray50") +
    geom_vline(xintercept = beta_choose, linetype = 2, color = "#e66101") +
    geom_point(data = filter(sft_df, Power == beta_choose),
               aes(y = signed_R2), color = "#e66101", size = 6) +
    annotate("text", x = beta_choose + 0.8, y = 0.2,
             label = as.expression(bquote(beta == .(beta_choose))),
             color = "#e66101", size = 6) +
    annotate("text", x = 12, y = rsquare.cut,
             label = as.expression(bquote(R^2 == .(rsquare.cut))),
             color = "gray40", size = 5) +
    labs(x = "Soft-thresholding power (beta)",
         y = "Scale-free topology fit, signed R2") +
    scale_x_continuous(breaks = seq(0, max(powers), 5),
                       limits = c(0, max(powers))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0, size = 18),
          panel.grid = element_blank()) +
    ggtitle("(A)")

pB <- ggplot(sft_df, aes(x = Power)) +
    geom_line(aes(y = mean.k.), color = "#92c5de", linewidth = 1.2) +
    geom_point(aes(y = mean.k.), color = "#1f78b4", size = 3) +
    geom_vline(xintercept = beta_choose, linetype = 2, color = "#e66101") +
    geom_point(data = filter(sft_df, Power == beta_choose),
               aes(y = mean.k.), color = "#e66101", size = 6) +
    annotate("text", x = beta_choose + 0.8, y = k_val * 1.8,
             label = as.expression(bquote(beta == .(beta_choose))),
             color = "#e66101", size = 6) +
    labs(x = "Soft-thresholding power (beta)", y = "Mean connectivity") +
    scale_x_continuous(breaks = seq(0, max(powers), 5),
                       limits = c(0, max(powers))) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0, size = 18),
          panel.grid = element_blank()) +
    ggtitle("(B)")

ggsave(file.path(main_dir, "sft_left.pdf"),  pA, width = 5, height = 5)
ggsave(file.path(main_dir, "sft_right.pdf"), pB, width = 5, height = 5)

# ============================================================
# 3. Adjacency and TOM (signed + bicor; dimnames explicitly set)
# ============================================================
adjacency_mat <- adjacency(
    datExpr0,
    power  = softPower,
    type   = networkType,
    corFnc = corType,
    corOptions = list(maxPOutliers = maxPOutliers)
)
dimnames(adjacency_mat) <- list(colnames(datExpr0), colnames(datExpr0))

TOM <- TOMsimilarity(adjacency_mat, TOMType = TOMType)
dimnames(TOM) <- list(colnames(datExpr0), colnames(datExpr0))
dissTOM <- 1 - TOM

stopifnot(identical(rownames(TOM), colnames(datExpr0)),
          identical(colnames(TOM), colnames(datExpr0)),
          identical(rownames(adjacency_mat), colnames(datExpr0)),
          identical(colnames(adjacency_mat), colnames(datExpr0)))

# ============================================================
# 4. Gene clustering and module detection
# ============================================================
geneTree <- hclust(as.dist(dissTOM), method = "average")
cat("Genes in tree:", length(geneTree$order), "\n")

set.seed(12345)
dynamicMods <- cutreeDynamic(
    dendro = geneTree,
    distM  = dissTOM,
    deepSplit = deepSplit,
    pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize,
    verbose = 3
)
dynamicColors <- labels2colors(dynamicMods)
cat("Initial modules (incl. grey):", length(unique(dynamicColors)), "\n")
print(table(dynamicColors))

# ============================================================
# 5. Module eigengenes and merging
# ============================================================
MEList <- moduleEigengenes(datExpr0, colors = dynamicColors, impute = TRUE)
MEs    <- MEList$eigengenes

mergeRes <- mergeCloseModules(
    datExpr0,
    dynamicColors,
    cutHeight = mergeCutHeight,
    corFnc = corType,
    corOptions = list(maxPOutliers = maxPOutliers),
    verbose = 3
)
mergedColors <- mergeRes$colors
mergedMEs    <- mergeRes$newMEs

MEDissMat <- 1 - bicor(MEs, maxPOutliers = maxPOutliers)
METree    <- hclust(as.dist(MEDissMat), method = "average")

pdf(file.path(main_dir, "step_ME_clustering.pdf"), width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = mergeCutHeight, col = "red")
dev.off()

cat("Merged modules (incl. grey):", length(unique(mergedColors)), "\n")
print(table(mergedColors))

pdf(file.path(main_dir, "step_geneTree_before_after.pdf"), width = 12, height = 8)
plotDendroAndColors(
    geneTree,
    cbind(dynamicColors, mergedColors),
    c("Dynamic Tree Cut", "Merged dynamic"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = "Gene dendrogram: before and after module merging"
)
dev.off()

# ============================================================
# 6. Gene-module table
# ============================================================
moduleLabels <- match(mergedColors, c("grey", standardColors(50))) - 1
gene_module_color <- data.frame(
    Gene = colnames(datExpr0),
    ModuleLabel = moduleLabels,
    ModuleColor = mergedColors,
    stringsAsFactors = FALSE
)
write.table(gene_module_color,
            file.path(main_dir, "step_gene_module_color.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# ============================================================
# 6b. Module size barplot
# ============================================================
module_size_df <- data.frame(
    Module  = names(table(mergedColors)),
    N_genes = as.integer(table(mergedColors)),
    stringsAsFactors = FALSE
)
write.table(module_size_df,
            file.path(main_dir, "module_size.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

grey_n <- module_size_df$N_genes[module_size_df$Module == "grey"]
if (length(grey_n) == 0) grey_n <- 0L

module_size_main <- module_size_df %>%
    filter(Module != "grey") %>%
    arrange(desc(N_genes))
module_size_main$Module <- factor(module_size_main$Module,
                                  levels = module_size_main$Module)

p_bar <- ggplot(module_size_main,
                aes(x = Module, y = N_genes,
                    fill = as.character(Module))) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    geom_text(aes(label = N_genes), vjust = -0.4, size = 3.5) +
    scale_fill_identity() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Module", y = "Number of genes",
         title = "Gene counts per module",
         subtitle = paste0("Grey module (unassigned): ", grey_n, " genes")) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.x        = element_text(angle = 45, hjust = 1),
        plot.title         = element_text(face = "bold", hjust = 0)
    )

ggsave(file.path(main_dir, "module_size_barplot.pdf"),
       p_bar, width = 7, height = 5)

# ============================================================
# 7. kME (bicor; two-sided P values)
# ============================================================
MEs_ordered <- orderMEs(mergedMEs)
bp_res      <- bicorAndPvalue(datExpr0, MEs_ordered, maxPOutliers = maxPOutliers)
kme_all_mat <- bp_res$bicor
p_all_mat   <- bp_res$p

gene_names <- rownames(kme_all_mat)
stopifnot(identical(gene_names, colnames(datExpr0)))

kme_df <- data.frame(
    Gene = gene_names,
    ModuleColor = mergedColors,
    kME = NA,
    kME_p.value = NA,
    stringsAsFactors = FALSE
)
for (i in seq_along(gene_names)) {
    mod_col <- mergedColors[i]
    me_name <- paste0("ME", mod_col)
    if (me_name %in% colnames(kme_all_mat)) {
        kme_df$kME[i]         <- kme_all_mat[i, me_name]
        kme_df$kME_p.value[i] <- p_all_mat[i, me_name]
    }
}
write.table(kme_df,
            file.path(main_dir, "step_All_Gene_KME.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

geneModuleMembership <- as.data.frame(kme_all_mat)
write.table(
    data.frame(Gene = rownames(geneModuleMembership),
               geneModuleMembership, check.names = FALSE),
    file.path(main_dir, "step_gene_kME_all_modules.txt"),
    sep = "\t", row.names = FALSE, quote = FALSE
)

# ============================================================
# 8. Module-trait associations (bicor + BH; ggplot heatmap)
# ============================================================
trait_file <- file.path(main_dir, "traits.xlsx")
if (file.exists(trait_file)) {
    traits <- read.xlsx(trait_file, rowNames = TRUE)

    if (any(sapply(traits, is.character))) {
        char_cols <- names(traits)[sapply(traits, is.character)]
        for (cc in char_cols) {
            if (all(rownames(traits) == traits[[cc]])) {
                rownames(traits) <- traits[[cc]]
                traits[[cc]] <- NULL
                break
            }
        }
    }

    traits[] <- lapply(traits, function(x) as.numeric(as.character(x)))
    if (anyNA(traits)) stop("traits.xlsx contains NA after numeric conversion.")

    traits <- traits[rownames(datExpr0), , drop = FALSE]
    if (anyNA(traits)) stop("Sample names in traits.xlsx do not match datExpr0.")

    str(traits)

    useMEs <- MEs_ordered[, !colnames(MEs_ordered) %in% "MEgrey", drop = FALSE]
    stopifnot(all(sapply(useMEs, is.numeric)))

    bp_trait        <- bicorAndPvalue(useMEs, traits, maxPOutliers = maxPOutliers)
    moduleTraitCor  <- bp_trait$bicor
    moduleTraitP    <- bp_trait$p

    moduleTraitPadj <- matrix(
        p.adjust(moduleTraitP, method = "BH"),
        nrow = nrow(moduleTraitP),
        dimnames = dimnames(moduleTraitP)
    )

    write.table(moduleTraitCor,
                file.path(main_dir, "moduleTrait_bicor.txt"),
                sep = "\t", quote = FALSE, col.names = NA)
    write.table(moduleTraitP,
                file.path(main_dir, "moduleTrait_p.txt"),
                sep = "\t", quote = FALSE, col.names = NA)
    write.table(moduleTraitPadj,
                file.path(main_dir, "moduleTrait_BHpadj.txt"),
                sep = "\t", quote = FALSE, col.names = NA)

    cor_long <- as.data.frame(moduleTraitCor) %>%
        rownames_to_column("Module") %>%
        mutate(Module = sub("^ME", "", Module)) %>%
        pivot_longer(-Module, names_to = "Trait", values_to = "r")

    p_long <- as.data.frame(moduleTraitPadj) %>%
        rownames_to_column("Module") %>%
        mutate(Module = sub("^ME", "", Module)) %>%
        pivot_longer(-Module, names_to = "Trait", values_to = "p")

    cor_long <- left_join(cor_long, p_long, by = c("Module", "Trait"))

    module_levels <- sub("^ME", "", colnames(useMEs))
    trait_levels  <- colnames(traits)

    cor_long$Module <- factor(cor_long$Module, levels = module_levels)
    cor_long$Trait  <- factor(cor_long$Trait,  levels = rev(trait_levels))

    cor_long <- cor_long %>%
        mutate(
            r_txt = sprintf("%.2f", r),
            p_txt = case_when(
                is.na(p)   ~ "",
                p >= 0.001 ~ sprintf("%.3f", p),
                TRUE       ~ sprintf("%.1e", p)
            ),
            label = paste0(r_txt, "\n(", p_txt, ")"),
            face  = ifelse(!is.na(p) & p < 0.05, "bold", "plain")
        )

    color_luminance <- function(col) {
        rgb_mat <- grDevices::col2rgb(col) / 255
        0.2126 * rgb_mat[1, ] + 0.7152 * rgb_mat[2, ] + 0.0722 * rgb_mat[3, ]
    }
    lum      <- color_luminance(module_levels)
    is_light <- lum > 0.75

    module_html <- ifelse(
        is_light,
        paste0("<span style='color:", module_levels,
               "; text-shadow: 0.8px 0.8px 0.8px #666; font-size:12pt;'>\u25A0</span> ",
               module_levels),
        paste0("<span style='color:", module_levels,
               "; font-size:12pt;'>\u25A0</span> ",
               module_levels)
    )
    html_map <- setNames(module_html, module_levels)

    grad_col <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                  "#F7F7F7",
                  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")

    p_heat <- ggplot(cor_long, aes(x = Module, y = Trait, fill = r)) +
        geom_tile(color = "white", linewidth = 0.6) +
        geom_text(aes(label = label, fontface = face),
                  size = 3.0, color = "black", lineheight = 0.9) +
        scale_x_discrete(labels = html_map) +
        scale_fill_gradientn(
            colours = grad_col,
            limits  = c(-1, 1),
            breaks  = seq(-1, 1, 0.5),
            name    = "bicor",
            guide   = guide_colorbar(
                barheight = unit(3, "cm"),
                barwidth  = unit(0.4, "cm"),
                ticks.colour = "black",
                frame.colour = "black"
            )
        ) +
        labs(x = "Module", y = NULL,
             title = "Module-trait relationships (bicor, BH-adjusted P)") +
        theme_bw(base_size = 11, base_family = "sans") +
        theme(
            panel.grid       = element_blank(),
            panel.border     = element_rect(color = "black", linewidth = 0.6),
            axis.text.x      = element_markdown(angle = 45, hjust = 1, vjust = 1,
                                                size = 10),
            axis.text.y      = element_text(size = 10, color = "black"),
            axis.ticks       = element_line(color = "black"),
            plot.title       = element_text(hjust = 0, face = "bold", size = 12),
            legend.position  = "right",
            legend.title     = element_text(size = 9, face = "bold"),
            legend.text      = element_text(size = 8),
            plot.margin      = margin(8, 8, 8, 8)
        )

    ggsave(file.path(main_dir, "moduleTrait_heatmap.pdf"),
           p_heat, width = 10, height = 8)
} else {
    warning("traits.xlsx not found; skipping module-trait association.")
}

# ============================================================
# 9. Hub and edge selection
# ============================================================
hub_list  <- list()
edge_list <- list()

for (mod in setdiff(unique(mergedColors), "grey")) {
    col <- paste0("ME", mod)
    if (!col %in% colnames(kme_all_mat)) next

    genes <- colnames(datExpr0)[mergedColors == mod]
    kme   <- kme_all_mat[genes, col]

    sel <- kme >= kME_hub_cutoff
    if (sum(sel) > 0) {
        hub_list[[mod]] <- data.frame(
            module = mod,
            gene   = genes[sel],
            kME    = kme[sel]
        )
    }

    subTOM <- TOM[genes, genes]
    subTOM[lower.tri(subTOM, diag = TRUE)] <- NA
    idx <- which(subTOM >= edge_tom_cutoff & !is.na(subTOM), arr.ind = TRUE)
    if (nrow(idx) > 0) {
        edge_list[[mod]] <- data.frame(
            module = mod,
            from   = genes[idx[, 1]],
            to     = genes[idx[, 2]],
            TOM    = subTOM[idx]
        )
    }
}

if (length(hub_list) > 0) {
    hub_df <- do.call(rbind, hub_list)
    write.table(hub_df,
                file.path(main_dir, "hub_genes.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}
if (length(edge_list) > 0) {
    edge_df <- do.call(rbind, edge_list)
    write.table(edge_df,
                file.path(main_dir, "module_edges.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

# ============================================================
# 9b. Export Cytoscape input files (one pair per module)
#     Edge threshold = edge_tom_cutoff (0.1)
# ============================================================
TOM_mat <- as.matrix(TOM)

for (module in sub("^ME", "", colnames(mergedMEs))) {
    if (module == "grey") next

    inModule  <- mergedColors == module
    modProbes <- colnames(datExpr0)[inModule]
    if (length(modProbes) < 2) next

    modTOM <- TOM_mat[inModule, inModule, drop = FALSE]
    dimnames(modTOM) <- list(modProbes, modProbes)

    exportNetworkToCytoscape(
        modTOM,
        edgeFile  = paste("CytoscapeInput-edges-", module, ".txt", sep = ""),
        nodeFile  = paste("CytoscapeInput-nodes-", module, ".txt", sep = ""),
        weighted  = TRUE,
        threshold = edge_tom_cutoff,
        nodeNames = modProbes,
        nodeAttr  = mergedColors[inModule]
    )
    cat("Cytoscape files exported for module:", module,
        "(", length(modProbes), "genes )\n")
}

# ============================================================
# 10. Parameter sensitivity (OFAT)
#     minModuleSize: 20 / 30 / 50
# ============================================================
if (run_sensitivity) {
    sens_results <- list()

    r2_sensitivity_cut <- 0.80
    powers_even <- seq(2, max(powers), by = 2)
    power_grid  <- sft_table$Power[sft_table$Power %in% powers_even &
                                   sft_table$Signed_R2 >= r2_sensitivity_cut]
    power_grid  <- sort(unique(power_grid))

    if (length(power_grid) == 0) {
        stop("No even power satisfies signed R2 >= ", r2_sensitivity_cut)
    }
    if (!(softPower %in% power_grid)) {
        warning("Selected softPower (", softPower,
                ") is not in the even-power grid; adding it explicitly.")
        power_grid <- sort(unique(c(power_grid, softPower)))
    }

    cat("Power grid (even powers, signed R2 >= ", r2_sensitivity_cut, "): ",
        paste(power_grid, collapse = ", "), "\n", sep = "")

    fixed_for_power <- paste0("deepSplit=", deepSplit,
                              "; minModuleSize=", minModuleSize,
                              "; mergeCutHeight=", mergeCutHeight)

    for (p in power_grid) {
        adj_p  <- adjacency(datExpr0, power = p, type = networkType,
                            corFnc = corType,
                            corOptions = list(maxPOutliers = maxPOutliers))
        TOM_p  <- TOMsimilarity(adj_p, TOMType = TOMType)
        tree_p <- hclust(as.dist(1 - TOM_p), method = "average")
        mods_p <- cutreeDynamic(tree_p, distM = 1 - TOM_p,
                                deepSplit = deepSplit,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
        col_p  <- labels2colors(mods_p)
        merged_p <- mergeCloseModules(datExpr0, col_p,
                                      cutHeight = mergeCutHeight,
                                      corFnc = corType,
                                      corOptions = list(maxPOutliers = maxPOutliers),
                                      verbose = 0)$colors
        sens_results[[paste0("power_", p)]] <- data.frame(
            analysis     = "power",
            value        = p,
            n_modules    = length(setdiff(unique(merged_p), "grey")),
            grey_pct     = mean(merged_p == "grey"),
            fixed_params = fixed_for_power,
            stringsAsFactors = FALSE
        )
        cat("  power", p, "done:",
            sens_results[[paste0("power_", p)]]$n_modules, "modules\n")
    }

    fixed_for_ds <- paste0("power=", softPower,
                           "; minModuleSize=", minModuleSize,
                           "; mergeCutHeight=", mergeCutHeight)
    for (ds in 1:4) {
        mods_ds <- cutreeDynamic(geneTree, distM = dissTOM,
                                 deepSplit = ds,
                                 pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize)
        col_ds  <- labels2colors(mods_ds)
        merged_ds <- mergeCloseModules(datExpr0, col_ds,
                                       cutHeight = mergeCutHeight,
                                       corFnc = corType,
                                       corOptions = list(maxPOutliers = maxPOutliers),
                                       verbose = 0)$colors
        sens_results[[paste0("deepSplit_", ds)]] <- data.frame(
            analysis     = "deepSplit",
            value        = ds,
            n_modules    = length(setdiff(unique(merged_ds), "grey")),
            grey_pct     = mean(merged_ds == "grey"),
            fixed_params = fixed_for_ds,
            stringsAsFactors = FALSE
        )
    }

    fixed_for_mms <- paste0("power=", softPower,
                            "; deepSplit=", deepSplit,
                            "; mergeCutHeight=", mergeCutHeight)
    for (m in c(20, 30, 50)) {
        mods_m <- cutreeDynamic(geneTree, distM = dissTOM,
                                deepSplit = deepSplit,
                                pamRespectsDendro = FALSE,
                                minClusterSize = m)
        col_m  <- labels2colors(mods_m)
        merged_m <- mergeCloseModules(datExpr0, col_m,
                                      cutHeight = mergeCutHeight,
                                      corFnc = corType,
                                      corOptions = list(maxPOutliers = maxPOutliers),
                                      verbose = 0)$colors
        sens_results[[paste0("minModuleSize_", m)]] <- data.frame(
            analysis     = "minModuleSize",
            value        = m,
            n_modules    = length(setdiff(unique(merged_m), "grey")),
            grey_pct     = mean(merged_m == "grey"),
            fixed_params = fixed_for_mms,
            stringsAsFactors = FALSE
        )
    }

    fixed_for_mch <- paste0("power=", softPower,
                            "; deepSplit=", deepSplit,
                            "; minModuleSize=", minModuleSize)
    for (h in c(0.10, 0.15, 0.20, 0.25, 0.30)) {
        merged_h <- mergeCloseModules(datExpr0, dynamicColors,
                                      cutHeight = h,
                                      corFnc = corType,
                                      corOptions = list(maxPOutliers = maxPOutliers),
                                      verbose = 0)$colors
        sens_results[[paste0("mergeCutHeight_", h)]] <- data.frame(
            analysis     = "mergeCutHeight",
            value        = h,
            n_modules    = length(setdiff(unique(merged_h), "grey")),
            grey_pct     = mean(merged_h == "grey"),
            fixed_params = fixed_for_mch,
            stringsAsFactors = FALSE
        )
    }

    sens_df <- do.call(rbind, sens_results)
    write.table(sens_df,
                file.path(main_dir, "sensitivity_parameters.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

# ============================================================
# 11. Stability: bootstrap Jaccard overlap of modules
# ============================================================
if (run_stability) {
    set.seed(123)
    ref_modules <- setdiff(unique(mergedColors), "grey")
    ref_genes   <- lapply(ref_modules, function(m)
        colnames(datExpr0)[mergedColors == m])
    names(ref_genes) <- ref_modules

    stab_jaccard <- matrix(NA, nrow = n_bootstrap, ncol = length(ref_modules),
                           dimnames = list(NULL, ref_modules))

    for (b in seq_len(n_bootstrap)) {
        idx <- sample(nSamples, round(0.8 * nSamples))
        subExpr <- datExpr0[idx, , drop = FALSE]

        adj_b  <- adjacency(subExpr, power = softPower, type = networkType,
                            corFnc = corType,
                            corOptions = list(maxPOutliers = maxPOutliers))
        TOM_b  <- TOMsimilarity(adj_b, TOMType = TOMType)
        tree_b <- hclust(as.dist(1 - TOM_b), method = "average")
        mods_b <- cutreeDynamic(tree_b, distM = 1 - TOM_b,
                                deepSplit = deepSplit,
                                pamRespectsDendro = FALSE,
                                minClusterSize = minModuleSize)
        col_b  <- labels2colors(mods_b)
        merged_b <- mergeCloseModules(subExpr, col_b,
                                      cutHeight = mergeCutHeight,
                                      corFnc = corType,
                                      corOptions = list(maxPOutliers = maxPOutliers),
                                      verbose = 0)$colors
        names(merged_b) <- colnames(subExpr)

        for (m in ref_modules) {
            ref_set <- intersect(ref_genes[[m]], colnames(subExpr))
            if (length(ref_set) < 5) next

            sub_mods <- setdiff(unique(merged_b), "grey")
            jac <- sapply(sub_mods, function(sm) {
                sub_set <- names(merged_b)[merged_b == sm]
                inter <- length(intersect(ref_set, sub_set))
                union <- length(union(ref_set, sub_set))
                if (union == 0) 0 else inter / union
            })

            if (length(jac) == 0 || all(is.na(jac))) {
                stab_jaccard[b, m] <- NA
            } else {
                stab_jaccard[b, m] <- max(jac, na.rm = TRUE)
            }
        }
        if (b %% 10 == 0) cat("Bootstrap", b, "/", n_bootstrap, "done\n")
    }

    stab_summary <- data.frame(
        module = ref_modules,
        mean_jaccard = colMeans(stab_jaccard, na.rm = TRUE),
        sd_jaccard   = apply(stab_jaccard, 2, sd, na.rm = TRUE),
        n_bootstrap  = n_bootstrap
    )
    write.table(stab_summary,
                file.path(main_dir, "stability_bootstrap.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as.data.frame(stab_jaccard),
                file.path(main_dir, "stability_bootstrap_raw.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("All done.\n")