# ============================================================
# Freeze-Thaw Cytokine Analysis Pipeline
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(emmeans)
  library(cowplot)
  library(grid)
  library(readr)
  library(ggforce)   # for rounded, labeled hulls on PCA
})

# --------------------------
# Output folders
# --------------------------
cat("Creating output directories...\n")
dir.create("freeze_thaw_analysis", recursive = TRUE, showWarnings = FALSE)
dir.create("freeze_thaw_analysis/individual_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("freeze_thaw_analysis/pca_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("freeze_thaw_analysis/heatmaps", recursive = TRUE, showWarnings = FALSE)
dir.create("freeze_thaw_analysis/summary_tables", recursive = TRUE, showWarnings = FALSE)

# Basic path check
cat("Checking directories:\n")
cat("  freeze_thaw_analysis/ exists:", dir.exists("freeze_thaw_analysis"), "\n")
cat("  Current working directory:", getwd(), "\n")

# --------------------------
# 1) Load and combine panel data
# --------------------------
first_existing <- function(paths) {
  for (p in paths) if (file.exists(p)) return(p)
  stop("File not found in any of these locations:\n", paste(paths, collapse = "\n"))
}

load_panel <- function(path){
  df <- read_csv(path, show_col_types = FALSE) %>%
    rename_with(~str_replace_all(., "\\s+", "_")) %>%
    rename(
      Freeze_Thaw_Cycle      = any_of(c("Freeze_Thaw_Cycle","Freeze_Thaw")),
      Aliquot                = any_of(c("Aliquot","Aliquot_1")),
      Protein_Name           = any_of(c("Protein_Name","Analyte")),
      Cytokine_Concentration = any_of(c("Cytokine_Concentration","Concentration","Value"))
    ) %>%
    mutate(
      Matrix = trimws(as.character(Matrix)),
      Matrix = recode(Matrix, "PL-EDTA"="PL-EDTA", "PL-HE"="PL-HE", "SE"="SE", .default = Matrix),
      Freeze_Thaw_Cycle = as.numeric(Freeze_Thaw_Cycle),
      Donor_ID   = as.numeric(Donor_ID),
      Protein_Name = as.character(Protein_Name),
      Cytokine_Concentration = as.numeric(Cytokine_Concentration),
      Log_Protein_Amount = log10(Cytokine_Concentration + 1)
    )
  
  # If Aliquot is missing/NA, create a stable sample id
  if (!"Aliquot" %in% names(df) || all(is.na(df$Aliquot))) {
    df <- df %>%
      mutate(Aliquot = paste0("S_", Donor_ID, "_", Matrix, "_C", Freeze_Thaw_Cycle, "_", row_number()))
  }
  
  df %>%
    select(Protein_Name, Log_Protein_Amount, Matrix, Freeze_Thaw_Cycle, Donor_ID, Aliquot) %>%
    drop_na(Protein_Name, Log_Protein_Amount, Matrix, Freeze_Thaw_Cycle, Donor_ID)
}

panel1 <- load_panel(first_existing(c(
  "data/Panel1_vert_clean.csv","Panel1_vert_clean.csv","/mnt/data/Panel1_vert_clean.csv")))
panel2 <- load_panel(first_existing(c(
  "data/Panel2_vert_clean.csv","Panel2_vert_clean.csv","/mnt/data/Panel2_vert_clean.csv")))
panel3 <- load_panel(first_existing(c(
  "data/Panel3_vert_clean.csv","Panel3_vert_clean.csv","/mnt/data/Panel3_vert_clean.csv")))

all_data <- bind_rows(panel1, panel2, panel3) %>%
  mutate(
    Matrix = factor(Matrix, levels = c("PL-EDTA","PL-HE","SE")),
    Freeze_Thaw_Cycle = as.numeric(Freeze_Thaw_Cycle)
  )

# --------------------------
# 2) Matrix-by-cycle model for each cytokine
# --------------------------
analyze_one_cytokine <- function(protein_name, data, show_details = FALSE) {
  d <- data %>% filter(Protein_Name == protein_name) %>% droplevels()
  if (nrow(d) < 12) {
    return(tibble(Protein_Name = protein_name, Status = "Not_Enough_Data"))
  }
  
  tryCatch({
    # Compare models with and without the matrix-by-cycle interaction
    main_ml  <- lme(Log_Protein_Amount ~ Matrix + Freeze_Thaw_Cycle,
                    random = ~1|Donor_ID, data = d, method = "ML")
    inter_ml <- lme(Log_Protein_Amount ~ Matrix * Freeze_Thaw_Cycle,
                    random = ~1|Donor_ID, data = d, method = "ML")
    a <- anova(main_ml, inter_ml)
    interaction_p <- a$`p-value`[2]
    
    # Refit the interaction model with REML for slope estimates
    inter_reml <- update(inter_ml, method = "REML")
    
    # Estimate freeze-thaw slope within each matrix
    tr <- emtrends(inter_reml, ~ Matrix, var = "Freeze_Thaw_Cycle")
    tr_sum <- summary(tr, infer = c(TRUE, TRUE)) %>%
      as_tibble() %>%
      mutate(Matrix = trimws(as.character(Matrix))) %>%
      transmute(
        Matrix,
        Slope = `Freeze_Thaw_Cycle.trend`,
        SE    = SE,
        df    = df,
        t     = `t.ratio`,
        p     = `p.value`
      ) %>%
      group_by(Matrix) %>% slice(1) %>% ungroup() %>%
      mutate(Matrix = factor(Matrix, levels = c("PL-EDTA","PL-HE","SE"))) %>%
      arrange(Matrix)
    
    out <- tr_sum %>%
      pivot_wider(
        names_from  = Matrix,
        values_from = c(Slope, SE, df, t, p),
        names_sep   = "_",
        values_fn   = dplyr::first
      )
    
    # Overall change from cycle 1 to cycle 4
    c1 <- mean(d$Log_Protein_Amount[d$Freeze_Thaw_Cycle == 1], na.rm = TRUE)
    c4 <- mean(d$Log_Protein_Amount[d$Freeze_Thaw_Cycle == 4], na.rm = TRUE)
    pc <- { o1 <- 10^c1 - 1; o4 <- 10^c4 - 1; if (is.finite(o1) && o1 > 0) 100*(o4-o1)/o1 else NA_real_ }
    
    out %>%
      mutate(
        Protein_Name        = protein_name,
        Interaction_P_Value = interaction_p,
        Percent_Change      = pc,
        Status              = "Success"
      ) %>%
      relocate(Protein_Name, Interaction_P_Value, Percent_Change, Status)
    
  }, error = function(e) {
    tibble(Protein_Name = protein_name, Status = "Analysis_Failed", Error = as.character(e$message))
  })
}

# --------------------------
# 3) Run model across cytokines
# --------------------------
proteins <- sort(unique(all_data$Protein_Name))
cat("Analyzing", length(proteins), "cytokines...\n")
res_list <- lapply(seq_along(proteins), function(i){
  if (i %% 10 == 0) cat("Progress:", i, "of", length(proteins), "\n")
  analyze_one_cytokine(proteins[i], all_data, show_details = FALSE)
})
uniform_results <- bind_rows(res_list) %>% filter(Status == "Success")
cat("Successfully analyzed:", nrow(uniform_results), "cytokines\n")

# --------------------------
# 4) Summary table
# --------------------------
make_mentor_summary <- function(uniform_results) {
  uniform_results %>%
    select(-starts_with("SE_"), -starts_with("df_"), -starts_with("t_"),
           -starts_with("FDR_"), -any_of("Panel")) %>%
    mutate(
      `Matrix analysis` = ifelse(
        is.na(Interaction_P_Value), "NA",
        paste0(ifelse(Interaction_P_Value < 0.001, "<0.001", sprintf("%.4f", Interaction_P_Value)),
               " (", ifelse(Interaction_P_Value < 0.05, "Significant", "Not Significant"), ")")
      )
    ) %>%
    rowwise() %>%
    mutate(
      `Freeze Thaw Effect on Cytokine` = {
        if (is.na(Interaction_P_Value)) {
          "N/A"
        } else if (Interaction_P_Value < 0.05) {
          "N/A: Matrix Difference"
        } else {
          s  <- c(`Slope_PL-EDTA`,`Slope_PL-HE`,`Slope_SE`)
          pv <- c(`p_PL-EDTA`,`p_PL-HE`,`p_SE`)
          if (all(s < 0, na.rm = TRUE) && all(pv < 0.05, na.rm = TRUE)) "Decreasing"
          else if (all(s > 0, na.rm = TRUE) && all(pv < 0.05, na.rm = TRUE)) "Increasing"
          else "Stable"
        }
      }
    ) %>%
    ungroup() %>%
    relocate(Protein_Name, `Freeze Thaw Effect on Cytokine`, `Matrix analysis`,
             Interaction_P_Value, Percent_Change, everything())
}

mentor_summary <- make_mentor_summary(uniform_results)
cat("Writing mentor summary to:", file.path(getwd(), "freeze_thaw_analysis/summary_tables/Uniform_Summary_ForMentor.csv"), "\n")
write_csv(mentor_summary, "freeze_thaw_analysis/summary_tables/Uniform_Summary_ForMentor.csv")

# --------------------------
# 5) Identify matrix with strongest interaction signal
# --------------------------
make_matrix_diff_color_table <- function(uniform_results) {
  uniform_results %>%
    select(Protein_Name, Interaction_P_Value, `Slope_PL-EDTA`, `Slope_PL-HE`, `Slope_SE`) %>%
    rowwise() %>%
    mutate(
      Different_Matrix = if (is.na(Interaction_P_Value) || Interaction_P_Value >= 0.05) {
        NA_character_
      } else {
        s <- c(`Slope_PL-EDTA`,`Slope_PL-HE`,`Slope_SE`)
        labs <- c("PL-EDTA","PL-HE","SE")
        labs[ which.max(abs(s - mean(s, na.rm = TRUE))) ]
      },
      Color = case_when(
        Different_Matrix == "PL-EDTA" ~ "yellow",
        Different_Matrix == "PL-HE"   ~ "red",
        Different_Matrix == "SE"      ~ "green",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    select(Protein_Name, Different_Matrix, Color)
}
matrix_diff_colors <- make_matrix_diff_color_table(uniform_results)
cat("Writing matrix difference colors to:", file.path(getwd(), "freeze_thaw_analysis/summary_tables/Matrix_Difference_ColorTable.csv"), "\n")
write_csv(matrix_diff_colors, "freeze_thaw_analysis/summary_tables/Matrix_Difference_ColorTable.csv")

# --------------------------
# 6) Per-cytokine plots
# --------------------------
create_uniform_plot <- function(protein_name, data, results_df) {
  d <- data %>% filter(Protein_Name == protein_name) %>% droplevels()
  r <- results_df %>% filter(Protein_Name == protein_name)
  if (nrow(d) == 0 || nrow(r) == 0) return(NULL)
  
  mats   <- c("PL-EDTA","PL-HE","SE")
  colors <- c("PL-EDTA"="#E31A1C","PL-HE"="#1F78B4","SE"="#33A02C")
  
  get_lab <- function(m) {
    s <- r[[paste0("Slope_", m)]]
    p <- r[[paste0("p_", m)]]
    s_txt <- ifelse(is.na(s), "N/A", sprintf("%.4f", s))
    p_txt <- ifelse(is.na(p), "N/A", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
    paste0("Slope = ", s_txt, ", p = ", p_txt)
  }
  
  pct <- r$Percent_Change[1]
  top_title <- protein_name
  top_subtitle <- paste0(
    "Overall change (Cycle 1 → 4): ",
    ifelse(is.na(pct), "N/A", sprintf("%.1f%%", pct)),
    "   |   Matrix Interaction p = ",
    ifelse(r$Interaction_P_Value[1] < 0.001, "<0.001", sprintf("%.4f", r$Interaction_P_Value[1]))
  )
  
  y_rng <- range(d$Log_Protein_Amount, na.rm = TRUE); buf <- diff(y_rng) * 0.15
  y_min <- y_rng[1] - buf; y_max <- y_rng[2] + buf
  
  jar_plots <- lapply(mats, function(m) {
    dm <- d %>% filter(Matrix == m)
    ggplot(dm, aes(Freeze_Thaw_Cycle, Log_Protein_Amount)) +
      geom_smooth(method = "lm", se = FALSE, color = colors[m],
                  linewidth = 3.0, alpha = 0.95) +
      geom_point(aes(color = factor(Donor_ID)), size = 4, alpha = 0.9) +
      geom_line(aes(group = Donor_ID, color = factor(Donor_ID)),
                linewidth = 1.2, linetype = "dotted", alpha = 0.7) +
      scale_color_brewer(palette = "Set1", name = "Donor") +
      scale_x_continuous(breaks = 1:4, labels = 1:4) +
      scale_y_continuous(limits = c(y_min, y_max)) +
      labs(title = m, subtitle = get_lab(m),
           x = "Freeze-Thaw Cycles", y = "log10 MFI") +
      theme_bw(base_size = 18) +
      theme(
        plot.title = element_text(face="bold", size = 20),
        plot.subtitle = element_text(color = colors[m], face="bold", size = 16),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        axis.text.y = element_text(size = 16, face = "bold"),
        legend.position = "bottom",
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        panel.grid.minor = element_blank()
      )
  })
  
  # Figure title band
  title_grob <- ggdraw() +
    draw_label(top_title, fontface="bold", size=20, color="white", y=0.7) +
    draw_label(top_subtitle, y = 0.3, size=14, color="white") +
    theme(plot.background = element_rect(fill = "gray10", color = NA))
  
  # Matrix-specific panels
  grid <- plot_grid(plotlist = jar_plots, ncol = length(jar_plots))
  
  # Combine title and panels
  final_plot <- plot_grid(title_grob, grid, ncol = 1, rel_heights = c(0.12, 0.88))
  
  return(final_plot)
}

cat("Creating plots...\n")
for (i in seq_len(nrow(uniform_results))) {
  pr <- uniform_results$Protein_Name[i]
  p  <- create_uniform_plot(pr, all_data, uniform_results)
  if (!is.null(p)) {
    ggsave(file.path("freeze_thaw_analysis/individual_plots/",
                     paste0(gsub("[^A-Za-z0-9_-]", "_", pr), "_UNIFORM.png")),
           p, width = 15, height = 8.5, dpi = 300)
  }
  if (i %% 20 == 0) cat("  ...", i, "of", nrow(uniform_results), "plots saved\n")
}

# --------------------------
# 7) PCA of cytokine profiles
# --------------------------
prepare_pca_data <- function(data) {
  cat("\nPreparing PCA data...\n")
  
  # Convert long cytokine measurements to a sample-by-cytokine matrix
  wide_data <- data %>%
    select(Aliquot, Donor_ID, Matrix, Freeze_Thaw_Cycle, Protein_Name, Log_Protein_Amount) %>%
    # Average duplicate sample-cytokine measurements
    group_by(Aliquot, Donor_ID, Matrix, Freeze_Thaw_Cycle, Protein_Name) %>%
    summarise(Log_Protein_Amount = mean(Log_Protein_Amount, na.rm = TRUE), .groups = "drop") %>%
    # Wide format: one row per sample, one column per cytokine
    pivot_wider(names_from = Protein_Name, values_from = Log_Protein_Amount, values_fill = NA)
  
  cat("Wide data dimensions:", nrow(wide_data), "samples x", ncol(wide_data)-4, "proteins\n")
  
  # Separate sample metadata from cytokine values
  meta <- wide_data %>% select(Aliquot, Donor_ID, Matrix, Freeze_Thaw_Cycle)
  protein_data <- wide_data %>% select(-Aliquot, -Donor_ID, -Matrix, -Freeze_Thaw_Cycle)
  
  # Remove cytokines with high missingness
  missing_prop <- colMeans(is.na(protein_data))
  keep_proteins <- missing_prop <= 0.75
  protein_data <- protein_data[, keep_proteins, drop = FALSE]
  
  cat("After removing proteins with >75% missing:", ncol(protein_data), "proteins remain\n")
  
  if (ncol(protein_data) == 0) {
    stop("ERROR: No proteins have sufficient data (all have >75% missing values)")
  }
  
  # Remove samples with high missingness
  sample_missing_prop <- rowMeans(is.na(protein_data))
  keep_samples <- sample_missing_prop <= 0.5
  protein_data <- protein_data[keep_samples, , drop = FALSE]
  meta <- meta[keep_samples, , drop = FALSE]
  
  cat("After removing samples with >50% missing:", nrow(protein_data), "samples remain\n")
  
  if (nrow(protein_data) == 0) {
    stop("ERROR: No samples have sufficient data")
  }
  
  # Mean-impute remaining missing values by cytokine
  for (i in 1:ncol(protein_data)) {
    col_mean <- mean(protein_data[[i]], na.rm = TRUE)
    if (!is.na(col_mean)) {
      protein_data[[i]][is.na(protein_data[[i]])] <- col_mean
    }
  }
  
  # Remove cytokines with near-zero variance
  protein_vars <- apply(protein_data, 2, var, na.rm = TRUE)
  keep_variable <- protein_vars > 1e-10 & !is.na(protein_vars)
  protein_data <- protein_data[, keep_variable, drop = FALSE]
  
  cat("After removing constant proteins:", ncol(protein_data), "proteins remain\n")
  
  if (ncol(protein_data) < 2) {
    stop("ERROR: Need at least 2 variable proteins for PCA")
  }
  
  # Scale cytokine features before PCA
  scaled_data <- scale(protein_data)
  
  # Replace any non-finite scaled values
  if (any(!is.finite(scaled_data))) {
    cat("WARNING: Non-finite values after scaling, replacing with 0\n")
    scaled_data[!is.finite(scaled_data)] <- 0
  }
  
  cat("Final PCA matrix:", nrow(scaled_data), "samples x", ncol(scaled_data), "proteins\n")
  
  return(list(X = scaled_data, meta = meta, protein_names = colnames(scaled_data)))
}

cat("Creating PCA plots...\n")
pca_data <- prepare_pca_data(all_data)
pca_fit <- prcomp(pca_data$X, center = FALSE, scale. = FALSE)
pve <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
pc_lab <- function(k){ paste0("PC", k, " (", sprintf("%.1f", 100*pve[k]), "%)") }

# PCA scores with sample metadata
pca_scores <- as_tibble(pca_fit$x[, 1:10]) %>%  keep first 10 PCs
  bind_cols(pca_data$meta) %>%
  mutate(
    Freeze_Thaw_Cycle = factor(Freeze_Thaw_Cycle, levels = 1:4, labels = paste("Cycle", 1:4)),
    Donor_ID = factor(Donor_ID),
    Matrix = factor(Matrix, levels = c("PL-EDTA", "PL-HE", "SE"))
  )

# Add group counts for hull labels
pca_scores_with_counts <- pca_scores %>%
  group_by(Freeze_Thaw_Cycle, Matrix) %>%
  mutate(n_group = n(),
         hull_label = paste0(Matrix, " (n=", n_group, ")")) %>%
  ungroup()

# Multi-panel PCA plot
create_multipanel_pca <- function(pca_scores, pve) {
  
  # Axis labels with variance explained
  pc_label <- function(pc_num) {
    paste0("PC", pc_num, " (", sprintf("%.1f%%", 100 * pve[pc_num]), " variance)")
  }
  
  # Panel 1: By Donor
  p1 <- ggplot(pca_scores, aes(PC1, PC2)) +
    geom_point(aes(color = Donor_ID), size = 3, alpha = 0.8) +
    stat_ellipse(aes(color = Donor_ID), level = 0.68, type = "norm",
                 linetype = "dashed", size = 1, alpha = 0.8) +
    scale_color_brewer(palette = "Set1", name = "Donor") +
    labs(
      title = "Principal Component Analysis: By Donor",
      x = pc_label(1),
      y = pc_label(2)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Panel 2: By Freeze-Thaw Cycle
  p2 <- ggplot(pca_scores, aes(PC1, PC2)) +
    geom_point(aes(color = Freeze_Thaw_Cycle), size = 3, alpha = 0.8) +
    stat_ellipse(aes(color = Freeze_Thaw_Cycle), level = 0.68, type = "norm",
                 linetype = "dashed", size = 1, alpha = 0.8) +
    scale_color_viridis_d(name = "Freeze-Thaw\nCycle", option = "plasma") +
    labs(
      title = "Principal Component Analysis: By Freeze-Thaw Cycle",
      x = pc_label(1),
      y = pc_label(2)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Panel 3: By Matrix Type
  p3 <- ggplot(pca_scores, aes(PC1, PC2)) +
    geom_point(aes(color = Matrix), size = 3, alpha = 0.8) +
    stat_ellipse(aes(color = Matrix), level = 0.68, type = "norm",
                 linetype = "dashed", size = 1, alpha = 0.8) +
    scale_color_manual(
      values = c("PL-EDTA" = "#2E8B57", "PL-HE" = "#FF6347", "SE" = "#4169E1"),
      name = "Matrix"
    ) +
    labs(
      title = "Principal Component Analysis: By Matrix Type",
      x = pc_label(1),
      y = pc_label(2)
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 11)
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))
  
  # Combine panels
  plot_grid(p1, p2, p3, ncol = 1, align = "v", axis = "lr")
}

# Save PCA plots
multipanel_pca <- create_multipanel_pca(pca_scores, pve)

# Save the clean multi-panel plot
ggsave("freeze_thaw_analysis/pca_plots/PCA_MultiPanel_Analysis.png",
       multipanel_pca, width = 8, height = 16, dpi = 300)

# Faceted PCA with group outlines
p_pca <- ggplot(pca_scores_with_counts, aes(PC1, PC2)) +
  ggforce::geom_mark_hull(
    aes(group = interaction(Freeze_Thaw_Cycle, Matrix),
        fill = Matrix, label = hull_label),
    alpha = 0.15, expand = unit(3, "pt"), concavity = 2, radius = unit(4, "pt"),
    label.fontsize = 8, label.fill = "white", label.colour = "black",
    label.buffer = unit(4, "pt")
  ) +
  geom_point(aes(color = Matrix, shape = Donor_ID), size = 3, alpha = 0.8) +
  facet_wrap(~ Freeze_Thaw_Cycle, nrow = 1) +
  scale_color_brewer(palette = "Set1", name = "Matrix") +
  scale_fill_brewer(palette = "Set1", name = "Matrix") +
  scale_shape_discrete(name = "Donor") +
  labs(
    title = "PCA of Cytokine Profiles: Matrix Effects Across Freeze-Thaw Cycles",
    subtitle = paste("Based on", ncol(pca_data$X), "cytokines across", nrow(pca_data$X), "samples"),
    x = pc_lab(1), y = pc_lab(2)
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom"
  )

ggsave("freeze_thaw_analysis/pca_plots/PCA_Hull_Faceted.png",
       p_pca, width = 16, height = 5, dpi = 300)

cat("PCA plots saved successfully.\n")

# --------------------------
# 8) Heatmap of matrix-by-cycle cytokine means
# --------------------------
make_matrix_cycle_heat <- function(data){
  cat("Creating heatmap data...\n")
  
  avg <- data %>%
    select(Protein_Name, Matrix, Freeze_Thaw_Cycle, Log_Protein_Amount) %>%
    group_by(Protein_Name, Matrix, Freeze_Thaw_Cycle) %>%
    summarise(mean_log = mean(Log_Protein_Amount, na.rm = TRUE), .groups = "drop") %>%
    mutate(Col = paste(Matrix, paste0("C", Freeze_Thaw_Cycle), sep = "_")) %>%
    select(-Matrix, -Freeze_Thaw_Cycle) %>%
    pivot_wider(names_from = Col, values_from = mean_log)
  
  cat("Heatmap dimensions:", nrow(avg), "proteins x", ncol(avg)-1, "conditions\n")
  
  mat <- avg %>%
    select(-Protein_Name) %>%
    as.matrix()
  
  rownames(mat) <- avg$Protein_Name
  
  complete_rows <- rowSums(!is.na(mat)) > 0
  mat <- mat[complete_rows, , drop = FALSE]
  
  cat("After removing empty rows:", nrow(mat), "proteins\n")
  
  mat_z <- t(scale(t(mat)))
  mat_z[is.nan(mat_z)] <- 0
  
  cat("Z-score matrix created successfully\n")
  return(mat_z)
}

cat("Creating heatmap...\n")
mat_z <- make_matrix_cycle_heat(all_data)

heat_df <- as_tibble(mat_z, rownames = "Protein_Name") %>%
  pivot_longer(-Protein_Name, names_to = "Column", values_to = "Z") %>%
  separate(Column, into = c("Matrix","Cycle"), sep = "_", remove = FALSE)

p_heat <- ggplot(heat_df, aes(x = Column, y = Protein_Name, fill = Z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                       midpoint = 0, name = "Z-score") +
  labs(title = "Cytokine Heatmap (Row Z-scored)",
       subtitle = "Columns = Matrix × Cycle (mean across donors)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

ggsave("freeze_thaw_analysis/heatmaps/Cytokines_Heatmap_MatrixCycle.png",
       p_heat, width = 10, height = 14, dpi = 300)

# --------------------------
# 9) Create tables by freeze-thaw effect category
# --------------------------
create_effect_category_tables <- function(mentor_summary) {
  cat("Creating effect category tables...\n")
  
  decreasing_cytokines <- mentor_summary %>%
    filter(`Freeze Thaw Effect on Cytokine` == "Decreasing") %>%
    select(Protein_Name, `Matrix analysis`, Interaction_P_Value, Percent_Change,
           `Slope_PL-EDTA`, `p_PL-EDTA`, `Slope_PL-HE`, `p_PL-HE`, `Slope_SE`, `p_SE`) %>%
    arrange(Protein_Name)
  
  increasing_cytokines <- mentor_summary %>%
    filter(`Freeze Thaw Effect on Cytokine` == "Increasing") %>%
    select(Protein_Name, `Matrix analysis`, Interaction_P_Value, Percent_Change,
           `Slope_PL-EDTA`, `p_PL-EDTA`, `Slope_PL-HE`, `p_PL-HE`, `Slope_SE`, `p_SE`) %>%
    arrange(Protein_Name)
  
  stable_cytokines <- mentor_summary %>%
    filter(`Freeze Thaw Effect on Cytokine` == "Stable") %>%
    select(Protein_Name, `Matrix analysis`, Interaction_P_Value, Percent_Change,
           `Slope_PL-EDTA`, `p_PL-EDTA`, `Slope_PL-HE`, `p_PL-HE`, `Slope_SE`, `p_SE`) %>%
    arrange(Protein_Name)
  
  matrix_effect_cytokines <- mentor_summary %>%
    filter(`Freeze Thaw Effect on Cytokine` == "N/A: Matrix Difference") %>%
    select(Protein_Name, `Matrix analysis`, Interaction_P_Value, Percent_Change,
           `Slope_PL-EDTA`, `p_PL-EDTA`, `Slope_PL-HE`, `p_PL-HE`, `Slope_SE`, `p_SE`) %>%
    arrange(Protein_Name)
  
  cat("Writing effect category tables...\n")
  write_csv(decreasing_cytokines, "freeze_thaw_analysis/summary_tables/Decreasing_Cytokines_Table.csv")
  cat("  Saved:", file.path(getwd(), "freeze_thaw_analysis/summary_tables/Decreasing_Cytokines_Table.csv"), "\n")
  
  write_csv(increasing_cytokines, "freeze_thaw_analysis/summary_tables/Increasing_Cytokines_Table.csv")
  cat("  Saved:", file.path(getwd(), "freeze_thaw_analysis/summary_tables/Increasing_Cytokines_Table.csv"), "\n")
  
  write_csv(stable_cytokines, "freeze_thaw_analysis/summary_tables/Stable_Cytokines_Table.csv")
  cat("  Saved:", file.path(getwd(), "freeze_thaw_analysis/summary_tables/Stable_Cytokines_Table.csv"), "\n")
  
  write_csv(matrix_effect_cytokines, "freeze_thaw_analysis/summary_tables/Matrix_Effect_Cytokines_Table.csv")
  cat("  Saved:", file.path(getwd(), "freeze_thaw_analysis/summary_tables/Matrix_Effect_Cytokines_Table.csv"), "\n")
  
  cat("\n=== CYTOKINE EFFECT SUMMARY ===\n")
  cat("Decreasing cytokines:", nrow(decreasing_cytokines), "\n")
  cat("Increasing cytokines:", nrow(increasing_cytokines), "\n")
  cat("Stable cytokines:", nrow(stable_cytokines), "\n")
  cat("Matrix effect cytokines:", nrow(matrix_effect_cytokines), "\n")
  cat("Total analyzed:", nrow(mentor_summary), "\n\n")
  
  if (nrow(decreasing_cytokines) > 0) {
    cat("DECREASING CYTOKINES:\n")
    cat(paste(decreasing_cytokines$Protein_Name, collapse = ", "), "\n\n")
  }
  
  if (nrow(increasing_cytokines) > 0) {
    cat("INCREASING CYTOKINES:\n")
    cat(paste(increasing_cytokines$Protein_Name, collapse = ", "), "\n\n")
  }
  
  if (nrow(stable_cytokines) > 0) {
    cat("STABLE CYTOKINES:\n")
    cat(paste(stable_cytokines$Protein_Name, collapse = ", "), "\n\n")
  }
  
  if (nrow(matrix_effect_cytokines) > 0) {
    cat("MATRIX EFFECT CYTOKINES:\n")
    cat(paste(matrix_effect_cytokines$Protein_Name, collapse = ", "), "\n\n")
  }
  
  return(list(
    decreasing = decreasing_cytokines,
    increasing = increasing_cytokines,
    stable = stable_cytokines,
    matrix_effect = matrix_effect_cytokines
  ))
}

effect_tables <- create_effect_category_tables(mentor_summary)

# --------------------------
# 10) Save full results table
# --------------------------
write_csv(uniform_results, "freeze_thaw_analysis/summary_tables/Uniform_Interaction_Results.csv")
cat("Writing full results to:", file.path(getwd(), "freeze_thaw_analysis/summary_tables/Uniform_Interaction_Results.csv"), "\n")

# Print output files
cat("\n=== FILES CREATED ===\n")
result_files <- list.files("freeze_thaw_analysis", recursive = TRUE, full.names = TRUE)
for(f in result_files) {
  cat("-", f, "\n")
}

cat("\nAnalysis complete.\n",
    "Files written to freeze_thaw_analysis/ ...\n",
    " • summary_tables/Uniform_Interaction_Results.csv\n",
    " • summary_tables/Uniform_Summary_ForMentor.csv\n",
    " • summary_tables/Matrix_Difference_ColorTable.csv\n",
    " • summary_tables/Decreasing_Cytokines_Table.csv\n",
    " • summary_tables/Increasing_Cytokines_Table.csv\n",
    " • summary_tables/Stable_Cytokines_Table.csv\n",
    " • summary_tables/Matrix_Effect_Cytokines_Table.csv\n",
    " • individual_plots/ (one PNG per protein)\n",
    " • pca_plots/PCA_MultiPanel_Analysis.png\n",
    " • pca_plots/PCA_Hull_Faceted.png\n",
    " • heatmaps/Cytokines_Heatmap_MatrixCycle.png\n", sep = "")

