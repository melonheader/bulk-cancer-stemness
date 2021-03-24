### INFO: Plotting functions for Stem-TCGA project
### DATE: 23.12.2019
### AUTHOR: Artem Baranovskii


# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# Plotting domain
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Main theme for ggplot
pal_comp <- list(c("#7cad3e", "#236477"), 
                 c("#82b21e", "#006beb"), 
                 c("#59821e", "#3c8d93"))


theme_tcga <- function(base_family = "sans", legend.position = "right", legend.direction = "vertical", aspect.ratio = 1, angle.t.x = 45, angle.t.y = 45, ...){
  theme_bw(base_family = base_family, ...) + 
    theme(text = element_text(colour = "#333333", hjust = 1),
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          legend.background = element_rect(fill = alpha('white', 0.1)),
          axis.text.y = element_text(angle = angle.t.y, margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
          axis.text.x = element_text(angle = angle.t.x, margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
          axis.line = element_line(colour = "black"),
          axis.ticks.length = unit(-0.05, "in"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          panel.background = element_rect(fill = "white", 
                                          colour = "white",
                                          size = 0.5, 
                                          linetype = "solid"))
}



# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Legend collector
get_legend <- function(ggobj){ 
  tmp <- ggplot_gtable(ggplot_build(ggobj)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  ggleg <- ggpubr::as_ggplot(legend)
  return(ggleg)
} 

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## PCA-plot drawer
plot_nice_pca <- function(pca_list, pcs = c('PC1', 'PC2'), density = FALSE, ...) {
  #
  pca_eigvals <- pca_list[[2]]
  ## project-wise nice correlation plots PC1 vs Purity
  plot_tab <- pca_list[[1]] %>%
    mutate(project_id = str_remove_all(project_id, "TCGA-"), 
           col_var = ifelse(project_id == c("Differentiating iPSC"), "iPSC", "tmp"), 
           col_var = ifelse(project_id != c("Differentiating iPSC") & sample_type == "Solid Tissue Normal", "Normal tissue", col_var), 
           col_var = ifelse(project_id != c("Differentiating iPSC") & sample_type != "Solid Tissue Normal", "Tumor", col_var),
           col_var = factor(col_var, levels = c("iPSC", "Tumor", "Normal tissue"))) 
  #
  if (density == TRUE) {
    ggobj <- ggplot(data = plot_tab, 
                    aes(x = !! sym(pcs[1]), 
                        y = !! sym(pcs[2]))) +
      stat_density_2d(data = plot_tab %>% filter(col_var == "Tumor"), 
                      aes(fill = stat(density)), 
                      geom = "raster", 
                      contour = FALSE) + 
      geom_point(data = plot_tab %>% filter(col_var != "Tumor"), 
                 aes(color = col_var, 
                     shape = col_var), 
                 alpha = 0.4, 
                 size = 0.4) +
      scale_fill_viridis_c(guide = FALSE) + 
      scale_shape_manual(values = c(16, 15), guide = FALSE) +
      scale_color_manual(values = c("#FF1493", "grey70")) +
      guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
      labs(x = paste0("PC", as.numeric(str_sub(pcs[1], 3)), " ", round((pca_eigvals[as.numeric(str_sub(pcs[1], 3))] / sum(pca_eigvals))*100, 0), "%"), 
           y = paste0("PC", as.numeric(str_sub(pcs[2], 3)), " ", round((pca_eigvals[as.numeric(str_sub(pcs[2], 3))] / sum(pca_eigvals))*100, 0), "%")) + 
      theme_tcga(...)
    return(ggobj)
  } else {
    ggobj <- ggplot(data = plot_tab, 
                    aes(x = !! sym(pcs[1]), 
                        y = !! sym(pcs[2]))) +
      geom_point(data = plot_tab %>% filter(col_var != "iPSC"), 
                 aes(color = col_var, 
                     shape = col_var), 
                 alpha = 0.2, 
                 size = 1.0) +
      geom_point(data = plot_tab %>% filter(col_var == "iPSC"), 
                 aes(fill = rep(c(1, 2, 3, 4), each = 8)[2:32]),
                 shape = 25,
                 alpha = 0.8, 
                 size = 2, 
                 color = "#666666") + 
      scale_shape_manual(values = c(16, 15), guide = FALSE) +
      scale_color_manual(values = c("#f88e38", "#83a9c1")) +
      scale_fill_gradientn("iPSC differentiation stage: ", 
                           colours = c("red", "blue")) +
      guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
      labs(x = paste0("PC", as.numeric(str_sub(pcs[1], 3)), " ", round((pca_eigvals[as.numeric(str_sub(pcs[1], 3))] / sum(pca_eigvals))*100, 0), "%"), 
           y = paste0("PC", as.numeric(str_sub(pcs[2], 3)), " ", round((pca_eigvals[as.numeric(str_sub(pcs[2], 3))] / sum(pca_eigvals))*100, 0), "%"), 
           color = "Sample type:") + 
      theme_tcga(...)
    return(ggobj)
  }
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Jitter-plot drawer
plot_nice_jitter <- function(pca_table, cluster_table, pc = 1, ...) {
  #Correct project names
  pca_table %<>% mutate(project_id = str_remove_all(project_id, "TCGA-"))
  cluster_table %<>% mutate(project_id = str_remove_all(project_id, "TCGA-"))
  
  #Join
  table_clust <- pca_table %>% 
    dplyr::select(gene_set, cases, project_id, sample_type, site_of_resection_or_biopsy, primary_diagnosis, gender, !! sym(paste0("PC", pc))) %>% 
    full_join(cluster_table, 
              ., 
              by = c("cases", "project_id", "sample_type")) %>% 
    mutate(!! sym(pca_table[[1]][1]) := ifelse(sample_type == "Solid Tissue Normal", "Solid Tissue Normal", 
                                               ifelse(project_id == "Differentiating iPSC", "Differentiating iPSC", !! sym(pca_table[[1]][1]))), 
           project_id = ifelse(project_id == "Differentiating iPSC" | sample_type == "Solid Tissue Normal", sample_type, project_id))
  #Median
  if (pca_table$gene_set[[1]] == "EMT_MET") {
    median_tum <- -5
  } else {
    median_tum <- table_clust %>% 
      filter(!sample_type %in% c("iPSC", "Primitive Streak", "Endoderm Progenitors", "Definitive Endoderm", "Solid Tissue Normal")) %>% 
      pull(!! sym(paste0("PC", pc))) %>% 
      median(., na.rm = TRUE)
  }
  #Sort
  table_clust_order <- table_clust %>% 
    filter(!sample_type %in% c("iPSC", "Primitive Streak", "Endoderm Progenitors", "Definitive Endoderm", "Solid Tissue Normal")) %>%
    group_by(project_id, !! sym(table_clust[[5]][1])) %>% 
    summarise(count = n(), .groups = "drop") %>% 
    spread(., 
           key = !! sym(table_clust[[5]][1]), 
           value = count)
  table_clust_order[is.na(table_clust_order)] <- 0
  table_clust_order <- table_clust_order %>% 
    mutate(intensity = (ESC - Normal) / (ESC + Normal)) %>% 
    arrange(desc(intensity))
  table_clust_order <- c('iPSC',
                         'Primitive Streak', 
                         'Endoderm Progenitors', 
                         'Definitive Endoderm', 
                         table_clust_order$project_id, 
                         'Solid Tissue Normal')
  table_clust <- table_clust %>% 
    mutate(project_id = factor(project_id, levels = table_clust_order))
  #Plot
  ggobj <- ggplot() + 
    geom_jitter(data = table_clust, 
                aes(x = project_id, y = !! sym(paste0("PC", pc)), group = project_id, color = !! sym(table_clust[[5]][1])), 
                size = 0.7, 
                alpha = 0.3) + 
    scale_color_manual(values = c("#CD5C5C", "#fc4a1a", "#F8AA38", "#83a9c1", "#f7b733")) + 
    geom_hline(yintercept = median_tum, size = 1, color = "grey30") +  
    coord_flip() +
    #scale_alpha(guise = FALSE) +
    #scale_size(guide = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    theme_tcga(...) +
    labs(subtitle = paste0(table_clust[[5]][1], "\nsignature"),
         color = "Legend:", 
         y = paste0("PC", pc), 
         x = "")
  return(ggobj)
  
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Scatter-plot drawer
plot_nice_scatter <- function(cluster_stem, 
                              cluster_prolif, 
                              purity, 
                              meta,
                              col_var = "purity",
                              shape_var = NA, 
                              comp_medians = FALSE, 
                              pca_tables_list = NULL, 
                              ...) {
  if (comp_medians == TRUE) {
    if (purrr::is_null(pca_tables_list) == TRUE) {
      stop("Please, provide 'pca_tables_list'")
    }
    medians_list <- list()
    for (i in 1:length(pca_tables_list)) {
      medians_list[[i]] <- pca_tables_list[[i]][[1]] %>% 
        dplyr::select(cases, project_id, sample_type, PC1) %>% 
        mutate(PC1 = norm_to_one(-1 * PC1)) %>%
        filter(!sample_type %in% c('Solid Tissue Normal', 'Differentiating iPSC'), 
               !project_id %in% c('Solid Tissue Normal', 'Differentiating iPSC')) %>%
        group_by(project_id) %>% 
        summarise(!! sym(paste0(str_split_fixed(pca_tables_list[[i]][[1]][1, 1], "_", 2)[, 1], "_Intensity")) := median(PC1), .groups = "drop")
    }
    medians_table <- purrr::reduce(medians_list, inner_join, by = "project_id") %>% 
      mutate(project_id = str_split_fixed(project_id, "-", 2)[, 2])
    purity_proj <- purity %>% 
      group_by(project_id) %>% 
      summarise(median_purity = median(purity, na.rm = TRUE), .groups = "drop") %>% 
      mutate(project_id = str_split_fixed(project_id, "-", 2)[, 2])
    cluster_freq <- medians_table %>% 
      inner_join(., purity_proj, 
                 by = "project_id") %>% 
      mutate(purity_factorial = ifelse(median_purity < 0.72, "Low", 
                                       ifelse(median_purity >= 0.72 & median_purity < 0.83, "Medium", 
                                              ifelse(median_purity >= 0.83, "High", median_purity)))) %>% 
      mutate(purity_factorial = factor(purity_factorial, levels = c("Low", "Medium", "High")))
  } else {
    clusters <- list(cluster_stem, 
                     cluster_prolif) %>% 
      map(., 
          ~ dplyr::filter(., 
                          project_id != "Differentiating iPSC", 
                          sample_type != "Solid Tissue Normal")) %>% 
      purrr::reduce(., inner_join, by = c("cases", "project_id", "sample_type")) %>% 
      dplyr::select(project_id, Stemness, Proliferation)
    stem_freq <- clusters %>% 
      dplyr::count(project_id, Stemness, sort = TRUE) %>% 
      mutate(Stemness = str_c(colnames(.[, 2]), Stemness, sep = "_")) %>% 
      spread(., Stemness, n) %>% 
      mutate(Stemness_ESC = ifelse(is.na(Stemness_ESC), 0, Stemness_ESC), 
             Stemness_Normal = ifelse(is.na(Stemness_Normal), 0, Stemness_Normal), 
             Stemness_Intensity = Stemness_ESC / (Stemness_Normal + Stemness_ESC))
    prolif_freq <- clusters %>% 
      dplyr::count(project_id, Proliferation, sort = TRUE) %>% 
      mutate(Proliferation = str_c(colnames(.[, 2]), Proliferation, sep = "_")) %>% 
      spread(., Proliferation, n) %>% 
      mutate(Proliferation_ESC = ifelse(is.na(Proliferation_ESC), 0, Proliferation_ESC), 
             Proliferation_Normal = ifelse(is.na(Proliferation_Normal), 0, Proliferation_Normal),  
             Proliferation_Intensity = Proliferation_ESC / (Proliferation_Normal + Proliferation_ESC))
    cluster_freq <- list(stem_freq, 
                         prolif_freq) %>% 
      purrr::reduce(., inner_join, by = "project_id") %>% 
      mutate(project_id = str_split_fixed(project_id, "-", 2)[, 2]) %>% 
      dplyr::select(project_id, Stemness_Intensity, Proliferation_Intensity)
    # Update 29/08/19 - compute median Purity for each project and attach to freq_table
    purity_proj <- purity %>% 
      group_by(project_id) %>% 
      summarise(median_purity = median(purity, na.rm = TRUE), .groups = "drop") %>% 
      mutate(project_id = str_split_fixed(project_id, "-", 2)[, 2]) %>% 
      mutate(purity_factorial = ifelse(median_purity < 0.72, "Low", 
                                       ifelse(median_purity >= 0.72 & median_purity < 0.83, "Medium", 
                                              ifelse(median_purity >= 0.83, "High", median_purity)))) %>% 
      mutate(purity_factorial = factor(purity_factorial, levels = c("Low", "Medium", "High")))
    ## Join
    cluster_freq %<>% 
      inner_join(., purity_proj, 
                 by = "project_id")
    # Update 06/11/19 - compute average tumor_stage for each project and attach to freq_tablev
    avg_stage <- meta %>% filter(grepl("TCGA", project_id), 
                                 !tumor_stage %in% c("-", "--", "i/ii nos", "is", "not reported", "stage 0", "stage x")) %>% 
      mutate(tumor_stage_simp = str_count(tumor_stage, "i"), 
             tumor_stage_simp = ifelse(grepl("iv", tumor_stage), 4, tumor_stage_simp)) %>% 
      tabyl(project_id, tumor_stage_simp) %>% 
      mutate(project_id = str_remove_all(project_id, "TCGA-")) %>% 
      gather(., 
             -contains("project"), 
             key = "stage_simp", 
             value = "count") %>% 
      mutate(stage_simp = as.numeric(stage_simp)) %>% 
      group_by(project_id) %>% 
      summarise(stage_avg = weighted.mean(x = stage_simp, w = (count / sum(count))), .groups = "drop") %>% 
      ungroup() %>% 
      mutate(stage_avg_factorial = ifelse(stage_avg < 2, "Low (1-2)", 
                                          ifelse(stage_avg >= 2 & stage_avg < 3, "Medium (2-3)", 
                                                 ifelse(stage_avg >= 3, "High (3-4)", stage_avg))), 
             stage_avg_factorial = factor(stage_avg_factorial, levels = c("Low (1-2)", "Medium (2-3)", "High (3-4)")))
    ## Joim 
    cluster_freq %<>% 
      left_join(., avg_stage, 
                by = "project_id")
    
  }
  
  #Nice scatter
  sc1 <- ggplot(data = cluster_freq, 
                aes(x = Proliferation_Intensity, 
                    y = Stemness_Intensity, 
                    label = project_id, 
                    color = !!sym(paste0(col_var, "_factorial")), 
                    shape = !!sym(paste0(shape_var, "_factorial")))) +
    geom_hline(yintercept = 0.5) + 
    geom_vline(xintercept = 0.5) +
    geom_point(size = 5) + 
    ggrepel::geom_label_repel(force = 40, seed = 42, size = 2.7) + 
    xlim(c(0, 1)) + 
    ylim(c(0, 1)) + 
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1)), 
           shape = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    theme_tcga(...) +
    labs(color = ifelse(col_var == "purity", "Median\nproject-wise\npurity", 
                        ifelse(col_var == "stage_avg", "Average\nproject-wise\ntumor stage")), 
         shape = ifelse(shape_var == "purity", "Median\nproject-wise\npurity", 
                        ifelse(shape_var == "stage_avg", "Average\nproject-wise\ntumor stage", 
                               ifelse(is.na(shape_var), NA))), 
         y = "Stemness Intensity", x = "Proliferation Intensity")
  return(sc1)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Heatmap drawer
plot_nice_heatmap <- function(cluster_stem, cluster_prolif,
                              legend.position = "right", 
                              legend.direction = "vertical",
                              comp_medians = FALSE, 
                              pca_tables_list = NULL, 
                              ...) {
  
  if (comp_medians == TRUE) {
    if (purrr::is_null(pca_tables_list) == TRUE) {
      stop("Please, provide 'pca_tables_list'")
    }
    medians_list <- list()
    for (i in 1:length(pca_tables_list)) {
      medians_list[[i]] <- pca_tables_list[[i]][[1]] %>% 
        dplyr::select(cases, project_id, sample_type, PC1) %>% 
        mutate(PC1 = norm_to_one(-1 * PC1)) %>%
        filter(!sample_type %in% c('Solid Tissue Normal', 'Differentiating iPSC'), 
               !project_id %in% c('Solid Tissue Normal', 'Differentiating iPSC')) %>%
        group_by(project_id) %>% 
        summarise(!! sym(paste0(str_split_fixed(pca_tables_list[[i]][[1]][1, 1], "_", 2)[, 1], "_Intensity")) := median(PC1), 
                  .groups = "drop")
    }
    cluster_freq <- purrr::reduce(medians_list, inner_join, by = "project_id") %>% 
      mutate(project_id = str_split_fixed(project_id, "-", 2)[, 2])
  } else {
    clusters <- list(cluster_stem, 
                     cluster_prolif) %>% 
      map(., 
          ~ filter(., 
                   project_id != "Differentiating iPSC", 
                   sample_type != "Solid Tissue Normal")) %>% 
      purrr::reduce(., inner_join, by = c("cases", "project_id", "sample_type")) %>% 
      dplyr::select(project_id, Stemness, Proliferation)
    stem_freq <- clusters %>% 
      dplyr::count(project_id, Stemness, sort = TRUE) %>% 
      mutate(Stemness = str_c(colnames(.[, 2]), Stemness, sep = "_")) %>% 
      spread(., Stemness, n) %>% 
      mutate(Stemness_ESC = ifelse(is.na(Stemness_ESC), 0, Stemness_ESC), 
             Stemness_Normal = ifelse(is.na(Stemness_Normal), 0, Stemness_Normal), 
             Stemness_Intensity = Stemness_ESC / (Stemness_Normal + Stemness_ESC))
    prolif_freq <- clusters %>% 
      dplyr::count(project_id, Proliferation, sort = TRUE) %>% 
      mutate(Proliferation = str_c(colnames(.[, 2]), Proliferation, sep = "_")) %>% 
      spread(., Proliferation, n) %>% 
      mutate(Proliferation_ESC = ifelse(is.na(Proliferation_ESC), 0, Proliferation_ESC), 
             Proliferation_Normal = ifelse(is.na(Proliferation_Normal), 0, Proliferation_Normal),  
             Proliferation_Intensity = Proliferation_ESC / (Proliferation_Normal + Proliferation_ESC))
    cluster_freq <- list(stem_freq, 
                         prolif_freq) %>% 
      purrr::reduce(., inner_join, by = "project_id") %>% 
      mutate(project_id = str_split_fixed(project_id, "-", 2)[, 2]) %>% 
      dplyr::select(project_id, Stemness_Intensity, Proliferation_Intensity)
  }
  
  
  #Hierarchical order
  # Dissimilarity matrix
  d <- cluster_freq %>% 
    dplyr::select(-project_id) %>% 
    as.matrix(.) %>% 
    `rownames<-`(cluster_freq$project_id) %>%
    dist(., method = "euclidean")
  # Hierarchical clustering using Complete Linkage
  hc1 <- hclust(d, method = "complete")
  #Make tall
  cluster_freq_tall <- cluster_freq %>%
    gather(.,
           contains("Intensity"), 
           key = "Signature", 
           value = "Intensity"
    ) %>% 
    mutate(project_id = factor(project_id, levels = hc1$labels[hc1$order]), 
           Signature = str_remove_all(Signature, "_Intensity"))
  #Nice heatmap
  hm1 <- ggplot(data = cluster_freq_tall, 
                aes(x = project_id, 
                    y = Signature)) + 
    geom_tile(aes(fill = Intensity)) + 
    scale_fill_viridis_c() +
    labs(y = "", x = "") +
    theme_tcga(...)
  return(hm1)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Cox-reg boxplot drawer 
plot_surv_boxplots <- function(cluster_table, meta, ...) {
  #change project names
  cluster_table %<>% mutate(project_id = str_remove_all(project_id, "TCGA-"))
  meta %<>% mutate(project_id = str_remove_all(project_id, "TCGA-"))
  #Clustering as an independent predictor of survival
  cluster <- names(cluster_table[4])
  names(cluster_table)[4] <- "Cluster"
  
  #Prepare data
  projects <- cluster_table %>% 
    filter(!project_id %in% c("Differentiating iPSC", "KICH", "PRAD")) %>% 
    distinct(project_id) %>%
    pull(project_id)
  results.list <- list()
  for (i in 1:length(projects)) {
    data.surv <- cluster_table %>% 
      filter(project_id != "Differentiating iPSC") %>% 
      filter(sample_type != "Solid Tissue Normal") %>% 
      filter(project_id == projects[i]) %>%
      inner_join(., meta, 
                 by = c("cases", "project_id", "sample_type")) %>%
      dplyr::rename(Site = site_of_resection_or_biopsy, 
                    Diagnosis = primary_diagnosis) %>%
      mutate(days_to_event = as.numeric(ifelse(days_to_death == "--", days_to_last_follow_up, days_to_death)), 
             event = ifelse(vital_status == "dead", 1, 0)) %>%
      filter(days_to_event > 0) %>%
      #na.omit() %>% 
      mutate(Cluster = factor(Cluster, levels = c("Normal", "ESC")))
    
    surv.obj <- survival::Surv(time = data.surv$days_to_event, event = data.surv$event, type = "right")
    fit.coxph <- survival::coxph(surv.obj ~ Cluster, 
                                 data = data.surv)
    #Prepare table of results
    results <- unclass(summary(fit.coxph))
    data.to.write <- cbind(results$coefficients, results$conf.int)
    data.to.write <- data.to.write %>%
      as.tibble() %>%
      dplyr::select(-V6) %>%
      mutate(Covariates = rownames(data.to.write),
             TCGA_project = projects[i],
             Concordance = results$concordance[[1]],
             Conc.SE = results$concordance[[2]],
             Rsq = results$rsq[[1]],
             loglike.test = results$logtest[[1]], 
             loglike.pval = results$logtest[[3]]) %>%
      dplyr::select(TCGA_project, Covariates, everything())
    #Attach to list
    results.list[[projects[i]]] <- data.to.write
  }
  
  #Write when fully ready
  results_table <- purrr::reduce(results.list, rbind)
  
  #Attach stats for plotting
  #Attach plotting variables
  res.surv <- results_table %>%
    mutate(p_adj = p.adjust(`Pr(>|z|)`, method = "fdr"), 
           Significance = ifelse(p_adj < 5e-2, 
                                 "*", 
                                 "NS"), 
           Significance = ifelse(p_adj < 1e-2, 
                                 "**", 
                                 Significance), 
           Significance = ifelse(p_adj < 1e-3, 
                                 "***", 
                                 Significance)) %>% 
    mutate(`Risk of Adverse Event` = ifelse(`lower .95` > 1 & Significance != "NS", 
                                            "Higher risk", 
                                            ifelse(`upper .95` < 1 & Significance != "NS", 
                                                   "Lower risk", 
                                                   "Not significant")),
           `Risk of Adverse Event` = factor(`Risk of Adverse Event`, levels = c("Higher risk", "Lower risk", "Not significant")))
  #Prepare palette
  # palette_check <- res.surv %>% 
  #   pull(`Risk of Adverse Event`) %>% 
  #   unique(.)
  # 
  # if (length(palette_check) == 2 & any(str_detect(palette_check, "Higher Risk")) & any(str_detect(palette_check, "Lower Risk", negate = TRUE))) {
  #   palette <- c("#F8AA38", "grey80")
  # } else if (length(palette_check) == 2 & any(str_detect(palette_check, "Lower Risk")) & any(str_detect(palette_check, "Higher Risk", negate = TRUE))) { 
  #   palette <- c("#5A92B7", "grey80")
  # } else if (length(palette_check) == 2 & any(str_detect(palette_check, "Not significant", negate = TRUE))) {
  #   palette <- c("#F8AA38", "#5A92B7")
  # } else if (length(palette_check) == 3) {
  #   palette <- c("#F8AA38","#5A92B7", "grey80")
  # } else if (length(palette_check) == 1 & any(str_detect(palette_check, "Not significant"))) {
  #   palette <- c("grey80")
  # }
  palette <- c("#F8AA38","#5A92B7", "grey80")
  
  
  #Plot
  ggobj <- ggplot(res.surv, aes(y = `exp(coef)`, x = TCGA_project, color = `Risk of Adverse Event`)) + 
    geom_point(shape = 15, size =2) +
    geom_errorbar(aes(ymin = `lower .95`, ymax = `upper .95`), 
                  width = .2, 
                  position = position_dodge(0.05)) + 
    scale_color_manual(values = palette) +
    geom_text(aes(label = Significance, x = TCGA_project, y = `upper .95`), 
              hjust = 2, nudge_x = 0.1,
              size = 2, 
              color = "black") +
    geom_hline(yintercept = 1, color = "red", size = 0.5) + 
    coord_flip(ylim = c(0, 10)) +
    theme_tcga(...) +
    xlab("TCGA project id") +
    ylab("Hazard Ratio") +
    labs(title = paste0(cluster, " signature"), 
         subtitle = paste0("Survival analysis between ESC and Normal clusters", "\n* pval < 0.05, ** pval < 0.005, *** pval < 0.001"))
  return(ggobj)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Cox-reg curves drawer
plot_nice_curves <- function(cluster_table, meta, project, file_path, save = FALSE, show = TRUE, ...) {
  #Clustering as an independent predictor of survival
  cluster <- names(cluster_table[4])
  names(cluster_table)[4] <- "Cluster"
  #Prepare data
  data.surv <- cluster_table %>% 
    filter(project_id != "Differentiating iPSC") %>% 
    filter(sample_type != "Solid Tissue Normal") %>% 
    filter(grepl(project, project_id)) %>%
    inner_join(., meta, 
               by = c("cases", "project_id", "sample_type")) %>%
    dplyr::rename(Site = site_of_resection_or_biopsy, 
                  Diagnosis = primary_diagnosis) %>%
    mutate(days_to_event = as.numeric(ifelse(days_to_death == "--", days_to_last_follow_up, days_to_death)), 
           event = ifelse(vital_status == "dead", 1, 0)) %>%
    filter(days_to_event > 0) %>%
    na.omit() %>% 
    mutate(Cluster = factor(Cluster, levels = c("Normal", "ESC")))
  
  #surv.obj <- survival::Surv(time = data.surv$days_to_event, event = data.surv$event, type = "right")
  fit.surv <- do.call(survival::survfit, 
                      list(formula = survival::Surv(time = days_to_event, event = event, type = "right") ~ Cluster, data = data.surv))
  
  #survival::survfit(surv.obj ~ Cluster, 
  #data = data.surv)
  
  
  surv.curve <- survminer::ggsurvplot(fit = fit.surv, 
                                      data = data.surv, pval = TRUE, pval.size = 3.5,
                                      censor = TRUE, 
                                      palette = c("#5A92B7", "#F8AA38"),  
                                      legend.title = "Cluster:", 
                                      legend.labs = c("Normal", 
                                                      "ESC"),  
                                      xlab = "Time, hours",
                                      subtitle = paste0(project, " ", cluster, "\ncluster")
  )
  ggobj <- surv.curve$plot + theme_tcga(...)
  #Branch
  if (show == TRUE) {
    return(ggobj)
  } else if (save == TRUE) {
    #Write plot
    svg(filename = file.path(file_path, paste0(project, "_SCurve", ".svg")), width = 6, height = 6)
    ggobj
    dev.off()
    #Talk
    return(print(paste0("Survival curve for ", cluster, " cluster of ", project)))
  }
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## PC1 vs Purity cor-plots drawer
plot_nice_corr <- function(pca_table, meta, purity, pc = 1, project = "global") {
  if (project == "global") {
    data_tmp <- pca_table %>% 
      dplyr::filter(grepl("Tumor", sample_type)) %>%
      dplyr::select(cases, project_id, sample_type, !! sym(paste0("PC", pc))) %>% 
      inner_join(., purity[, !grepl("project_id", names(purity))], 
                 by = "cases")
    cor_metric <- cor(data_tmp[[paste0("PC", pc)]], data_tmp$purity, method = "spearman")
    cor_pval <- cor.test(data_tmp[[paste0("PC", pc)]], data_tmp$purity, method = "spearman", exact = FALSE)$p.value
    ggobj <- ggplot(data = data_tmp, 
                    aes(x = purity, 
                        y = !! sym(paste0("PC", pc)))) + 
      geom_point(color = "grey50", 
                 size = 0.5) + 
      geom_smooth(method = "lm", se = FALSE) +
      theme_bw(base_family = "", 
               base_size = base_size,
               base_line_size = base_size/22, 
               base_rect_size = base_size/22) + 
      theme(text = element_text(colour = "#49463C", hjust = 1), 
            plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(face = "italic", colour = "grey20"),
            axis.text.y = element_text(angle = 0),
            axis.text.x = element_text(angle = 45),
            legend.position = "right",
            legend.direction = "vertical",
            #legend.justification = c(1, 1),
            legend.background = element_rect(fill = alpha('white', 0.1)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line = element_line(),
            panel.background = element_rect(fill = "white", 
                                            colour = "white",
                                            size = 0.5, 
                                            linetype = "solid")) +
      xlab("purity estimate") +
      ylab(paste0("PC", pc)) +
      labs(subtitle = paste0(project, " correlation between PC", pc, " and purity", "\nSpearman correlation = ", round(cor_metric, 3), "; pvalue = ", round(cor_pval, 6))) 
    return(ggobj)
  } else {
    data_tmp <- pca_table %>% 
      dplyr::filter(grepl("Tumor", sample_type), 
                    grepl(project, project_id)) %>%
      dplyr::select(cases, project_id, sample_type, !! sym(paste0("PC", pc))) %>% 
      inner_join(., purity[, !grepl("project_id", names(purity))], 
                 by = "cases")
    cor_metric <- cor(data_tmp[[paste0("PC", pc)]], data_tmp$purity, method = "spearman")
    cor_pval <- cor.test(data_tmp[[paste0("PC", pc)]], data_tmp$purity, method = "spearman", exact = FALSE)$p.value
    ggobj <- ggplot(data = data_tmp, 
                    aes(x = purity, 
                        y = !! sym(paste0("PC", pc)))) + 
      geom_point(color = "grey50", 
                 size = 0.5) + 
      geom_smooth(method = "lm", se = FALSE) +
      theme_bw(base_family = "", 
               base_size = base_size,
               base_line_size = base_size/22, 
               base_rect_size = base_size/22) + 
      theme(text = element_text(colour = "#49463C", hjust = 1), 
            plot.title = element_text(face = "bold"),
            plot.subtitle = element_text(face = "italic", colour = "grey20"),
            axis.text.y = element_text(angle = 0),
            axis.text.x = element_text(angle = 45),
            legend.position = "right",
            legend.direction = "vertical",
            #legend.justification = c(1, 1),
            legend.background = element_rect(fill = alpha('white', 0.1)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), 
            axis.line = element_line(),
            panel.background = element_rect(fill = "white", 
                                            colour = "white",
                                            size = 0.5, 
                                            linetype = "solid")) +
      xlab("purity estimate") +
      ylab(paste0("PC", pc)) +
      labs(subtitle = paste0(project, " correlation between PC", pc, " and purity", "\nSpearman correlation = ", round(cor_metric, 3), "; pvalue = ", round(cor_pval, 6)))
    return(ggobj)
  }
}


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Dot-plot for purity drawer
dotplot_purity <- function(data_purity, ...) {
  # remove TCGA-
  data_purity <- data_purity %>% 
    mutate(project_id = str_remove_all(project_id, "TCGA-"))
  # arrange
  order <- data_purity %>% 
    filter(!is.na(purity)) %>% 
    group_by(project_id) %>% 
    summarise(median_purity = median(purity), .groups = "drop_last") %>% 
    arrange(desc(median_purity)) %>% 
    pull(project_id)
  data_purity %<>% 
    mutate(project_id = factor(project_id, levels = order))
  ggplot(data = data_purity, 
         aes(x = project_id, 
             y = purity)) + 
    geom_boxplot(width = 0.4, outlier.shape = NA) +
    #geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.003) +
    theme_tcga(...) + 
    xlab("") + ylab("Conssensus purity estimate, %") 
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Hist. drawer to show effect of regressing tumor purity
hist_medians <- function(medians_table, ..., c_pal) {
  medians_table <- medians_table %>%  mutate(Group = ifelse(Group == "Median", "Raw (1 + cpm)\ncounts", "CPE-corrected\ncounts"))
  ggplot(data = medians_table, 
         aes(x = Log2CPM, 
             fill = Group)) + 
    geom_histogram(position = "identity",
                   color = "grey30",
                   alpha = 0.2, 
                   binwidth = 0.5) +
    scale_fill_manual(values = c_pal) + 
    theme_tcga(...) + 
    xlab("log2(1 + cpm)") + 
    ylab("Frequency")
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Scatter plot drawer (median counts)
scatter_medians <- function(medians_table, ...) {
  
  medians_table_wide <- medians_table %>% 
    spread(., key = "Group", value = "Log2CPM")
  # ---------------------------------------------
  ggplot(data = medians_table_wide) + 
    geom_point(aes(x = Median, 
                   y = Median_Cor), 
               size = 0.4, 
               alpha = 0.1,
               color = "grey20") +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    theme_tcga(...) +
    ylab("Median CPE-corrected\ncounts") + xlab("Median log2(1 + cpm)")
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## PCA-plot drawer (color by gradient of selected gene)
gradient_pca <- function(data_total, mset, gene, pcs = c("PC1", "PC2"), correct = TRUE, cor_v = 1, ...) {
  #
  if (grepl("stem", mset)) {
    mset <- "stem"
  } 
  if (grepl("prolif", mset)) {
    mset <- "prolif"
  }
  if (grepl("emt", mset)) {
    mset <- "emt"
  }
  #
  exp_table <- data_total$Norm_Counts[[paste0(mset, "_counts_table")]]
  purity <- data_total$Purity
  # get ens id
  if (!grepl("ENS", gene)) {
    gns <- data_total[["GNames"]] %>% 
      dplyr::filter(gene_name %in% gene)
    if (purrr::is_null(data_total$Total_table)) {
      # check
      gns <- gns %>% filter(ensembl_gene_id %in% colnames(exp_table))
      if (dim(gns)[1] == 0) {
        stop("No genes detected")
      }
    } else {
      gns <- gns %>% filter(ensembl_gene_id %in% colnames(data_total$Total_table))
      if (dim(gns)[1] == 0) {
        stop("No genes detected")
      }
    }
  }
  # prepare table
  if (correct == TRUE) {
    exp_table <- correct_cpe(exp_table, purity, version = cor_v)
    ## get expression levels of genes of interest
    if (purrr::is_null(data_total$Total_table)) {
      gns_table <- exp_table[, c("cases", gns$ensembl_gene_id)]
    } else {
      gns_table <- data_total$Total_table[, c("cases", "project_id", "sample_type", "site_of_resection_or_biopsy", "primary_diagnosis", "gender", gns$ensembl_gene_id)]
      gns_table <- correct_cpe(gns_table, purity, version = cor_v)
      gns_table <- gns_table %>% dplyr::select(cases, contains("ENS"))
    }
  } else {
    if (purrr::is_null(data_total$Total_table)) {
      gns_table <- exp_table[, c("cases", gns$ensembl_gene_id)]
    } else {
      gns_table <- data_total$Total_table[, c("cases", gns$ensembl_gene_id)]
    }
  }
  pca_obj <- do_pca(exp_table)
  pca_obj[[1]] <- pca_obj[[1]] %>% inner_join(., gns_table, by = "cases")  
  pca_df <- pca_obj[[1]]
  pca_eigvals <- pca_obj[[2]]
  
  # organise for plotting
  ggobj <- ggplot(data = pca_df, 
                  aes(x = !! sym(pcs[1]), 
                      y = !! sym(pcs[2]))) +
    geom_point(data = pca_df[pca_df$sample_type == "Solid Tissue Normal", ], 
               color = "grey30",
               alpha = 0.4, 
               size = 1) +
    geom_point(data = pca_df[pca_df$sample_type %in% c("Primary solid Tumor", "Metastatic", "Recurrent Solid Tumor", 
                                                       "Additional Metastatic", "Additional - New Primary"), ], 
               aes(color = !! sym(gns$ensembl_gene_id)),
               alpha = 0.7, 
               size = 0.5) + 
    scale_color_gradient(name = paste0(gns$gene_name, ", log2(1 + cpm)")) +
    geom_point(data = pca_df[pca_df$project_id == "Differentiating iPSC", ], 
               color = "#ba370b",
               alpha = 1, 
               size = 1.2, 
               shape = 17) + 
    labs(x = paste("PC", as.numeric(str_sub(pcs[1], 3)), round((pca_eigvals[as.numeric(str_sub(pcs[1], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " "), 
         y = paste("PC", as.numeric(str_sub(pcs[2], 3)), round((pca_eigvals[as.numeric(str_sub(pcs[2], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " "), 
         color = gene) + 
    theme_tcga(...)
  return(ggobj)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Qunatile plot drawer (show change in expression of selected gene per quantile of PC1)
quantile_pca <- function(data_total, mset, gene, correct = TRUE, cor_v = 1, f_scales = "fixed", ...) {
  #
  if (grepl("stem", mset)) {
    mset <- "stem"
  } 
  if (grepl("prolif", mset)) {
    mset <- "prolif"
  }
  if (grepl("emt", mset)) {
    mset <- "emt"
  }
  #
  exp_table <- data_total$Norm_Counts[[paste0(mset, "_counts_table")]]
  purity <- data_total$Purity
  # get ens id
  if (!any(grepl("ENS", gene))) {
    gns <- data_total[["GNames"]] %>% 
      dplyr::filter(gene_name %in% gene)
    if (purrr::is_null(data_total$Total_table)) {
      # check
      gns <- gns %>% filter(ensembl_gene_id %in% colnames(exp_table))
      if (dim(gns)[1] == 0) {
        stop("No genes detected")
      }
    } else {
      gns <- gns %>% filter(ensembl_gene_id %in% colnames(data_total$Total_table))
      if (dim(gns)[1] == 0) {
        stop("No genes detected")
      }
    }
  }
  # prepare table
  if (correct == TRUE) {
    exp_table <- correct_cpe(exp_table, purity, version = cor_v)
    ## get expression levels of genes of interest
    if (purrr::is_null(data_total$Total_table)) {
      gns_table <- exp_table[, c("cases", gns$ensembl_gene_id)]
    } else {
      gns_table <- data_total$Total_table[, c("cases", "project_id", "sample_type", "site_of_resection_or_biopsy", "primary_diagnosis", "gender", gns$ensembl_gene_id)]
      gns_table <- correct_cpe(gns_table, purity, version = cor_v)
      gns_table <- gns_table %>% dplyr::select(cases, contains("ENS"))
    }
  } else {
    if (purrr::is_null(data_total$Total_table)) {
      gns_table <- exp_table[, c("cases", gns$ensembl_gene_id)]
    } else {
      gns_table <- data_total$Total_table[, c("cases", gns$ensembl_gene_id)]
    }
  }
  pca_obj <- do_pca(exp_table)
  pca_obj[[1]] <- pca_obj[[1]] %>% inner_join(., gns_table, by = "cases")  
  pca_df <- pca_obj[[1]] %>% filter(sample_type %in% c("Primary solid Tumor", "Metastatic", "Recurrent Solid Tumor", 
                                                       "Additional Metastatic", "Additional - New Primary"))
  # compute quantiles
  quantiles <- quantile(pca_df[["PC1"]])
  pca_df %<>% mutate(Quartile = "I", 
                     Quartile = ifelse(PC1 <= quantiles[[2]], "I", 
                                       ifelse(PC1 > quantiles[[2]] & PC1 <= quantiles[[3]], "II", 
                                              ifelse(PC1 > quantiles[[3]] & PC1 <= quantiles[[4]], "III", 
                                                     ifelse(PC1 > quantiles[[4]], "IV", Quartile)))))
  # box plot
  ## vectorized
  plot_tab <- pca_df %>% 
    `[`(c("Quartile", gns$ensembl_gene_id)) %>% 
    tidyr::gather(contains("ENS"), key = "gene_id", value = "val") %>% 
    inner_join(., gns, 
               by = c("gene_id" = "ensembl_gene_id"))
    
  
  ggobj <- ggplot(data = plot_tab, 
                  aes(x = Quartile, 
                      y = val)) + 
    geom_boxplot(width = 0.6) + 
    labs(x = "PC1, quartiles", y = "log2(1+ cpm)") +
    ggpubr::stat_compare_means(method = "anova", 
                               label.x = 1.5,
                               label.y = 13,
                               label = "p.format") +
    theme_tcga(...) + 
    facet_grid(. ~ gene_name, scales = f_scales)
  
  return(ggobj)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Box plot drawer (show expression of selected gene project-wise)
gene_boxppr <- function(data_total, gene, correct = T, cor_v = 1, aspect.ratio = 0.25, base_size = base.size) {
  purity <- data_total$Purity
  total_exp_table <- data_total$Total_table
  # get ens id
  if (!any(grepl("ENS", gene))) {
    gns <- data_total[["GNames"]] %>% 
      dplyr::filter(gene_name %in% gene) 
  }
  gns <- gns %>% filter(ensembl_gene_id %in% colnames(data_total$Total_table))
  if (dim(gns)[1] == 0) {
    stop("No genes detected")
  }
  # prepare table
  exp_table <- total_exp_table[, c(names(total_exp_table)[1:6], gns$ensembl_gene_id)] %>% 
    mutate(project_id = str_remove_all(project_id, "TCGA-")) %>%
    filter(!project_id %in% c("PAAD", "SARC", "PCPG", "THYM", "ESCA", "STAD")) # - too little samples - remove
  proj_ord <- names(table(exp_table$project_id)) %>% grep("iPSC", ., value = TRUE, invert = TRUE)
  if (correct == T) {
    exp_table <- correct_cpe(exp_table, purity, version = cor_v)
  }
  names(exp_table)[7:length(names(exp_table))] <- paste0(gns$gene_name)
  exp_table <- exp_table %>% 
    mutate(sample_type = ifelse(project_id != "Differentiating iPSC" & sample_type != "Solid Tissue Normal", "Tumor", sample_type), 
           project_id = ifelse(project_id == "Differentiating iPSC", sample_type, project_id), 
           sample_type = ifelse(sample_type %in% c("Definitive Endoderm", "Endoderm Progenitors", "iPSC", "Primitive Streak"), "Differentiating iPSC", sample_type), 
           project_id = factor(project_id, levels = c("iPSC", "Primitive Streak", "Endoderm Progenitors", "Definitive Endoderm", proj_ord))) %>%
    pivot_longer(., matches(paste(gns$gene_name, collapse = "|")), names_to = "gene", values_to = "log2CPM")
  #
  palette <- c("#ba370b", "#83a9c1", "#f88e38")
  #
  ggobj <- ggplot(data = exp_table, 
                  aes(x = project_id, 
                      y = `log2CPM`)) + 
    geom_boxplot(aes(fill = sample_type), 
                 outlier.shape = NA) + 
    #geom_hline(yintercept = 0, color = "black", size = 0.5) +
    scale_fill_manual(values = palette) +
    facet_grid(gene ~ .) + 
    theme_bw(base_family = "sans") + 
    theme(text = element_text(colour = "#333333", hjust = 1),
          aspect.ratio = aspect.ratio,
          legend.position = "right",
          #legend.direction = legend.direction,
          legend.background = element_rect(fill = alpha('white', 0.1)),
          axis.text.y = element_text(angle = 45, margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
          axis.text.x = element_text(angle = 45, margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
          axis.line = element_line(colour = "black"),
          axis.ticks.length = unit(-0.05, "in"),
          panel.grid.major.x = element_line(size = 0.25, color = "grey40"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(), 
          panel.background = element_rect(fill = "white", 
                                          colour = "white",
                                          size = 0.5, 
                                          linetype = "solid"), 
          strip.background = element_rect(colour = "grey30", fill = "white")) + 
    ylab("log2(1 + cpm)") + xlab("") + labs(fill = "Sample type:")
  return(ggobj)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Loadings table drawer
plot_loadings <- function(exp_table, purity, pc = "PC1", n = 10, clust = "esc") {
  
  require(ggpubr)
  # prepare table
  exp_table_cor <- correct_cpe(exp_table, purity, version = 1)
  pca_obj <- do_pca(exp_table_cor)
  
  # conditionally revert PC1
  if (!is_empty(grep("ENSG00000004897", names(exp_table_cor)))) {
    pca_obj[[3]][, 1] <-  pca_obj[[3]][, 1] * -1 #--------- Revert PC1 if prolif
  } 
  
  # select direction
  cf <- ifelse(clust == "esc", 1, -1)
  pca_obj[[3]][, 1] <-  pca_obj[[3]][, 1] * cf
  
  
  # sort loadings
  loadings <- pca_obj[[3]] %>% 
    as_tibble() %>% 
    mutate(GeneID = rownames(pca_obj[[3]])) %>% 
    dplyr::select(GeneID, !! sym(pc)) %>% 
    arrange(!! sym(pc)) %>% 
    dplyr::slice(1:n)
  
  names(loadings) <- c("GeneID", paste0(pc, "_loadings"))
  
  # get gene names
  require(org.Hs.eg.db, quietly = TRUE, warn.conflicts = FALSE)
  gene_names <- invisible(mapIds(org.Hs.eg.db,
                                 keys = loadings$GeneID,
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first"))
  gene_names <- tibble(GeneID = names(gene_names), GeneName = unname(gene_names))
  loadings %<>% left_join(., gene_names, by = "GeneID") %>% dplyr::select(GeneID, GeneName, !! sym(paste0(pc, "_loadings")))
  # output
  tab <- ggtexttable(x = loadings, 
                     rows = NULL, 
                     theme = ttheme("minimal"))
  
  return(tab)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Plot results of Tukey HSD
plot_tuk <- function(tab, signature, ...) {
  order <- tab %>% arrange(diff) %>% pull(pair)
  tuk_tab <- tab %>% mutate(pair = factor(pair, levels = order))
  ggplot(data = tuk_tab, aes(x = pair, y = diff, fill = sf)) + 
    geom_errorbar(aes(ymin = lwr, ymax = upr)) +
    geom_point(shape = 21, size = 3) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    theme_bw(...) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(colour = "#49463C", hjust = 1, size = 12), 
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(face = "italic", size = 12, colour = "grey40"),
          axis.text.y = element_text(angle = 0),
          axis.text.x = element_blank(),
          legend.position = "right") +
    xlab("Comparison pair") +
    ylab("Between-tumor difference in mean PC1") + labs(fill = "Significance:", 
                                                        title = paste0(signature, " signature"), 
                                                        subtitle = "Post-hoc Tukey comparison of means")
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## PCA plot for lung tumors
plot_pair_main <- function(exp_table, pcs = c('PC1', 'PC2'), pair, ...) {
  #
  cmp <- ifelse(grepl(pattern = "brain", x = pair, ignore.case = T), c("LGG", "GBM"), c("LUAD", "LUSC"))
  #
  pca_list <- do_pca(exp_table)
  pca_eigvals <- pca_list[[2]]
  #
  plot_tab <- pca_list[[1]] %>%
    mutate(project_id = str_remove_all(project_id, "TCGA-"), 
           col_var = ifelse(project_id %in% c("LUAD", "LUSC", "Differentiating iPSC"), project_id, "else"), 
           col_var = ifelse(project_id %in% c("LUAD", "LUSC") & sample_type == "Solid Tissue Normal", "Normal tissue", col_var), 
           col_var = factor(col_var, levels = c("Differentiating iPSC", "LUSC", "LUAD", "Normal tissue", "else")))
  #
  ggplot(data = plot_tab, 
         aes(x = PC1, y = PC2)) + 
    geom_point(data = plot_tab %>% filter(col_var == "else"), 
               alpha = 0.2, 
               color = "#999999", 
               size = 0.5) + 
    geom_point(data = plot_tab %>% filter(!col_var %in% c("else", "Differentiating iPSC")), 
               aes(color = col_var, 
                   shape = col_var), 
               alpha = 0.8, 
               size = 1) +
    scale_shape_manual(values = c(16, 16, 15), guide = FALSE) +
    scale_color_manual(values = c("#fc4a1a", "#F8AA38", "#83a9c1")) + 
    geom_point(data = plot_tab %>% filter(col_var == "Differentiating iPSC"), 
               aes(fill = rep(c(1, 2, 3, 4), each = 8)[2:32]),
               color = "#666666",
               shape = 25,
               alpha = 1, 
               size = 2) + 
    scale_fill_gradientn("iPSC differentiation day: ", 
                         colours = c("red", "blue")) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    labs(x = paste("PC", as.numeric(str_sub(pcs[1], 3)), round((pca_eigvals[as.numeric(str_sub(pcs[1], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " "), 
         y = paste("PC", as.numeric(str_sub(pcs[2], 3)), round((pca_eigvals[as.numeric(str_sub(pcs[2], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " "), 
         color = "Sample type") + 
    theme_tcga(...)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## gradient PCA plot for lung tumors
plot_lung_gradient <- function(data_total, mset, gene, correct = T, cor_v = 1, pcs = c("PC1", "PC2"), ...) {
  #
  if (grepl("stem", mset)) {
    mset <- "stem"
  } 
  if (grepl("prolif", mset)) {
    mset <- "prolif"
  }
  if (grepl("emt", mset)) {
    mset <- "emt"
  }
  #
  exp_table <- data_total$Norm_Counts[[paste0(mset, "_counts_table")]]
  purity <- data_total$Purity
  # get ens id
  if (!any(grepl("ENS", gene))) {
    gns <- data_total[["GNames"]] %>% 
      dplyr::filter(gene_name %in% gene)
    if (purrr::is_null(data_total$Total_table)) {
      # check
      gns <- gns %>% filter(ensembl_gene_id %in% colnames(exp_table))
      if (dim(gns)[1] == 0) {
        stop("No genes detected")
      }
    } else {
      gns <- gns %>% filter(ensembl_gene_id %in% colnames(data_total$Total_table))
      if (dim(gns)[1] == 0) {
        stop("No genes detected")
      }
    }
  }
  # prepare table
  if (correct == TRUE) {
    exp_table <- correct_cpe(exp_table, purity, version = cor_v)
    ## get expression levels of genes of interest
    if (purrr::is_null(data_total$Total_table)) {
      gns_table <- exp_table[, c("cases", gns$ensembl_gene_id)]
    } else {
      gns_table <- data_total$Total_table[, c("cases", "project_id", "sample_type", "site_of_resection_or_biopsy", "primary_diagnosis", "gender", gns$ensembl_gene_id)]
      gns_table <- correct_cpe(gns_table, purity, version = cor_v)
      gns_table <- gns_table %>% dplyr::select(cases, contains("ENS"))
    }
  } else {
    if (purrr::is_null(data_total$Total_table)) {
      gns_table <- exp_table[, c("cases", gns$ensembl_gene_id)]
    } else {
      gns_table <- data_total$Total_table[, c("cases", gns$ensembl_gene_id)]
    }
  }
  pca_obj <- do_pca(exp_table)
  pca_obj[[1]] %<>% inner_join(., gns_table, by = "cases") %>% dplyr::rename(!!sym(gns$gene_name) := !!sym(gns$ensembl_gene_id))
  pca_df <- pca_obj[[1]]
  pca_eigvals <- pca_obj[[2]]
  
  pca_df <- pca_df %>%
    mutate(project_id = str_remove_all(project_id, "TCGA-"), 
           col_var = ifelse(project_id %in% c("LUAD", "LUSC", "Differentiating iPSC"), project_id, "else"), 
           col_var = ifelse(project_id %in% c("LUAD", "LUSC") & sample_type == "Solid Tissue Normal", "Lung normal tissue", col_var), 
           col_var = factor(col_var, levels = c("Differentiating iPSC", "LUSC", "LUAD", "Lung normal tissue", "else")))
  
  # organise for plotting
  ggobj <- ggplot(data = pca_df, 
                  aes(x = !! sym(pcs[1]), 
                      y = !! sym(pcs[2]))) +
    geom_point(data = pca_df %>% filter(col_var == "else"), 
               alpha = 0.2, 
               color = "#999999", 
               size = 0.5) +
    geom_point(data = pca_df %>% filter(col_var %in% c("Lung normal tissue", "Differentiating iPSC")), 
               aes(fill = col_var,
                   shape = col_var),
               #color = "#333333",
               alpha = 0.8, 
               size = 2) +
    scale_shape_manual(values = c(25, 22), guide = FALSE) +
    scale_fill_manual(values = c("red", "#83a9c1")) +
    geom_point(data = pca_df %>% filter(col_var %in% c("LUSC", "LUAD")), 
               aes(color = !!sym(gns$gene_name)),
               alpha = 0.7, 
               size = 1) + 
    #scale_color_gradientn(paste0("Log2CPM of ", names(enid)), 
    #                      colours = c("blue", "red")) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 3.5, alpha = 1))) + 
    labs(x = paste("PC", as.numeric(str_sub(pcs[1], 3)), round((pca_eigvals[as.numeric(str_sub(pcs[1], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " "), 
         y = paste("PC", as.numeric(str_sub(pcs[2], 3)), round((pca_eigvals[as.numeric(str_sub(pcs[2], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " "), 
         fill = "Sample type: ", 
         color = paste0(gene, ", log2(CPM)")) + 
    theme_tcga(...)
  return(ggobj)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Plot tukey hsd for tumor stages as barplots (deprecated)
plot_diff_bar <- function(emm_out, aspect.ratio = 0.7, ...) {
  means <- emm_out$`emmeans of project_id, tumor_stage_simp` %>% as_tibble()
  diffs <- emm_out$`pairwise differences of project_id, tumor_stage_simp` %>% as_tibble()
  diffs %<>% mutate(pair_l = str_split_fixed(contrast, " - ", 2)[, 1], 
                    pair_r = str_split_fixed(contrast, " - ", 2)[, 2], 
                    stage_l = str_split_fixed(pair_l, ",", 2)[, 2], 
                    stage_r = str_split_fixed(pair_r, ",", 2)[, 2], 
                    pair_l = str_split_fixed(pair_l, ",", 2)[, 1], 
                    pair_r = str_split_fixed(pair_r, ",", 2)[, 1], 
                    ssign = ifelse(p.value < 0.05, TRUE, FALSE)) %>% 
    filter(pair_l == pair_r, 
           stage_l == 1) %>% 
    dplyr::select(-contrast)
  
  ggplot(data = diffs, 
         aes(x = pair_l, 
             y = estimate)) + 
    geom_col(aes(fill = stage_r, 
                 color = ssign), 
             position = "dodge2", 
             size = 0.5) + 
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30", size = 1.2) +
    geom_errorbar(aes(ymin = estimate - SE, 
                      ymax = estimate + SE, 
                      group = stage_r), 
                  position = "dodge") + 
    scale_color_manual(values = c("#999999", "red")) +
    scale_fill_grey() + 
    theme_tcga(aspect.ratio = aspect.ratio, base_size = 14) + 
    ylab("Difference from Stage 1, -1 * PC1 coordinate") + xlab("") +
    labs(fill = "Tumor stage", color = "Significant:")
}

