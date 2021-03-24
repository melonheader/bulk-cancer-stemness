### INFO: Functions for Stem-TCGA project
### DATE: 20.09.2019
### AUTHOR: Artem Baranovskii



# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# Plotting domain
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Main theme for ggplot
theme_tcga <- function(base_family = "sans", legend.position = "right", legend.direction = "vertical", aspect.ratio = 1, ...){
  theme_bw(base_family = base_family, ...) + 
    theme(text = element_text(colour = "#333333", hjust = 1),
          aspect.ratio = aspect.ratio,
          legend.position = legend.position,
          legend.direction = legend.direction,
          legend.background = element_rect(fill = alpha('white', 0.1)),
          axis.text.y = element_text(angle = 45, margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
          axis.text.x = element_text(angle = 45, margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")),
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
  #
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
                 alpha = 0.8, 
                 size = 0.4) +
      scale_fill_viridis_c(guide = FALSE) + 
      scale_shape_manual("Sample type: ", 
                         values = c(16, 17)) +
      scale_color_manual("Sample type: ", 
                         values = c("#FF1493", "grey70"),) +
      labs(x = paste("PC", as.numeric(str_sub(pcs[1], 3)), round((pca_eigvals[as.numeric(str_sub(pcs[1], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " "), 
           y = paste("PC", as.numeric(str_sub(pcs[2], 3)), round((pca_eigvals[as.numeric(str_sub(pcs[2], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " ")) + 
      theme_tcga(...)
    return(ggobj)
  } else {
    ggobj <- ggplot(data = plot_tab, 
                    aes(x = !! sym(pcs[1]), 
                        y = !! sym(pcs[2]))) +
      geom_point(data = plot_tab %>% filter(col_var != "iPSC"), 
                 aes(color = col_var), 
                 alpha = 0.4, 
                 size = 1.0) +
      geom_point(data = plot_tab %>% filter(col_var == "iPSC"), 
                 aes(fill = rep(c(1, 2, 3, 4), each = 8)[2:32]),
                 shape = 24,
                 alpha = 0.8, 
                 size = 2, 
                 color = "#666666") + 
      scale_color_manual("Sample type: ", 
                         values = c("#f88e38", "#83a9c1")) +
      scale_fill_gradientn("iPSC differentiation stage: ", 
                           colours = c("red", "blue")) +
      labs(x = paste(pcs[1], round((pca_eigvals[as.numeric(str_sub(pcs[1], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " "), 
           y = paste(pcs[2], round((pca_eigvals[as.numeric(str_sub(pcs[2], 3))] / sum(pca_eigvals))*100, 0), "%", sep = " ")) + 
      theme_tcga(...)
    return(ggobj)
  }
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Jitter-plot drawer
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
    summarise(count = n()) %>% 
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
                aes(x = project_id, y = !! sym(paste0("PC", pc)), group = project_id, col = !! sym(table_clust[[5]][1])), 
                size = 1, 
                alpha = 0.6) + 
    scale_color_manual(values = c("#CD5C5C", "#fc4a1a", "#F8AA38", "#83a9c1", "#f7b733")) + 
    geom_hline(yintercept = median_tum, size = 1, color = "grey30") +  
    coord_flip() +
    theme_tcga(...) +
    labs(subtitle = paste0(table_clust[[5]][1], "\nsignature"),
         color = "Legend:", 
         y = paste0("PC", pc), 
         x = "Group name")
  return(ggobj)
  
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Scatter-plot drawer
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
        summarise(!! sym(paste0(str_split_fixed(pca_tables_list[[i]][[1]][1, 1], "_", 2)[, 1], "_Intensity")) := median(PC1))
    }
    medians_table <- purrr::reduce(medians_list, inner_join, by = "project_id") %>% 
      mutate(project_id = str_split_fixed(project_id, "-", 2)[, 2])
    purity_proj <- purity %>% 
      group_by(project_id) %>% 
      summarise(median_purity = median(purity, na.rm = TRUE)) %>% 
      mutate(project_id = str_split_fixed(project_id, "-", 2)[, 2])
    cluster_freq <- medians_table %>% 
      inner_join(., purity_proj, 
                 by = "project_id") %>% 
      mutate(purity_factorial = ifelse(median_purity < 0.72, "Low", 
                                       ifelse(median_purity >= 0.72 & median_purity < 0.83, "Medium", 
                                              ifelse(median_purity >= 0.83, "High", median_purity)))) %>% 
      mutate(purity_factorial = factor(purity_factorial, levels = c("Low", "Medium", "High")))
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
      summarise(stage_avg = weighted.mean(x = stage_simp, w = (count / sum(count)))) %>% 
      ungroup() %>% 
      mutate(stage_avg_factorial = ifelse(stage_avg < 2, "Low (1-2)", 
                                          ifelse(stage_avg >= 2 & stage_avg < 3, "Medium (2-3)", 
                                                 ifelse(stage_avg >= 3, "High (3-4)", stage_avg))), 
             stage_avg_factorial = factor(stage_avg_factorial, levels = c("Low (1-2)", "Medium (2-3)", "High (3-4)")))
    ## Joim 
    cluster_freq %<>% 
      left_join(., avg_stage, 
                by = "project_id")
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
      summarise(median_purity = median(purity, na.rm = TRUE)) %>% 
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
      summarise(stage_avg = weighted.mean(x = stage_simp, w = (count / sum(count)))) %>% 
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
    ggrepel::geom_label_repel(force = 20) + 
    xlim(c(0, 1)) + ylim(c(0, 1)) + 
    theme_tcga(...) +
  labs(color = ifelse(col_var == "purity", "Median\nproject-wise\npurity", 
                      ifelse(col_var == "stage_avg", "Average\nproject-wise\ntumor stage")), 
       shape = ifelse(shape_var == "purity", "Median\nproject-wise\npurity", 
                      ifelse(shape_var == "stage_avg", "Average\nproject-wise\ntumor stage", 
                             ifelse(is.na(shape_var), NA))))
  return(sc1)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Heatmap drawer
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
        summarise(!! sym(paste0(str_split_fixed(pca_tables_list[[i]][[1]][1, 1], "_", 2)[, 1], "_Intensity")) := median(PC1))
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
    mutate(project_id = factor(project_id, levels = hc1$labels[hc1$order]))
  #Nice heatmap
  hm1 <- ggplot(data = cluster_freq_tall, 
                aes(x = project_id, 
                    y = Signature)) + 
    geom_tile(aes(fill = Intensity)) + 
    scale_fill_viridis_c() +
    labs(fill = "Signature\nintensity") +
    theme_tcga(...)
  return(hm1)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Cox-reg boxplot drawer 
plot_surv_boxplots <- function(cluster_table, meta, ...) {
  #change project names
  cluster_table %<>% mutate(project_id = str_remove_all(project_id, "TCGA-"))
  meta %<>% mutate(project_id = str_remove_all(project_id, "TCGA-"))
  #Clustering as an independent predictor of survival
  cluster <- names(cluster_table[4])
  names(cluster_table)[4] <- "Cluster"
  # Save quadrant in analysis
  quadrant <- cluster_table[["quadrant"]] %>% unique()
  
  
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
    mutate(`Risk of Adverse Event` = ifelse(`lower .95` > 1, 
                                            "Higher Risk", 
                                            ifelse(`upper .95` < 1, 
                                                   "Lower Risk", 
                                                   "Not significant")), 
           `Risk of Adverse Event` = factor(`Risk of Adverse Event`, 
                                            levels = c("Higher Risk", "Lower Risk", "Not significant"))) %>%
    mutate(Significance = ifelse(`Pr(>|z|)` < 0.05, 
                                 "*", 
                                 "NS"), 
           Significance = ifelse(`Pr(>|z|)` < 0.005, 
                                 "**", 
                                 Significance), 
           Significance = ifelse(`Pr(>|z|)` < 0.001, 
                                 "***", 
                                 Significance))
  #Prepare palette
  palette <- c("#F8AA38","#5A92B7", "grey80")
  
  #Plot
  ggobj <- ggplot(res.surv, 
                  aes(y = `exp(coef)`, 
                      x = TCGA_project, 
                      color = `Risk of Adverse Event`)
                  ) + 
    geom_point(shape = 15, size =3) +
    geom_errorbar(aes(ymin = `lower .95`, ymax = `upper .95`), 
                  width = .2, 
                  position = position_dodge(0.05)
                  ) + 
    scale_color_manual(values = palette, drop = FALSE) +
    geom_text(aes(label = Significance, x = TCGA_project, y = `upper .95`), 
              hjust = 0.5, 
              nudge_x = 0,
              nudge_y = 0.25,
              size = 3, 
              color = "black") +
    geom_hline(yintercept = 1, color = "red", size = 0.5) + 
    ylim(c(0, 10)) +
    #coord_flip(ylim = c(0, 10)) +
    theme_tcga(...) +
    xlab(NULL) +
    ylab("Hazard Ratio") #+
    #labs(subtitle = paste0(cluster, " signature; ", quadrant))
  return(ggobj)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Cox-reg curves drawer
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
                                      data = data.surv, 
                                      pval = TRUE, 
                                      censor = TRUE, 
                                      palette = c("#5A92B7", "#F8AA38"),  
                                      legend.title = "Cluster:", 
                                      legend.labs = c("Normal", 
                                                      "ESC"),  
                                      xlab = "Time, hours",
                                      subtitle = paste0(project, " ", cluster, " cluster")
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



# project-wise nice correlation plots PC1 vs Purity
# ---------------------------------------------------------------------------- #
plot_nice_corr <- function(pca_table, meta, purity, pc = 1, project = "global", base_size = 11) {
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

