### INFO: Functions for Stem-TCGA project
### DATE: 23.12.2019
### AUTHOR: Artem Baranovskii




# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# Misc. functions
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Spread table keeping key in new colnames
myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Normalize to 0-1 scale
norm_to_one <- function(x) {
  x_norm <- (x - min(x)) / (max(x) - min(x))
  return(x_norm)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Compute anova and extract differences from Tukey HSD
anova_tuk <- function(pca_tabs_list, signature = "stemness", aspect.ratio = 0.2, ...) {
  # select signature
  choice <- names(pca_tabs_list) %>% grep(str_sub(signature, 1, 4), ., ignore.case = TRUE, value = TRUE)
  ## one-way ANOVA for PC1 and projects
  tab <- pca_tabs_list[[choice]][[1]] %>% dplyr::select(2:9)
  # prs to keep
  prs_keep <- janitor::tabyl(tab, project_id) %>% filter(n > 60) %>% pull(project_id)
  # filter tab
  tab %<>% filter(project_id %in% prs_keep, 
                  sample_type != "Solid Tissue Normal")
  # anova
  fit <- aov(PC1 ~ project_id, data = tab)
  
  # tukey post-hoc
  tuk <- TukeyHSD(fit, "project_id", ordered = TRUE)
  tuk_tab <- as_tibble(tuk$project_id) %>% 
    mutate(pair = rownames(tuk$project_id), 
           sf = ifelse(`p adj` < 0.05, "Significant", "NS"), 
           pair = str_remove_all(pair, "TCGA-")) %>% 
    dplyr::select(pair, sf, everything())

  
  # output
  out <- list(fit, tuk_tab)
  
  return(out)
}


# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# Data processing
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Load preprocessed data
# Get annotation
load_data <- function(data_path = file.path("/local", "artem", "Projects", "Draft", "Data"), short = TRUE) {
  
  # Load normalized counts  -------------------------------------------
  if (short == FALSE) {
    message("loading total table.....")
    total_table <- read_tsv(file.path(data_path, "Counts", "total_norm_counts_t.tsv"), 
                            col_names = TRUE, 
                            col_types = cols(.default = col_double(), 
                                             cases = col_character(), 
                                             project_id = col_character(), 
                                             sample_type = col_character(), 
                                             site_of_resection_or_biopsy = col_character(), 
                                             primary_diagnosis = col_character(), 
                                             gender = col_character()))
    message("loading subsetted tables.....")
    t_paths <- c("stem_counts_table.tsv", 
                 "prolif_counts_table.tsv",
                 "emt_counts_table.tsv")
    exp_tables <- purrr::map(t_paths, function(x) read_tsv(file = file.path(data_path, "Counts", x), 
                                                           col_names = TRUE, 
                                                           col_types = cols(.default = col_double(), 
                                                                            cases = col_character(), 
                                                                            project_id = col_character(), 
                                                                            sample_type = col_character(), 
                                                                            site_of_resection_or_biopsy = col_character(), 
                                                                            primary_diagnosis = col_character(), 
                                                                            gender = col_character()))) %>% 
      set_names(., str_remove_all(t_paths, ".tsv"))
  } else {
    if (short == TRUE) {
      total_table <- NULL
      message("loading subsetted tables.....")
      t_paths <- c("stem_counts_table.tsv", 
                   "prolif_counts_table.tsv",
                   "emt_counts_table.tsv")
      exp_tables <- purrr::map(t_paths, function(x) read_tsv(file = file.path(data_path, "Counts", x), 
                                                             col_names = TRUE, 
                                                             col_types = cols(.default = col_double(), 
                                                                              cases = col_character(), 
                                                                              project_id = col_character(), 
                                                                              sample_type = col_character(), 
                                                                              site_of_resection_or_biopsy = col_character(), 
                                                                              primary_diagnosis = col_character(), 
                                                                              gender = col_character()))) %>% 
        set_names(., str_remove_all(t_paths, ".tsv"))
    }
  }
  
  
  # Meta ---------------------------------------------------------------
  message("loading metadata.....")
  meta.tcga <- read_tsv(file = file.path(data_path, "Meta", "Meta.tsv"), 
                        col_names = TRUE, 
                        col_types = cols(.default = col_character()))
  meta.de <- read_tsv(file = file.path(data_path, "Meta", "Meta_iPSC.tsv"),
                      col_names = TRUE, 
                      col_types = cols(.default = col_character()))
  meta <- rbind(meta.de, meta.tcga)
  
  # TCGA auxiliary -----------------------------------------------------
  pairs <- read_tsv(file = file.path(data_path, "Meta", "Paired_Samples.tsv"),  
                    col_names = TRUE, 
                    col_types = cols(.default = col_character()))
  purity <- read_tsv(file = file.path(data_path, "Meta", "Purity.tsv"),
                     col_names = TRUE, 
                     col_types = cols(.default = col_character(), 
                                      purity = col_double())) 
  markers <- read_tsv(file = file.path(data_path, "Meta", "Markers.tsv"),
                      col_names = TRUE, 
                      col_types = cols(.default = col_integer(), 
                                       gene_id = col_character(), 
                                       gene_name = col_character()))
  # H-sapiens gene names ------------------------------------------------
  hs_genes <- read_tsv(file = file.path(data_path, "Meta", "Annotation", "hs_genes.tsv"), 
                       col_names = TRUE, 
                       col_types = cols(.default = col_character()))
  # output
  message("done!")
  data_output <- list(Total_table = total_table,
                      Norm_Counts = exp_tables, 
                      Markers = markers, 
                      Meta = meta, 
                      Purity = purity, 
                      Pairs = pairs,
                      GNames = hs_genes)
  return(data_output)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Process counts
prepare_counts_pr <- function(projects, meta, cut_bot = -2, cut_up = 9) {
  # paths
  # ---------------------------------------------------------------------------- #
  f_path <- file.path("~/Scripts/Scripts/Draft")
  counts_path <- file.path("/local", "artem", "Projects", "Draft", "Data", "Counts")
  meta_path <- file.path("/local", "artem", "Projects", "Draft", "Data", "Meta")
  gdc_path <- file.path(counts_path, "GDCprocessed")
  # libraries
  # ---------------------------------------------------------------------------- #
  require(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
  require(magrittr, quietly = TRUE, warn.conflicts = FALSE)
  require(edgeR, quietly = TRUE, warn.conflicts = FALSE)
  #
  if (isSingleString(projects)) {
    projects %<>% str_squish(.) %>% str_split_fixed(., ",", Inf) %>% as.character() %>% paste(., collapse = "|")
  } 
  
  # load GDC
  # ---------------------------------------------------------------------------- #
  message("loading GDC counts.....")
  gdc_raw_counts <- purrr::map(list.files(path = gdc_path, pattern = projects, full.names = TRUE), 
                               function(x) read_tsv(x, 
                                                    col_names = TRUE, 
                                                    col_types = cols(.default = col_double(), 
                                                                     GeneID = col_character()))) %>% 
    purrr::reduce(., full_join, by = "GeneID")
  # total
  # ---------------------------------------------------------------------------- #
  # prepare matrix for normalisation
  counts_matrix <- gdc_raw_counts %>% dplyr::select(-GeneID) %>% as.data.frame()
  rownames(counts_matrix) <- gdc_raw_counts %>% pull(GeneID)
  # quality
  message("quality check.....")
  cpm.Log <- cpm(gdc_raw_counts %>% dplyr::select(-GeneID), 
                 log = TRUE)
  median.log2.cpm <- apply(cpm.Log, 1, median) %>%
    tibble::enframe(name = NULL) %>% 
    mutate(GeneID = gdc_raw_counts$GeneID)
  # normalize with edgeR
  message("normalizing.....")
  y <- DGEList(counts_matrix, genes = rownames(counts_matrix))
  y <- calcNormFactors(y, method = "TMM")
  counts_matrix_norm <- cpm(y, 
                            normalized.lib.sizes = TRUE, 
                            log = TRUE, 
                            prior.count = 1) %>% 
    as_tibble() %>%
    mutate(GeneID = y$genes$genes) %>%
    dplyr::select(GeneID, everything()) 
  # subset signature tables
  message("subsetting tables.....")
  meta <- purrr::map(c("Meta.tsv", "Meta_iPSC.tsv"), 
                     function(x) read_tsv(file = file.path(meta_path, x), 
                                          col_names = TRUE, 
                                          col_types = cols(.default = col_character()))) %>% 
    purrr::reduce(., base::rbind)
  ## read in markers
  markers <- read_tsv(file = file.path(meta_path, "Markers.tsv"), 
                      col_names = TRUE, 
                      col_types = cols(gene_id = col_character(), 
                                       gene_name = col_character(),
                                       Stem_Markers_DP = col_double(),
                                       Prolif_Markers = col_double(),
                                       Emt_Met_Markers = col_double()
                      ))
  exp_tables <- purrr::map(c("Stem", "Prolif", "Emt"), 
                           function(x) subset_exp(counts_matrix_norm, 
                                                  markers = markers, 
                                                  meta = meta, 
                                                  gene_set = x)) 
  return(exp_tables)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Normalize counts table
normalize_counts <- function(counts_table) {
  #Prepare groups
  Group_1 <- character()
  for (i in 1:length(names(counts_table[, -1]))) {
    Group_1[i] <- meta %>% 
      filter(cases == names(counts_table[, -1][i])) %>%
      pull(project_id)
  }
  #factorize
  Group_1 <- factor(Group_1)
  #Prepare DGE
  counts.matrix <- as.matrix(counts_table[, -1]) %>% 
    `rownames<-`(counts_table$Name)
  y <- DGEList(counts = counts.matrix, 
               genes = rownames(counts.matrix),
               group = Group_1)
  #Filter out lowly expressed genes
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  #Normalize libraries 
  y <- calcNormFactors(y)
  #Get output matrix
  output_table<- cpm(y, 
                     normalized.lib.sizes = TRUE, 
                     log = TRUE, 
                     prior.count = 1) %>%
    as.tibble() %>%
    mutate(GeneID = y$genes$genes) %>%
    dplyr::select(GeneID, everything())
  #Return
  return(output_table)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Subset normalized counts table by set of genes
subset_exp <- function(exp_table, meta, markers, gene_set = 'Stem') {
  if (is.data.frame(markers)) {
    #Get set name
    name <- rlang::sym(names(markers[, grepl(gene_set, names(markers))]))
    #Get genes
    genes <- markers %>% 
      dplyr::filter(!! name == 1) %>% 
      pull(gene_id)
  } else {
    if (!any(stringr::str_detect("ENS", markers))) {
      # get gene ann tab
      genes <- read_tsv("/local/artem/Projects/Draft/Data/Meta/Annotation/hs_genes.tsv", 
                        col_types = cols(
                          ensembl_gene_id = col_character(),
                          gene_name = col_character(),
                          gene_biotype = col_character()
                        )) %>% 
        dplyr::filter(gene_name %in% markers) %>% 
       pull(ensembl_gene_id)
    } else {
      genes <- markers
    }
  }
  # Get shortcut for metadata columns
  names_meta <- rlang::syms(names(meta))
  #Filter genes and transpose
  exp_table <- exp_table %>% 
    dplyr::filter(GeneID %in% genes) %>% 
    gather(., 
           2:ncol(.), 
           key = 'cases', 
           value = 'Counts') %>%
    spread(GeneID, Counts) 
  #Annotate
  exp_table <- exp_table %>% 
    inner_join(., meta, 
               by = c('cases')) %>% 
    dplyr::select(cases, !! names_meta[[2]], !! names_meta[[3]], !! names_meta[[4]], !! names_meta[[5]], !! names_meta[[8]], contains('ENS'))
  return(exp_table)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Regress effect of tumor purity from normlaized counts
correct_cpe <- function(subsetted_exp, purity, correct = TRUE, version = 1, keep_all = FALSE) { 
  #______________________________________APPROVED_____________________________________#
  #--------------------------------------------------#
  #Regress the effect of purity                      #
  #Compute global lm for all TCGA-Projects           #
  #using gene expression                             #
  #--------------------------------------------------#--------------------------------#
  # Output: First three PCs computed from gene expressino matrix corrected for purity #
  #-----------------------------------------------------------------------------------#
  if (version == 1) { 
    
    data.tum.list <- list()
    #Load data and subset
    if (keep_all == TRUE) {
      data.tum.reg <- subsetted_exp %>% 
        filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor")) %>%
        right_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>% 
        distinct(cases, .keep_all = TRUE)
    } else {
      data.tum.reg <- subsetted_exp %>% 
        filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor")) %>%
        inner_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>% 
        distinct(cases, .keep_all = TRUE)
    }
    
    # Artificialy assign purity = 0.05 to normal tissues and correct them as well
    data.norm.reg <- subsetted_exp %>% 
      filter(sample_type %in% c("Solid Tissue Normal")) %>%
      mutate(purity = rnorm(length(.$cases), 
                            quantile(purity$purity, na.rm = TRUE, probs = 0.000005), 
                            quantile(purity$purity, na.rm = TRUE, probs = 0.000005) / 3)
             ) %>%
      distinct(cases, .keep_all = TRUE)
    
    # Artificialy assign purity = 0.05 to iPSCs and correct them as well
    data.ipsc.reg <- subsetted_exp %>% 
      filter(project_id == "Differentiating iPSC") %>%
      mutate(purity = rnorm(length(.$cases), 
                            quantile(purity$purity, na.rm = TRUE, probs = 0.005), 
                            quantile(purity$purity, na.rm = TRUE, probs = 0.005) / 3)
             ) %>%
      distinct(cases, .keep_all = TRUE)
    
    #-------------------#
    #    Correct CPE    #
    #-------------------#
    if (correct == TRUE) {
      #Correct gene expression
      for (x in 1:length(names(data.tum.reg %>% dplyr::select(starts_with("ENSG"))))) {
        lm.fit <- lm(data.tum.reg[[x + 7]] ~ purity, 
                     data = data.tum.reg, 
                     na.action = "na.exclude"
        )
        # Correct tumor samples
        data.tum.reg[[(x + 7)]] <- data.tum.reg[[(x + 7)]] - data.tum.reg[['purity']] * lm.fit$coefficients[[2]]     #------changed +5 to +7 because of new meta
        # Artificially correct normal tissues
        data.norm.reg[[(x + 7)]] <- data.norm.reg[[(x + 7)]] - data.norm.reg[['purity']] * lm.fit$coefficients[[2]]
        # Artificailly correct iPSCs
        data.ipsc.reg[[(x + 7)]] <- data.ipsc.reg[[(x + 7)]] - data.ipsc.reg[['purity']] * lm.fit$coefficients[[2]]
      }
      rm(x)
    }
    
    
    
    #Output
    data.reg.full <- rbind(data.tum.reg, 
                           data.norm.reg, 
                           data.ipsc.reg) %>% 
      dplyr::select(-purity) %>% 
      distinct(cases, .keep_all = TRUE)
  }
  
  #______________________________________UNSURE__________________________________________#
  #--------------------------------------------------#
  #Regress the effect of purity                      #
  #Compute lm projectwise for every TCGA-Project     #
  #using gene expression                             #
  #use an error variable as a corrected estimator    #
  #of gene expression                                #
  #--------------------------------------------------#--------#
  #Output: Gene expression matrix with CPM corrected for CPE  #
  #-----------------------------------------------------------#
  if (version == 5) { 
    
    #Extract vector projects
    projects <- purity %>%
      distinct(project_id) %>%
      pull(project_id)
    
    #------------------#
    #Correct tumors    #
    #------------------#
    #Correct tumor samples for CPE 
    data.tum.list <- list()
    
    #Load data
    data <- subsetted_exp %>% 
      filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor"))
    for (i in 1:length(projects)) { 
      #Prepare data
      data.tum.reg <- data %>%  
        inner_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        filter(project_id == projects[i]) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>%
        distinct(cases, .keep_all = TRUE)
      
      #Skip if empty
      if(dim(data.tum.reg)[1] == 0) { 
        next
      }
      
      #Correct gene expression
      for (x in 1:length(names(data.tum.reg %>% dplyr::select(starts_with("ENSG"))))) {
        lm.fit <- lm(data.tum.reg[[x + 7]] ~ purity, 
                     data = data.tum.reg, 
                     na.action = "na.exclude"
        )
        data.tum.reg[[(x + 7)]] <- data.tum.reg[[(x + 7)]] - data.tum.reg[['purity']] * lm.fit$coefficients[[2]] 
      }
      rm(x)
      #Store
      data.tum.list[[projects[i]]] <- data.tum.reg
    }
    
    #------------------#
    #Correct normal    #
    #------------------#
    #As we substracted only B*x correction for normal tissues and stem cells is not required
    #Load data
    data.norm.reg <- subsetted_exp %>% 
      filter(sample_type %in% c("Solid Tissue Normal")) %>%
      mutate(purity = runif(length(.$cases), 0.05, 0.02)) %>%
      distinct(cases, .keep_all = TRUE)
    
    #------------------#
    #Correct iPSC      #
    #------------------#
    #Artificialy assign purity = 1 for all stem cell samples and correct them as well
    #Load data
    data.ipsc.reg <- subsetted_exp %>% 
      filter(project_id == "Differentiating iPSC") %>%
      mutate(purity = runif(runif(length(.$cases), 0.05, 0.02))) %>%
      distinct(cases, .keep_all = TRUE)
    
    
    #Get things together
    data.tum.full <- purrr::reduce(data.tum.list, rbind)
    data.reg.full <- rbind(data.tum.full, 
                           data.norm.reg, 
                           data.ipsc.reg) %>% 
      distinct(cases, .keep_all = TRUE) %>%
      dplyr::select(-purity)
  }
  
  #______________________________________UNSURE__________________________________________#
  #--------------------------------------------------#
  #Regress the effect of purity                      #
  #Compute global lm for all TCGA-Projects           #
  #using gene expression (correct normal)            #
  #--------------------------------------------------#--------#
  #Output: Gene expression matrix with CPM corrected for CPE  #
  #-----------------------------------------------------------#
  if (version == 3) { 
    #------------------#
    #Prepare tumors    #
    #------------------#
    #Correct tumor samples for purity 
    data.tum.list <- list()
    
    #Load data
    if (keep_all == TRUE) {
      data.tum.reg <- subsetted_exp %>% 
        filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor")) %>%
        right_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>% 
        distinct(cases, .keep_all = TRUE)
    } else {
      data.tum.reg <- subsetted_exp %>% 
        filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor")) %>%
        inner_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>% 
        distinct(cases, .keep_all = TRUE)
    }
    
    #------------------#
    #Prepare normal    #
    #------------------#
    #As we substracted only B*x correction for normal tissues and stem cells is not required
    #Load data
    data.norm.reg <- subsetted_exp %>% 
      filter(sample_type %in% c("Solid Tissue Normal")) %>%
      mutate(purity = rnorm(length(.$cases), 0, 0.1)) %>%
      distinct(cases, .keep_all = TRUE) %>% 
      dplyr::select(cases, purity, everything())
    
    #------------------#
    #Prepare iPSC      #
    #------------------#
    #Artificialy assign purity = 1 for all stem cell samples and correct them as well
    #Load data
    data.ipsc.reg <- subsetted_exp %>% 
      filter(project_id == "Differentiating iPSC") %>%
      mutate(purity = rnorm(length(.$cases), 0.05, 0.02)) %>%
      distinct(cases, .keep_all = TRUE) %>% 
      dplyr::select(cases, purity, everything())
    #
    data.reg.full <- rbind(data.tum.reg, 
                           data.norm.reg, 
                           data.ipsc.reg) %>% 
      distinct(cases, .keep_all = TRUE)
    
    if (correct == TRUE) {
      #Correct gene expression
      for (x in 1:length(names(data.reg.full %>% dplyr::select(starts_with("ENSG"))))) {
        lm.fit <- lm(data.reg.full[[x + 7]] ~ purity, 
                     data = data.reg.full, 
                     na.action = "na.exclude"
        )
        data.reg.full[[(x + 7)]] <- data.reg.full[[(x + 7)]] - data.reg.full[['purity']] * lm.fit$coefficients[[2]]
        #data.norm.reg[[(x + 7)]] <- data.reg.full[[(x + 7)]] - lm.fit$coefficients[[1]]
      }
      rm(x)
    }
  }  
  
  
  #______________________________________UNSURE__________________________________________#
  #--------------------------------------------------#
  #Regress the effect of purity                      #
  #Compute global lm for all TCGA-Projects           #
  #using gene expression (correct normal)            #
  #--------------------------------------------------#--------#
  #Output: Gene expression matrix with CPM corrected for CPE  #
  #-----------------------------------------------------------#
  if (version == 2) { 
    #------------------#
    #Prepare tumors    #
    #------------------#
    #Correct tumor samples for purity 
    data.tum.list <- list()
    
    #Load data
    if (keep_all == TRUE) {
      data.tum.reg <- subsetted_exp %>% 
        filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor")) %>%
        right_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>% 
        distinct(cases, .keep_all = TRUE)
    } else {
      data.tum.reg <- subsetted_exp %>% 
        filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor")) %>%
        inner_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>% 
        distinct(cases, .keep_all = TRUE)
    }
    
    #------------------#
    #Prepare normal    #
    #------------------#
    #As we substracted only B*x correction for normal tissues and stem cells is not required
    #Load data
    data.norm.reg <- subsetted_exp %>% 
      filter(sample_type %in% c("Solid Tissue Normal")) %>%
      mutate(purity = rnorm(length(.$cases), 0.05, 0.02)) %>%
      distinct(cases, .keep_all = TRUE) %>% 
      dplyr::select(cases, purity, everything())
    
    #------------------#
    #Prepare iPSC      #
    #------------------#
    #Artificialy assign purity = 1 for all stem cell samples and correct them as well
    #Load data
    data.ipsc.reg <- subsetted_exp %>% 
      filter(project_id == "Differentiating iPSC") %>%
      mutate(purity = rnorm(length(.$cases), 0.05, 0.02)) %>%
      distinct(cases, .keep_all = TRUE) %>% 
      dplyr::select(cases, purity, everything())
    #
    data.reg.full <- rbind(data.tum.reg, 
                           data.norm.reg, 
                           data.ipsc.reg) %>% 
      distinct(cases, .keep_all = TRUE)
    
    if (correct == TRUE) {
      #Correct gene expression
      for (x in 1:length(names(data.reg.full %>% dplyr::select(starts_with("ENSG"))))) {
        lm.fit <- lm(data.tum.reg[[x + 7]] ~ purity, 
                     data = data.tum.reg, 
                     na.action = "na.exclude"
        )
        data.reg.full[[(x + 7)]] <- data.reg.full[[(x + 7)]] - data.reg.full[['purity']] * lm.fit$coefficients[[2]]
        #data.norm.reg[[(x + 7)]] <- data.reg.full[[(x + 7)]] - lm.fit$coefficients[[1]]
      }
      rm(x)
    }
  }
  
  
  #______________________________________UNSURE__________________________________________#
  #--------------------------------------------------#
  #Regress the effect of purity                      #
  #Compute global lm for all TCGA-Projects           #
  #using gene expression (correct normal)            #
  #--------------------------------------------------#--------#
  #Output: Gene expression matrix with CPM corrected for CPE  #
  #-----------------------------------------------------------#
  if (version == 4) { 
    #------------------#
    #Prepare tumors    #
    #------------------#
    #Correct tumor samples for purity 
    data.tum.list <- list()
    
    #Load data
    if (keep_all == TRUE) {
      data.tum.reg <- subsetted_exp %>% 
        filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor")) %>%
        right_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>% 
        distinct(cases, .keep_all = TRUE)
    } else {
      data.tum.reg <- subsetted_exp %>% 
        filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor")) %>%
        inner_join(purity[, c("cases", "purity")], ., by = c("cases")) %>%
        mutate(purity = as.numeric(purity)) %>% 
        filter(purity != "NaN") %>% 
        distinct(cases, .keep_all = TRUE)
    }
    
    #------------------#
    #Prepare normal    #
    #------------------#
    #As we substracted only B*x correction for normal tissues and stem cells is not required
    #Load data
    data.norm.reg <- subsetted_exp %>% 
      filter(sample_type %in% c("Solid Tissue Normal")) %>%
      mutate(purity = rnorm(length(.$cases), 0.3, 0.1)) %>%
      distinct(cases, .keep_all = TRUE) %>% 
      dplyr::select(cases, purity, everything())
    
    #------------------#
    #Prepare iPSC      #
    #------------------#
    #Artificialy assign purity = 1 for all stem cell samples and correct them as well
    #Load data
    data.ipsc.reg <- subsetted_exp %>% 
      filter(project_id == "Differentiating iPSC") %>%
      mutate(purity = rnorm(length(.$cases), 0.3, 0.1)) %>%
      distinct(cases, .keep_all = TRUE) %>% 
      dplyr::select(cases, purity, everything())
    
    
    if (correct == TRUE) {
      #Correct gene expression
      for (x in 1:length(names(data.tum.reg %>% dplyr::select(starts_with("ENSG"))))) {
        lm.fit <- lm(data.tum.reg[[x + 7]] ~ purity, 
                     data = data.tum.reg, 
                     na.action = "na.exclude"
        )
        data.tum.reg[[(x + 7)]] <- data.tum.reg[[(x + 7)]] - data.tum.reg[['purity']] * lm.fit$coefficients[[2]]
        data.norm.reg[[(x + 7)]] <- data.norm.reg[[(x + 7)]] - lm.fit$coefficients[[1]]
        data.ipsc.reg[[(x + 7)]] <- data.ipsc.reg[[(x + 7)]] - lm.fit$coefficients[[1]]
      }
      rm(x)
    }
    #
    data.reg.full <- rbind(data.tum.reg, 
                           data.norm.reg, 
                           data.ipsc.reg) %>% 
      distinct(cases, .keep_all = TRUE)
  }
  return(data.reg.full)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## PCA
do_pca <- function(ann_exp_table) {
  set.seed(42)
  #Compute prcomp
  res.pca <- prcomp(ann_exp_table[, grepl('ENS', names(ann_exp_table))], 
                    center = TRUE)
  eig.vals <- res.pca$sdev^2
  loadings <- res.pca$rotation
  
  #Prepare table
  table.pca <- as.data.frame(res.pca$x)
  table.pca <- table.pca %>% 
    as_tibble()
  data.output <- cbind(ann_exp_table[, !grepl('ENS', names(ann_exp_table))], table.pca)
  
  #Specify gene set used for cluster computation
  if (!is_empty(grep("ENSG00000004897", names(ann_exp_table)))) {
    data.output <- data.output %>%
      mutate(gene_set = 'Proliferation') %>%
      dplyr::select(gene_set, everything()) %>%
      mutate(PC1 = -1 * PC1) %>% #---------Revert PC1 if prolif to minimize confusion
      as_tibble()
  }  
  if (!is_empty(grep("ENSG00000204531", names(ann_exp_table)))) {
    data.output <- data.output %>%
      mutate(gene_set = 'Stemness') %>% 
      dplyr::select(gene_set, everything()) %>% 
      as_tibble()
  } 
  if (!is_empty(grep("ENSG00000177606", names(ann_exp_table))) | !is_empty(grep("ENSG00000019549", names(ann_exp_table)))) {
    data.output <- data.output %>%
      mutate(gene_set = 'EMT_MET') %>% 
      dplyr::select(gene_set, everything()) %>% 
      mutate(PC1 = -1 * PC1) %>% #---------Revert PC1 if emt to minimize confusion
      as_tibble()
  } 
  #Enlist
  list.output <- list(data.output, eig.vals, loadings)
  
  return(list.output)
}

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
## Process counts
get_clusters <- function(pca_table, individual = FALSE, reverse = FALSE) {
  if (reverse == TRUE) {
    pca_table <- pca_table %>% 
      mutate(PC1 = -1 * PC1)
  }
  if (individual == FALSE) {
    if (pca_table[[1]][1] == "EMT_MET") {
      median.pc1 <- -5
      clusters <- pca_table %>% 
        dplyr::select(cases, project_id, sample_type, PC1) %>% 
        filter(!sample_type %in% c('Solid Tissue Normal', 'Differentiating iPSC')) %>% 
        mutate(!! pca_table[[1]][1] := ifelse(PC1 < median.pc1, 'ESC', 'Normal')) %>%
        dplyr::select(-PC1)
    } else {
      median.pc1 <- pca_table %>%  
        dplyr::select(cases, project_id, sample_type, PC1) %>% 
        filter(!sample_type %in% c('Solid Tissue Normal', 'Differentiating iPSC')) %>% 
        pull(PC1) %>% 
        median(.)
      clusters <- pca_table %>% 
        dplyr::select(cases, project_id, sample_type, PC1) %>% 
        filter(!sample_type %in% c('Solid Tissue Normal', 'Differentiating iPSC')) %>% 
        mutate(!! pca_table[[1]][1] := ifelse(PC1 > median.pc1, 'ESC', 'Normal')) %>%
        dplyr::select(-PC1)
    }
    return(clusters)
  } else {
    if (individual == TRUE) {
      medians <- pca_table %>% 
        dplyr::select(cases, project_id, sample_type, PC1) %>% 
        filter(!sample_type %in% c('Solid Tissue Normal', 'Differentiating iPSC')) %>%
        group_by(project_id) %>% 
        summarise(median.pc1 = median(PC1), .groups = "drop")
      clusters <- pca_table %>% 
        dplyr::select(cases, project_id, sample_type, PC1) %>% 
        filter(!sample_type %in% c("Solid Tissue Normal"), 
               !project_id %in% c("Differentiating iPSC")) %>% 
        mutate(Cluster = '-')
      for (i in 1:length(medians$project_id)) { 
        #It is a mess - anyone who reads this, please forgive me
        project <- medians$project_id[i]
        median.pc1 <- medians[medians$project_id == project, ]$median.pc1
        if (pca_table[[1]][1] == "EMT_MET") {
          clusters[clusters$project_id == project & clusters$PC1 <= median.pc1, ]$Cluster <- "ESC"
          clusters[clusters$project_id == project & clusters$PC1 > median.pc1, ]$Cluster <- "Normal"
        } else {
          clusters[clusters$project_id == project & clusters$PC1 >= median.pc1, ]$Cluster <- "ESC"
          clusters[clusters$project_id == project & clusters$PC1 < median.pc1, ]$Cluster <- "Normal"
        }
      }
      clusters <- clusters %>% 
        mutate(!! sym(pca_table[[1]][1]) := Cluster) %>% 
        dplyr::select(-Cluster, -PC1)
      return(clusters)
    } else {
      stop('Specify clustreing type')
    }
  } 
}


# data processing functions
# ---------------------------------------------------------------------------- #
# Medians for QC and visualization of CPE correction
prepare_medians <- function(total_exp_norm, meta, purity) {
  
  # keep only tumors
  samps_to_keep <- meta %>% 
    filter(sample_type %in% c('Additional - New Primary', 'Additional Metastatic', 'Metastatic', 'Recurrent Solid Tumor', "Primary solid Tumor"), 
           cases %in% purity$cases) %>% 
    pull(cases)
  gene_ids <- total_exp_norm$GeneID
  total_exp_norm %<>% dplyr::select(one_of(c("GeneID", samps_to_keep)))
  
  # transpose
  message("transposing.....")
  total_exp_norm <- total_exp_norm %>% 
    gather(., 
           2:ncol(.), 
           key = 'cases', 
           value = 'Counts') %>%
    spread(GeneID, Counts) %>% 
    right_join(., purity[, c("cases", "purity")], by = "cases") %>% 
    dplyr::select(cases, purity, everything())
  
  # compute medians for uncorrected expression values
  message("computing medians for original counts.....")
  medians_total_exp_norm <- sapply(total_exp_norm %>% dplyr::select(-cases, -purity), median, na.rm = TRUE) %>% 
    as.tibble() %>% 
    mutate(GeneID = colnames(total_exp_norm %>% dplyr::select(-cases, -purity)))
  
  # correct gene expression
  message("regressing out the effect of purity.....")
  for (x in 1:(length(colnames(total_exp_norm)) - 2)) {
    lm.fit <- lm(total_exp_norm[[x + 2]] ~ purity, 
                 data = total_exp_norm, 
                 na.action = "na.exclude"
    )
    total_exp_norm[[(x + 2)]] <- total_exp_norm[[(x + 2)]] - total_exp_norm[['purity']] * lm.fit$coefficients[[2]]
  }
  rm(x)
  
  # compute medians for CPE corrected expression values
  message("computing medians for corrected counts,,,,,")
  medians_total_exp_norm_cor <- sapply(total_exp_norm %>% dplyr::select(-cases, -purity), median, na.rm = TRUE) %>% 
    as.tibble() %>% 
    mutate(GeneID = colnames(total_exp_norm %>% dplyr::select(-cases, -purity)))
  
  # turn to table
  colnames(medians_total_exp_norm) <- c("Median", "GeneID")
  colnames(medians_total_exp_norm_cor) <- c("Median_Cor", "GeneID")
  medians_table <- inner_join(medians_total_exp_norm, medians_total_exp_norm_cor, by = "GeneID") %>% 
    dplyr::select(GeneID, everything()) %>% 
    gather(., 
           contains("Median"), 
           key = "Group", 
           value = "Log2CPM")
  return(medians_table)
}
# Control over random genes
pull_rgenes <- function(total_exp_ann, markers, mset, size = 1000) {
  mset <- markers %>% dplyr::select(gene_id, gene_name, contains(mset)) %>% filter_at(., vars(contains(mset)), any_vars(. > 0))
  samp <- sample(1:dim(total_exp_ann)[[1]], size)
  tab <- total_exp_ann %>% 
    dplyr::slice(samp) %>% 
    pivot_longer(., starts_with("ENS"), names_to = "gene_id", values_to = "count") %>% 
    group_by(gene_id) %>% 
    summarise(mean_count = mean(count), .groups = "drop_last") %>% 
    arrange(desc(mean_count))
  inds <- purrr::map(mset$gene_id, 
                     ~ (which(tab$gene_id == .x))) %>% 
    empt_omit_list() %>%
    purrr::map(., 
               function(x) if (x < 50) {
                 x + sample(0:50, 1)
               } else {
                 if ((x - dim(tab)[[1]]) < 50) {
                   x + sample(-50:0, 1) 
                 } else {
                   x + sample(-50:50, 1)
                 }
               }) %>% 
    purrr::reduce(., rbind) %>%
    as.numeric() %>% `+`(., 6) # +6 for correct subsetting
  tab_out <- total_exp_ann %>% dplyr::select(c(1:6, inds))
  return(tab_out)
}
# Compute intensities
get_intensity <- function(cluster_stem, cluster_prolif) {
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
  
  return(cluster_freq)
}
