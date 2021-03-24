### INFO: Script to quantify gene expressiom
### DATE: 15.10.2019
### AUTHOR: Artem Baranovskii

# arguments
# ---------------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
cut_bot <- as.numeric(args[1])
cut_up <- as.numeric(args[2])
out_path <- as.character(args[3])

# libraries
# ---------------------------------------------------------------------------- #
require(tidyverse, quietly = TRUE, warn.conflicts = FALSE)
require(magrittr, quietly = TRUE, warn.conflicts = FALSE)
require(edgeR, quietly = TRUE, warn.conflicts = FALSE)

# paths
# ---------------------------------------------------------------------------- #
f_path <- file.path("~/Scripts/Scripts/Draft")
counts_path <- file.path("/local", "artem", "Projects", "Draft", "Data", "Counts")
meta_path <- file.path("/local", "artem", "Projects", "Draft", "Data", "Meta")
gdc_path <- file.path(counts_path, "GDCprocessed")

# functions
# ---------------------------------------------------------------------------- #
source(file.path(f_path, "F_stemtcga.R"))

# variables
# ---------------------------------------------------------------------------- #
if (is.na(cut_bot)) {
  cut_bot <- -2
}

if (is.na(cut_up)) {
  cut_up <- 9
}

if (is.na(out_path)) {
  out_path <- counts_path
}

# speak
message(paste0("selected cutoffs are ", cut_bot, " and ", cut_up, " log2(CPM)"))
message(paste0("outpath defined as ", out_path))


# load GSE98411 counts
# ---------------------------------------------------------------------------- #
message("loading GSE98411 counts.....")
ipsc_raw_counts <- read_tsv(file = file.path(counts_path, "GSE98411_featureCounts.tsv"), 
                            col_names = TRUE, 
                            col_types = cols(.default = col_double(), 
                                             Name = col_character())) %>% 
  dplyr::select(-SRR5491216) # drop outlier SRR5491216 (from original publication)

# load GDC
# ---------------------------------------------------------------------------- #
message("loading GDC counts.....")
gdc_raw_counts <- purrr::map(list.files(path = gdc_path, full.names = TRUE), 
                             function(x) read_tsv(x, 
                                                  col_names = TRUE, 
                                                  col_types = cols(.default = col_double(), 
                                                                   GeneID = col_character()))) %>% 
  purrr::reduce(., full_join, by = "GeneID")

# total
# ---------------------------------------------------------------------------- #
message("joining tables.....")
total_raw_counts <- inner_join(gdc_raw_counts, ipsc_raw_counts, by = c("GeneID" = "Name")) %>% distinct(GeneID, .keep_all = TRUE)

# prepare matrix for normalisation
counts_matrix <- total_raw_counts %>% dplyr::select(-GeneID) %>% as.data.frame()
rownames(counts_matrix) <- total_raw_counts %>% pull(GeneID)

# quality
message("quality check.....")
cpm.Log <- cpm(total_raw_counts %>% dplyr::select(-GeneID), 
               log = TRUE)
median.log2.cpm <- apply(cpm.Log, 1, median) %>%
  as_tibble() %>% 
  mutate(GeneID = total_raw_counts$GeneID)

# plot hist
cpm_hist <- ggplot(median.log2.cpm, aes(value)) + 
  geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = c(cut_bot, cut_up), colour = c("blue", "red"), size = 1) +
  theme_minimal() + 
  xlab("Median log2(CPM)")
ggsave(file.path(out_path, "Median_log2_CPM_hist.pdf"), plot = cpm_hist, device = "pdf", dpi = "print")

# filter
counts_matrix <- counts_matrix[median.log2.cpm > cut_bot & median.log2.cpm < cut_up, ]

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

# sample sheet
sample_sheet <- tibble(sample = rownames(y$samples), 
                       lib_size = y$samples$lib.size, 
                       norm_factor = y$samples$norm.factors)


# subsetted tables
message("subsetting tables.....")
meta <- purrr::map(c("Meta.tsv", "Meta_iPSC.tsv"), 
                   function(x) read_tsv(file = file.path(meta_path, x), 
                                        col_names = TRUE, 
                                        col_types = cols(.default = col_character()))) %>% 
  purrr::reduce(., base::rbind)

markers <- read_tsv(file = file.path(meta_path, "Markers.tsv"))

exp_tables <- purrr::map(c("Stem", "Prolif", "Emt"), 
                         function(x) subset_exp(counts_matrix_norm, 
                                                markers = markers, 
                                                meta = meta, 
                                                gene_set = x))
## paths to write subsets
t_paths <- c("stem_counts_table.tsv", 
             "prolif_counts_table.tsv",
             "emt_counts_table.tsv")

# transposing
message("transposing.....")
chunks <- list()
for (i in 1:round(dim(counts_matrix_norm)[[1]] / 2000, 0)) {
  # cut & transpose
  if (i != round(dim(counts_matrix_norm)[[1]] / 2000, 0)) {
    chunks[[i]] <- counts_matrix_norm %>% 
      dplyr::slice(., ((2000 * (i - 1)) + 1):(i * 2000)) %>% 
      pivot_longer(., 
                   2:ncol(.), 
                   names_to = 'cases', 
                   values_to = 'Counts') %>%
      pivot_wider(names_from = GeneID, values_from = Counts) 
    message(paste0(i * 2000, " rows processed"))
  } else {
    chunks[[i]] <- counts_matrix_norm %>% 
      dplyr::slice(., ((2000 * (i - 1)) + 1):dim(counts_matrix_norm)[[1]]) %>% 
      gather(., 
             2:ncol(.), 
             key = 'cases', 
             value = 'Counts') %>%
      spread(GeneID, Counts) 
    message(paste0("total of ", dim(counts_matrix_norm)[[1]], " rows processed"))
  }
}

# annotate
message("annotating.....")
names_meta <- rlang::syms(names(meta))
counts_matrix_norm_t <- chunks %>% 
  purrr::reduce(., inner_join, by = "cases") %>% 
  inner_join(., meta, 
             by = c('cases')) %>% 
  dplyr::select(cases, !! names_meta[[2]], !! names_meta[[3]], !! names_meta[[4]], !! names_meta[[5]], !! names_meta[[8]], contains('ENS'))


# write results
message("writing down results.....")
write_tsv(path = file.path(out_path, "total_counts_medians.tsv"), median.log2.cpm, col_names = TRUE)
write_tsv(path = file.path(out_path, "total_norm_counts.tsv"), counts_matrix_norm, col_names = TRUE)
write_tsv(path = file.path(out_path, "total_norm_counts_t.tsv"), counts_matrix_norm_t, col_names = TRUE)
write_tsv(path = file.path(out_path, "total_sample_sheet.tsv"), sample_sheet, col_names = TRUE)
exp_tables <- purrr::map2(.x = exp_tables, 
                          .y = t_paths, 
                          function(x, y) write_tsv(x, path = file.path(counts_path, y), col_names = TRUE))
message("done!")