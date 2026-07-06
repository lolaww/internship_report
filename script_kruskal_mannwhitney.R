#  Kruskal–Wallis + pairwise Mann–Whitney U
#  + summary table of mean/SD per population

library(tidyverse)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)

input_dir    <- "."
file_pattern <- "CRISPResso_quantification_of_editing_frequency_.*\\.(txt|tsv)$"

read_crispr_file <- function(filepath) {
  df      <- read.delim(filepath, sep = "\t", header = TRUE,
                        check.names = FALSE, stringsAsFactors = FALSE)
  mod_col <- grep("^Modified%$|^Modified_pct$|^Modified %$", names(df), value = TRUE)
  if (length(mod_col) == 0) stop("Coluna Modified% não encontrada: ", filepath)
  mod_val <- df[[mod_col[1]]]
  if (is.character(mod_val)) mod_val <- as.numeric(gsub(",", ".", mod_val))
  parts   <- strsplit(tools::file_path_sans_ext(basename(filepath)), "_")[[1]]
  data.frame(clone = parts[length(parts)-1], population = tail(parts,1),
             modified_pct = mod_val, stringsAsFactors = FALSE)
}

files <- list.files(input_dir, pattern = file_pattern, full.names = TRUE)
if (length(files) == 0) stop("Nenhum ficheiro encontrado.")

dat <- map_dfr(files, read_crispr_file) %>%
  mutate(population = factor(population, levels = c("MM","M","P","PP"))) %>%
  filter(!is.na(population))

clones <- sort(unique(dat$clone))
clone_colors <- setNames(
  c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00",
    "#A65628","#F781BF","#999999","#00CED1")[seq_len(length(clones))],
  clones
)

# TABELA: média, desvio-padrão, mediana por população
summary_tbl <- dat %>%
  group_by(population) %>%
  summarise(
    n      = sum(!is.na(modified_pct)),
    mean   = round(mean(modified_pct, na.rm = TRUE), 1),
    sd     = round(sd(modified_pct,   na.rm = TRUE), 1),
    median = round(median(modified_pct, na.rm = TRUE), 1),
    .groups = "drop"
  )
cat("\n── Mean / SD / median of Modified% per population ──\n")
print(summary_tbl, n = Inf)
write.csv(summary_tbl, file.path(input_dir, "editing_efficiency_summary_per_population.csv"),
          row.names = FALSE)
cat("Guardado: editing_efficiency_summary_per_population.csv\n")

# Kruskal–Wallis (não-paramétrico) 
kw_res <- dat %>% kruskal_test(modified_pct ~ population)
cat("\n── Kruskal–Wallis (global) ──\n"); print(kw_res)

#Mann–Whitney U 2 a 2, Bonferroni 
pairwise_res <- dat %>%
  wilcox_test(modified_pct ~ population,
              comparisons = list(c("MM","M"), c("MM","P"), c("MM","PP")),
              p.adjust.method = "bonferroni",
              paired = FALSE) %>%
  add_significance(p.col = "p.adj",
                   cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                   symbols   = c("***", "**", "*", "ns")) %>%
  add_xy_position(x = "population")

cat("\n── Mann–Whitney U (MM vs cada), Bonferroni ──\n")
print(pairwise_res[, c("group1","group2","n1","n2","statistic","p","p.adj","p.adj.signif")])

write.csv(pairwise_res[, c("group1","group2","n1","n2","statistic","p","p.adj","p.adj.signif")],
          file.path(input_dir, "editing_efficiency_mannwhitney_results.csv"), row.names = FALSE)

y_data_max <- max(dat$modified_pct, na.rm = TRUE)
pairwise_res$y.position <- c(y_data_max + 5, y_data_max + 12, y_data_max + 19)

#GRÁFICO 
p <- ggplot(dat, aes(x = population, y = modified_pct)) +
  geom_boxplot(outlier.shape = NA, fill = "grey95",
               color = "grey35", width = 0.5, linewidth = 0.45) +
  geom_point(aes(fill = clone), shape = 21, size = 3.5,
             color = "white", stroke = 0.5,
             position = position_jitter(width = 0.15, height = 0)) +
  stat_pvalue_manual(pairwise_res, label = "p.adj.signif",
                     tip.length = 0.008, bracket.size = 0.45, hide.ns = FALSE) +
  scale_fill_manual(values = clone_colors, name = "Clone") +
  scale_y_continuous(limits = c(0, y_data_max + 25),
                     breaks = seq(0, 120, by = 20),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title    = "CRISPR/Cas9 editing frequency across GFP-sorted populations",
    subtitle = paste0("n = ", length(clones),
                      " clones  ·  Kruskal–Wallis p = ",
                      signif(kw_res$p, 3),
                      "  ·  pairwise Mann–Whitney U, Bonferroni"),
    x = "GFP population", y = "Modified reads (%)"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4, color = "grey30", stroke = 0.6))) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey45", size = 9.5),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black"),
    legend.title  = element_text(face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

print(p)
ggsave(file.path(input_dir, "editing_efficiency_kruskal_mannwhitney.pdf"), p, width = 7.5, height = 5.5)
ggsave(file.path(input_dir, "editing_efficiency_kruskal_mannwhitney.png"), p, width = 7.5, height = 5.5, dpi = 300)
cat("\nGuardado: editing_efficiency_kruskal_mannwhitney.pdf / .png\n")
