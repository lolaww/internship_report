#  Enhancer chr8 (1) — deletion hotspots (% reads, MM+M+P+PP pooled)
#  Ficheiros: Modification_count_vectors_<CLONE>_<POP>.txt (POP = MM, M, P, PP) e sgRNA_positions_chr8.xlsx

library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)

data_dir <- "."     # <- pasta com os Modification_count_vectors_*.txt

clones <- c("E22P1A3", "E25B2", "E25B3", "E25E10")
pops   <- c("MM", "M", "P", "PP")

sgRNA_file  <- "sgRNA_positions_chr8.xlsx"
sgRNA_sheet <- 1

y_mode          <- "percent"
enhancer_length <- 200

read_deletions_counts <- function(clone, pop) {
  fpath <- file.path(data_dir,
                     paste0("Modification_count_vectors_", clone, "_", pop, ".txt"))
  if (!file.exists(fpath)) stop("Ficheiro não encontrado: ", fpath)

  df <- read.delim(fpath, header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
  
  labels <- df[[1]]

  del_row   <- which(labels == "Deletions")
  total_row <- which(labels == "Total")
  if (length(del_row) == 0) stop("Sem linha 'Deletions' em: ", fpath)

  del_counts <- as.numeric(df[del_row, 2:ncol(df)])
  total      <- if (length(total_row) > 0)
                  as.numeric(df[total_row, 2]) else NA_real_

  list(counts = del_counts, total = total)
}

#PARA CADA CLONE: pool das 4 populações, % independente
build_clone <- function(clone) {
  per_pop <- lapply(pops, function(p) read_deletions_counts(clone, p))

  mat           <- do.call(rbind, lapply(per_pop, function(x) x$counts))
  summed_counts <- colSums(mat)                                  
  total_reads   <- sum(sapply(per_pop, function(x) x$total))     

  value <- if (y_mode == "percent") summed_counts / total_reads * 100 else summed_counts

  data.frame(
    Position = seq_along(value),
    Value    = value,
    Clone    = clone,
    stringsAsFactors = FALSE
  )
}

del <- bind_rows(lapply(clones, build_clone))
del$Clone <- factor(del$Clone, levels = clones)

cat("Clones:", paste(clones, collapse = ", "), "\n")
cat("Eixo Y:", y_mode, "\n")
cat("Pico:", format(round(max(del$Value), 2), big.mark = ","),
    "na posição", del$Position[which.max(del$Value)], "\n")

# POSIÇÕES DOS sgRNAs
sg_pos <- read_excel(sgRNA_file, sheet = sgRNA_sheet) %>%
  pull(`position sgRNA`) %>%
  na.omit() %>% as.integer() %>% unique() %>% sort()

# GRÁFICO 
cores_clones <- c(
  "E22P1A3" = "#2E75B6",
  "E25B2"   = "#548235",
  "E25B3"   = "#C55A11",
  "E25E10"  = "#7030A0"
)

y_label <- if (y_mode == "percent")
  "% reads with deletion (MM+M+P+PP pooled)" else
  "Reads with deletion (sum of MM+M+P+PP)"

grafico_final <- ggplot(del, aes(x = Position, y = Value, color = Clone)) +
  # cut-sites como linhas verticais cinza, por trás das curvas
  geom_vline(xintercept = sg_pos, color = "grey80",
             linewidth = 0.3, linetype = "solid") +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = cores_clones) +
  scale_x_continuous(limits = c(1, enhancer_length),
                     breaks = seq(0, enhancer_length, by = 25),
                     expand = c(0.005, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Position in enhancer (bp)", y = y_label, color = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    legend.key.width = unit(1.0, "cm"),
    axis.text        = element_text(color = "black"),
    plot.margin      = margin(t = 6, r = 10, b = 6, l = 10)
  )


out <- paste0("chr8_deletions_countss_", y_mode)
ggsave(paste0(out, ".pdf"), grafico_final, width = 22, height = 14, units = "cm")
ggsave(paste0(out, ".png"), grafico_final, width = 22, height = 14, units = "cm", dpi = 300)

message("✓ Guardado: ", out, ".pdf / .png")
