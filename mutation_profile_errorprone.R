library(tidyverse)

files <- c(
  "Control" = "CRISPResso_quantification_of_editing_frequency_ctrl.txt",
  "Round 1" = "CRISPResso_quantification_of_editing_frequency_1.txt",
  "Round 2" = "CRISPResso_quantification_of_editing_frequency_2.txt",
  "Round 3" = "CRISPResso_quantification_of_editing_frequency_3.txt"
)
input_dir <- "."

# LER E EXTRAIR AS CLASSES DE MUTAÇÃO 
read_one <- function(label, fname) {
  df <- read.delim(file.path(input_dir, fname), sep = "\t",
                   header = TRUE, check.names = FALSE)
  al <- as.numeric(df[["Reads_aligned"]][1])
  data.frame(
    Sample        = label,
    Substitutions = as.numeric(df[["Only Substitutions"]][1]) / al * 100,
    Deletions     = as.numeric(df[["Only Deletions"]][1])     / al * 100,
    Insertions    = as.numeric(df[["Only Insertions"]][1])    / al * 100,
    stringsAsFactors = FALSE
  )
}

dat <- imap_dfr(files, ~ read_one(.y, .x)) %>%
  pivot_longer(c(Substitutions, Deletions, Insertions),
               names_to = "Class", values_to = "Pct") %>%
  mutate(
    Sample = factor(Sample, levels = names(files)),
    Class  = factor(Class, levels = c("Substitutions", "Deletions", "Insertions"))
  )

# GRÁFICO
class_colors <- c(
  "Substitutions" = "#C0392B",   # vermelho
  "Deletions"     = "#2E75B6",   # azul
  "Insertions"    = "#E9A23B"    # amarelo
)

p <- ggplot(dat, aes(x = Sample, y = Pct, fill = Class)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(aes(label = sprintf("%.1f", Pct)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3, color = "grey25") +
  scale_fill_manual(values = class_colors, name = "Mutation class\n(reads with only this type)") +
  scale_y_continuous(limits = c(0, 100),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Mutation profile of error-prone PCR (chr8:41511631-41511831)",
    subtitle = "Control vs increasing rounds of error-prone amplification",
    x = NULL, y = "Reads (% of aligned)"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey45", size = 10),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black"),
    legend.title  = element_text(face = "bold", size = 9),
    legend.position = "right",
    panel.grid    = element_blank()
  )

print(p)
ggsave(file.path(input_dir, "errorprone_mutation_profile.pdf"), p, width = 8, height = 5)
ggsave(file.path(input_dir, "errorprone_mutation_profile.png"), p, width = 8, height = 5, dpi = 300)
cat("Guardado: errorprone_mutation_profile.pdf / .png\n")
