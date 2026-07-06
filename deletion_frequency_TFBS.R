library(ggplot2)
library(dplyr)

data_dir <- "."   # pasta com Effect_vector_deletion_MM.txt e _PP.txt

read_effect <- function(fname, pop) {
  df <- read.delim(file.path(data_dir, fname), header = TRUE,
                   sep = "\t", comment.char = "", check.names = FALSE)
  names(df) <- c("Position", "Pct")
  df$Pct <- as.numeric(gsub(",", ".", df$Pct))
  df$Population <- pop
  df
}

del <- bind_rows(
  read_effect("Effect_vector_deletion_MM.txt", "MM (loss of activity)"),
  read_effect("Effect_vector_deletion_PP.txt", "PP (gain of activity)")
)

ymax <- max(del$Pct, na.rm = TRUE)

# REGIÕES DE TFBS (JASPAR) 
# MM hotspot: SP1 (67-77), KLF4 (67-78), GATA2 (80-84) -> banda 67-84
# PP hotspot: NKX2-5 (122-132) -> banda 122-132
tfbs <- data.frame(
  xmin  = c(67, 122),
  xmax  = c(84, 132),
  label = c("SP1 · KLF4 · GATA2", "NKX2-5"),
  fill  = c("#C0392B", "#2E75B6"),  
  lx    = c(75.5, 127),             
  stringsAsFactors = FALSE
)


pop_colors <- c("MM (loss of activity)" = "#C0392B",
                "PP (gain of activity)" = "#2E75B6")

#GRÁFICO 
p <- ggplot() +
  
  geom_rect(data = tfbs,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax * 1.08),
            fill = tfbs$fill, alpha = 0.12) +
  
  geom_text(data = tfbs,
            aes(x = lx, y = ymax * 1.05, label = label),
            angle = 0, size = 3, fontface = "bold",
            color = tfbs$fill, hjust = 0.5) +
 
  geom_line(data = del, aes(x = Position, y = Pct, color = Population),
            linewidth = 0.8) +
  scale_color_manual(values = pop_colors, name = NULL) +
  scale_x_continuous(limits = c(1, 200), breaks = seq(0, 200, 25),
                     expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, ymax * 1.12),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title    = "Deletion frequency and predicted TF binding sites",
    subtitle = "(clone E25E10) · shaded bands = predicted TFBS (JASPAR)",
    x = "Position in enhancer (bp)",
    y = "Reads with deletion (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = "grey45", size = 9.5),
    axis.title    = element_text(face = "bold"),
    axis.text     = element_text(color = "black"),
    legend.position = c(0.98, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3),
    panel.grid    = element_blank()
  )

print(p)
ggsave(file.path(data_dir, "deletion_TFBS_overlay.pdf"), p, width = 9, height = 5.5)
ggsave(file.path(data_dir, "deletion_TFBS_overlay.png"), p, width = 9, height = 5.5, dpi = 300)
cat("Guardado: deletion_TFBS_overlay.pdf / .png\n")
cat("Pico MM na posição", del$Position[which.max(ifelse(del$Population=="MM (loss of activity)", del$Pct, -1))], "\n")
