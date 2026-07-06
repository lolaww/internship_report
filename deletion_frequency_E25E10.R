library(tidyverse)

clone_name   <- "E25E10"
enhancer     <- "chr8:41511631-41511831"
input_dir    <- "."          # pasta com os ficheiros Effect_vector_deletion_*.txt
populations  <- c("MM", "M", "P", "PP")

pop_labels <- c(
  "MM" = "MM (GFP very low)",
  "M"  = "M (GFP low)",
  "P"  = "P (GFP high)",
  "PP" = "PP (GFP very high)"
)
pop_colours <- c(
  "MM" = "#A32D2D", 
  "M"  = "#E07B54",  
  "P"  = "#4CAF50",  
  "PP" = "#185FA5"   
)


#LEITURA
read_effect <- function(pop, dir) {
  fname <- file.path(dir, paste0("Effect_vector_deletion_", pop, ".txt"))
  if (!file.exists(fname)) { warning("Não encontrado: ", fname); return(NULL) }
  read_tsv(fname, comment = "#", col_names = c("position", "effect"),
           col_types = cols(position = col_integer(), effect = col_double())) %>%
    mutate(population = pop)
}

df <- map(populations, read_effect, dir = input_dir) %>%
  bind_rows() %>%
  mutate(
    population = factor(population, levels = populations),
    pop_label  = recode(population, !!!pop_labels)
  )


# GRÁFICO — 4 populações sobrepostas

p_line <- ggplot(df, aes(x = position, y = effect,
                         colour = population, fill = population)) +
  geom_area(alpha = 0.12, position = "identity") +
  geom_line(linewidth = 0.75) +
  scale_colour_manual(values = pop_colours, labels = pop_labels) +
  scale_fill_manual(  values = pop_colours, labels = pop_labels) +
  scale_x_continuous(breaks = seq(0, 200, 25), expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title    = paste0("Deletion frequency per amplicon position — ", clone_name),
    subtitle = paste0("Enhancer ", enhancer,
                      "  ·  Values = % of total aligned reads with a deletion at each position"),
    x        = "Position in amplicon (bp)",
    y        = "Reads with deletion (%)",
    colour   = NULL, fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 9, colour = "grey45", margin = margin(b = 10)),
    legend.position = "bottom",
    legend.text     = element_text(size = 10),
    axis.line       = element_line(colour = "grey65", linewidth = 0.4),
    axis.ticks      = element_line(colour = "grey65", linewidth = 0.4),
    panel.grid      = element_blank()
  )

ggsave(
  filename = paste0(clone_name, "_deletion_frequency_by_position.png"),
  plot     = p_line,
  width    = 9, height = 5, dpi = 300
)


# TABELA — top 10 hotspots 
hotspots <- df %>%
  group_by(population) %>%
  slice_max(order_by = effect, n = 10) %>%
  arrange(population, desc(effect)) %>%
  mutate(
    rank       = row_number(),
    effect_pct = sprintf("%.2f%%", effect)
  ) %>%
  select(population, rank, position, effect_pct)


hotspots_wide <- hotspots %>%
  pivot_wider(
    names_from  = population,
    values_from = c(position, effect_pct),
    names_glue  = "{population}_{.value}"
  ) %>%
  select(rank,
         MM_position, MM_effect_pct,
         M_position,  M_effect_pct,
         P_position,  P_effect_pct,
         PP_position, PP_effect_pct)

print(hotspots_wide)
write_csv(hotspots_wide, paste0(clone_name, "_deletion_hotspots_top10.csv"))


#TABELA COMO FIGURA 
library(gridExtra)
library(grid)

colnames(hotspots_wide) <- c(
  "Rank",
  "MM pos", "MM %",
  "M pos",  "M %",
  "P pos",  "P %",
  "PP pos", "PP %"
)

tg <- tableGrob(
  hotspots_wide,
  rows = NULL,
  theme = ttheme_minimal(
    colhead = list(
      bg_params = list(fill = c("grey90", "#A32D2D","#A32D2D",
                                "#E07B54","#E07B54",
                                "#4CAF50","#4CAF50",
                                "#185FA5","#185FA5")),
      fg_params = list(col = c("grey20", rep("white", 8)),
                       fontface = "bold", fontsize = 9)
    ),
    core = list(
      bg_params = list(fill = rep(c("#FFF5F5","#FFFFFF"), length.out = 10)),
      fg_params = list(fontsize = 9)
    )
  )
)

ggsave(
  filename = paste0(clone_name, "_deletion_hotspots_table.png"),
  plot     = tg,
  width    = 10, height = 3.5, dpi = 300
)

cat("\nDone! Ficheiros gerados:\n")
cat(" -", paste0(clone_name, "_deletion_frequency_by_position.png\n"))
cat(" -", paste0(clone_name, "_deletion_hotspots_top10.csv\n"))
cat(" -", paste0(clone_name, "_deletion_hotspots_table.png\n"))
