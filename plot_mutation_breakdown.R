library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# mudar estas duas linhas a depender do enhancer/clone
clone_name <- "E25E10"
input_dir  <- "."

#
populations <- c("MM", "M", "P", "PP")

read_crispresso <- function(clone, pop, dir) {
  fname <- file.path(dir,
    paste0("CRISPResso_quantification_of_editing_frequency_", clone, "_", pop, ".txt"))

  if (!file.exists(fname)) {
    warning("Ficheiro não encontrado: ", fname)
    return(NULL)
  }

  df  <- read_tsv(fname, show_col_types = FALSE)
  row <- df[df$Amplicon == "Reference", ]

  data.frame(
    population = pop,
    aligned    = row$Reads_aligned,
    unmod      = row$Unmodified,
    only_ins   = row$`Only Insertions`,
    only_del   = row$`Only Deletions`,
    only_sub   = row$`Only Substitutions`,
    del_containing = row$`Insertions and Deletions` +
                     row$`Deletions and Substitutions` +
                     row$`Insertions Deletions and Substitutions`,
    ins_sub    = row$`Insertions and Substitutions`
  )
}

data <- lapply(populations, read_crispresso,
               clone = clone_name, dir = input_dir) |>
  bind_rows()


# Converter para % de reads alinhadas
data <- data %>%
  mutate(across(c(unmod, only_ins, only_del, only_sub, del_containing, ins_sub),
                ~ . / aligned * 100))

long <- data %>%
  select(population, unmod, only_del, del_containing, only_sub, only_ins, ins_sub) %>%
  pivot_longer(-population, names_to = "category", values_to = "pct")

long$population <- factor(long$population, levels = c("MM", "M", "P", "PP"))

long$category <- factor(long$category,
  levels = c("unmod", "only_del", "del_containing", "only_sub", "only_ins", "ins_sub"),
  labels = c(
    "Unmodified",
    "Deletion only",
    "Deletion-containing (+ ins and/or sub)",
    "Substitution only",
    "Insertion only",
    "Insertion + substitution"
  )
)

colours <- c(
  "Unmodified"                          = "#B4B2A9",
  "Deletion only"                       = "#2171B5",
  "Deletion-containing (+ ins and/or sub)" = "#6BAED6",
  "Substitution only"                   = "#7F77DD",
  "Insertion only"                      = "#5DCAA5",
  "Insertion + substitution"            = "#EF9F27"
)


# GRÁFICO
p <- ggplot(long, aes(x = population, y = pct, fill = category)) +
  geom_col(width = 0.45) +
  scale_fill_manual(
    values = colours,
    guide = guide_legend(
      ncol      = 1,        
      keywidth  = unit(0.5, "cm"),
      keyheight = unit(0.5, "cm")   
    )
  ) +
  scale_y_continuous(
    labels = scales::label_percent(scale = 1),
    expand = c(0, 0),
    limits = c(0, 101)
  ) +
  labs(
    title    = paste0(clone_name, " \u2014 Mutation breakdown by FACS population"),
    subtitle = "Each bar represents 100% of aligned reads, decomposed by mutation class",
    x        = "FACS population",
    y        = "% of aligned reads",
    fill     = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title       = element_text(size = 10, face = "bold",
                                    margin = margin(b = 3)),
    plot.subtitle    = element_text(size = 9, colour = "grey45",
                                    margin = margin(b = 12)),
    plot.margin      = margin(t = 14, r = 14, b = 10, l = 10),
    axis.line        = element_line(colour = "grey65", linewidth = 0.4),
    axis.ticks       = element_line(colour = "grey65", linewidth = 0.4),
    axis.text        = element_text(colour = "grey30", size = 11),
    axis.title.x     = element_text(colour = "grey20", size = 11,
                                    margin = margin(t = 8)),
    axis.title.y     = element_text(colour = "grey20", size = 11,
                                    margin = margin(r = 8)),
    legend.position  = "right",
    legend.text      = element_text(size = 8, colour = "grey25", lineheight = 1.3),
    legend.key.size  = unit(0.5, "cm"),     # ← todos iguais
    legend.spacing.y = unit(0.3, "cm"),
    panel.grid       = element_blank()
  )

ggsave(
  filename = paste0(clone_name, "_mutation_breakdown.png"),
  plot     = p,
  width    = 6.5,    # ← mais largo para acomodar legenda vertical
  height   = 3,
  dpi      = 300
)


library(gridExtra)
library(grid)

# Tabela  
tab <- long %>%
  mutate(pct_fmt = sprintf("%.2f%%", pct)) %>%
  select(population, category, pct_fmt) %>%
  pivot_wider(names_from = category, values_from = pct_fmt) %>%
  mutate(
    `SUM check` = long %>%
      group_by(population) %>%
      summarise(s = sprintf("%.2f%%", sum(pct))) %>%
      arrange(factor(population, levels = c("MM","M","P","PP"))) %>%
      pull(s)
  ) %>%
  rename(Population = population)

tg <- tableGrob(
  tab,
  rows = NULL,
  theme = ttheme_minimal(
    colhead = list(
      bg_params = list(fill = "#2171B5"),
      fg_params = list(col = "white", fontface = "bold", fontsize = 9)
    ),
    core = list(
      bg_params = list(fill = rep(c("#EEF4FB","#FFFFFF"), length.out = nrow(tab))),
      fg_params = list(fontsize = 9)
    )
  )
)

ggsave(
  filename = paste0(clone_name, "_percentage_table.png"),
  plot     = tg,
  width    = 11,
  height   = 1.8,
  dpi      = 300
)
