library(tidyverse)

# Carregar os ficheiros - ajustar os caminhos conforme a pasta
del_MM <- read_tsv("Effect_vector_deletion_MM.txt", comment = "#",
                   col_names = c("position", "effect_MM"))

del_PP <- read_tsv("Effect_vector_deletion_PP.txt", comment = "#",
                   col_names = c("position", "effect_PP"))


df <- left_join(del_MM, del_PP, by = "position")

# Calcular diferença MM - PP
df <- df %>%
  mutate(diff_MM_PP = effect_MM - effect_PP)

# Gráfico: deleções MM x PP por posição
df_long <- df %>%
  pivot_longer(cols = c(effect_MM, effect_PP),
               names_to = "population",
               values_to = "effect") %>%
  mutate(population = recode(population,
                             "effect_MM" = "MM (perda de atividade)",
                             "effect_PP" = "PP (ganho de atividade)"))

ggplot(df_long, aes(x = position, y = effect, color = population, fill = population)) +
  geom_line(linewidth = 0.8) +
  geom_area(alpha = 0.15, position = "identity") +
  scale_color_manual(values = c("MM (perda de atividade)" = "#A32D2D",
                                "PP (ganho de atividade)"  = "#185FA5")) +
  scale_fill_manual(values  = c("MM (perda de atividade)" = "#A32D2D",
                                "PP (ganho de atividade)"  = "#185FA5")) +
  labs(
    title    = "Frequência de deleções por posição no amplicon",
    subtitle = "Clone E25E10 · enhancer chr8:41511631–41511831",
    x        = "Posição no amplicon (pb)",
    y        = "Frequência de reads com deleção (%)",
    color    = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom",
        plot.title      = element_text(face = "bold"))

ggsave("effect_vector_del_MM_PP.pdf", width = 10, height = 5)

#Tabela de hotspots

hotspots <- df %>%
  mutate(
    max_effect   = pmax(effect_MM, effect_PP),
    dominant_pop = if_else(effect_MM >= effect_PP,
                           "MM (perda de atividade)",
                           "PP (ganho de atividade)")
  ) %>%
  arrange(desc(max_effect)) %>%
  slice_head(n = 20) %>%         
  select(position, effect_MM, effect_PP, diff_MM_PP, dominant_pop)

print(hotspots)

write_csv(hotspots, "hotspots_delecao_MM_PP.csv")

# (Opcional) Gráfico da diferença MM - PP
ggplot(df, aes(x = position, y = diff_MM_PP)) +
  geom_line(color = "#7F77DD", linewidth = 0.8) +
  geom_area(fill = "#7F77DD", alpha = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  annotate("text", x = 5, y = max(df$diff_MM_PP) * 0.9,
           label = "MM tem mais deleções", hjust = 0, size = 3.5, color = "#A32D2D") +
  annotate("text", x = 5, y = min(df$diff_MM_PP) * 0.9,
           label = "PP tem mais deleções", hjust = 0, size = 3.5, color = "#185FA5") +
  labs(
    title    = "Diferença de deleções (MM − PP) por posição",
    subtitle = "Clone E25E10 · enhancer chr8:41511631–41511831",
    x        = "Posição no amplicon (pb)",
    y        = "Diferença MM − PP (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("effect_vector_diff_MM_PP.pdf", width = 10, height = 5)