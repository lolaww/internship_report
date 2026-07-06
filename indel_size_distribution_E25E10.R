library(tidyverse)

# Carregar ficheiros
mm <- read_tsv("Indel_histogram_MM.txt") %>% mutate(population = "MM")
pp <- read_tsv("Indel_histogram_PP.txt") %>% mutate(population = "PP")

# Normalizar para % (dividir pelo total de reads de cada populacao)
mm <- mm %>% mutate(pct = fq / sum(fq) * 100)
pp <- pp %>% mutate(pct = fq / sum(fq) * 100)

df <- bind_rows(mm, pp)

# Filtrar range de interesse (ex: -30 a +15, exclui o pico de 0) -> remover o 0 para ver melhor os indels 
df_plot <- df %>% filter(indel_size != 0, indel_size >= -30, indel_size <= 15)

# Grafico MM vs PP sobrepostos
ggplot(df_plot, aes(x = indel_size, y = pct, fill = population)) +
  geom_col(position = "dodge", width = 0.8) +
  scale_fill_manual(
    values = c("MM" = "#A32D2D", "PP" = "#185FA5"),
    labels = c("MM (perda de atividade)", "PP (ganho de atividade)")
  ) +
  scale_x_continuous(breaks = seq(-30, 15, by = 5)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1, suffix = "%")) +
  labs(
    title    = "Distribuicao de indels por tamanho",
    subtitle = "Clone E25E10 - enhancer chr8:41511631-41511831",
    x        = "Tamanho do indel (pb)  |  negativo = delecao  |  positivo = insercao",
    y        = "Frequencia de reads (%)",
    fill     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    plot.title      = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("indel_histogram_MM_PP.pdf", width = 10, height = 5)

#Versao com o pico 0 incluido (ver proporcao real)
df_full <- df %>% filter(indel_size >= -30, indel_size <= 15)

ggplot(df_full, aes(x = indel_size, y = pct, fill = population)) +
  geom_col(position = "dodge", width = 0.8) +
  scale_fill_manual(
    values = c("MM" = "#A32D2D", "PP" = "#185FA5"),
    labels = c("MM (perda de atividade)", "PP (ganho de atividade)")
  ) +
  scale_x_continuous(breaks = seq(-30, 15, by = 5)) +
  scale_y_continuous(labels = scales::percent_format(scale = 1, suffix = "%")) +
  labs(
    title    = "Indel size distribution",
    subtitle = "Clone E25E10 - enhancer chr8:41511631-41511831",
    x        = "indel size (bp)  |  0 = no indel  |  negative = deletion  |  positive = insertion",
    y        = "read frequency (%)",
    fill     = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("indel_histogram_MM_PP_com_zero.pdf", width = 10, height = 5)

# Tabela: indels mais frequentes em cada populacao
top_mm <- mm %>%
  filter(indel_size != 0) %>%
  arrange(desc(pct)) %>%
  slice_head(n = 10) %>%
  select(indel_size, fq, pct) %>%
  mutate(population = "MM")

top_pp <- pp %>%
  filter(indel_size != 0) %>%
  arrange(desc(pct)) %>%
  slice_head(n = 10) %>%
  select(indel_size, fq, pct) %>%
  mutate(population = "PP")

tabela <- bind_rows(top_mm, top_pp) %>%
  arrange(population, desc(pct))

print(tabela)
write_csv(tabela, "top_indels_MM_PP.csv")
