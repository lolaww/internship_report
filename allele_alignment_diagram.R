library(tidyverse)
library(magick)

clone_id  <- "E25E10"
enhancer  <- "chr8:41511631\u201341511831"
file_MM   <- "Alleles_frequency_table_MM_E25E10.txt"
file_PP   <- "Alleles_frequency_table_PP_E25E10.txt"

# posições
targets <- list(
  list(pop = "MM", size = 22),
  list(pop = "MM", size = 14),
  list(pop = "PP", size = 18)
)

context_nt <- 18                      
out_base   <- paste0("allele_diagrams_", clone_id)
png_scale  <- 3                       

# SVG
CH        <- 8.6    
FS        <- 13     
MARGIN_L  <- 78     
PANEL_H   <- 92     
TOP       <- 62     
ROW_REF_DY<- 30     
ROW_SEQ_DY<- 52     
POS_DY    <- 12     
HEAD_DY   <- -8     

style <- list(
  MM = list(title = "#C0392B", fill = "#FADBD8", border = "#C0392B", dash = "#C0392B"),
  PP = list(title = "#2471A3", fill = "#D6EAF8", border = "#2471A3", dash = "#2471A3")
)

# Tabela
tbls <- list(
  MM = read_tsv(file_MM, show_col_types = FALSE),
  PP = read_tsv(file_PP, show_col_types = FALSE)
)

# Filtra deleção PURA (nada inserido/substituído)
get_top_allele <- function(tbl, size_bp) {
  tbl %>%
    filter(n_inserted == 0,       
           n_mutated  == 0,       
           n_deleted  == size_bp,  
           Read_Status == "MODIFIED") %>%
    arrange(desc(`#Reads`)) %>%    
    slice_head(n = 1)              
}

# Parse + numeraçāo
parse_allele <- function(pop, size_bp) {
  row <- get_top_allele(tbls[[pop]], size_bp)
  if (nrow(row) == 0) stop("Alelo nao encontrado: ", pop, " del-", size_bp)

  ref <- row$Reference_Sequence
  aln <- row$Aligned_Sequence

  # definir posicao 1
  ref_chars <- strsplit(ref, "")[[1]]
  nz <- which(ref_chars != "-")
  rs <- min(nz); re_ <- max(nz)     
  ref_core <- substr(ref, rs, re_)   
  aln_core <- substr(aln, rs, re_)  

  # Encontrar posiçāo da deleçāo
  del_pos <- which(strsplit(aln_core, "")[[1]] == "-")
  ds <- min(del_pos)   # posição (no enhancer, sem primers) da 1ª base deletada 
  de <- max(del_pos)   

  cs <- max(1, ds - context_nt)
  ce <- min(nchar(ref_core), de + context_nt)

  list(
    pop = pop, size = size_bp,
    reads = row$`#Reads`, pct = row$`%Reads`,
    ref_win = substr(ref_core, cs, ce),   
    aln_win = substr(aln_core, cs, ce),   
    del_s = ds - cs + 1, del_e = de - cs + 1,   
    amp_s = ds, amp_e = de                       
  )
}

# SVG
esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}

nt_x <- function(i) MARGIN_L + (i - 0.5) * CH

# Painel
panel_svg <- function(d, y0) {
  st <- style[[d$pop]]
  p <- character(0)

  title <- sprintf("%s population \u2014 deletion \u2212%d bp", d$pop, d$size)
  stats <- sprintf("(%s reads = %.3f%% of all reads)",
                   formatC(d$reads, format = "d", big.mark = ","), d$pct)
  p <- c(p, sprintf('<text x="%d" y="%d" font-size="12.5" font-weight="700" fill="%s" font-family="Arial">%s</text>',
                    MARGIN_L, y0 + HEAD_DY, st$title, esc(title)))
  p <- c(p, sprintf('<text x="%d" y="%d" font-size="11.5" fill="#333" font-family="Arial">%s</text>',
                    MARGIN_L + 232, y0 + HEAD_DY, esc(stats)))

  box_x <- nt_x(d$del_s) - CH/2
  box_w <- (d$del_e - d$del_s + 1) * CH
  box_y <- y0 + ROW_REF_DY - FS + 1
  p <- c(p, sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="3" fill="%s" stroke="%s" stroke-width="1.4"/>',
                    box_x, box_y, box_w, FS + 5, st$fill, st$border))

  for (k in 1:2) {
    pos_i <- if (k == 1) d$del_s else d$del_e   # ONDE (posição na janela -> x)
    amp   <- if (k == 1) d$amp_s else d$amp_e   # O QUE (nº no enhancer: 67 / 88)
    p <- c(p, sprintf('<text x="%.1f" y="%d" font-size="10" font-style="italic" fill="#7F8C8D" text-anchor="middle" font-family="Arial">%d</text>',
                      nt_x(pos_i), y0 + POS_DY, amp))
  }

  # "5'-Ref:" e "5'-Seq:" a esquerda
  for (k in 1:2) {
    dy  <- if (k == 1) ROW_REF_DY else ROW_SEQ_DY
    lab <- if (k == 1) "Ref:" else "Seq:"
    p <- c(p, sprintf('<text x="%d" y="%d" font-size="12" fill="#333" text-anchor="end" font-family="Courier New" font-weight="700">5\u2032&#8211;%s</text>',
                      MARGIN_L - 8, y0 + dy, lab))
  }

  ref_c <- strsplit(d$ref_win, "")[[1]]
  for (i in seq_along(ref_c)) {
    p <- c(p, sprintf('<text x="%.1f" y="%d" font-size="%d" fill="#2C3E50" text-anchor="middle" font-family="Courier New">%s</text>',
                      nt_x(i), y0 + ROW_REF_DY, FS, ref_c[i]))
  }
  
  aln_c <- strsplit(d$aln_win, "")[[1]]
  for (i in seq_along(aln_c)) {
    if (aln_c[i] == "-") next
    p <- c(p, sprintf('<text x="%.1f" y="%d" font-size="%d" fill="#2C3E50" text-anchor="middle" font-family="Courier New">%s</text>',
                      nt_x(i), y0 + ROW_SEQ_DY, FS, aln_c[i]))
  }
  # traço horizontal - representar a deleção na linha Seq
  lx0 <- nt_x(d$del_s) - CH/2; lx1 <- nt_x(d$del_e) + CH/2
  ly  <- y0 + ROW_SEQ_DY - FS/2 + 2
  p <- c(p, sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.2"/>',
                    lx0, ly, lx1, ly, st$dash))

  end_x <- MARGIN_L + length(ref_c) * CH + 6
  for (k in 1:2) {
    dy <- if (k == 1) ROW_REF_DY else ROW_SEQ_DY
    p <- c(p, sprintf('<text x="%.1f" y="%d" font-size="11" fill="#888" font-family="Courier New">&#8211;3\u2032</text>',
                      end_x, y0 + dy))
  }

  list(svg = paste(p, collapse = "\n"), end_x = end_x)
}

# SVG completo
panels <- lapply(targets, function(t) parse_allele(t$pop, t$size))

body <- character(0); maxx <- 0
for (k in seq_along(panels)) {
  y0 <- TOP + (k - 1) * PANEL_H       
  r  <- panel_svg(panels[[k]], y0)
  body <- c(body, r$svg); maxx <- max(maxx, r$end_x)
}
W <- ceiling(maxx + 46)                
H <- ceiling(TOP + length(panels) * PANEL_H)

svg <- sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">
<rect width="%d" height="%d" fill="white"/>
<text x="26" y="26" font-size="15" font-weight="700" fill="#1a1a1a" font-family="Arial">Most frequent deletion alleles</text>
<text x="26" y="42" font-size="10.5" fill="#777" font-family="Arial">Clone %s &#183; enhancer %s</text>
%s
</svg>', W, H, W, H, W, H, clone_id, esc(enhancer), paste(body, collapse = "\n"))

svg_file <- paste0(out_base, ".svg")
writeLines(svg, svg_file)             

# Converter para PNG/PDF
img <- image_read_svg(svg_file, width = W * png_scale)   
image_write(img, paste0(out_base, ".png"), format = "png")
image_write(img, paste0(out_base, ".pdf"), format = "pdf")

cat("Guardado:", out_base, ".svg / .png / .pdf\n\n")
cat("Alelos mostrados (posições em coordenadas do enhancer 1-200):\n")
for (p in panels) {
  cat(sprintf("  %s  del -%d bp  |  pos %d-%d  |  %s reads (%.3f%%)\n",
              p$pop, p$size, p$amp_s, p$amp_e,
              formatC(p$reads, format = "d", big.mark = ","), p$pct))
}
