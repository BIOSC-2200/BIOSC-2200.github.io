ggplot2::theme_set(ggplot2::theme_light())

base_font <- 10

mytheme <- theme(axis.title = element_text(size = base_font + 2,
                                           face = "bold"),
                 plot.title = element_text(size = base_font + 6,
                                           face = "bold"),
                 plot.subtitle = element_text(size = base_font + 4,
                                              face = "italic"),
                 axis.text = element_text(size = base_font),
                 panel.border = element_blank(), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line = element_line(colour = "black"),
                 strip.background = element_rect(fill = "white",
                                                 color = "black"),
                 strip.text = element_text(color = "black",
                                           size = base_font + 2))

genotypes_plot <- function(n_genes = 2) {
  n_genes <- round(n_genes)
  n_phenos <- 2 * n_genes
  
  DD <- data.frame(ways = choose(n_phenos, 0:n_phenos))
  DD$Pct <- DD$ways / sum(DD$ways) * 100
  DD$ID <- seq_len(nrow(DD))
  DD$ways_str <- format(DD$ways, big.mark=",")
  
  if (n_genes == 1) {
    title_str = "1 Gene"
  } else {
    title_str = paste0(n_genes, " Genes")
  }
  
  ways_str <- format(sum(DD$ways), big.mark=",")
  
  P <- ggplot(DD) +
    geom_bar(aes(x = ID, y = Pct), stat = "identity",
             fill = "gray70") +
    scale_y_continuous(limits = c(0, 1.1 * max(DD$Pct))) +
    labs(x = "Number of A alleles",
         y = "Percent",
         title = title_str,
         subtitle = paste0("Total combinations: ",
                           ways_str)) +
    mytheme
  
  if (n_genes < 15) {
    P <- P +
      scale_x_continuous(breaks = DD$ID,
                         labels = DD$ID - 1)
  } else {
    P <- P +
      scale_x_continuous(breaks = DD$ID,
                         labels = c("0",
                                    rep("", times = n_phenos / 2 - 1),
                                    (max(DD$ID) - 1) / 2,
                                    rep("", times = n_phenos / 2 - 1),
                                    max(DD$ID) - 1))
  }
  
  if (n_genes < 9) {
    P <- P + 
      geom_text(aes(x = ID, y = Pct,
                    label = ways_str, vjust = -0.2),
                color = "firebrick",
                fontface = "bold",
                size = 5)
  }
  
  return(P)
}
