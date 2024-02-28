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


# Plot genotypes

genotypes_plot <- function(n_genes = 2) {
  n_genes <- round(n_genes)
  n_alleles <- 2 * n_genes
  
  DD <- data.frame(ways = choose(n_alleles, 0:n_alleles))
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
                                    rep("", times = n_alleles / 2 - 1),
                                    (max(DD$ID) - 1) / 2,
                                    rep("", times = n_alleles / 2 - 1),
                                    max(DD$ID) - 1))
  }
  
  if (n_genes < 9) {
    P <- P + 
      geom_text(aes(x = ID, y = Pct,
                    label = ways_str, vjust = -0.2),
                color = "firebrick",
                fontface = "bold",
                size = 4)
  }
  
  return(P)
}


## Normal distribution by simulation

simulate_heights <- function(n_genes = 25) {
  set.seed(3242343)
  
  download.file("https://raw.githubusercontent.com/BIOSC-2200/BIOSC-2200.github.io/main/QT/NHANES.csv",
                "NHANES.csv")
  XX <- read.csv("NHANES.csv") |> 
    filter(Genotype == "XX" & Age > 20)
  
  n_individuals <- nrow(XX)
  
  n_genes <- 25
  n_alleles <- n_genes * 2
  
  h_bar <- mean(XX$Height)
  h_sd <- sd(XX$Height)
  h_q <- quantile(XX$Height, c(0.025, 0.975))
  h_range <- range(XX$Height)
  
  q_range <- as.numeric(h_q[2] - h_q[1])
  
  n_A <- rbinom(n = n_individuals, size = n_alleles, prob = 0.5)
  n_a <- n_alleles - n_A
  
  correction <- (q_range / max(n_A - n_a))
  
  HTcomp <- bind_rows(
    tibble(Height = XX$Height,
           Set = "NHANES"),
    tibble(Height = h_bar + 
             (n_A - n_a) * correction + 
             rnorm(n = n_individuals, 0, 0.5),
           Set = "Simulation")
  )
  
  ggplot(HTcomp, aes(Height, fill = Set)) +
    geom_histogram(bins = 30) +
    facet_grid(Set ~ ., scales = "free_y") +
    scale_fill_manual(values = c("goldenrod", "darkblue"), guide = "none") +
    labs(x = "Height (cm)", y = "Count",
         title = paste0(n_alleles, " Alleles"),
         subtitle = paste0("Each allele = ±", round(correction, 2), " cm")) +
    mytheme
  
}