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


# Plot genotypes #############################################################

genotypes_plot <- function(n_genes = 2) {
  
  # Use only integers
  n_genes <- round(n_genes) |> as.integer()
  
  # Two alleles for each gene
  n_alleles <- 2 * n_genes
  
  # Calculate the ways to get different numbers of major allele
  DD <- data.frame(ways = choose(n_alleles, 0:n_alleles))
  DD$Pct <- DD$ways / sum(DD$ways) * 100
  
  # Add a column with just 1:n_alleles
  DD$ID <- seq_len(nrow(DD))
  
  # Format a string to use for labeling the bars
  DD$ways_str <- format(DD$ways, big.mark = ",")
  
  # Format title string
  if (n_genes == 1) {
    title_str = "1 Gene"
  } else {
    title_str = paste0(n_genes, " Genes")
  }
  
  # Format the "ways" string
  ways_str <- format(sum(DD$ways), big.mark = ",")
  
  # Create the base plot. Add a little vertical space to accommodate
  # the number of ways at the top of the bar.
  P <- ggplot(DD) +
    geom_bar(aes(x = ID, y = Pct), stat = "identity",
             fill = "gray70") +
    scale_y_continuous(limits = c(0, 1.1 * max(DD$Pct))) +
    labs(x = "Allele Combination",
         y = "Percent",
         title = title_str,
         subtitle = paste0("Total combinations: ",
                           ways_str)) +
    mytheme
  
  if (n_genes < 9) {
    P <- P + 
      geom_text(aes(x = ID, y = Pct,
                    label = ways_str, vjust = -0.2, hjust = 0.5),
                color = "firebrick",
                fontface = "bold",
                size = 4)
  }

  if (n_genes < 12) {
    AA <- DD$ID - 1
    TT <- n_alleles - AA
    x_labs <- paste0(AA, " A\n", TT, " T")
    
    P <- P +
      scale_x_continuous(breaks = DD$ID,
                         labels = x_labs)
  } else {
    P <- P +
      scale_x_continuous(breaks = DD$ID,
                         labels = c(paste0("0 A\n", n_alleles, " T"),
                                    rep("", times = n_alleles / 2 - 1),
                                    paste0(n_alleles / 2, " A\n",
                                           n_alleles / 2, " T"),
                                    rep("", times = n_alleles / 2 - 1),
                                    paste0(n_alleles, " A\n0 T")))
  }
  
  return(P)
}


## Normal distribution by simulation #########################################

simulate_heights <- function(n_genes = 25, XX = NULL) {
  set.seed(3242343)
  
  # If XX is not supplied, then we are working in webr.
  # Load locally from the downloaded file from the setup chunk.
  # Otherwise, XX will be loaded locally in an {r} chunk.
  if (is.null(XX)) {
    XX <- read.csv("NHANES.csv") |> 
      filter(Genotype == "XX" & Age > 20)
  }
  
  n_individuals <- nrow(XX)
  
  n_alleles <- n_genes * 2
  
  h_bar <- mean(XX$Height)
  h_q <- quantile(XX$Height, c(0.025, 0.975))

  q_range <- as.numeric(h_q[2] - h_q[1])
  
  n_A <- rbinom(n = n_individuals, size = n_alleles, prob = 0.5)
  n_a <- n_alleles - n_A
  
  # This is the correction factor to rescale the allele effect size to
  # properly generate a distribution with the same phenotypic variance
  # as the observed data
  correction <- (q_range / max(n_A - n_a))
  
  # Combine the NHANES data with the simulated data
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
         title = paste0(n_genes, " Genes, ", n_alleles, " Alleles"),
         subtitle = paste0("Each allele = ±", round(correction, 2), " cm")) +
    mytheme
  
}
