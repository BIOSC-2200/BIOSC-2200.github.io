# Adapted from https://github.com/cjbattey/driftR

Simulate_Population <- function(
    n_Populations = 20,
    n_Generations = 100,
    
    Population_Size = 100,
    Initial_Frequency = 0.5,
    
    Fitness_AA = 1,
    Fitness_AB = 1,
    Fitness_BB = 1,
    
    Migration = 0,
    
    Mutation_AB = 0,
    Mutation_BA = 0) {
  
  # Limits
  if (n_Populations > 1000) {
    n_Populations <- 1000
  }
  
  if (Population_Size > 100000) {
    Population_Size <- 100000
  }
  
  if (n_Generations > 1000) {
    n_Generations <- 1000
  }
  
  n <- Population_Size
  p <- Initial_Frequency
  
  stats <- "p"
  infinitePop <- FALSE
  
  allele.freq <- data.frame(matrix(ncol = 4 * n_Populations)) # initialize summary stat matrix
  if(length(p)>1){
    allele.freq[1,(1:n_Populations)] <- p
  } else {
    allele.freq[1,(1:n_Populations)] <- rep(p, n_Populations) #starting allele freqs
  }
  
  for(i in 1:n_Generations){ 
    if(length(n) > 1){ #weight mean allele freq (used for migration) by relative population size if different
      ps <- allele.freq[i,(1:n_Populations)] |> unlist()
      mean.p <- weighted.mean(ps,n)
    } else {
      mean.p <- as.numeric(rowMeans(data.frame(allele.freq[i,(1:n_Populations)])))
    }
    for(j in 1:n_Populations){
      p <- allele.freq[i,j]
      if(length(n)>1){ #get population-specific pop size if multiple listed
        n2 <- n[j]
      } else {
        n2 <- n
      }
      p <- (1 - Mutation_AB) * p + Mutation_BA * (1 - p) # mutation
      p <- p * (1 - Migration) + Migration * mean.p # migration
      q <- 1 - p
      
      if(p > 0 && p < 1){ #if alleles are not fixed
        w <- p*p*Fitness_AA+2*p*q*Fitness_AB+q*q*Fitness_BB #population average fitness
        freq.aa <- (p*p*Fitness_AA)/w #post-selection genotype frequencies (weighted by relative fitness)
        freq.ab <- (2*p*q*Fitness_AB)/w
        if(infinitePop == FALSE){ 
          Naa <- rbinom(1,n2,freq.aa)
          if(freq.aa<1){ 
            Nab <- rbinom(1,(n2-Naa),(freq.ab/(1-freq.aa)))
          }
          else {
            Nab <- 0
          }
          p <- ((2 * Naa)+Nab) / (2*n2)
          q <- 1 - p
          allele.freq[(i+1),j] <- p #new p after drift in columns 1:n_Populations
          allele.freq[(i+1),(j+n_Populations)] <- Nab/n2 #Ho in columns (n_Populations+1):(n_Populations*2)
          allele.freq[(i+1),(j+2*n_Populations)] <- 2*p*q #He in columns (n_Populations*2+1):n_Populations*3
          allele.freq[(i+1),(j+3*n_Populations)] <- w #pop mean fitness in last columns
        } 
        else { #no drift (infinite population) conditions
          p <- freq.aa+(freq.ab/2)
          q <- 1-p
          allele.freq[(i+1),j] <- p
          allele.freq[(i+1),(j+n_Populations)] <- freq.ab 
          allele.freq[(i+1),(j+2*n_Populations)] <- 2*p*q
          allele.freq[(i+1),(j+3*n_Populations)] <- w
        }
      } else { #if alleles are fixed
        if (p <= 0){
          p <- 0
          w <- p * p * Fitness_AA + 2 * p * q * Fitness_AB + q * q * Fitness_BB
        } else {
          p <- 1
          w <- p*p*Fitness_AA + 2 * p * q * Fitness_AB + q * q * Fitness_BB
        }
        allele.freq[(i + 1), j] <- p
        allele.freq[(i + 1), (j + n_Populations)] <- 0
        allele.freq[(i + 1), (j + 2 * n_Populations)] <- 0
        allele.freq[(i + 1), (j + 3 * n_Populations)] <- w
      }
    } #end populations loop
  } #end generations loop
  
  # summary stats
  names <- c()
  for(i in 1:n_Populations) {
    names[i] <- paste0("p", i)
  }
  for(i in (n_Populations + 1):(2 * n_Populations)) {
    names[i] <- paste0("Ho", i - n_Populations)
  }
  for(i in (n_Populations * 2 + 1):(3 * n_Populations)) {
    names[i] <- paste0("He",i - 2 * n_Populations)
  }
  for(i in (n_Populations * 3 + 1):(4 * n_Populations)) {
    names[i] <- paste0("W", i - 3 * n_Populations)
  }
  colnames(allele.freq) <- names
  
  allele.freq$meanHo <- 
    rowMeans(allele.freq[(n_Populations+1):(n_Populations*2)])
  allele.freq$meanHe <- 
    rowMeans(allele.freq[(n_Populations*2+1):(n_Populations*3)])
  allele.freq$Fis <- 
    1 - (allele.freq$meanHo/allele.freq$meanHe)
  allele.freq$mean.p <- rowMeans(allele.freq[1:n_Populations]) 
  allele.freq$Hs <- 
    rowMeans(allele.freq[(n_Populations * 2 + 1):(n_Populations * 3)])
  allele.freq$Ht <- 2 * allele.freq$mean.p * (1 - allele.freq$mean.p)
  allele.freq$Fst <- (allele.freq$Ht - allele.freq$Hs) / allele.freq$Ht
  allele.freq$Fst[allele.freq$Fst < 0] <- 0
  allele.freq$n_Generations <- 0:n_Generations
  allele.freq$Fst[allele.freq$n_Generations == 0] <- NA
  
  allele.freq_melt <- meltPlotData(allele.freq,
                                   n_Populations,
                                   n_Generations)
  plotSingleRun(allele.freq_melt,
                n_Populations = n_Populations,
                n_Generations = n_Generations)
}

meltPlotData <- function(allele.freq.df,
                         n_Populations = n_Populations,
                         n_Generations = n_Generations){
  stats <- "p"
  df <- reshape::melt(allele.freq.df, 
                      id.vars = "n_Generations")
  
  df$dataType <- c(rep("p", (n_Populations*(n_Generations + 1))),
                   rep("Ho", n_Populations*(n_Generations + 1)),
                   rep("He", n_Populations*(n_Generations + 1)),
                   rep("W", n_Populations*(n_Generations + 1)),
                   rep("meanHo", (n_Generations + 1)),
                   rep("meanHe", (n_Generations + 1)),
                   rep("Fis", (n_Generations + 1)),
                   rep("mean.p", (n_Generations + 1)),
                   rep("Hs", (n_Generations + 1)),
                   rep("Ht", (n_Generations + 1)),
                   rep("Fst", (n_Generations + 1)))
  df <- subset(df, dataType %in% stats)
  return(df)
}

plotSingleRun <- function(sim,
                          n_Populations,
                          n_Generations,
                          legend = FALSE,
                          scales = "fixed"){
  P <- ggplot(sim,
              aes(x = n_Generations,
                  y = value,
                  color = variable)) + 
    theme_bw() +
    ylim(0, 1) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold"),
          strip.background = element_blank())+
    scale_color_viridis(discrete = TRUE)+
    labs(x = "Generation", y = "Allele frequency") +
    geom_line()
  return(P)
}
