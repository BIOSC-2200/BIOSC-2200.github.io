library(tidyverse)

NH_demo <- foreign::read.xport("P_DEMO.XPT") |> 
  select(SEQN, RIAGENDR, RIDAGEYR) |> 
  mutate(Genotype = if_else(RIAGENDR == 1, "XY", "XX")) |> 
  rename(Age = RIDAGEYR)
  
NH_data <- foreign::read.xport("P_BMX.XPT") |> 
  select(SEQN, BMXWT, BMXHT) |> 
  rename(Weight = BMXWT,
         Height = BMXHT)

NH <- left_join(NH_demo, NH_data, by = join_by(SEQN)) |> 
  drop_na() |> 
  select(-SEQN, -RIAGENDR) |> 
  select(-Weight) |> 
  as_tibble()

write_csv(NH, file = "../NHANES.csv")

