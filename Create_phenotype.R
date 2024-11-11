
setwd("/home/midus_1114420/edwa0506")

library(tidyverse)
library(haven)
library(fixest)

fam_file <- read_table("data/filtered_binary_files/M2MRMJ_EUR.fam", col_names = FALSE)

princ <- read_table("data/pca.eigenvec", 
                    col_names = c("MIDUSID", "FAMID", paste0("PC", 1:20)))

midus <- read_sav("data/Midus_phenotypes.sav", user_na = FALSE) |>
         as_tibble() |>
         zap_labels() |>
         filter(MIDUSID %in% fam_file$X1) |>
         left_join(princ, by = "MIDUSID")

vote_names <- c("RA4Q11S",#MR1P4 RA1PRAGE
                "B4Q11S", # M2P4 B1PRAGE_2019
                "C4Q11S", # M3P4  C1PRAGE
                "A1SK7L") # M1P1 A1PRAGE_2019
ages <- c("RA1PRAGE",
          "B1PRAGE_2019",
          "C1PRAGE",
          "A1PRAGE_2019")

# Japanese only names
#"K2Q47B",
#"J2Q47B"

midus |> select(contains("AGE"), contains("RAGE")) |> names()

midus |> select(contains("AGE"), contains("RAGE")) |> View()

# Which age variables correspond to each vote variable? 


map(vote_names, ~ table(midus[[.x]]))

# Quite a few people have multiple measurements. Need to aggregate them. 
n_vote <- rowSums(!is.na(midus[vote_names])) 
table(n_vote)
hist(n_vote)

residual_df <- map(1:length(vote_names), \(x) as.formula(paste(vote_names[x]," ~ ", ages[x]," + ",
                                                paste("PC", 1:20, sep ="", collapse =" + "), sep=""))) |> 
  map(\(x) lm(x, data = midus, na.action = na.exclude)) |>
  map(\(x) residuals(x)) |>
  map(\(x) scale(x)) |>
  bind_cols() 

colnames(residual_df) <- paste(vote_names, "_residuals", sep = "")

midus <- midus |> cbind(residual_df) %>%
  mutate(vote_phenotype = rowMeans(select(., ends_with("_residuals")), 
                                   na.rm = TRUE))

# Solve this, There shouldn't be any NA..
# sum(is.na(midus$vote_phenotype))
# [1] 2

list_rbind(test) |> head()


results <- map(1:length(vote_names), \(x) {
  as.formula(paste(vote_names[x], " ~ ", ages[x], " + ",
                   paste("PC", 1:20, sep ="", collapse =" + "), sep=""))
}) |> 
  map(\(x) lm(x, data = midus, na.action = na.exclude)) |>
  map(\(x) residuals(x)) |>
  map(\(x) scale(x)) |>

# Combine the residuals into a tibble with the appropriate column names
residuals_tibble <- tibble(
  `resid_` = results
)

# Rename the columns of the tibble
colnames(residuals_tibble) <- paste0("resid_", vote_names)


scale <-function(x) {(x - mean(x, na.rm = TRUE))/sd(x, na.rm=TRUE)}
