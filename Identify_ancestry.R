setwd("/home/midus_1114420/edwa0506")

# install.packages("rgl")
library(tidyverse)
library(rgl)

#######################################
#
# Create PCs and look at self-identified Whites on PC plot
#
#######################################

system2("./plink",
        args = c("--bfile", "data/Imputed_SNPs/Imputed/M2MRMJ_Imputed",
                 "--chr", "1-22",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imputed_autosomes"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_autosomes",
                 "--indep-pairwise", "50", "10", "0.1",
                 "--out", "data/all_races_ld_pruned_snps_for_pca"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile",  "data/filtered_binary_files/M2MRMJ_Imputed_autosomes",
                 "--extract","data/all_races_ld_pruned_snps_for_pca.prune.in",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imputed_autosomes_all_races_ld_pruned"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile",  "data/filtered_binary_files/M2MRMJ_Imputed_autosomes_all_races_ld_pruned",
                 "--pca",
                 "--out", "data/all_races_pca"),
        stdout = TRUE,
        stderr = TRUE)

PC <- read_table("data/all_races_pca.eigenvec", col_names = c("ID", "FAMID", paste0("PC",1:20)))


MJ <-  read_csv("data/Imputed_SNPs/MIDJA Phenotypic Data/MIDJAGeneticsMetadata_n=328_07-12-22.csv")
M2 <-   read_csv("data/Imputed_SNPs/MIDUS 2 Metadata/M2GeneticsMetadata_n=980_07-12-22.csv")
MR <-  read_csv("data/Imputed_SNPs/MIDUS Refresher Metadata/MRGeneticsMetadata_n=809_07-12-22.csv")
# 2117 people altogether

MIDUS <- rbind.data.frame(MJ, M2, MR)

PC <- PC |>
      mutate(ID = as.character(ID)) |>
      left_join(MIDUS, by = c("ID" = "M2MRMJ_ID_S")) 

ggplot(PC, aes(x = PC1, y = PC2, color = RACE)) +
  geom_point() +
  labs(x = "PC 1", y = "PC 2", color = "Race") +
  theme_minimal()

PC |> filter(RACE %in% c("WHITE",
                         "REFUSED",
                         "OTHER (SPECIFY)")) |>
  ggplot(aes(x = PC1, y = PC2, color = RACE)) +
  geom_point() +
  labs(x = "PC 1", y = "PC 2", color = "Race") +
  theme_minimal()

PC_plot <- function(data, PCX, PCY) {
  data |>
    ggplot(aes_string(x = PCX, y = PCY, color = "RACE")) +
    geom_point() +
    labs(x = PCX, y = PCY, color = "Race") +
    theme_minimal()
}

PC_plot_EUR <- function(data, PCX, PCY) {
  data |>
    filter(RACE %in% c("WHITE", "REFUSED", "OTHER (SPECIFY)")) |>
    ggplot(aes_string(x = PCX, y = PCY, color = "RACE")) +
    geom_point() +
    labs(x = PCX, y = PCY, color = "Race") +
    theme_minimal()
}

PC_plot(PC, "PC1", "PC2")
PC_plot(PC, "PC2", "PC3")
PC_plot(PC, "PC3", "PC4")
PC_plot(PC, "PC4", "PC5")

PC_plot_EUR(PC, "PC1", "PC2")
PC_plot_EUR(PC, "PC2", "PC3")
PC_plot_EUR(PC, "PC3", "PC4")
PC_plot_EUR(PC, "PC4", "PC5")

PCfilt <- PC |>
  filter(PC1 < - 0.010,
         PC2 < 0.02,
         PC2 > 0.009,
         PC3 < 0.003,
         PC4 < 0.001,
         PC4 > -0.0015,
         PC5 < 0.004,
         PC5 > -0.25,
         PC6 < 0.005,
         PC6 > -0.013,
         PC7 < 0.0025,
         PC7 > -0.003) 
  
PC_plot(PCfilt, "PC1", "PC2")
PC_plot(PCfilt, "PC2", "PC3")
PC_plot(PCfilt, "PC3", "PC4")
PC_plot(PCfilt, "PC4", "PC5")
PC_plot(PCfilt, "PC5", "PC6") #
PC_plot(PCfilt, "PC6", "PC7")
PC_plot(PCfilt, "PC7", "PC8")
PC_plot(PCfilt, "PC8", "PC9")
PC_plot(PCfilt, "PC9", "PC10")



PCfilt |> select(ID, FAMID) |>
  write.table("data/ID_PC_outlier_removed.txt", quote = FALSE, sep = " ", row.names = FALSE, col.name = FALSE)


#######################################
#
# QC data and check new PCs
#
#######################################


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_autosomes",
                 "--keep", "data/ID_PC_outlier_removed.txt",
                 "--MAF", "0.01",
                 "--"
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imp_EUR_QC"),
        stdout = TRUE, 
        stderr = TRUE)



#######################################
#
# Keep only self-ID white and create PCs
#
#######################################
# Select only self-identified Europeans
# Create PCs
# Examine PC plot for outliers

MJ <-  read_csv("data/Imputed_SNPs/MIDJA Phenotypic Data/MIDJAGeneticsMetadata_n=328_07-12-22.csv")
M2 <-   read_csv("data/Imputed_SNPs/MIDUS 2 Metadata/M2GeneticsMetadata_n=980_07-12-22.csv")
MR <-  read_csv("data/Imputed_SNPs/MIDUS Refresher Metadata/MRGeneticsMetadata_n=809_07-12-22.csv")
# 2117 people altogether

rbind.data.frame(MJ, M2, MR) |>
  filter(RACE == "WHITE") |>
  select(M2MRMJ_ID_S) |>
  mutate(famid = M2MRMJ_ID_S) |>
      write.table("data/self_idwhites.txt",
                  quote = FALSE, sep = " ", row.names =FALSE, col.names = FALSE)
# 1347

iteration <- "1"
system2("./plink",
        args = c("--bfile", "data/Imputed_SNPs/Imputed/M2MRMJ_Imputed",
                 "--keep", "data/self_idwhites.txt",
                 "--chr", "1-22",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white"),
        stdout = TRUE, 
        stderr = TRUE)

        
system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white",
                 "--indep-pairwise", "50", "10", "0.1",
                 "--out", paste0("data/ld_pruned_snps_for_pca_",iteration)),
        stdout = TRUE,
        stderr = TRUE)
        
        
system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white",
                 "--extract",paste0("data/ld_pruned_snps_for_pca_",iteration),
                 "--make-bed",
                 "--out", paste("data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_ld_pruned_", iteration)),
        stdout = TRUE,
        stderr = TRUE)
      
        
system2("wc",
        args = c("-l", 
                 "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_ld_pruned.bim"),
        stdout = TRUE,
        stderr = TRUE)
# Leaves 5 million variants for PCA. Probably far too many tbh.
# But not really a problem given the small sample size. 
        
system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_ld_pruned",
                 "--pca",
                 "--out", "data/principal_components"),
        stdout = TRUE,
        stderr = TRUE)

PC <- read_table("data/principal_components.eigenvec", col_names = c("ID", "FAMID", paste0("PC",1:20)))

plot(PC$PC1, PC$PC2)
plot(PC$PC2, PC$PC3)
plot(PC$PC3, PC$PC4)
plot(PC$PC4, PC$PC5)
plot(PC$PC5, PC$PC6)
plot(PC$PC7, PC$PC8)
plot(PC$PC9, PC$PC10)


# Remove, PC1 > 0.4, PC2 < -0.4, PC3 > 0.4, 

PCfilt <- filter(PC, PC1 < 0.4, 
                 PC2 > -0.4,
                 PC3 < 0.4,
                 PC3 > -0.2,
                 PC4 < 0.175,
                 PC4 > -0.1,
                 PC5 < 0.4,
                 PC5 > -0.2,
                 PC6 < 0.1,
                 PC6 > 0.1,
                 PC7 < 0.1,
                 PC7 > -0.1,
                 PC8 < 0.1,
                 PC8 > 0.1)
nrow(PCfilt)
#1341 versys 1347 before filter

PCfilt |> select(ID, FAMID) |>
  write.table("data/ID_PC_outlier_removed.txt", quote = FALSE, sep = " ", row.names = FALSE, col.name = FALSE)


#######################################
#
# Rereate PCs without outliers
#
#######################################

system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white",
                 "--keep", "data/ID_PC_outlier_removed.txt",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers"),
        stdout = TRUE, 
        stderr = TRUE)


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers",
                 "--indep-pairwise", "50", "10", "0.1",
                 "--out", "data/ld_pruned_snps_for_pca_after_outlier_removed"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers",
                 "--extract","data/ld_pruned_snps_for_pca_after_outlier_removed.prune.in",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers_ld_pruned"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers_ld_pruned",
                 "--pca",
                 "--out", "data/principal_components_outliers_removed"),
        stdout = TRUE,
        stderr = TRUE)


PC <- read_table("data/principal_components_outliers_removed.eigenvec", col_names = c("ID", "FAMID", paste0("PC",1:20)))

plot(PC$PC1, PC$PC2)
plot(PC$PC2, PC$PC3)
plot(PC$PC3, PC$PC4)
plot(PC$PC4, PC$PC5)
plot(PC$PC5, PC$PC6)

######################################################
#
# Second iteration: Rereate PCs for without outliers
#
######################################################

PCfilt <- PC |> filter(PC1 > -0.05, 
                       PC2 > -0.2,
                       PC2 < 0.1,
                       PC3 < 0.1,
                       PC3 > -0.2,
                       PC4 < 0.175,
                       PC4 > -0.1,
                       PC5 < 0.1,
                       PC5 > -0.2)

nrow(PCfilt)
#1319


PCfilt |> select(ID, FAMID) |>
  write.table("data/ID_PC_outlier_removed.txt", quote = FALSE, sep = " ", row.names = FALSE, col.name = FALSE)

system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers",
                 "--keep", "data/ID_PC_outlier_removed.txt",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers"),
        stdout = TRUE, 
        stderr = TRUE)


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers",
                 "--indep-pairwise", "50", "10", "0.1",
                 "--out", "data/ld_pruned_snps_for_pca_after_outlier_removed"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers",
                 "--extract","data/ld_pruned_snps_for_pca_after_outlier_removed.prune.in",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers_ld_pruned"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_selfid_white_no_PC_outliers_ld_pruned",
                 "--pca",
                 "--out", "data/principal_components_outliers_removed"),
        stdout = TRUE,
        stderr = TRUE)


PC <- read_table("data/principal_components_outliers_removed.eigenvec", col_names = c("ID", "FAMID", paste0("PC",1:20)))

plot(PC$PC1, PC$PC2)
plot(PC$PC2, PC$PC3)
plot(PC$PC3, PC$PC4)
plot(PC$PC4, PC$PC5)
plot(PC$PC5, PC$PC6)

