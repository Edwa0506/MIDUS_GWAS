#######################################
#
# Setup and useful functions
#
#######################################

setwd("/home/midus_1114420/edwa0506")

library(tidyverse)


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

#######################################
#
# Create PCs and look at self-identified Whites on PC plot
#
#######################################



# First remove non-autosomes
system2("./plink",
        args = c("--bfile", "data/Imputed_SNPs/Imputed/M2MRMJ_Imputed",
                 "--chr", "1-22",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imputed_autosomes"),
        stdout = TRUE,
        stderr = TRUE)


# Now identify related individuals so we can exclude them from PC calculation
# we will project PCs onto them 

system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_autosomes",
                 "--maf", "0.10",
                 "--rel-cutoff", "0.025",
                 "--out", "data/related_individuals"),
        stdout = TRUE,
        stderr = TRUE)

fam_file <- read_table("data/filtered_binary_files/M2MRMJ_Imputed_autosomes.fam",
                       col_names =c("FID", "IID", "M", "P", "S", "pheno"))

related <- read_table("data/related_individuals.rel.id", # Check ending
                      col_names = c("IID"))

fam_file$group <- ifelse(fam_file$IID %in% related$IID, "g2", "g1")

fam_file |> select(FID, IID, group) |> 
            write.table("data/related_cluster_file.txt",
                         sep = " ",
                         quote = FALSE,
                         col.names = FALSE,
                         row.names = FALSE) 

#FID, IID, cluster name 


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_autosomes",
                 "--maf", "0.01",
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
                 "--within", "data/related_cluster_file.txt",
                 "--pca",
                 "--pca-cluster-names", "g1",
                 "--out", "data/all_races_pca"),
        stdout = TRUE,
        stderr = TRUE)

PC <- read_table("data/all_races_pca.eigenvec", 
                 col_names = c("ID", "FAMID", paste0("PC",1:20)))


MJ <-  read_csv("data/Imputed_SNPs/MIDJA Phenotypic Data/MIDJAGeneticsMetadata_n=328_07-12-22.csv")
M2 <-   read_csv("data/Imputed_SNPs/MIDUS 2 Metadata/M2GeneticsMetadata_n=980_07-12-22.csv")
MR <-  read_csv("data/Imputed_SNPs/MIDUS Refresher Metadata/MRGeneticsMetadata_n=809_07-12-22.csv")
# 2117 people altogether

MIDUS <- rbind.data.frame(MJ, M2, MR)

PC <- PC |>
      mutate(ID = as.character(ID)) |>
      left_join(MIDUS, by = c("ID" = "M2MRMJ_ID_S")) 

PC_plot(PC, "PC1", "PC2")
PC_plot(PC, "PC2", "PC3")
PC_plot(PC, "PC3", "PC4")
PC_plot(PC, "PC4", "PC5")

PC_plot_EUR(PC, "PC1", "PC2")
PC_plot_EUR(PC, "PC2", "PC3")
PC_plot_EUR(PC, "PC3", "PC4")
PC_plot_EUR(PC, "PC4", "PC5")

PCfilt <- PC |>
  filter(PC1 > 0.005,
         PC2 < -0.015,
         PC3 < 0.10,
         PC4 > -0.01,
         PC4 < 0.01,
         PC5 < 0.004,
         PC5 > -0.015,
         PC6 > -0.02,
         PC6 < 0.02,
         PC7 < 0.025,
         PC8 < 0.025,
         PC8 > -0.025,
         PC9 < 0.02,
         PC10 < 0.1,
         PC10 > -0.1,
         PC11 < 0.02,
         PC12 < 0.070,
         PC14 < 0.02)
  
PC_plot(PCfilt, "PC1", "PC2")
PC_plot(PCfilt, "PC3", "PC4")
PC_plot(PCfilt, "PC4", "PC5")
PC_plot(PCfilt, "PC5", "PC6") #
PC_plot(PCfilt, "PC6", "PC7")
PC_plot(PCfilt, "PC7", "PC8")
PC_plot(PCfilt, "PC8", "PC9")
PC_plot(PCfilt, "PC9", "PC10")
PC_plot(PCfilt, "PC11", "PC12")
PC_plot(PCfilt, "PC13", "PC14")



PCfilt |> select(ID, FAMID) |>
  write.table("data/ID_PC_outlier_removed.txt", quote = FALSE, sep = " ", 
              row.names = FALSE, col.name = FALSE)


#######################################
#
# Create PCs
#
#######################################


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_Imputed_autosomes",
                 "--keep", "data/ID_PC_outlier_removed.txt",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_EUR"),
        stdout = TRUE, 
        stderr = TRUE)


system2("./plink",
        args = c("--bfile", "data/filtered_binary_files/M2MRMJ_EUR",
                 "--maf", "0.05",
                 "--indep-pairwise", "50", "10", "0.1",
                 "--out", "data/ld_pruned_snps_for_pca"),
        stdout = TRUE,
        stderr = TRUE)


system2("./plink",
        args = c("--bfile",  "data/filtered_binary_files/M2MRMJ_EUR",
                 "--extract","data/ld_pruned_snps_for_pca.prune.in",
                 "--make-bed",
                 "--out", "data/filtered_binary_files/M2MRMJ_Imp_outlier_removed_ld_pruned"),
        stdout = TRUE,
        stderr = TRUE)

system2("./plink",
        args = c("--bfile",  "data/filtered_binary_files/M2MRMJ_Imp_outlier_removed_ld_pruned",
                 "--within", "data/related_cluster_file.txt",
                 "--pca",
                 "--pca-cluster-names", "g1",
                 "--out", "data/pca"),
        stdout = TRUE,
        stderr = TRUE)



PC <- read_table("data/pca.eigenvec", col_names = c("ID", "FAMID", paste0("PC",1:20)))


MJ <-  read_csv("data/Imputed_SNPs/MIDJA Phenotypic Data/MIDJAGeneticsMetadata_n=328_07-12-22.csv")
M2 <-   read_csv("data/Imputed_SNPs/MIDUS 2 Metadata/M2GeneticsMetadata_n=980_07-12-22.csv")
MR <-  read_csv("data/Imputed_SNPs/MIDUS Refresher Metadata/MRGeneticsMetadata_n=809_07-12-22.csv")

MIDUS <- rbind.data.frame(MJ, M2, MR)

PC <- PC |>
  mutate(ID = as.character(ID)) |>
  left_join(MIDUS, by = c("ID" = "M2MRMJ_ID_S")) 
dim(PC)


PC_plot(PC, "PC1", "PC2")
PC_plot(PC, "PC2", "PC3")
PC_plot(PC, "PC3", "PC4")
PC_plot(PC, "PC5", "PC6")
PC_plot(PC, "PC7", "PC8")
PC_plot(PC, "PC9", "PC10")
