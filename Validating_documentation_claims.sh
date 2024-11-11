# Documentation doesn't clearly match the files.
# e.g. documentations claims QC was performed and 8 million SNPs were remaining
# Below code checks that info cutoff has been applied and what happens
# After QC has been applied


# working directory
cd /home/midus_1114420/edwa0506

module load plink 

wc -l data/Imputed_SNPs/Imputed/M2MRMJ_Imputed.bim
# 15249697
wc -l data/Imputed_SNPs/Imputed/WellImputedSNPs.txt
# 15249720

plink --bfile data/Imputed_SNPs/Imputed/M2MRMJ_Imputed \
      --extract data/Imputed_SNPs/Imputed/WellImputedSNPs.txt \
      --make-bed \
      --out data/M2MRMJ_WellImputed
      

# Warning: 307245 het. haploid genotypes present (see data/M2MRMJ_WellImputed.hh
# ); many commands treat these as missing.
# Surely just choose autosomal only SNPs? Should remove this warning
      
      
wc -l data/M2MRMJ_WellImputed.bim
# 15249697
# Ok, we've confirmed the file M2MRMJ_Imputed already filters by R2
# According to the documentation this should be  R2 < 0.5

rm data/M2MRMJ_WellImputed*
# remove results of above test


plink --bfile data/Imputed_SNPs/Imputed/M2MRMJ_Imputed \
      --maf 0.01 \
      --geno 0.1 \
      --make-bed \
      --out data/M2MRMJ_filtered
      
wc -l data/M2MRMJ_filtered.bim
# 9124155
# Documentation says 8 million SNPs remaining
# Perhaps they removed autosomes too? 

rm data/M2MRMJ_filtered*


plink --bfile data/Imputed_SNPs/Imputed/M2MRMJ_Imputed \
      --maf 0.01 \
      --geno 0.1 \
      --chr 1-22 \
      --make-bed \
      --out data/M2MRMJ_filtered
      
wc -l data/M2MRMJ_filtered.bim
# 8872449 data/M2MRMJ_filtered.bim
# Closer to 8 million...


# have duplicated BP been deleted? 

plink --bfile data/Imputed_SNPs/Imputed/M2MRMJ_Imputed \
      --list-duplicate-vars 
    
head plink.dupvar 
# Yep, no duplicates

rm plink.dupvar*

      


