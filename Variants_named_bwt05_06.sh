### The Variants have to be named before the LD pruning and PCA base steps:

plink \
  --bfile /scratch/project_2000886/Hoedric/GWAS_2025/PreGWAS/result/GWAS_qc_maf001_geno005_mind005 \
  --allow-extra-chr \
  --set-missing-var-ids @:#:\$1:\$2 \
  --make-bed \
  --out /scratch/project_2000886/Hoedric/GWAS_2025/PreGWAS/result/GWAS_qc_named
