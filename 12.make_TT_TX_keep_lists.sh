#!/bin/bash

BASE=/scratch/project_2000886/Hoedric/GWAS_2025

mkdir -p ${BASE}/Genetics_Analysis/Diversity

awk '$3=="TT"{print $1,$2}' OFS="\t" \
${BASE}/Genetics_Analysis/Population_structure/TT_TX_clusters.txt \
> ${BASE}/Genetics_Analysis/Diversity/TT.keep

awk '$3=="TX"{print $1,$2}' OFS="\t" \
${BASE}/Genetics_Analysis/Population_structure/TT_TX_clusters.txt \
> ${BASE}/Genetics_Analysis/Diversity/TX.keep

echo TT:
wc -l ${BASE}/Genetics_Analysis/Diversity/TT.keep

echo TX:
wc -l ${BASE}/Genetics_Analysis/Diversity/TX.keep
