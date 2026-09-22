suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  stop("Package VennDiagram not available. The previous Venn script probably used it; load the same R environment/module.")
}

BASE <- "/scratch/project_2000886/Hoedric/GWAS_2025"
CAMP <- file.path(BASE, "Genetics_Analysis/Campbell_exact_reproduction")
INDIR <- file.path(CAMP, "strict_HetDev_venn")
OUTDIR <- file.path(CAMP, "strict_HetDev_venn_figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

get_counts <- function(dt) {
  A <- as_bool(dt$FST)
  B <- as_bool(dt$DXY)
  C <- as_bool(dt$muLD)
  D <- as_bool(dt$HetDev_TX)

  list(
    area1 = sum(A),
    area2 = sum(B),
    area3 = sum(C),
    area4 = sum(D),
    n12 = sum(A & B),
    n13 = sum(A & C),
    n14 = sum(A & D),
    n23 = sum(B & C),
    n24 = sum(B & D),
    n34 = sum(C & D),
    n123 = sum(A & B & C),
    n124 = sum(A & B & D),
    n134 = sum(A & C & D),
    n234 = sum(B & C & D),
    n1234 = sum(A & B & C & D)
  )
}

plot_one <- function(label, mode) {
  file <- file.path(INDIR, paste0("membership_", label, "_", mode, ".tsv"))
  dt <- fread(file)

  cnt <- get_counts(dt)

  title <- paste0("Campbell SNP-level Venn ", label, " - ", mode)

  png(
    file.path(OUTDIR, paste0("Campbell_SNP_level_Venn_", label, "_", mode, ".png")),
    width = 2200, height = 1900, res = 250
  )

  grid.newpage()

  venn <- VennDiagram::draw.quad.venn(
    area1 = cnt$area1,
    area2 = cnt$area2,
    area3 = cnt$area3,
    area4 = cnt$area4,
    n12 = cnt$n12,
    n13 = cnt$n13,
    n14 = cnt$n14,
    n23 = cnt$n23,
    n24 = cnt$n24,
    n34 = cnt$n34,
    n123 = cnt$n123,
    n124 = cnt$n124,
    n134 = cnt$n134,
    n234 = cnt$n234,
    n1234 = cnt$n1234,
    category = c("FST", "DXY", "muLD", "HetDev_TX"),
    fill = c("#c6dbef", "#fdd0a2", "#c7e9c0", "#fcbba1"),
    alpha = 0.45,
    lwd = 1.2,
    cex = 0.75,
    cat.cex = 0.9,
    margin = 0.08,
    ind = FALSE
  )

  grid.draw(venn)
  grid.text(title, y = unit(0.97, "npc"), gp = gpar(fontsize = 13, fontface = "bold"))

  dev.off()

  cat("Wrote Venn:", label, mode, "\n")
  print(as.data.table(cnt))
}

for (label in c("1pct", "5pct")) {
  for (mode in c("projected", "strictHet")) {
    plot_one(label, mode)
  }
}
