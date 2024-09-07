#--------1. Setup--------#
#working_dir <- "L:/basic/divg/KVI/HIV_RNA_seq/Thijs/set1"
working_dir<-'C:/Users/tcdor/OneDrive/Desktop/Set1Scripts'

# Set up libraries
library(dplyr)
library(readxl)
library(writexl) 
library(ggvenn)



setwd(working_dir)

# Set conditions
cond1 <- "MDM_IL4_Mock"
cond2 <- "MDM_IL4_R5"
cond3 <- "MDM_IL10_Mock"
cond4 <- "MDM_IL10_R5"
cond5 <- "MDM_IFNgTNFa_Mock"
cond6 <- "MDM_IFNgTNFa_R5"
cond7 <- "MDM_Mock"  # For comparison
cond8 <- "MDM_R5"   # For comparison

# Read the Excel files
MDM_IL4_Mock_vs_IL4_R5 <- read_excel(paste0("set1_analysis/", cond1, "_vs_", cond2, "/", cond1, "_vs_", cond2, "_sigLog_inf.xlsx"))
MDM_IL10_Mock_vs_IL10_R5 <- read_excel(paste0("set1_analysis/", cond3, "_vs_", cond4, "/", cond3, "_vs_", cond4, "_sigLog_inf.xlsx"))
MDM_IFNgTNFa_Mock_vs_IFNgTNFa_R5 <- read_excel(paste0("set1_analysis/", cond5, "_vs_", cond6, "/", cond5, "_vs_", cond6, "_sigLog_inf.xlsx"))
MDM_Mock_vs_R5 <- read_excel(paste0("set1_analysis/", cond7, "_vs_", cond8, "/", cond7, "_vs_", cond8, "_sigLog_inf.xlsx"))

# SET UP THE VENN

# Extract gene names and regulation direction from each data frame
# Consolidate gene names and directions
genes_IL4 <- unique(c(MDM_IL4_Mock_vs_IL4_R5$ILMN_GENE))
direction_IL4 <- ifelse(MDM_IL4_Mock_vs_IL4_R5$logFC > 0, "Upregulated", "Downregulated")

genes_IL10 <- unique(c(MDM_IL10_Mock_vs_IL10_R5$ILMN_GENE))
direction_IL10 <- ifelse(MDM_IL10_Mock_vs_IL10_R5$logFC > 0, "Upregulated", "Downregulated")

genes_IFNgTNFa <- unique(c(MDM_IFNgTNFa_Mock_vs_IFNgTNFa_R5$ILMN_GENE))
direction_IFNgTNFa <- ifelse(MDM_IFNgTNFa_Mock_vs_IFNgTNFa_R5$logFC > 0, "Upregulated", "Downregulated")

genes_Mock_vs_R5 <- unique(c(MDM_Mock_vs_R5$ILMN_GENE))
direction_Mock_vs_R5 <- ifelse(MDM_Mock_vs_R5$logFC > 0, "Upregulated", "Downregulated")

# Create a list of gene sets
gene_sets1 <- list(
  IL4 = genes_IL4,
  IL10 = genes_IL10,
  IFNgTNFa = genes_IFNgTNFa
)



# Define colors for the Venn diagram
light_colors1 <- c("#ADD8E6", "#FFCCCB","#90EE90")
darker_colors1 <- c("#5F9CC4","#ff9594", "#56A564")
light_colors2 <- c("#ADD8E6", "#FFCCCB",  "#90EE90","#D8BFD8" )
darker_colors2 <- c("#5F9CC4","#ff9594", "#56A564","#6A0D91")


venn_plot_inf <- ggvenn(
  gene_sets1,
  fill_color = light_colors1,
  stroke_size = 0,
  text_size = 4,
  set_name_size = 5,
  set_name_color = darker_colors1,
  show_percentage = TRUE)

ggsave(filename = "venn_plot_inf.pdf", plot = venn_plot_inf, width = 8, height = 6, dpi = 300)

