#--------INTRODUCTION-----#
#This script is written to use DESeq2 and WGCNA to analyze Langerhans cell data obtained from RNA-seq data. 
#The initial input is the count table.
#Initially, ,DESeq2-based transcriptomic analysis is performed
#Next, WGCNA is performed
#Lastly, the intersect of these 2 is obtained
#please check out the subsetting section a bit down below for setting specific analysis subsets



#--------1. Setup--------#

working_dir<- "C:/Users/tcdor/OneDrive/Desktop/DEG_Analysis"
#working_dir<- "L:/basic/divg/KVI/FR/Theseus/Analyse_MaartjeLC/Analysis_salmon/DEG_Analysis_LC"
#working_dir<- "L:/basic/divg/KVI/FR/Theseus/Analyse_MaartjeLC/Analysis_STAR/"

# set up libraries
library(dplyr)
library(VennDiagram)
library(ggvenn)
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(here)
library(readr)
library(readxl)
library(writexl)
library(mixOmics)
library(limma)
library(tidyr)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(DOSE)
library(fgsea)
library(annotate)
library(pheatmap)
library(viridis)
library(enrichplot)
library(CorLevelPlot)
library(DT)
library(DOSE)
library(optparse)
library(biomaRt)
library(DESeq2)
library(edgeR)
library(apeglm)
library(WGCNA)
library(impute)
library(gridExtra)
library(EnsDb.Hsapiens.v86)
library(glue)
library(dendextend)

setwd(working_dir)


#setting initial Conditions

# choose conditions of interest
comp <- "maturity" # if you change the comparison you will need to edit the design and contrast further down below!

# create essential directories
dir.create(comp)
dir.create(paste0(comp, "/excel_exports"))
dir.create(paste0(comp, "/plots"))

# choose significance threshold & FC thres
pthres<- 0.05
adjpthres<-0.05
log_FC<- 1

# set up file selection
inf_extract <- "interferome_genes.xlsx" # this is optional in case you want to focus in inteferon stimulated genes specifically, for example
#normalized_data <- "salmon.merged.gene_counts_length_scaled.tsv" # this is your input count table for the salmon only data
normalized_data <- "tx_salmon_counts_length_scaled.tsv" # this is your input count table for the star-salmon data

#--------2. Data Load--------#

# load interferome list (OPTIONAL)
inf_genes_table <- read_xlsx(inf_extract, sheet =4)
inf_genes_names <- unique(inf_genes_table$Gene_Name) # create gene names vector
inf_genes_gene_id<-unique(inf_genes_table$`Ensembl Id`)

# load expression data
data_import <- read_tsv(normalized_data)

original_colnames <- colnames(data_import)[3:26]

# create dataframe with integers suitable for deseq analysis
data_d <- data_import[3:26] #slelect only the samples
data_d <- as.data.frame(lapply(data_d, function(x) round(x)))
rownames(data_d) <-data_import$gene_id
# reapply the original column names (with hyphens)
colnames(data_d) <- original_colnames
metadata_df <- data.frame(read_name = colnames(data_import)[3:26])#set up metadata info

#add donor numbers
donor_sequence <- rep(1:6, each = 4)
metadata_df$donor_num <- paste0("donor ", donor_sequence)

#add cell type
metadata_df$cell_type <- rep(c("iLC", "miLC"), each=2, length.out = 24)

#add cd_status
metadata_df$cd_status <- rep(c("CD1A+_CD3-", "CD1A+_CD3+"), each=1, length.out = 24)


#------------SUBSETTTING FOR ANALYSIS-------------------#
# we want to focus only on the CD1A+_CD3- samples
metadata_df<- metadata_df[metadata_df$cd_status == "CD1A+_CD3-", ]
head(data_d)
head(metadata_df)
data_d <- data_d[, colnames(data_d) %in% metadata_df$read_name]
head(data_d)
#------------END OF SUBSETTING
# plot PCA
pca_d<-mixOmics::pca(t(data_d))
dir.create(paste0(comp, "/plots/PCA"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/PCA/", comp, "_PCA.pdf")))
plotIndiv(pca_d, group=metadata_df$cell_type, pch=16,legend=TRUE, title = '', )
dev.off()

#----------------2. DESEQ analysis-----------------------
#run deseq2
design<- ~donor_num + cell_type
dds <- DESeqDataSetFromMatrix(countData = data_d, colData = metadata_df,design=design)

#filter for low read counts
smallestGroupSize <- 6
keep <- rowSums(counts(dds) >= 15) >= smallestGroupSize # filter on whether gene is expressed in at least 6 samples with a count of 15 or higher
dds <- dds[keep,]

#run & extract results
dds<- DESeq(dds)
results_d<- results(dds, contrast = c("cell_type", "iLC", "miLC"), alpha= adjpthres)

#apply logFC shrinkage
resLFC <- lfcShrink(dds, coef = "cell_type_miLC_vs_iLC", type="apeglm")

#perform initial annotation
results_df<- as.data.frame(resLFC)
results_df$gene_id <- rownames(results_df)

#annotate usoing ensembldb
ensemble_ids <- results_df$gene_id
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v86, keys= ensemble_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))


#merge step
results_df_annot<- merge(results_df, geneIDs1, by.x = "gene_id", by.y = "GENEID", all.x = TRUE)
results_df_annot_inf<- results_df_annot[results_df_annot$SYMBOL %in% inf_genes_names, ]

#export the full rersults
write_xlsx(results_df_annot, paste0(comp, "/excel_exports/unfiltered_annot_results.xlsx"))

#set external_gene_name column to gene_name for downstream functionality
colnames(results_df_annot)[7]<- "gene_name"

#inital adj-p val filtering
results_df_annot_sig<- select <- results_df_annot[(results_df_annot$padj < adjpthres & !is.na( results_df$padj ) ), ]

# log fc filtering
results_df_annot_sig_log <- results_df_annot_sig[abs(results_df_annot_sig$log2FoldChange)> log_FC,]


#create tables for inferome genes only
results_df_annot_inf_sig <- results_df_annot_sig  [results_df_annot_sig$gene_name %in% inf_genes_names, ]
results_df_annot_inf_sig_log <- results_df_annot_sig_log[results_df_annot_sig_log$gene_name %in% inf_genes_names, ]


#write all additional excels 
write_xlsx(results_df_annot_sig, paste0(comp, "/excel_exports/annot_", comp, "_sig.xlsx"))
write_xlsx(results_df_annot_sig_log, paste0(comp, "/excel_exports/annot_", comp, "_sig_log.xlsx"))
write_xlsx(results_df_annot_inf_sig, paste0(comp, "/excel_exports/annot_", comp, "_inf_sig.xlsx"))
write_xlsx(results_df_annot_inf_sig_log, paste0(comp, "/excel_exports/annot_", comp, "_inf_sig_log.xlsx"))

# DESEQ PLOTS-------------#
# barplot up/dowm
res_df_for_bar <- results_df_annot_sig_log[order(results_df_annot_sig_log$log2FoldChange, decreasing = TRUE), ]

#define the number of top genes 
top_n <- 10

#select top upregulated genes
top_upregulated <- head(res_df_for_bar, top_n)

#select top downregulated genes
top_downregulated <- tail(res_df_for_bar, top_n)

#combine the top upregulated and downregulated genes
top_genes <- rbind(top_upregulated, top_downregulated)

# Create a bar plot of top upregulated and downregulated genes
dir.create(paste0(comp, "/plots/up_down"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/up_down/", comp, "_up_down_DESEQ.pdf")))
nasty_test<-ggplot(top_genes, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", color = "NA") +
  coord_flip() +
  scale_fill_manual(values = c("lightblue", "darkred"), labels = c("Downregulated", "Upregulated")) +
  labs(title = "miLC vs iLC gene expression", x = "Gene", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(legend.title = element_blank())
plot(nasty_test)
dev.off()

# HEATMAPS
selected_data <- results_df_annot_sig_log #variable depending on what you want to plot
selected_genes <- selected_data$gene_name
selection_plot <- as.data.frame(data_import[data_import$gene_name %in% selected_genes, ])

# Filter selection_plot to include only the columns present in data_d
data_d_columns <- colnames(data_d)
selection_plot <- selection_plot %>%
  dplyr::select(gene_id, gene_name, all_of(data_d_columns))

gene_expr <- selection_plot %>%
  dplyr::select(-gene_id, -gene_name) %>%  # Remove gene_id and gene_name columns
  as.matrix()

metadata_df_plot<- metadata_df # optional, depends on the samples of interest

#add rownames and remove non-numeric columns
rownames(metadata_df_plot)<- metadata_df$read_name
metadata_df_plot$read_name=NULL
metadata_df_plot$donor_num=NULL
colnames(metadata_df_plot)<- c("Cell Type", "CD-status")

# remove the ones that have NAs after scaling to avoid heatmap error
scaled_gene_expr <- t(scale(t(gene_expr)))
rows_with_na <- apply(scaled_gene_expr, 1, function(x) any(is.na(x)))
summary(rows_with_na)
cleaned_gene_expr <- gene_expr[!rows_with_na, ]

#plot the heatmap
dir.create(paste0(comp, "/plots/heatmaps"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/heatmaps/", comp, "_heatmap.pdf")))
pheatmap(
  cleaned_gene_expr,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "",
  fontsize = 8,
  annotation_col = metadata_df_plot,  # Provide the vector of cell types directly
  color=inferno(50),
  annotation_names_col = FALSE,
  show_colnames = FALSE
)
dev.off()

#Vulcano plots

#----VULCANO1: Plot of all DEGs
selected_data <- results_df_annot_sig

dir.create(paste0(comp, "/plots/vulcano"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/vulcano/", comp, "_vulcano_all.pdf")),width = 30, height = 10)

plot(EnhancedVolcano(selected_data,
                lab = selected_data$gene_name,
                ylab= bquote(-log[10]~'Padj'),
                title="",
                subtitle = "",
                x="log2FoldChange",
                y = "padj",
                FCcutoff = log_FC,
                legendPosition = "none",
                labSize = 4,
                pCutoff = adjpthres,
))
dev.off()


#----VULCANO2: Plot of inf DEGs

selected_data2 <- results_df_annot_inf
head(selected_data2)
pdf(file.path(paste0(working_dir, "/", comp, "/plots/vulcano/", comp, "_vulcano_inf.pdf")),width = 30, height = 10)

plot(EnhancedVolcano(selected_data2,
                lab = selected_data2$SYMBOL,
                ylab= bquote(-log[10]~'Padj'),
                title="",
                subtitle = "",
                x="log2FoldChange",
                y = "padj",
                FCcutoff = log_FC,
                legendPosition = "none",
                labSize = 4,
                pCutoff = adjpthres,
))
dev.off()


#-----------Pathway analysis

#------------1: GSEA analysis

## load pathway halmarks
url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt"
local_file <- "h.all.v2023.2.Hs.symbols.gmt"
if (!file.exists(local_file)) {
  # Create the directory if it does not exist
  dir.create(dirname(local_file), recursive = TRUE, showWarnings = FALSE)
  
  # Download the file
  download.file(url, destfile = local_file, method = "auto")
  
  cat("File downloaded successfully.\n")
} else {
  cat("File already exists.\n")
}
#annotation for GSEA (no signficance filtering)
results_df_annot <- merge(results_df, data_import[, c("gene_id", "gene_name")], by = "gene_id", all.x = TRUE)


# set up signed ranking for GSEA. Upregulated genes will have positive rankings, downregulated genes will have negative rankings.
GSEA_rank<- sign(results_df_annot$log2FoldChange)*(-log10(results_df_annot$pvalue))
names(GSEA_rank)<-results_df_annot$gene_name

# we need to exclude duplicate genes
# ceate a data frame to handle duplicates
gsea_df <- data.frame(
  Gene = names(GSEA_rank),
  Rank = GSEA_rank,
  stringsAsFactors = FALSE
)

# remove duplicates by taking the maximal rank for duplicate genes
gsea_df_unique <- gsea_df %>%
  group_by(Gene) %>%
  summarise(Rank = max(Rank)) %>%
  ungroup()


# sort the genes based on the calculated ranks
gene_list_GSEA <- gsea_df_unique %>%
  arrange(desc(Rank))

# convert to a named vector and sort
gene_list_GSEA <- setNames(gene_list_GSEA$Rank, gene_list_GSEA$Gene)
gene_list_GSEA <- sort(GSEA_rank, decreasing = TRUE)

#load pathways
pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")

# rung fsgsea
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=gene_list_GSEA,nPermSimple=1000)


# tidy the results:
# Sort by Normalized enrichment scores
GSEA_table <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))


# Plot in nice table
GSEA_table <- GSEA_table %>%
  dplyr::mutate(
    pathway = str_replace_all(pathway, "HALLMARK_", ""), #remove HALLMARK_ prefix
    pathway = str_replace_all(pathway, "_", " ")        # replace underscore
  ) %>%
  dplyr::arrange(padj)


# write to excel
write_xlsx(GSEA_table, paste0(comp, "/excel_exports/", comp, "_GSEA.xlsx"))

#write to plot
dir.create(paste0(comp, "/plots/GSEA"))

pdf(file.path(paste0(working_dir, "/", comp, "/plots/GSEA/", comp, "_GSEA.pdf")),width = 10, height = 10)
plot(ggplot(GSEA_table, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="") + 
  theme_minimal())

dev.off()


#---------Gene Ontology
# save the enrich dotplot
dir.create(paste0(comp, "/plots/dotplot"))

#GO1: GO of all siglog filtered genes but not inf filtered
pdf(file.path(paste0(working_dir, "/", comp, "/plots/dotplot/", comp, "_dotplotGO_allsiglog.pdf")),width = 10, height = 10)
ego2 <- enrichGO(gene         = unique(results_df_annot_sig_log$gene_name),
                 OrgDb         = org.Hs.eg.db,
                 universe      = unique(results_df_annot$gene_name),
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = pthres,
                 qvalueCutoff  = adjpthres,
                 readable=T
)
dotplot(ego2)

dev.off()

#------GO2: GO of top hits only 
top_hits_number <- as.numeric(100)
# Sort by logFC to get top 100 upregulated and downregulated genes
top_upregulated <- head(results_df_annot_sig_log[order(results_df_annot_sig_log$log2FoldChange, decreasing = TRUE), ], top_hits_number)
top_downregulated <- head(results_df_annot_sig_log[order(results_df_annot_sig_log$log2FoldChange), ], top_hits_number)
head(top_upregulated)

# Combine the two sets of genes
combined_genes <- c(top_upregulated$gene_name, top_downregulated$gene_name)


pdf(file.path(paste0(working_dir, "/", comp, "/plots/dotplot/", comp, "_dotplotGO_top100hits.pdf")),width = 10, height = 10)
ego_combined <- enrichGO(
  gene         = unique(combined_genes),
  OrgDb        = org.Hs.eg.db,
  universe     = unique(results_df_annot$gene_name),
  keyType      = 'SYMBOL',
  ont          = "BP",  # Biological Process
  pAdjustMethod = "fdr",
  pvalueCutoff = pthres,
  qvalueCutoff = adjpthres,
  readable     = TRUE
)
dotplot(ego_combined)
dev.off()

#------GO2: GO of inf genes only 
pdf(file.path(paste0(working_dir, "/", comp, "/plots/dotplot/", comp, "_dotplot_inf.pdf")),width = 10, height = 10)
ego_inf <- enrichGO(gene         = unique(results_df_annot_inf_sig_log$gene_name),
                 OrgDb         = org.Hs.eg.db,
                 universe      = unique(results_df_annot$gene_name),
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = pthres,
                 qvalueCutoff  = adjpthres,
                 readable=T
)
dotplot(ego_inf)
dev.off()


# HERE STARTS THE WGCNA ANALYSIS!

#-------------WGCNA ANALYSIS----------------#
dir.create(paste0(comp, "/plots/WGCNA"))

wgcna_data<- data_d

#QC outliers (this is a bit redundant as similar processes are done in DESeq2 later on)
gsg<-goodSamplesGenes(t(wgcna_data))

table(gsg$goodGenes)
table(gsg$goodSamples)

#exclude outlier genes
wgcna_data<-wgcna_data[gsg$goodGenes==TRUE,]

wgcna_data2<-wgcna_data
colnames(wgcna_data2)<- metadata_df$cell_type
#detect outlier samples through hierarchical clustering
htree<-hclust(dist(t(wgcna_data2)), method="average")

dir.create(paste0(comp, "/plots/WGCNA/clusterplot"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/clusterplot/", comp, "_clusterplot.pdf")),width = 10, height = 10)
plot(htree, main="")
dev.off()


#detect outlier samples through pca
pca<-prcomp(t(wgcna_data))
pca.dat<-pca$x
pca.var<-pca$sdev^2
pca.var.percent<-round(pca.var/sum(pca.var)*100, digits=2)
ggplot(pca.dat,aes(PC1,PC2)) +
  geom_point()+
  geom_text(label=rownames(pca.dat))+
  labs(x= paste0('PC1:', pca.var.percent[1],'%'),
    y= paste0('PC2:',pca.var.percent[2],  '%'))

# in conclusion we can say that no samples need to be excluded


# normalisation through deseq2
rownames(metadata_df)<- metadata_df$read_name
all(rownames(metadata_df) %in% colnames(wgcna_data))# check
all(rownames(metadata_df) == colnames(wgcna_data))# check

#create dds2
dds2<- DESeqDataSetFromMatrix(countData=wgcna_data, colData=metadata_df, design = ~1) #dont specify metrix for WGCNA

dds75<-dds2[rowSums(counts(dds2) >=15)>= 8,] #filter for genes with minimum 15 counts in more than 75% of samples

#variance stabilisation
dds_norm <- vst(dds75)

#get normalized counts
norm_counts<- assay(dds_norm) %>% t()

# Network construction
#choose powers 
powers <- c(c(1:10), seq(from = 12, to = 50, by =2))

#call network topolgy analysis function
sft<-pickSoftThreshold(norm_counts,
                  powerVector = powers,
                  networkType = 'signed',
                  verbose = 5)

# visualisation to choose the right power.
sft.data<- sft$fitIndices
names(sft.data)

# we use R2 and mean connektivity for selecting power. We want a high R2 with minimum konnektivity. Let's visualize
r2<-ggplot(sft.data, aes(Power, SFT.R.sq, label= Power))+
  geom_point()+
  geom_text(nudge_y= 0.1)+
  geom_hline(yintercept = 0.8, color ='red')+
  labs(x='Power', y='Scale free topology model fit, signed R2')+
  theme_classic()

k2 <- ggplot(sft.data, aes(Power, mean.k., label=Power)) +
  geom_point() +
  geom_text(nudge_y= 0.1)+
  labs(x='Power', y='Mean Connectivity')+
  theme_minimal()

dir.create(paste0(comp, "/plots/WGCNA/R2K2"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/R2K2/", comp, "_r2_kon.pdf")),width = 10, height = 10)
grid.arrange(r2, k2, nrow=2)
dev.off()
#10 seems suitable

#convert matrix to numeric
norm_counts[] <- sapply(norm_counts, as.numeric)

soft_power <-10
temp_cor<-cor # set this so we can set it back to normal cor function later on
cor <-WGCNA::cor #force that this opne is used

#memory estimated
bwnet <- blockwiseModules(norm_counts,
                maxBlockSize = 14000,
                minModuleSize = 30,
                TOMtype = 'signed',
                power = soft_power,
                mergeCutHeight = 0.25,
                numericLabels = FALSE,
                randomSeed= 1234, 
                verbose =3)

cor<-temp_cor #reassign to standard cor function

# extract modules Eigengenes and colors
module_eigengenes <- bwnet$MEs
moduleColors <- bwnet$colors 

# no of genes per module
table(bwnet$colors) 

#Dendogram plot
dir.create(paste0(comp, "/plots/WGCNA/dendogram"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/dendogram/", comp, "_dendogram.pdf")),width = 10, height = 10)

plotDendroAndColors(bwnet$dendrograms[[1]],
               cbind(bwnet$unmergedColors, bwnet$colors),c("unmerged","merged"),
               dendroLabels=FALSE,
               addGuide = TRUE,
               hang= 0.03,
               guideHang = 0.05)
dev.off()



#-------- Relate modules to traits

#creat traits file for cell_type
traits <- metadata_df %>%
  dplyr::select(cell_type)

# split the columsn for clearer plot into a design matrix
cell_type_matrix <- model.matrix(~ cell_type + 0, data = metadata_df)
traits <- as.data.frame(cell_type_matrix)
head(traits)
# rename columns
traits <- traits %>%
  dplyr::rename(iLC = cell_typeiLC, miLC = cell_typemiLC)

#define no of genes and samples5
nSamples <- nrow(norm_counts)
nGenes <- ncol(norm_counts)
names(traits)

module.trait.corr <-cor(module_eigengenes, traits, use ='p')

module.trait.corr.Pvalue <-corPvalueStudent(module.trait.corr,nSamples)


#visualize association
heatmap.data <-merge(module_eigengenes, traits, by='row.names')
head(heatmap.data)

#move rownames to columns
heatmap.data<- heatmap.data %>%
  column_to_rownames(var='Row.names')

names(heatmap.data)
#correlation plot
dir.create(paste0(comp, "/plots/WGCNA/corplot"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/corplot/", comp, "_corplot.pdf")),width = 10, height = 10)
CorLevelPlot(heatmap.data,
              x = names(heatmap.data)[17:18],
              y = names(heatmap.data)[1:16],
              col=c('#5db6da',"white", "#ca2626"))
dev.off()

# LETS START ANALYZING
#creat binary traits file for cell_type
traits2 <- metadata_df %>%
  mutate(cell_type_bin = ifelse(grepl('miLC', cell_type),0,1)) %>% dplyr::select(5)

#measure scores (MM)
module.membership.measures <-cor(module_eigengenes, norm_counts, use='p')# scores
module.membership.measures.pvals<-corPvalueStudent(module.membership.measures, nSamples) #pvalues
mm_scores <-t(module.membership.measures.pvals) #transpose

#Gene significance scores (GS)
gene.signf.corr<-cor(norm_counts, traits2$cell_type, use='p') #GS
colnames(gene.signf.corr)[1] <- "GS"

# Gene significance score pvals
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)#Pvals
colnames(gene.signf.corr.pvals)[1] <- "Pval"

# Set up the full dataframe
gs_scores<- data.frame(
  gene_id = rownames(gene.signf.corr),
  row.names = rownames(gene.signf.corr),
  gs_score = as.vector(gene.signf.corr),
  pval = as.vector(gene.signf.corr.pvals)
  )


#function to perform filtering by color and obtain module membership
filter_module_by_color <- function(color) {

  full_color <- paste0("ME", color)

  selected_genes<-as.data.frame(bwnet$colors) %>%
  dplyr::filter(bwnet$colors == color) %>%
  rownames()

  filtered_mm_scores <- module.membership.measures[full_color, selected_genes]
  filtered_mm_pvals <- module.membership.measures.pvals[full_color, selected_genes ]
  
  filtered_data <- data.frame(
    gene_id = names(filtered_mm_scores),
    mm_score = filtered_mm_scores,
    mm_pvals = filtered_mm_pvals,
    row.names = names(filtered_mm_scores)
  )
  return(filtered_data)
}

#-----------------MM GS Scatterplots
MM_GS_scatter <- function(color) {
  mod_col <- filter_module_by_color(color)
  gs_col <- gs_scores[row.names(mod_col),2]
  
  plot<- verboseScatterplot(abs(mod_col[,2]),
                   abs(gs_col), 
                   xlab = "Module Membership",
                   ylab = "Gene significance",
                   main = "",
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = color)

return(plot)

}


# MM AND GS Filtering for hub selection!

hub_filtering <- function(gs_scores, color) {

  mm_filtered_frame <- filter_module_by_color(color)

  #filter GS - This function filters GS >0.5 and MM > 0.8 and returns the intersect
  gs_scores_0.5 <- gs_scores %>%
    dplyr::filter(pval < 0.05) %>%
    dplyr::filter(abs(gs_score) > 0.5)

  #filter MM
  mm_filtered_0.8 <- mm_filtered_frame %>%
    dplyr::filter(mm_pvals < 0.05) %>%
    dplyr::filter(abs(mm_score) > 0.8)

  #find overlap
  intersection <- gs_scores_0.5 %>%
  inner_join(mm_filtered_0.8, by = "gene_id")

  
  return(intersection)
}

# now find the overlap with deseq2 data

mm_gs_deseq_intersector <- function(mm_gs_intersect, deseq2data) {
filtered_turq_deseq2 <- merge(mm_gs_intersect, deseq2data, by = "gene_id")
filtered_turq_deseq2 <- filtered_turq_deseq2 %>% arrange(log2FoldChange)
return(filtered_turq_deseq2)
}



#venn intersector
venn_intersector <- function(mm_gs_intersect, deseq2data, plotcolor) {
 
   # Create a named list of gene sets for Venn diagram
  gene_sets <- list(
    `Hub Genes` = mm_gs_intersect$gene_id,
    `DEGs` = deseq2data$gene_id,
    `ISGs` = inf_genes_gene_id
  )
  
  # Create a Venn diagram using ggvenn
   venn_plot <- ggvenn(
    data = gene_sets,
    show_percentage = TRUE,
    fill_color = c(plotcolor, "#E7B800", "#FF474C"),
    set_name_color = c(plotcolor, "#E7B800", "#FF474C"),
    stroke_size = 0.5,  # Thinner border around sections
    set_name_size = 6,  # Increase the set label size
    text_size = 5       # Increase the size of the percentages
    ) 
  
  # Return the ggvenn plot
  return(venn_plot)
}

#--------for turquoise
#scatter
dir.create(paste0(comp, "/plots/WGCNA/scatter_MM_GS"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/scatter_MM_GS/", comp, "_scatter_turq.pdf")),width = 10, height = 10)
scat_turq <- MM_GS_scatter("turquoise") #scatter of gs-mm
dev.off()

# hub filtered
hub_turquoise <-hub_filtering(gs_scores, "turquoise" ) # hub filtering
write_xlsx(hub_turquoise, paste0(comp, "/excel_exports/hub_turquoise", comp, "_sig.xlsx"))

# hub deseq intersect
turq_intersect_MM_GS_deseq <- mm_gs_deseq_intersector(hub_turquoise,results_df_annot_sig_log ) # intersect of hubfilter with deseq results
write_xlsx(turq_intersect_MM_GS_deseq, paste0(comp, "/excel_exports/hub_deseq_turquoise", comp, "_sig.xlsx"))

# hub deseq isg intersect
intersect_inf_turq <- turq_intersect_MM_GS_deseq[turq_intersect_MM_GS_deseq$gene_name %in% inf_genes_names, ] # intersect of hubfilter with deseq results filtered for ISGs
write_xlsx(intersect_inf_turq, paste0(comp, "/excel_exports/hub_deseq_ISG_turquoise", comp, "_sig.xlsx"))

# venn plot
turq_intersect_venn<- venn_intersector(hub_turquoise,results_df_annot_sig_log, "#4ed8ca") #venn diagram of hubfilter with deseq with inf
dir.create(paste0(comp, "/plots/WGCNA/venn"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/venn/", comp, "_venn_turquoise.pdf")),width = 10, height = 10)
turq_intersect_venn
dev.off()

#--------for blue
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/scatter_MM_GS/", comp, "_scatter_blue.pdf")),width = 10, height = 10)
scat_blue <- MM_GS_scatter("blue") #scatter of gs-mm
dev.off()
# hub filtered
hub_blue <-hub_filtering(gs_scores, "blue" ) # hub filtering
write_xlsx(hub_blue, paste0(comp, "/excel_exports/hub_blue", comp, "_sig.xlsx"))

# hub deseq intersect
blue_intersect_MM_GS_deseq <- mm_gs_deseq_intersector(hub_blue,results_df_annot_sig_log ) # intersect of hubfilter with deseq results
write_xlsx(blue_intersect_MM_GS_deseq, paste0(comp, "/excel_exports/hub_deseq_blue", comp, "_sig.xlsx"))

# hub deseq isg intersect
intersect_inf_blue<- blue_intersect_MM_GS_deseq[blue_intersect_MM_GS_deseq$gene_name %in% inf_genes_names, ] # intersect of hubfilter with deseq results filtered for ISGs
write_xlsx(intersect_inf_blue, paste0(comp, "/excel_exports/hub_deseq_ISG_blue", comp, "_sig.xlsx"))

# venn plot
blue_intersect_venn<-venn_intersector(hub_blue,results_df_annot_sig_log, "blue") #venn diagram of hubfilter with deseq with inf
dir.create(paste0(comp, "/plots/WGCNA/venn"))
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/venn/", comp, "_venn_blue.pdf")),width = 10, height = 10)
blue_intersect_venn
dev.off()

# ------ UP DOWN------#
dir.create(paste0(comp, "/plots/WGCNA/up_down"))

#first we create the figures for up down plots of each module

#turquoise
top_up_turq<- tail(turq_intersect_MM_GS_deseq,n=10)
top_up_turq["1994", "gene_name"] <- "ZNF775-AS1" #manually look up ensg for NA value
top_down_turq<-head(turq_intersect_MM_GS_deseq,n=10)
combined_turq <- rbind(top_up_turq, top_down_turq)

pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/up_down/", comp, "_up_down_module_turq.pdf")))
turq_up_down_plot <- ggplot(combined_turq, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", color = NA) +
  coord_flip() +
  scale_fill_manual(values = c("#40E0D0", "#247e75"), labels = c("Downregulated", "Upregulated")) +
  labs(title = "", x = "Gene", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(legend.title = element_blank())
plot(turq_up_down_plot)
dev.off()

#blue
top_up_blue<- tail(blue_intersect_MM_GS_deseq,n=10)
top_up_blue["1674", "gene_name"] <- "CRLF2" #manually look up ensg for NA value
top_down_blue<-head(blue_intersect_MM_GS_deseq,n=10)
combined_blue <- rbind(top_up_blue, top_down_blue)

pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/up_down/", comp, "_up_down_module_blue.pdf")))

blue_up_down_plot <- ggplot(combined_blue, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", color = NA) +
  coord_flip() +
  scale_fill_manual(values = c("lightblue", "darkblue"), labels = c("Downregulated", "Upregulated")) +
  labs(title = "", x = "Gene", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(legend.title = element_blank())
plot(blue_up_down_plot)
dev.off()

# next: up down for the mod-deseq-inf intersects

#turquoise
top_up_turq2<- tail(intersect_inf_turq,n=10)
top_down_turq2<-head(intersect_inf_turq,n=10)
combined_turq2 <- rbind(top_up_turq2, top_down_turq2)

pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/up_down/", comp, "_up_down_modseqisg_turq.pdf")))
turq_up_down_modseqinf_plot <- ggplot(combined_turq2, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", color = NA) +
  coord_flip() +
  scale_fill_manual(values = c("#40E0D0", "#247e75"), labels = c("Downregulated", "Upregulated")) +
  labs(title = "", x = "Gene", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(legend.title = element_blank())
plot(turq_up_down_modseqinf_plot)
dev.off() 


#blue
top_up_blue2<- tail(intersect_inf_blue,n=10)
top_down_blue2<-head(intersect_inf_blue,n=10)
combined_blue2 <- rbind(top_up_blue2, top_down_blue2)

pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/up_down/", comp, "_up_down_modsequbf_blue.pdf")))
blue_up_down_modseqinf_plot <- ggplot(combined_blue2, aes(x = reorder(gene_name, log2FoldChange), y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", color = NA) +
  coord_flip() +
  scale_fill_manual(values = c("lightblue", "darkblue"), labels = c("Downregulated", "Upregulated")) +
  labs(title = "", x = "Gene", y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "right") +
  theme(legend.title = element_blank())
plot(blue_up_down_modseqinf_plot)
dev.off()

#-------GENE ONTOLOGY
#-------------plots----------------------
dir.create(paste0(comp, "/plots/WGCNA/GeneOntology"))

# PLOT 1 &2 : Module specific
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/GeneOntology/", comp, "_GO_turq_module.pdf")),width = 10, height = 10)
GO_turq_module<-enrichGO(
    gene = unique(hub_turquoise$gene_id),
    OrgDb = org.Hs.eg.db,
    keyType = 'ENSEMBL',
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = pthres,
    qvalueCutoff = adjpthres,
    readable = TRUE)

dotplot(GO_turq_module)
dev.off()

pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/GeneOntology/", comp, "_GO_blue_module.pdf")),width = 10, height = 10)

# now for blue
GO_blue_module<-enrichGO(
    gene = unique(hub_blue$gene_id),
    OrgDb = org.Hs.eg.db,
    keyType = 'ENSEMBL',
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = pthres,
    qvalueCutoff = adjpthres,
    readable = TRUE)

dotplot(GO_blue_module)
dev.off()



# Plot 3 &4 intrsect of module and deseq
# GO enrichment plot
  #turquoise
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/GeneOntology/", comp, "_GO_turq_GS_MM_DESEQ.pdf")),width = 10, height = 10)

GO_turq_deseq<-enrichGO(
    gene = unique(turq_intersect_MM_GS_deseq$gene_id),
    OrgDb = org.Hs.eg.db,
    keyType = 'ENSEMBL',
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = pthres,
    qvalueCutoff = adjpthres,
    readable = TRUE)
dotplot(GO_turq_deseq)
dev.off()


  #blue
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/GeneOntology/", comp, "_GO_blue_GS_MM_DESEQ.pdf")),width = 10, height = 10)
GO_blue_deseq<-enrichGO(
    gene = unique(blue_intersect_MM_GS_deseq$gene_id),
    OrgDb = org.Hs.eg.db,
    keyType = 'ENSEMBL',
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = pthres,
    qvalueCutoff = adjpthres,
    readable = TRUE)
dotplot(GO_blue_deseq)
dev.off()

# PLOT 5 & 6 : with filter for ISGs
pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/GeneOntology/", comp, "_GO_turq_GS_MM_DESEQ_INF.pdf")),width = 10, height = 10)
GO_turq_inf<-enrichGO(
    gene = unique(intersect_inf_turq$gene_id),
    OrgDb = org.Hs.eg.db,
    keyType = 'ENSEMBL',
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = pthres,
    qvalueCutoff = adjpthres,
    readable = TRUE)
dotplot(GO_turq_inf)
dev.off()

pdf(file.path(paste0(working_dir, "/", comp, "/plots/WGCNA/GeneOntology/", comp, "_GO_blue_GS_MM_DESEQ_INF.pdf")),width = 10, height = 10)
GO_blue_inf<-enrichGO(
    gene = unique(intersect_inf_blue$gene_id),
    OrgDb = org.Hs.eg.db,
    keyType = 'ENSEMBL',
    ont = "BP",
    pAdjustMethod = "fdr",
    pvalueCutoff = pthres,
    qvalueCutoff = adjpthres,
    readable = TRUE)
dotplot(GO_blue_inf)
dev.off()


 
sessionInfo()
