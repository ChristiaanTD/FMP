#!/C:/Program Files/R/R-4.3.3/bin/R.exe

#--------1. Setup--------#
#working_dir<- "L:/basic/divg/KVI/HIV_RNA_seq/Thijs/set1"
working_dir<-'C:/Users/tcdor/OneDrive/Desktop/Set1Scripts'

# set up libraries
library(dplyr)
library(stringr)
library(ggplot2)
library(here)
library(readxl)
library(writexl)
library(mixOmics)
library(tidyverse)
library(limma)
library(grid)
library(tidyr)
library(viridis)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GO.db)
library(DOSE)
library(fgsea)
library(annotate)
library(pheatmap)
library(seriation)
library(dendextend)
library(enrichplot)
library(org.Hs.eg.db)
library(DT)
library(optparse)
library(biomaRt)
setwd(working_dir)

# -----------optional parsing settings-------------
#parser <- OptionParser()

# define options for cond1 and cond2
#parser <- add_option(parser, "--cond1", dest="cond1", type="character", help="Condition 1")
#parser <- add_option(parser, "--cond2", dest="cond2", type="character", help="Condition 2")

# parse the command-line arguments
#options <- parse_args(parser)

# Access the values of cond1 and cond2
#cond1 <- options$cond1
#cond2 <- options$cond2



#--------1. Setting Conditions--------#

# choose conditions of interest
cond1 <- "MDM_IL10_Mock" 
cond2 <- "MDM_IL10_R5"

#conditions available to use : "MDM_IL4_R5","MDM_IFNgTNFa_R5","MDM_IL10_R5", compare them against their respective mock
#Alternatively one can compare "MDM_IL4_Mock",'MDM_IFNgTNFa_Mock','MDM_IL10_Mock' and compare them against 'MDM_Mock' if you want to compare in absence of HIV-1
# choose significance threshold & FC thres
pthres<- 0.05
adjpthres<-0.05
log_FC<- 1

# set up directories
inf_extract_dir <- "set1_analysis/Interferome_genes.csv" # set directory for list of interferome genes
normalized_data_dir <- "set1_raw/normalised_data_ThijsDormans.xlsx" # directory of normalized data xsls

#--------2. Data Load--------#

# load interferome list
inf_genes_table <- read.csv(inf_extract_dir, sep = ";", header=T)
inf_genes_names <- unique(inf_genes_table$Gene.Name) # create gene names vector

# load expression data
data_import <- read_excel(normalized_data_dir, sheet = "normalised_data.txt")

# set the rows and remove redundant columns
data<- as.data.frame(data_import[, c(2,4, 7:ncol(data_import))]
)

# change + to _ to avoid LIMMA error
colnames(data) <- gsub("\\+", "_", colnames(data))

# set rownames to probeID
rownames(data)<- data$ProbeID

#-------------------FOR LOOP FOR ANALYIS STARTS HERE-----------------------#
#for (cond2 in selected_con) {
  # Print the current conditions being processed
  #cat("Running analysis for", cond1, "vs", cond2, "\n")


# create directory name for result storing
dir_name <- paste0("set1_analysis/",cond1, "_vs_", cond2)

# Create the directory if it doesn't exist
if (!dir.exists(dir_name)) {
  dir.create(dir_name)
}


#--------3. data prep--------#

# create a df with probeID as rownames and only normalized expression data
prep_norm_gene<- data[, 3:ncol(data)]

#transpose data for PCA analysis format
t_norm_gene <- t(prep_norm_gene)

# set conditions vector
conditions <- c(
  rep("MDM_Mock", 3),
  rep("MDM_R5", 3),
  rep("MDM_X4", 3),
  rep("MDM_IL4_Mock", 3),
  rep("MDM_IL4_R5", 3),
  rep("MDM_IL4_X4", 3),
  rep("MDM_IFNgTNFa_Mock", 3),
  rep("MDM_IFNgTNFa_R5", 3),
  rep("MDM_IL10_Mock", 3),
  rep("MDM_IL10_R5", 3)
)
# Convert conditions to a factor with desired levels and order
conditions <- factor(conditions, levels = c("MDM_Mock", "MDM_R5", "MDM_X4", 
                                            "MDM_IL4_Mock", "MDM_IL4_R5", "MDM_IL4_X4", 
                                            "MDM_IFNgTNFa_Mock", "MDM_IFNgTNFa_R5", 
                                            "MDM_IL10_Mock", "MDM_IL10_R5"))

# we also set up Condtions (with captial C) for  more aesthetic plotting later on
Conditions <- c(
  rep("Mock", 3),
  rep("Mock + R5", 3),
  rep("Mock + X4", 3),
  rep("IL4 + Mock", 3),
  rep("IL4 + R5", 3),
  rep("IL4 + X4", 3),
  rep("IFNgTNFa + Mock", 3),
  rep("IFNgTNFa + R5", 3),
  rep("IL10 + Mock", 3),
  rep("IL10 + R5", 3)
)


# quick check if conditions vector match with actual column name
print(row.names(t_norm_gene))
print(conditions)


#---OPTIONAL-----4. PCA exploration---OPTIONAL----#
# for all genes
pca_result_all <- pca(t_norm_gene, scale = TRUE)

# plot for all samples
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_PCA.pdf")), width = 10, height = 10)
plotIndiv(pca_result_all, comp = c(1, 2), ind.names = T,
          group= (Conditions),
          legend = TRUE, title = '')
dev.off()


#--------5. LIMMA analysis--------#

# set up design matrix 
design <- model.matrix (~0+conditions)

# set up colnames
colnames(design) <- gsub("conditions", "", colnames(design))

# fit linear model
fit<-lmFit(prep_norm_gene, design)


# set up contrast matrix
contrast.matrix <- makeContrasts(
  paste(cond2, "-", cond1),  # contrast between cond2 and cond1
  levels = design
)

# apply the contrast fit 
fit2 <- contrasts.fit(fit, contrast.matrix)

# apply bayesian fitting
efit <- eBayes(fit2, trend=TRUE)


#--------6. DEGs--------#
top_table <- topTable(efit, coef = 1, number = Inf, adjust = "fdr") #extract all the values by setting number to inf

# add rownames and addd gene names and NM_ID
top_table$ProbeID <- rownames(top_table)
top_table_annot <- merge(top_table, data[, c("ProbeID", "ILMN_GENE")], by = "ProbeID", all.x = TRUE)
top_table_annot <- top_table_annot[, c("ILMN_GENE", names(top_table_annot)[!names(top_table_annot) %in% "ILMN_GENE"])] #add gene name
full_top_table <- merge(top_table_annot, data_import[, c("ProbeID", "SEARCH_KEY")], by = "ProbeID", all.x = TRUE)# add NM_ID


# adjusted p-value filtering, logFC filtering and sort
results_adj_p <- full_top_table[full_top_table$adj.P.Val<adjpthres,] #pval
results_adj_p <- results_adj_p[order(results_adj_p$logFC),] #logfc
results_adj_p <- results_adj_p[order(results_adj_p$logFC),] #logfc
results_adj_p_log <- results_adj_p[results_adj_p$logFC <log_FC,] #logfc



#functional annotation using biomaRt
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                     filters = "hgnc_symbol", 
                     values = results_adj_p$ILMN_GENE, 
                     mart = ensembl)

# remove duplicate annotations
annotations<- annotations[-which(duplicated(annotations$external_gene_name)),]

# update with newly generated columns
results_adj_p <- merge(results_adj_p, annotations, by.x = "ILMN_GENE", by.y = "external_gene_name", all.x = TRUE)

# extract only interferome genes
top_table_inf <- full_top_table[full_top_table$ILMN_GENE %in% inf_genes_names, ]

# extract only significant interferome genes
sig_top_table_inf <- results_adj_p[results_adj_p$ILMN_GENE %in% inf_genes_names, ]
sig_top_table_inf<- sig_top_table_inf[order(sig_top_table_inf$logFC),]

# sig but lofFC filtered
sigLOG_top_table_inf <- sig_top_table_inf[abs(sig_top_table_inf$logFC) >= log_FC, ]

# export both full results, significant results and interferome results to table
write_xlsx(full_top_table, file.path(dir_name,paste0(cond1, "_vs_", cond2, "_all.xlsx")))
write_xlsx(results_adj_p, file.path(dir_name, paste0(cond1, "_vs_", cond2, "_sig_all.xlsx")))
write_xlsx(results_adj_p_log, file.path(dir_name,paste0(cond1, "_vs_", cond2, "_sig_log_all.xlsx")))
write_xlsx(top_table_inf, file.path(dir_name, paste0(cond1, "_vs_", cond2, "_inf.xlsx")))
write_xlsx(sig_top_table_inf, file.path(dir_name, paste0(cond1, "_vs_", cond2, "_sig_inf.xlsx")))# for significant interferome genes
write_xlsx(sigLOG_top_table_inf, file.path(dir_name, paste0(cond1, "_vs_", cond2, "_sigLog_inf.xlsx")))# for significant interferome genes


#--------6. Plots--------#
#---- strongest DEGs plot for inf
top_n <- 20
deg_plot_strongest_inf <- sig_top_table_inf
deg_plot_strongest_inf$logFCabs <- abs(sig_top_table_inf$logFC)
deg_plot_strongest_inf <- deg_plot_strongest_inf[order(-deg_plot_strongest_inf$logFCabs), ]
top_strongest_deg_inf <- head(deg_plot_strongest_inf, top_n)


# ientify duplicates in the strongest DEGs
top_strongest_deg_inf$ILMN_GENE_DUPLICATE <- duplicated(top_strongest_deg_inf$ILMN_GENE) |
                                         duplicated(top_strongest_deg_inf$ILMN_GENE, fromLast = TRUE)

# set the gene label value based on presence in duplicates
top_strongest_deg_inf$GeneLabel <- ifelse(top_strongest_deg_inf$ILMN_GENE_DUPLICATE,
                                      paste0(top_strongest_deg_inf$ILMN_GENE, " (", top_strongest_deg_inf$ProbeID, ")"),
                                      top_strongest_deg_inf$ILMN_GENE)

#up down inf degs
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_infDEG_strongest.pdf")), width = 10, height = 10)
up_down_plot_strongest_inf <- ggplot(top_strongest_deg_inf, aes(x=reorder(GeneLabel, logFC), y=logFC, fill=logFC > 0)) +
  geom_bar(stat="identity", color='NA') +
  coord_flip() +
  scale_fill_manual(values=c("lightblue", "darkred"), labels=c("Downregulated", "Upregulated")) +
  labs(title="", x="Gene (Probe ID)", y="Log2 Fold Change") +
  theme_minimal() +
  theme(
    legend.position = 'right',
    legend.title = element_blank(),
    axis.text.x = element_text(size=16),        
    axis.text.y = element_text(size=16),        
    axis.title.x = element_text(size=16),       
    axis.title.y = element_text(size=16),       
    legend.text = element_text(size=18),         
  )
plot(up_down_plot_strongest_inf)
dev.off()


#----strongest all degs
all_deg_strongest <- results_adj_p
all_deg_strongest$logFCabs <- abs(results_adj_p$logFC)
all_deg_strongest <- all_deg_strongest[order(-all_deg_strongest$logFCabs), ]
top_strongest_deg_all <- head(all_deg_strongest, top_n)

# identify duplicates
top_strongest_deg_all$ILMN_GENE_DUPLICATE <- duplicated(top_strongest_deg_all$ILMN_GENE) | 
                                       duplicated(top_strongest_deg_all$ILMN_GENE, fromLast = TRUE)

#set the gene label value based on presence in duplicate
top_strongest_deg_all$GeneLabel <- ifelse(top_strongest_deg_all$ILMN_GENE_DUPLICATE,
                                    paste0(top_strongest_deg_all$ILMN_GENE, " (", top_strongest_deg_all$ProbeID, ")"),
                                    top_strongest_deg_all$ILMN_GENE)
# Create the plot
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_top_", top_n, "_up_down_all.pdf")), width = 10, height = 10)
up_down_plot_all <- ggplot(top_strongest_deg_all, aes(x=reorder(GeneLabel, logFC), y=logFC, fill=logFC > 0)) +
  geom_bar(stat="identity", color='NA') +
  coord_flip() +
  scale_fill_manual(values=c("lightblue", "darkred"), labels=c("Downregulated", "Upregulated")) +
  labs(title="", x="Gene", y="Log2 Fold Change") +
  theme_minimal() +
   theme(
    legend.position = 'right',
    legend.title = element_blank(),
    axis.text.x = element_text(size=16),        
    axis.text.y = element_text(size=16),        
    axis.title.x = element_text(size=16),       
    axis.title.y = element_text(size=16),       
    legend.text = element_text(size=16),         
  )
plot(up_down_plot_all)
dev.off()




#--------HEATMAPS-----------------

#----HEATMAPS SETUP-----

# set up column groups for col annotation of all conditions
column_groups_all <- data.frame(Conditions)
rownames(column_groups_all) <- colnames(data[,3:ncol(data)]) #ensure matching

# Remove the prefix from cond1
heat_cond1 <- sub("_.*$", "", gsub("^MDM_", "", cond1))

# Remove the prefix from cond2
heat_cond2 <- gsub("_", " + ", gsub("^MDM_", "", cond2))

# #set up clustering callback option to set logical order of the cluster presentation
callback = function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
}



#----------HEATMAP 1:  General heatmap of ALL DEGs for selected condition
selected_probeIDs1 <- results_adj_p$ProbeID

# set up data frame for only significant interferome probeIDs
deg_filtered_data <- data[data$ProbeID %in% selected_probeIDs1, ]

# select the sample columns corresponding to the variables
selected_cols <- grep(paste0("^", cond1, "_|^", cond2, "_"), colnames(deg_filtered_data))

# set up dataframe for adding the condition cluster coloring
heatmap1_deg_data<- deg_filtered_data[, selected_cols]

# set up column groups for col annotation of specific conditions
column_groups <- data.frame(
  "Condition" = factor(rep(c(heat_cond1, heat_cond2), each = 3), levels = c(heat_cond1, heat_cond2)))
rownames(column_groups) <- colnames(heatmap1_deg_data) #ensure matching


pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_allDEG_heat_comp.pdf")), width = 10, height = 10)
pheatmap(heatmap1_deg_data,
         scale = "row",             
         clustering_method = "ward.D2", 
         color = inferno(50), 
         border_color = "NA",
         labels_col = NULL,
         main = "",
         fontsize=15,
         clustering_callback = callback,
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_col = column_groups,
         annotation_names_col = FALSE
         
)
dev.off()

#----HEATMAP 2: All DEGs between the two selected condtions but shown for all conditions
heatmap2_deg_data<- deg_filtered_data[,3:ncol(deg_filtered_data)]

pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_allDEG_heat_all.pdf")), width = 10, height = 10)
pheatmap(heatmap2_deg_data,
         scale = "row",             
         clustering_method = "ward.D2", 
         color = inferno(50), 
         border_color = "NA",
         labels_col = NULL,
         main = "",
         fontsize=15,
         clustering_callback = callback,
         labels_row = deg_filtered_data$ILMN_GENE,
         annotation_col = column_groups_all,
         annotation_names_col = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE
)
dev.off()


#------ HEATMAP 3: Heatmap of Inf DEGs for selected conditions

# select probeIDs of DEGs
selected_probeIDs3 <- sig_top_table_inf$ProbeID

# set up data frame for only significant interferome probeIDs
inf_filtered_data <- data[data$ProbeID %in% selected_probeIDs3, ]

# select the sample columns corresponding to the variables
selected_cols <- grep(paste0("^", cond1, "_|^", cond2, "_"), colnames(inf_filtered_data))

# set up dataframe for adding the condition cluster coloring
heatmap_inf_data<- inf_filtered_data[, selected_cols]

# address duplicates
inf_filtered_data$ILMN_GENE_DUPLICATE <- duplicated(inf_filtered_data$ILMN_GENE) |
                                         duplicated(inf_filtered_data$ILMN_GENE, fromLast = TRUE)

# set the gene label value based on presence in duplicates
inf_filtered_data$GeneLabel <- ifelse(inf_filtered_data$ILMN_GENE_DUPLICATE,
                                      paste0(inf_filtered_data$ILMN_GENE, " (", inf_filtered_data$ProbeID, ")"),
                                      inf_filtered_data$ILMN_GENE)

# actual heatmap plot
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_infDEG_heat_comp.pdf")), width = 10, height = 10)

pheatmap(heatmap_inf_data,
         scale = "row",             
         clustering_method = "ward.D2", 
         color = inferno(50), 
         border_color = "NA",
         labels_col = NULL,
         main = "",
         fontsize=15,
         clustering_callback = callback,
         labels_row = inf_filtered_data$GeneLabel,
         annotation_col = column_groups,
         annotation_names_col = FALSE,
         show_colnames = FALSE,
         show_rownames = TRUE
)
dev.off()

#------ HEATMAP 4: Heatmap of Inf DEGs for all conditions
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_infDEG_heat_all.pdf")), width = 10, height = 10)
pheatmap(inf_filtered_data[,3:32],
         scale = "row",             
         clustering_method = "ward.D2", 
         color = inferno(50),
         border_color = NA,
         main = "",
         fontsize=15,
         labels_row = inf_filtered_data$GeneLabel,
         annotation_col = column_groups_all,
         annotation_names_col = FALSE,
         show_colnames=FALSE,
         show_rownames = TRUE
)
dev.off()

#remove the duplicate columns from inf_filtered_data for cleanliness

#------VULCANO PLOTS
# Vulcano of all DEGs
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_all_vulc.pdf")), width = 10, height = 10)
plot(EnhancedVolcano(full_top_table,
                title=NULL,
                subtitle=NULL,
                lab= full_top_table$ILMN_GENE,
                x="logFC",
                y = "adj.P.Val",
                FCcutoff = log_FC,
                labSize = 6,
                legendLabels = c('ns', expression(Log[2]~FC), 'adj. p-value', expression('adj. p-value & Log'[2]~FC)),
                ylab = expression(log[10]~adjusted~italic(P)),
                pCutoff = pthres))

dev.off()

# Vulcano of INF DEGs
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_infDEG_vulc.pdf")), width = 10, height = 10)
plot(EnhancedVolcano(top_table_inf,
                title=NULL,
                subtitle=NULL,
                lab= top_table_inf$ILMN_GENE,
                x="logFC",
                y = "adj.P.Val",
                FCcutoff = log_FC,
                labSize = 6,
                legendLabels = c('ns', expression(Log[2]~FC), 'adj. p-value', expression('adj. p-value & Log'[2]~FC)),
                ylab = expression(log[10]~adjusted~italic(P)),
                pCutoff = pthres))

dev.off()

# -----------FUNCTIONAL ANNOTATION & ENRICHMENT--------------#
#---------GSEA---------$
# set up signed ranking for GSEA. Upregulated genes will have positive rankings, downregulated genes will have negative rankings.
GSEA_rank<- sign(full_top_table$logFC)*(-log10(full_top_table$P.Value))
names(GSEA_rank)<-full_top_table$ILMN_GENE


# we need to exclude duplicate genes
# ceate a data frame to handle duplicates
gsea_df <- data.frame(
  Gene = names(GSEA_rank),
  Rank = GSEA_rank,
  stringsAsFactors = FALSE
)

# Count the number of occurrences of each gene
gene_counts <- gsea_df %>%
  count(Gene)

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



# load pathway hallmarks
pathways.hallmark <- gmtPathways("h.all.v2023.2.Hs.symbols.gmt")

# perform GSEA
fgseaRes <- fgsea(pathways=pathways.hallmark, stats=gene_list_GSEA)


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



# plot GSEA
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_GSEA_plot.pdf")), width = 10, height = 10)
ggplot(GSEA_table, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score")
       title= "" + 
  theme_minimal()
dev.off()



# write gsea results to excel
write_xlsx(GSEA_table, file.path(dir_name,paste0(cond1, "_vs_", cond2, "_all_GSEA.xlsx")))



# plot
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_sig_GO_PLOT_all.pdf")), width = 10, height = 15)
#------ GO enrichment---------#
GO_all <- enrichGO(gene        = unique(results_adj_p_log$ILMN_GENE),
                 OrgDb         = org.Hs.eg.db,
                 universe      = unique(data$ILMN_GENE),
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = pthres,
                 qvalueCutoff  = adjpthres,
                 readable=T
)
dotplot(GO_all,font.size=15)
#close plot
dev.off()


# plot
pdf(file.path(dir_name, paste0(cond1, "_vs_", cond2, "_sig_GO_PLOT_ISG.pdf")), width = 10, height = 15)
#------ GO-ISG enrichment---------#
GO_ISG <- enrichGO(gene         = unique(sigLOG_top_table_inf$ILMN_GENE),
                 OrgDb         = org.Hs.eg.db,
                 universe      = unique(data$ILMN_GENE),
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = pthres,
                 qvalueCutoff  = adjpthres,
                 readable=T
)
dotplot(GO_ISG, font.size=15)
#close plot
dev.off()



rm(list = ls())

