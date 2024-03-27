### SCPA analysis
### https://jackbibby1.github.io/SCPA/

setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/scRNA-seq/cell_ranger_count_output")

## instructions to install SCPA R package
install.packages("devtools")
devtools::install_version("crossmatch", version = "1.3.1", repos = "http://cran.us.r-project.org")
devtools::install_version("multicross", version = "2.1.0", repos = "http://cran.us.r-project.org")
devtools::install_github("jackbibby1/SCPA")

## load the required libraries
library(SCPA)
library(Seurat)
library(msigdbr)
library(magrittr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

## load the saved seurat integrated object of endometrium and peritoneal lesion samples

options(future.globals.maxSize = 4000 * 1024^2)

seurat_integrated <- readRDS(file = "clusters_seurat_integrated.rds")

#seurat_integrated@meta.data

# Assign identity of clusters 
Idents(object = seurat_integrated) <- "integrated_snn_res.0.6"

## rename idents by merged clusters names
new.cluster.ids <- c("Stromal",
                     "Endothelial",
                     "Stromal",
                     "Immune",
                     "Immune",
                     "Stromal",
                     "Perivascular",
                     "Perivascular",
                     "Macrophage",
                     "Perivascular",
                     "Unknown",
                     "Stromal",
                     "Bcell",
                     "Epithelial",
                     "Perivascular",
                     "Cyclingstromalcells",
                     "Endothelial",
                     "Immune")


levels(seurat_integrated)

names(new.cluster.ids) <- levels(seurat_integrated)

## change identity of clusters and rename it
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
levels(seurat_integrated)


levels(seurat_integrated@active.ident)
seurat_integrated@meta.data["celltype"] <- seurat_integrated@active.ident

head(seurat_integrated@meta.data)

saveRDS(seurat_integrated, file = "clusters_seurat_integrated_withcelltypes.rds")


###########################################################################################################################
#Extracting the expression matrices
#We then need the expression matrices for the populations we want to compare. 
#For this, we can use the seurat_extract function from the SCPA package, 
#which takes a Seurat object and subsets the data based on the metadata specified in the Seurat columns

#endometrium
endometrium <- seurat_extract(seurat_integrated,
                          meta1 = "sampleType", value_meta1 = "Endometrium",
                          )
dim(endometrium)

# peritoneal lesion
peritoneal <- seurat_extract(seurat_integrated,
                            meta1 = "sampleType", value_meta1 = "Peritonial_lesion")

dim(peritoneal)
endometrium[1:5,1:5]


### Defining metabolic pathways

#After we have our data, we can define the pathways to compare. 
#We curated a list of metabolic pathways.
#The file is in standard gmt file format with the pathway name in column 1 and genes of that pathway in subsequent columns. 
#Here we just converted a gmt file to csv, and manually curated a single file with all metabolic gene sets. 

## read the pathway file
pathways <- "KEGG_pathway14_genelist.csv"
pathways


###Running SCPA
##So now we have our samples and metabolic gene sets. 
#To compare these pathways, we can just run the compare_pathways function.

#########################################################################################################
############### compare pathways in peritoneal lesion vs eutopic endometrium ##########################
#########################################################################################################

endo_peri <- compare_pathways(samples = list(peritoneal, endometrium), 
                              pathways = pathways)

endo_peri

# Rank Plot
p2 <- plot_rank(endo_peri, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2



#Visualizing results
#Now we have the results stored in the endo_peri object, which we can then visualize.

endo_peri_viz <- endo_peri %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))


endo_peri_viz
aa_path <- endo_peri_viz %>% 
  filter(grepl(pattern = "a|e", ignore.case = T, x = Pathway))

aa_path


## enrichment plot
ggplot(endo_peri_viz, aes(-FC, qval, label = Pathway)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 1.6, shape = 21, fill = endo_peri_viz$color, stroke = 0.3) +
  geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  geom_text(hjust=0, vjust=0) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)

dev.off()



######################################################################################################################
################# To create rank plot of pathways in each cell type ################################################################

split_tissue <- SplitObject(seurat_integrated, split.by = "sampleType")
split_tissue

celltype <- levels(seurat_integrated@active.ident)
celltype

pathways <- "KEGG_pathway14_genelist.csv"
pathways


# We first create empty lists to store results from the for loop
endo_peri_celltype <- list()
for (i in celltype) {
  
  # We extract expression data using `seurat_extract` based on tissue, cell_type ("fine"), and stimulation ("none")
  endo <- seurat_extract(split_tissue$Endometrium, 
                       meta1 = "celltype", value_meta1 = i)
  
  peri<- seurat_extract(split_tissue$Peritonial_lesion, 
                       meta1 = "celltype", value_meta1 = i)
  
  
  # We then compare all tissues to the blood using `compare_pathways`
  print(paste("comparing", i))
  endo_peri_celltype[[i]] <- compare_pathways(list(peri, endo), pathways)
  
  
}

endo_peri_celltype

#endo_peri_celltype$Stromal
# Stromal cells
p2 <- plot_rank(endo_peri_celltype$Stromal, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2

# enothelial cells
p2 <- plot_rank(endo_peri_celltype$Endothelial, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2

# Immune cells
p2 <- plot_rank(endo_peri_celltype$Immune, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2

# Macrophage 
p2 <- plot_rank(endo_peri_celltype$Macrophage, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2

# PERIVASCULAR
p2 <- plot_rank(endo_peri_celltype$Perivascular, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2

# epithelial
p2 <- plot_rank(endo_peri_celltype$Epithelial, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2


# cycling stromal cells
p2 <- plot_rank(endo_peri_celltype$Cyclingstromalcells, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2

# Bcells
p2 <- plot_rank(endo_peri_celltype$Bcell, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2

# Unknown cells
p2 <- plot_rank(endo_peri_celltype$Unknown, pathway = "a|e",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")


p2

###################################################################################################################
########### create qvalue heatmap of celltypes in two populations are compared with each other ###########################

## Population 1: Peritoneal lesion samples
## Population 2: Eutopic endometrium samples

get_qvals <- function(scpa_out, name) {
  
  df <- list()
  for (i in names(scpa_out)) {
    df[[i]] <- scpa_out[[i]] %>%
      select(Pathway, qval)
  }
  
  col_names <- names(df)
  for (i in 1:length(df)) {
    df[[i]] <- set_colnames(df[[i]], c("pathway", paste(name, col_names[[i]], sep = "_")))
  }
  
  return(df)
  
}


cell_types <- unique(seurat_integrated$celltype)
cell_types

##################################################################################################################
######################## for heatmap of Qvalue ##################################################################
peri_vs_endo <- SplitObject(seurat_integrated, split.by = "sampleType")


scpa_out <- list()
for (i in cell_types) {
  
  Endometrium <- seurat_extract(peri_vs_endo$Endometrium, 
                            meta1 = "celltype", value_meta1 = i)
  
  peritoneal <- seurat_extract(peri_vs_endo$Peritonial_lesion, 
                          meta1 = "celltype", value_meta1 = i)
  
  print(paste("comparing", i))
  scpa_out[[i]] <- compare_pathways(list(peritoneal, Endometrium), pathways) %>%
    select(Pathway, qval) %>%
    set_colnames(c("Pathway", paste(i, "qval", sep = "_")))
  
}

scpa_out

scpa_out_ <- scpa_out %>% 
  reduce(full_join, by = "Pathway") %>% 
  set_colnames(gsub(colnames(.), pattern = " ", replacement = "_")) %>%
  select(c("Pathway", grep("_qval", colnames(.)))) %>%
  filter_all(any_vars(. > 2)) %>%
  column_to_rownames("Pathway")


scpa_out_

# rename the Adenine, Guanine and Glutamine metabolism to ADQ metabolism
rownames(scpa_out_)[8] <- "ADQ metabolism"

## create a object with pathway names
imp_pathways <- rownames(scpa_out_)
imp_pathways


### to plot heatmap 
#create row annotations for the heatmap to highlight these pathways, and colour scale
position <- which(rownames(scpa_out_) %in% imp_pathways)
position

# Row annotaions
row_an <- rowAnnotation(Genes = anno_mark(at = which(rownames(scpa_out_) %in% imp_pathways),
                                          labels = rownames(scpa_out_)[position],
                                          labels_gp = gpar(fontsize = 10),
                                          link_width = unit(2.5, "mm"),
                                          padding = unit(1, "mm"),
                                          link_gp = gpar(lwd = 0.5)))
row_an

## col annotations
col_hm <- colorRamp2(colors = c("blue", "white", "red"), breaks = c(0, 3, 6))

colnames(scpa_out_)

## plot the heatmap with the annotations
heatmap <- Heatmap(scpa_out_,
        col = col_hm,
        name = "Qval",
        show_row_names = T,
        #right_annotation = row_an,
        column_names_gp = gpar(fontsize = 8, fontface = "bold"),
        border = F,
        column_km = 3,
        row_km = 3,
        column_labels = c("Stromal", "Unknown",
                          "Perivascular", "Immune",             
                          "Endothelial", "Cyclingstromalcells",
                          "Macrophage", "Bcell",              
                          "Epithelial"))


## save heatmap
heatmap
save(heatmap, file = "heatmap_qvalue_scpa.RData")
scpa_out_
save(scpa_out_, file = "scpa_output.RData")
write.csv(scpa_out_, file = "scpa_out_peri_vs_endo.csv")
