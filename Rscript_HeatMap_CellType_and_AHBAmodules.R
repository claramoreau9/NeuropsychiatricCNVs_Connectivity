
#--------------------------------------------------------------------
# Script to generate Heatmaps for CellType and AHBA modules using ComplexHeatmap package (Figure 6)
#--------------------------------------------------------------------
# Figure 6 Panels: Clara et al. 
# Kuldeep Kumar
# Mar 2021
#
#--------------------------------------------------------------------
# Input: Pre-computed CorrPerGene, Pval-LabelShuffle and Pval-BrainSmash for FC maps
#--------------------------------------------------------------------
#
# NOTE: 1) Pre-computed Correlation, Pval-LabelShuffle and Pval-BrainSMASH pvalues are loaded
#         --> These can be computed using Label-Shuffle null-maps + BrainSMASH package in python
#         --> See Moreau et. al. Nat Commun 2020 (https://doi.org/10.1038/s41467-020-18997-2) 
#             for correlation between FC-maps and AHBA gene expression profiles
#       2) Cell-Type marker genes are explicitly coded in this script
#         ---> See Lake et. al. Nat Bio 2018 (https://doi.org/10.1038/nbt.4038); and 
#         ----> Anderson et al. Nat Comm. 2020 (https://doi.org/10.1038/s41467-020-16710-x)
#         ----> Anderson et al. PNAS 2020 (https://doi.org/10.1073/pnas.2008004117)
#              for the list of CellType marker genes
#       3) AHBA co-expression module eigen genes are explicitly coded in this script
#         ---> See Table 1 from Hawrylycz et al. Nat Neuro 2015 (https://doi.org/10.1038/nn.4171) 
#             for the list of AHBA co-expression marker genes
#
#--------------------------------------------------------------------
# Functional Connectivity Maps
#--------------------------------------------------------------------
# Note: 1. FC maps order and subset (out of 41) FC maps are defined manually
#        
#--------------------------------------------------------------------
# Output: Heatmaps + k-means clustering
#--------------------------------------------------------------------
#
#  ---> Heatmaps are made using CompleHeatmap package in R; 
#  ---> Each cell in heatmap is corresponds to Pearson correlation between 
#       FC-map and spatial pattern of Gene Expression (for the select gene)
#  ---> P-values for 1) Label Shuffle and 2) BrainSmash are shown in each plot
#  ---> Rows and Columns are clustered using  k-means clustering
#       (consensus clusters with 1000 iterations are reported)
# 
#--------------------------------------------------------------------



#---------------------------------------------------------------------
# set directories and paths
#---------------------------------------------------------------------

# Change this DIR/Folder 
in_wd <- 'D:/SSDN_KD/SSDN_Postdoc/Temp_Analysis/CrossCNV_fMRI_Clara/Shared_Code_for_Clara/FCmaps_GeneExpressionAnalysis_CellType_AHBAmodules_apr2021'
#in_wd <- getwd() 
setwd(in_wd)

data_dir <- in_wd   #'data'
plots_dir <- in_wd   #'plots'
#dir.create(plots_dir,showWarnings = FALSE)

#---------------------------------------------------------------------
# Libraries
#---------------------------------------------------------------------

library(readxl)
library(ggplot2)
library(circlize)      # needed for colormap definition
library(ComplexHeatmap)


#---------------------------------------------------------------------
# Parameters for plots
#---------------------------------------------------------------------


flag_legend <- 2        # 1: draw legend on left; ELSE: without legend
in_fontsize_stars <- 16
in_angle_rot_column_names <- 60
n_km_repeats <- 1000

# Color function for the HeatMap (uses circlize package)
col_fun = colorRamp2(c(-1, 0, 1), c("royalblue2","white", "red2"))

flag_dendogram <- FALSE  # => Show the dendograms in the heatmap
flag_transpose <- 1    # 1 => FC maps are ROWs; else FC-maps will be columns

flag_all_signatures <- 2  # 1: All 41 FCs; Else: 20 FC signatures from Fig 6 (Defined in 1.B.)

flag_select_AHBA_modules <- 1    # 1: 18 AHBA modules (Fig 6); Else: all 32 modules (Defined in 1.D.)

#---------------------------------------------------------------------
# 1.A. Data: Read Data -> Correlation + Pvalues -Per-Gene
#---------------------------------------------------------------------

flag_measure <- 'THAL'   # variable for ROIname

excel_filename <- paste0(data_dir,"/df_CorrPerGene_PvalLabelShuffle_PvalBrainSmash_MIST64_",flag_measure,".xlsx")
sheets <- readxl::excel_sheets(excel_filename)     # "Corr" "PvalLabelShuffle" "PvalBrainSmash"  

# 1. Correlation Per Gene csv files
df_CorrPerGene <- as.data.frame(readxl::read_excel(excel_filename, sheet = "Corr"))

# 2. Pvalue BrainSMASH
df_PvalCorrPerGene_BrainSmash <- as.data.frame(readxl::read_excel(excel_filename, sheet = "PvalBrainSmash"))

# 3. Pvalue Label-Shuffle
df_PvalCorrPerGene_LabShuffle <- as.data.frame(readxl::read_excel(excel_filename, sheet = "PvalLabelShuffle"))

#---------------------------------------------------------------------
# 1.B. FC-maps: select FC-maps for Heatmaps vs All-FC maps + modify FC map names
#---------------------------------------------------------------------

all_FC_signatures <- colnames(df_CorrPerGene)

# Modify as needed: currently selected on the basis of FC maps in Clara' et al.'s Figure 6

select_FC_signatures <- c("GeneSymbol",
                          "ASD","BIP","SZ",    
                          "PRS_ASD","PRS_MDD","PRS_SZ", "PRS_IQ","PRS_SA", 
                          "NT","CT","Gfactor","SA","FluidIntel",
                          "DEL1q21_1","DEL15q11_2","DEL16p11_2","DEL22q11_2", 
                          "DUP1q21_1","DUP16p11_2","DUP22q11_2")

# Simiplified FC map names for Plots
select_FC_signatures_simplified <- c("GeneSymbol",
                                     "ASD","BIP","SZ",
                                     "PGS ASD", "PGS MDD","PGS SZ", "PGS IQ", "PGS SA",
                                     "NT", "CT", "Gfactor", "SA","Fluid Intel",
                                     "DEL 1q21.1", "DEL 15q11.2","DEL 16p11.2","DEL 22q11.2",
                                     "DUP 1q21.1", "DUP 16p11.2","DUP 22q11.2")


if( flag_all_signatures == 1){
  select_FC_column_names <- all_FC_signatures
  suffix_FC_sigantures <- 'AllFCs'
} else {
  select_FC_column_names <- select_FC_signatures
  suffix_FC_sigantures <- 'SelectFCs'
}


#------------------------------------------------------------------------------
# 1.C. Cell-Type marker genes
#------------------------------------------------------------------------------ 

# Cell-Type marker genes are explicitly coded in this script
#         ---> See Lake et. al. Nat Bio 2018 (https://doi.org/10.1038/nbt.4038); and 
#         ----> Anderson et al. Nat Comm. 2020 (https://doi.org/10.1038/s41467-020-16710-x)
#         ----> Anderson et al. PNAS 2020 (https://doi.org/10.1073/pnas.2008004117)
#              for the list of CellType marker genes

# Define 3 lists (explicit):
#       1) List of genes; 2) Cell-Type (Marker-gene) for plot; 3) Cell-Type only

list_CellType_MarkerGenes = c('CBLN2','NEFM','TSHZ2','HS3ST2','NR4A2',
                              'CNR1','SHISA8','RELN','PVALB','SST',
                              'CLDN5','LAMA2','SLC4A4','MBP','LHFPL3','DOCK8')

list_namesCellType_and_MarkerGenes = c('Ex1 (CBLN2)','Ex3 (NEFM)','Ex4 (TSHZ2)','Ex5 (HS3ST2)','Ex8 (NR4A2)', #  # Ex (Excititory) Cell-Types
                             'In1 (CNR1)','In3 (SHISA8)','In4 (RELN)','In6 (PVALB)','SST', # # In (Inhibitory) Cell-Types
                             'End (CLDN5)','Per (LAMA2)','Ast (SLC4A4)','Oli (MBP)','OPC (LHFPL3)','Mic (DOCK8)') 

list_namesCellType = c('Ex1','Ex3','Ex4','Ex5','Ex8',
                        'In1','In3','In4','In6','SST',
                        'End','Per','Ast','Oli','OPC','Mic')


#--------------------------------------------------------------------------------
# 1.D. AHBA Co-expression modules eigen genes
#--------------------------------------------------------------------------------

# AHBA co-expression module eigen genes are explicitly coded in this script
# ---> See Table 1 from Hawrylycz et al. Nat Neuro 2015 (https://doi.org/10.1038/nn.4171) 
#      for the list of AHBA co-expression marker genes
#

if( flag_select_AHBA_modules > 1){
  # All modules
  list_AHBA_ModuleEigenGenes <- c("GABRB3","NEUROD6","KCNAB2","GABARAPL1",
                                  "MAFG","MEF2C","NGEF","BEX2",
                                  "PGAP1","PDE1B","NTNG1","SLC6A3",
                                  "KRT18","TLE6", "NEFH", "SLC47A1","PAXIP1","PTS",
                                  "VDAC2","B3GAT1","GBP4", "SYNE2","LZIC", "POGZ", "RGS10","C1orf194","NPC2",
                                  "SERPINA6","GAS5","VAMP3","BNIP2","SLC25A18")
  
  # Simplify Anatomy names: Remove WM, and DP; and max 3 Anatomy structures
  list_AHBA_Modules_and_EigenGenes_plusAnat <- c("M1-GABRB3 (Tel)", "M2-NEUROD6 (Hp,Amg)","M3-KCNAB2 (Hp,Th)", "M4-GABARAPL1 (Thc)",
                                                "M5-MAFG (Hp,Amg)","M6-MEF2C (Ncx,Cl)","M7-NGEF (Str,Ncx,Amg)", "M8-BEX2 (Cx,Amg,Hy)",
                                                "M9-PGAP1 (Hp,Amg,Hy)","M10-PDE1B (Str)","M11-NTNG1 (DT)","M12-SLC6A3 (SN,VTA)",
                                                "M13-KRT18 (Acn,RF)","M14-TLE6 (Hy)","M15-NEFH (DCbN,BS)","M16-SLC47A1 (DG)",
                                                "M17-PAXIP1 (CbCx)","M18-PTS (Hy,SN)","M19-VDAC2 (Th,CbN,BS)","M20-B3GAT1 (NCx,BG,VT)",
                                                "M21-GBP4 (SmN,Ch)","M22-SYNE2 (Epy,CB)","M23-LZIC (BS,GP)","M24-POGZ (CbCx,BG)",
                                                "M25-RGS10 (Epy,SN)","M26-C1orf194 (Epy)","M27-NPC2 (Ch,Epy)","M28-SERPINA6 (IB-HBn)",
                                                "M29-GAS5 (SN,GP)","M30-VAMP3 (VT,GP)","M31-BNIP2 (GP)","M32-SLC25A18 (Str,Amg,SN)" )
  
  list_AHBA_Modules_and_EigenGenes <- list_AHBA_Modules_and_EigenGenes_plusAnat
  
  list_AHBA_ModulesNames = c("M01", "M02","M03", "M04",
                             "M05","M06","M07", "M08",
                             "M09","M10", "M11", "M12",
                             "M13","M14", "M15", "M16",
                             "M17","M18", "M19", "M20",
                             "M21","M22", "M23","M24",
                             "M25",  "M26", "M27", "M28",
                             "M29", "M30", "M31","M32")
  
  
} else {
  # Select modules
  
  list_AHBA_ModuleEigenGenes <- c("GABRB3","KCNAB2","GABARAPL1",
                                  "MEF2C","NGEF",
                                  "PGAP1","PDE1B",
                                  "NTNG1","SLC6A3",
                                  "TLE6", "NEFH", 
                                  "PAXIP1","VDAC2",
                                  "B3GAT1","POGZ", 
                                  "GAS5","VAMP3","SLC25A18")
  
  # Simplify Anatomy names: Remove WM (WhiteMatter), BrainStem (BS),and DG; and max 3 Anatomy structures
  list_AHBA_Modules_and_EigenGenes_plusAnat <- c("M1-GABRB3 (Tel)","M3-KCNAB2 (Hp,Thal)","M4-GABARAPL1 (Thal-Cort)",
                                                "M6-MEF2C (Ncx,Cl)","M7-NGEF (Str,Ncx,Amg)",
                                                "M9-PGAP1 (Hp,Amg,Hy)","M10-PDE1B (Str)",
                                                "M11-NTNG1 (dThal)","M12-SLC6A3 (SN,VTA)",
                                                "M14-TLE6 (Hy)","M15-NEFH (DCbN,BS)",
                                                "M17-PAXIP1 (CbCx)","M19-VDAC2 (Thal,CbN)",
                                                "M20-B3GAT1 (NCx,BG,vThal)","M24-POGZ (CbCx,BG)",
                                                "M29-GAS5 (SN,GP)","M30-VAMP3 (vThal,GP)","M32-SLC25A18 (Str,Amg,SN)" )
  
  list_AHBA_Modules_and_EigenGenes <- list_AHBA_Modules_and_EigenGenes_plusAnat
  
  list_AHBA_ModulesNames = c("M01","M03","M04",
                             "M06","M07",
                             "M09","M10", 
                             "M11","M12",
                             "M14","M15", 
                             "M17","M19",
                             "M20","M24",
                             "M29","M30","M32")   # 18 modules
}


#--------------------------------------------------------------------------------
# 2. Function for Heatmap + k-means clustering
#--------------------------------------------------------------------------------

fComplexHeatmap_Corr_Pval_kmeans <- function(df_CorrPerGene,df_PvalCorrPerGene_BrainSmash,df_PvalCorrPerGene_LabShuffle,in_list_marker_genes,in_list_marker_genes_and_info,select_FC_column_names,col_fun,in_row_km,in_col_km,n_km_repeats,in_angle_rot_column_names,in_fontsize_stars,flag_transpose,flag_dendogram){
  
  
  #---- SubSet: FC-maps (Columns) and Gene set (Rows) ----------------------------------
  df_CorrPerGene_marker <- df_CorrPerGene[df_CorrPerGene$GeneSymbol %in% in_list_marker_genes,select_FC_column_names]
  df_PvalCorrPerGene_marker <- df_PvalCorrPerGene_BrainSmash[df_PvalCorrPerGene_BrainSmash$GeneSymbol %in% in_list_marker_genes,select_FC_column_names]
  df_PvalCorrPerGene_LabShuffle_marker <- df_PvalCorrPerGene_LabShuffle[ df_PvalCorrPerGene_LabShuffle$GeneSymbol  %in% in_list_marker_genes,select_FC_column_names]
  
  #------------ Create Matrices (drop GeneSymbol column)-------------------------------
  Corr_Mat <- as.matrix(df_CorrPerGene_marker[,c(2:ncol(df_CorrPerGene_marker))])
  pMat_BrainSmash <- as.matrix(df_PvalCorrPerGene_marker[,c(2:ncol(df_PvalCorrPerGene_marker))])
  pMat_Label <- as.matrix(df_PvalCorrPerGene_LabShuffle_marker[,c(2:ncol(df_PvalCorrPerGene_LabShuffle_marker))])
  
  rownames(Corr_Mat) <- df_CorrPerGene_marker[,1]
  rownames(pMat_BrainSmash) <- df_CorrPerGene_marker[,1]
  rownames(pMat_Label) <- df_CorrPerGene_marker[,1]
  
  colnames(Corr_Mat) <- select_FC_signatures_simplified[-1]
  colnames(pMat_BrainSmash) <- select_FC_signatures_simplified[-1]
  colnames(pMat_Label) <- select_FC_signatures_simplified[-1]
  
  #------- Check the order of GeneSet (explicit) --------------------
  Corr_Mat <- Corr_Mat[match(in_list_marker_genes,rownames(Corr_Mat)),]
  pMat_BrainSmash <- pMat_BrainSmash[match(in_list_marker_genes,rownames(pMat_BrainSmash)),]
  pMat_Label <- pMat_Label[match(in_list_marker_genes,rownames(pMat_Label)),]
  
  #------------------- Update Rownames: Gene Names with Cell-Type info --------------------------
  rownames(Corr_Mat) <- in_list_marker_genes_and_info
  rownames(pMat_BrainSmash) <- in_list_marker_genes_and_info
  rownames(pMat_Label) <- in_list_marker_genes_and_info
  
  #-------------------------------------------------------------------------------
  
  # Transpose to place Gene Markers as columns and FC-signatures as rows 
  # ELSE: FC-signatures are columns and Cell-Markers are rows
  if(flag_transpose == 1){
    Corr_Mat <- t(Corr_Mat)
    pMat_BrainSmash <- t(pMat_BrainSmash)       # Pvalue BrainSmash
    pMat_Label <- t(pMat_Label)  # P-value Label-Shuffle
  }
  
  # Heatmap
  ht = Heatmap(Corr_Mat, name = "Corr", row_km = in_row_km, column_km = in_col_km,
               row_km_repeats = n_km_repeats,column_km_repeats = n_km_repeats,
               column_names_rot = in_angle_rot_column_names,
               col = col_fun,
               # add *, + , x for p-value using cell_fun
               cell_fun = function(j, i, x, y, width, height, fill) {
                 
                 if( (pMat_BrainSmash[i, j] < 0.05) & (pMat_Label[i, j] < 0.05) ){
                   grid.text("*", x, y, gp = gpar(fontsize = in_fontsize_stars+1, fontface = "bold"))   # Both significant
                 } else if( (pMat_BrainSmash[i, j] < 0.05) & (pMat_Label[i, j] >= 0.05) ){
                   grid.text("x", x, y, gp = gpar(fontsize = in_fontsize_stars))   # BrainSMASH only
                 } else if( (pMat_BrainSmash[i, j] >= 0.05) & (pMat_Label[i, j] < 0.05) ){
                   grid.text("+", x, y, gp = gpar(fontsize = in_fontsize_stars))   # Label-Shuffle only
                 }
                 
               }, show_column_dend = flag_dendogram, show_row_dend = flag_dendogram,
               use_raster = TRUE,
               row_title = NULL, column_title = NULL
  )
  #draw(ht,heatmap_legend_side = "left")
  
  return(ht)  
}






#--------------------------------------------------------------------------------
# 3.A. HeatMap for Cell-Type marker genes
#--------------------------------------------------------------------------------

marker_gene_exp <- "CellTypeMarkerGenes"

#---- DEFINE: the lists of Marker Genes and Names etc.
in_list_marker_genes <- list_CellType_MarkerGenes  
in_list_marker_genes_and_info <- list_namesCellType_and_MarkerGenes   
in_list_markerGroups <- list_namesCellType 

#-------------------------------------------------------------------------------
nMarkerGenes <- length(in_list_marker_genes)
nFCmaps <- length(select_FC_column_names) -1  # remove GeneSymbol 

# Define number of Row and Column Clusters
in_row_km <- floor(nFCmaps/5)
in_col_km <- floor(nMarkerGenes/5)

# Define the width and height of the plot in inches
in_width <- 0.4*nMarkerGenes
in_height <- 0.4*nFCmaps

#------------------------------------------------------------------------------
# Make Heatmap

ht <- fComplexHeatmap_Corr_Pval_kmeans(df_CorrPerGene,df_PvalCorrPerGene_BrainSmash,df_PvalCorrPerGene_LabShuffle,in_list_marker_genes,in_list_marker_genes_and_info,select_FC_column_names,col_fun,in_row_km,in_col_km,n_km_repeats,in_angle_rot_column_names,in_fontsize_stars,flag_transpose,flag_dendogram)

#------------------------------------------------------------------------------
# Save plot

set.seed(1)  # Seed for k-means clustering
# save as SVG
plot_filename <- paste0(plots_dir,"/HeatMap_",marker_gene_exp,"_",flag_measure,"_",suffix_FC_sigantures,".svg")
svg(plot_filename,width = in_width, height = in_height)
if(flag_legend == 1){
  draw(ht,padding = unit(c(2, 12, 1, 1), "mm"),heatmap_legend_side = "left")
} else {
  draw(ht,padding = unit(c(2, 12, 1, 1), "mm"),show_heatmap_legend = FALSE)
}
dev.off()

# save as PDF
plot_filename <- paste0(plots_dir,"/HeatMap_",marker_gene_exp,"_",flag_measure,"_",suffix_FC_sigantures,".pdf")
pdf(plot_filename,width = in_width, height = in_height)
if(flag_legend == 1){
  draw(ht,padding = unit(c(2, 12, 1, 1), "mm"),heatmap_legend_side = "left")
} else {
  draw(ht,padding = unit(c(2, 12, 1, 1), "mm"),show_heatmap_legend = FALSE)
}
dev.off()




#--------------------------------------------------------------------------------
# 3.B. HeatMap for AHBA co-expression module eigen genes
#--------------------------------------------------------------------------------

marker_gene_exp <- "AHBAModuleEigenGenes"

#---- DEFINE: the lists of Marker Genes and Names etc.
in_list_marker_genes <- list_AHBA_ModuleEigenGenes
in_list_marker_genes_and_info <- list_AHBA_Modules_and_EigenGenes
in_list_markerGroups <- list_AHBA_ModulesNames

#-------------------------------------------------------------------------------
nMarkerGenes <- length(in_list_marker_genes)
nFCmaps <- length(select_FC_column_names) -1  # remove GeneSymbol 

# Define number of Row and Column Clusters
in_row_km <- floor(nFCmaps/5)
in_col_km <- floor(nMarkerGenes/5)

# Define the width and height of the plot in inches
in_width <- 0.4*nMarkerGenes
in_height <- 0.4*nFCmaps

#------------------------------------------------------------------------------
# Make Heatmap

ht <- fComplexHeatmap_Corr_Pval_kmeans(df_CorrPerGene,df_PvalCorrPerGene_BrainSmash,df_PvalCorrPerGene_LabShuffle,in_list_marker_genes,in_list_marker_genes_and_info,select_FC_column_names,col_fun,in_row_km,in_col_km,n_km_repeats,in_angle_rot_column_names,in_fontsize_stars,flag_transpose,flag_dendogram)
  
#------------------------------------------------------------------------------
# Save plot

set.seed(1)  # Seed for k-means clustering
# save as SVG
plot_filename <- paste0(plots_dir,"/HeatMap_",marker_gene_exp,"_",flag_measure,"_",suffix_FC_sigantures,".svg")
svg(plot_filename,width = in_width, height = in_height)
if(flag_legend == 1){
  draw(ht,padding = unit(c(2, 18, 1, 1), "mm"),heatmap_legend_side = "left")
} else {
  draw(ht,padding = unit(c(2, 18, 1, 1), "mm"),show_heatmap_legend = FALSE)
}
dev.off()

# save as PDF
plot_filename <- paste0(plots_dir,"/HeatMap_",marker_gene_exp,"_",flag_measure,"_",suffix_FC_sigantures,".pdf")
pdf(plot_filename,width = in_width, height = in_height)
if(flag_legend == 1){
  draw(ht,padding = unit(c(2, 18, 1, 1), "mm"),heatmap_legend_side = "left")
} else {
  draw(ht,padding = unit(c(2, 18, 1, 1), "mm"),show_heatmap_legend = FALSE)
}
dev.off()


#------------------------------------------------------------------------------
# End 
#------------------------------------------------------------------------------

