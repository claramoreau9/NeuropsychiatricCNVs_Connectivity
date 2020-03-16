
# Clara Moreau Feb 2020

library('circlize')
library('Matrix')
library('csv')
library('gdata')
library('RColorBrewer')
library('viridis')

setwd("your_repo_path/data/")


##load your networks names (in my case, MIST 64 and MIST 12)
label<- read.csv2 ("summary_label_64.csv",header = T, sep = ";")

#load your beta map, 64 * 64
eff<- read.table ("cnv_ukbb_del16p_vs_con_betas.tsv",header = F, sep = "\t")
eff = eff[-1,-1]
eff <- data.frame(apply(eff, 2, as.numeric)) 

#to threshold it with qval (64*64)
fdr<- read.csv2 ("cnv_ukbb_del16p_vs_con_qval.tsv",header = F, sep = "\t")
fdr = fdr[-1,-1]
fdr <- data.frame(apply(fdr, 2, as.numeric)) 

threshold = 0.8 #to plot 80%
thr_fdr = TRUE

# Make dataframes
eff <- as.matrix(eff)
rownames(eff) <- label$roi_label
colnames(eff) <- label$roi_label

fdr <- as.matrix(fdr)
rownames(fdr) <- label$roi_label
colnames(fdr) <- label$roi_label

# Make simple ord_names => 64 regions to 12 networks
ord_names <- label[label$s12_ID==1,]$roi_label
for (i in 2:12) {
  ord_names <- append(as.character(ord_names), as.character(label[label$s12_ID==i,]$roi_label))
}

# Mask the matrix
if (thr_fdr) {
  eff[fdr>0.05] <- 0
}

conn = upperTriangle(eff, diag=FALSE)
nz = abs(conn[conn!=0])
ind = sort(nz, index.return=TRUE, decreasing = TRUE)
slice = floor(length(ind$ix)*threshold)
thr = min(nz[ind$ix[1:slice]])

# Threshold the connections
eff[(eff)<thr] <- 0

# Find regions with no above threshold connections
fin <- which(rowSums(eff)==0)

# Give them a number which is the smallest connection sum and divide by the number of connections
numb <-min(rowSums(eff)[rowSums(eff)>0])/64

# Now give this number to the nodes without above threshold connections
eff[fin,] <- numb
eff[, fin] <- numb

tn <- colnames(eff)

n_neg=length(unique(eff[eff<0]))
col_neg=colorRampPalette(c("white", "blue"))(n_neg)

pos_values<-sort(unique(eff[eff>0]))
col_pos=colorRampPalette(c("white", "red"))(length(pos_values))

col_map = function(x) ifelse(x < 0, col_neg, ifelse(x==numb, "#00000000", ifelse(x==0,"00000000",col_pos))) #black background
colors <- c('#ffffb3','#8dd3c7','#ffed6f','#fccde5','#fb8072','#fdb462','#bc80bd','#bebada','#b3de69','#ccebc5','#80b1d3','#d9d9d9')

#check the number of connections to report
length(pos_values)
n_neg

# and plot with labels 
chordDiagram(eff,  col=col_map, grid.col = "grey25", order = ord_names, keep.diagonal = FALSE, preAllocateTracks = list(track.height = 0.05), symmetric=TRUE, annotationTrack = "grid")

net_names <- levels(label$s12_label)
for (i in seq_along(net_names)) {
  name <- net_names[i]
  rois <- label[label$s12_label==name,]$ROI
  highlight.sector(tn[rois], track.index = 1, col = colors[i],
                   text = name, cex = 0.9, text.col = "black", niceFacing = TRUE)
}

circos.track(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  xplot = get.cell.meta.data("xplot")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[2]+1, sector.name, facing = "clockwise",
              niceFacing = TRUE, cex = 0.85, adj = c(0.1, 0.5), col = "white") 
  #niceFacing = TRUE, cex = 0.7, adj = c(0, 0.5), col = "black") 
  
}, bg.border = NA)


#or without labels
chordDiagram(eff,  col=col_map, grid.col = "grey25", 
             order = ord_names, keep.diagonal = FALSE, 
             preAllocateTracks = list(track.height = 0.05), symmetric=TRUE, 
             annotationTrack = "grid")
net_names <- levels(label$s12_label)
for (i in seq_along(net_names)) {
  name <- net_names[i]
  rois <- label[label$s12_label==name,]$ROI
  highlight.sector(tn[rois], track.index = 1, col = colors[i], niceFacing = TRUE)
}
