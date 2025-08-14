# File name:      250708_BenedictAgee_CAEC_figures_eeb.R
# Author:         Emily Benedict, ebenedict@wustl.edu
# Created On:     2025-07-08
# Description:    This script generates Figures 1c, 3a-c, 4a-d and Supplemental Figures 1-4

# loading packages
library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(NatParksPalettes)
library(RColorBrewer)
library(phytools)
library(ggtree)
library(pheatmap)
library(tidytree)
library(ggnewscale)
library(lubridate)
library(kableExtra)
library(ggpubr)
library(qgraph)

# setting working directory
setwd('~/Box Sync/EEB/Dantas/aim2_ESBLe/manuscript/')

#### reading in data ####
## metadata
metadata = data.frame(read_excel('./supp_tables/250804_SupplementalTable2_Metadata.xlsx'), stringsAsFactors = F)

## MLST
mlst = data.frame(read_excel('./supp_tables/250605_SupplementalTable3_Genomics.xlsx'), stringsAsFactors = F)

## amrfinder data
amrfinder = data.frame(read_excel('./supp_tables/250605_SupplementalTable4_AMRFinder.xlsx'), stringsAsFactors = F)

## core genome SNP data
snps = data.frame(read_excel('./supp_tables/250813_SupplementalTable5_CoreGenomeSNPs.xlsx'), stringsAsFactors = F)

## core genome tree
tree = ape::read.tree('./tree_file/BenedictAgee_CAECtree_RAxML_bestTree.raxml_core_genome_tree')
# subsetting core genome tree to those in metadata
tree = drop.tip(tree, setdiff(tree$tip.label, metadata$Sample))
# midpoint-rooting tree
tree = midpoint.root(tree)

## plot colors
colors = data.frame(read_excel('./supp_tables/241216_PlotColors.xlsx'), stringsAsFactors = F)

## ordering MLST and metadata by tree$tip.label
ordered_mlst = data.frame(matrix(NA, nrow = nrow(mlst[mlst$Sample %in% tree$tip.label,]), ncol = ncol(mlst)), stringsAsFactors = F)
colnames(ordered_mlst) = colnames(mlst)
ordered_mlst$Sample = tree$tip.label
rownames(ordered_mlst) = tree$tip.label

ordered_metadata = data.frame(matrix(NA, nrow = nrow(metadata[metadata$Sample %in% tree$tip.label,]), ncol = ncol(metadata)), stringsAsFactors = F)
colnames(ordered_metadata) = colnames(metadata)
ordered_metadata$Sample = tree$tip.label
rownames(ordered_metadata) = tree$tip.label

for(i in 1:length(tree$tip.label)){
 ordered_mlst[i,] = mlst[mlst$Sample == ordered_mlst$Sample[i],]
 ordered_metadata[i,] = metadata[metadata$Sample == ordered_metadata$Sample[i],]
}

# renaming ordered_metadata media
ordered_metadata$media[ordered_metadata$media == 'MacConkey Cipro'] = 'CiproR only'
ordered_metadata$media[ordered_metadata$media == 'CHROM ESBL'] = 'ESBL only'

# joining MLST and metadata
ordered_metadata_mlst = full_join(ordered_metadata, ordered_mlst)

# setting MLST colors
mlst_colors = colors$color[colors$topic == 'MLST']
names(mlst_colors) = colors$name[colors$topic == 'MLST']

# setting media colors
media_colors = colors$color[colors$topic == 'Media']
names(media_colors) = colors$name[colors$topic == 'Media']

# setting ESBL colors
esbl_colors = colors$color[colors$topic == 'ESBL']
names(esbl_colors) = colors$name[colors$topic == 'ESBL']

# setting quinolone resistance colors
quinolone_colors = colors$color[colors$topic == 'Quinolone']
names(quinolone_colors) = colors$name[colors$topic == 'Quinolone']

# setting AMRFinder colors
amrfinder_colors = colors$color[colors$topic == 'AMRFinder']
names(amrfinder_colors) = colors$name[colors$topic == 'AMRFinder']

# setting AST colors
ast_colors = colors$color[colors$topic == 'AST']
names(ast_colors) = colors$name[colors$topic == 'AST']

# setting antibiotic colors
abx_colors = c(natparks.pals('Redwood',6)[1],'grey95')
names(abx_colors) = c('present','absent')

# setting year colors
year_colors = colors$color[colors$topic == 'Year']
names(year_colors) = colors$name[colors$topic == 'Year']

# subsetting AMRFinder data to AMR hits only
amrfinder = amrfinder[amrfinder$Element.type == 'AMR',]

## ordering AMRFinder data by tree$tip.label
ordered_amrfinder = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = (length(table(amrfinder$Gene.symbol))+1)), stringsAsFactors = F)
colnames(ordered_amrfinder) = c('Sample', names(table(amrfinder$Gene.symbol)))
ordered_amrfinder$Sample = tree$tip.label
rownames(ordered_amrfinder) = tree$tip.label

for(i in 1:nrow(ordered_amrfinder)){
  for(j in 2:ncol(ordered_amrfinder)){
    if(colnames(ordered_amrfinder)[j] %in% amrfinder$Gene.symbol[amrfinder$Name == ordered_amrfinder$Sample[i]]){
      ordered_amrfinder[i,j] = colnames(ordered_amrfinder)[j]
    }else{
      ordered_amrfinder[i,j] = 'not present'
    }
  }
}

## generating a class-based AMR dataframe
ordered_amr_classes = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = (length(table(amrfinder$Class)) + 1)), stringsAsFactors = F)
colnames(ordered_amr_classes) = c('Sample',names(table(amrfinder$Class)))
ordered_amr_classes$Sample = tree$tip.label
rownames(ordered_amr_classes) = tree$tip.label
for(i in 1:nrow(ordered_amr_classes)){
  for(j in 2:ncol(ordered_amr_classes)){
    if(colnames(ordered_amr_classes)[j] %in% amrfinder$Class[amrfinder$Name == ordered_amr_classes$Sample[i]]){
      ordered_amr_classes[i,j] = colnames(ordered_amr_classes)[j]
    }else{
      ordered_amr_classes[i,j] = 'not present'
    }
  }
}

## generating an ESBL-specific dataframe
# subsetting to those with 'extended-spectrum' in their name of their closest gene match
esbl_data = amrfinder[grep('extended-spectrum', amrfinder$Name.of.closest.sequence),]
# adding class of ESBL as a column
esbl_data$esbl_class = gsub('extended-spectrum ','', esbl_data$Name.of.closest.sequence)
esbl_data$esbl_class = gsub(' beta-lactamase .*','',esbl_data$esbl_class)
esbl_data$Name.of.closest.sequence = gsub('.*EC','blaEC', esbl_data$Name.of.closest.sequence)
esbl_data$Name.of.closest.sequence = gsub('.*CTX','blaCTX', esbl_data$Name.of.closest.sequence)
esbl_data$Name.of.closest.sequence = gsub('.*SHV','blaSHV', esbl_data$Name.of.closest.sequence)
esbl_data$Name.of.closest.sequence = gsub('.*TEM','blaTEM', esbl_data$Name.of.closest.sequence)
# subsetting esbl_data to only those isolates in this study
esbl_data = esbl_data[esbl_data$Name %in% tree$tip.label,]

# setting up data frame for esbl plotting
ordered_esbl = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = (length(table(esbl_data$Name.of.closest.sequence)) + length(table(esbl_data$esbl_class)) + 1)), stringsAsFactors = F)
colnames(ordered_esbl) = c('Sample', names(table(esbl_data$Name.of.closest.sequence)), names(table(esbl_data$esbl_class)))
ordered_esbl$Sample = tree$tip.label
rownames(ordered_esbl) = tree$tip.label
for(i in 1:nrow(ordered_esbl)){
  for(j in 2:ncol(ordered_esbl)){
    if(colnames(ordered_esbl)[j] %in% esbl_data$Name.of.closest.sequence[esbl_data$Name == ordered_esbl$Sample[i]] | colnames(ordered_esbl)[j] %in% esbl_data$esbl_class[esbl_data$Name == ordered_esbl$Sample[i]]){
      ordered_esbl[i,j] = colnames(ordered_esbl)[j]
    }else{
      ordered_esbl[i,j] = 'not present'
    }
  }
}

# making a simplified column so I can plot blaEC as one ring and the rest as another ring since blaEC is the only one overlapping with others
ordered_esbl$blaCTX = 'not present'
ordered_esbl$blaCTX[ordered_esbl$`blaCTX-M-14` != 'not present'] = 'blaCTX-M-14'
ordered_esbl$blaCTX[ordered_esbl$`blaCTX-M-15` != 'not present'] = 'blaCTX-M-15'
ordered_esbl$blaCTX[ordered_esbl$`blaCTX-M-27` != 'not present'] = 'blaCTX-M-27'
ordered_esbl$blaCTX[ordered_esbl$`blaCTX-M-3` != 'not present'] = 'blaCTX-M-3'
ordered_esbl$blaCTX[ordered_esbl$`blaCTX-M-55` != 'not present'] = 'blaCTX-M-55'
ordered_esbl$blaSHV = 'not present'
ordered_esbl$blaSHV[ordered_esbl$`blaSHV-12` != 'not present'] = 'blaSHV-12'
ordered_esbl$blaEC = 'not present'
ordered_esbl$blaEC[ordered_esbl$`blaEC-15` != 'not present'] = 'blaEC-15'
ordered_esbl$blaEC[ordered_esbl$`blaEC-19` != 'not present'] = 'blaEC-19'

## generating a quinolone resistance-specific dataframe
quinoloner_data = amrfinder[grep('QUINOLONE', amrfinder$Class),]

quinoloner_data$Quinolone = 'not present'
quinoloner_data$Quinolone[quinoloner_data$Class == 'AMINOGLYCOSIDE/QUINOLONE' | quinoloner_data$Class == 'QUINOLONE'] = quinoloner_data$Gene.symbol[quinoloner_data$Class == 'AMINOGLYCOSIDE/QUINOLONE' | quinoloner_data$Class == 'QUINOLONE']
quinoloner_data$simple_quinolone = gsub('qnr.*','qnr',quinoloner_data$Quinolone)
quinoloner_data$simple_quinolone = gsub('gyrA.*','gyrA', quinoloner_data$simple_quinolone)
quinoloner_data$simple_quinolone = gsub('parE.*','parE', quinoloner_data$simple_quinolone)
quinoloner_data$simple_quinolone = gsub('parC.*','parC',quinoloner_data$simple_quinolone)

# ordering ordered_quinolone by tree$tip.label
ordered_quinoloner = data.frame(matrix(NA, nrow = length(tree$tip.label), ncol = 6), stringsAsFactors = F)
colnames(ordered_quinoloner) = c('Isolate','aac','qnr','gyrA','parC','parE')
ordered_quinoloner$Isolate = tree$tip.label
ordered_quinoloner$aac = 'not present'
ordered_quinoloner$qnr = 'not present'
ordered_quinoloner$gyrA = 'not present'
ordered_quinoloner$parC = 'not present'
ordered_quinoloner$parE = 'not present'

for(i in 1:nrow(ordered_quinoloner)){
  if("aac(6')-Ib-cr5" %in% quinoloner_data$Gene.symbol[quinoloner_data$Name == ordered_quinoloner$Isolate[i]]){
    ordered_quinoloner$aac[i] = "aac(6')-Ib-cr5"
  }
  if('qnr' %in% quinoloner_data$simple_quinolone[quinoloner_data$Name == ordered_quinoloner$Isolate[i]]){
    ordered_quinoloner$qnr[i] = quinoloner_data$Gene.symbol[quinoloner_data$Name == ordered_quinoloner$Isolate[i] & quinoloner_data$simple_quinolone == 'qnr']
  }
  if('gyrA' %in% quinoloner_data$simple_quinolone[quinoloner_data$Name == ordered_quinoloner$Isolate[i]]){
    ordered_quinoloner$gyrA[i] = 'gyrA (D87N, D87Y, S83L)'
  }
  if('parC' %in% quinoloner_data$simple_quinolone[quinoloner_data$Name == ordered_quinoloner$Isolate[i]]){
    ordered_quinoloner$parC[i] = 'parC (A56T, E84G, E84V, S57T, S80I, S80R)'
  }
  if('parE' %in% quinoloner_data$simple_quinolone[quinoloner_data$Name == ordered_quinoloner$Isolate[i]]){
    ordered_quinoloner$parE[i] = 'parE (E460K, I529L, L416F, L445H, S458A)'
  }
}
rownames(ordered_quinoloner) = ordered_quinoloner$Isolate
ordered_quinoloner = ordered_quinoloner[,2:ncol(ordered_quinoloner), drop = F]


## Generating ordered antimicrobial susceptibility testing (AST) dataframe
ast_columns = setdiff(colnames(ordered_metadata), c('Sample','organism','media','Isolate_source','Study'))
ordered_ast = ordered_metadata[,colnames(ordered_metadata) %in% ast_columns]

# separating AST data into ESBL-related and quinolone-related
esbl_columns = c('Cefiderocol','Ceftazolin','Cefotetan','Ceftriaxone','Cefepime','Meropenem','Imipenem','Piperacillin.Tazo','Ceftolozane.Tazo','Ceftazidime.Avi','Ampicillin.Sul')
esbl_asts = ordered_ast[,colnames(ordered_ast) %in% esbl_columns]

quinolone_columns = c('Ciprofloxacin')
quinolone_asts = ordered_ast[,colnames(ordered_ast) %in% quinolone_columns, drop = F]

## Generating SNP matrix and then an adjacency matrix for SNP clusters
snp_matrix = matrix(NA, nrow = length(tree$tip.label), ncol = length(tree$tip.label))
colnames(snp_matrix) = tree$tip.label
rownames(snp_matrix) = tree$tip.label
snps = snps[snps$sample1 %in% tree$tip.label,]

for(i in 1:nrow(snp_matrix)){
  for(j in 1:ncol(snp_matrix)){
    if(rownames(snp_matrix)[i] == colnames(snp_matrix)[j]){
      snp_matrix[i,j] = 0
    }else{
      snp_matrix[i,j] = snps$core_snps[snps$sample1 == rownames(snp_matrix)[i] & snps$sample2 == colnames(snp_matrix)[j]]
    }
  }
}

adj_matrix = matrix(0, nrow = nrow(snp_matrix), ncol = nrow(snp_matrix))
colnames(adj_matrix) = rownames(snp_matrix)
rownames(adj_matrix) = rownames(snp_matrix)

for(i in 1:nrow(snps)){
  if(snps$core_snps[i] < 22){
    adj_matrix[snps$sample1[i], snps$sample2[i]] = 23 - snps$core_snps[i]
    adj_matrix[snps$sample2[i], snps$sample1[i]] = 23 - snps$core_snps[i]
  }
}

# reading in dataframe of AMRFinder hits by AMR class
ordered_allamr_class = read.csv('./manuscript/supp_tables/allARGs_by_class.csv', stringsAsFactors = F)
rownames(ordered_allamr_class) = ordered_allamr_class$Isolate
ordered_allamr_class = ordered_allamr_class[,-1]

#### Figure 1c ####
# plotting histogram of MLSTs in each patient ward and/or MLSTs for each culture condition on x axis
#pdf(paste0('./figures/',Sys.Date(),'_BenedictAgee_Fig1c.pdf'), width = 15, height = 11)
ggplot(ordered_metadata_mlst, aes(x = media, fill = MLST))+
  geom_bar(position = 'stack',stat = 'count')+
  labs(fill = 'MLST', y = 'Count of isolates', x = element_blank())+
  geom_text(
    aes(label = MLST),  # Label by gene
    stat = "count",
    position = position_stack(vjust = 0.5),  # Place text in the middle of each segment
    color = "white", size = 3  # Adjust label size and color
  )+
  theme_bw()+
  scale_fill_manual(values = mlst_colors, breaks = names(mlst_colors))+
  theme(axis.text.x = element_text(size = 14), axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12),
        strip.text = element_text(size = 16, face = 'bold'))
#dev.off()


#### Figure 3a-c ####
## 3a: SNP clusters; SNP numbers and sequence type added in InkScape
# making colors df
colors_df = data.frame(rownames(adj_matrix), stringsAsFactors = F)
colnames(colors_df) = 'Sample'
# Isolate-level year of collection will not be published; generating empty dataframe here so the following code will run for others
year_df = data.frame(tree$tip.label, rep(c('2017','2018','2019')))
colnames(year_df) = c('Sample','Year')
colors_df = full_join(colors_df, year_df)
year_color_df = data.frame(year_colors)
year_color_df$Year = rownames(year_color_df)
colors_df = full_join(colors_df, year_color_df)
node_colors = colors_df$year_colors

#pdf(file = paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig3a.pdf'), width = 20, height = 20)
qgraph(
  adj_matrix,
  layout = 'spring',
  color = node_colors,
  border.color = node_colors,
  shape = 'circle',
  border.width = 4,
  edge.color = 'grey30',
  cut = 0,
  vsize = 2,
  label.cex = 1,
  label.color = 'white',
  directed = F
)
#dev.off()


## 3b: Tree + AMR + AST + patient exposures; 'Beta-lactam' and 'Quinolone' labels added in InkScape
# adding mlst to treeplot
media_df = data.frame(ordered_metadata$Sample, ordered_metadata$media)
colnames(media_df) = c('Sample','Media')
treeplot = ggtree(tree) %<+% media_df
# generating mlst-tip treeplot
tip_treeplot = treeplot + geom_tiplab(align = TRUE, size = 0)+ geom_tippoint(aes(color = Media), size = 6) +
  scale_color_manual(values = media_colors, breaks = names(media_colors))
tip_treeplot

# adding MLST to treeplot
mlst_df = data.frame(ordered_mlst$MLST)
colnames(mlst_df) = 'MLST'
rownames(mlst_df) = ordered_mlst$Sample

mlst_media_tree = gheatmap(tip_treeplot, mlst_df, colnames_angle = 90, colnames_offset_y = -1,
                           offset = 0, width = 0.035, font.size = 6)+
  labs(fill = 'MLST')+
  scale_fill_manual(values = mlst_colors, breaks = names(mlst_colors))+ggtree::vexpand(0.2,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'),
        legend.box = 'horizontal')

mlst_media_tree2 = mlst_media_tree + new_scale_fill()

mlst_media_esbl_tree = gheatmap(mlst_media_tree2, ordered_esbl[,c('blaCTX','blaSHV','blaEC')], colnames_angle = 90, font.size = 6, colnames_offset_y = -1.5,
                                offset = 0.0015, width = 0.1)+
  labs(fill = 'ESBL gene')+
  scale_fill_manual(values = esbl_colors, breaks = names(esbl_colors))

mlst_media_esbl_tree2 = mlst_media_esbl_tree + new_scale_fill()

mlst_media_esbl_esblast_tree = gheatmap(mlst_media_esbl_tree2, esbl_asts, colnames_angle = 90, font.size = 6, colnames_offset_y = -2.5,
                                        offset = 0.005, width = 0.45)+
  labs(fill = 'Susceptibility')+
  scale_fill_manual(values = ast_colors, breaks = names(ast_colors))

mlst_media_esbl_esblast_tree2 = mlst_media_esbl_esblast_tree + new_scale_fill()

# Patient antibiotic exposure data will not be published; generating empty dataframe here so the following code will run for others
esbl_abx = data.frame(rep(NA, length(tree$tip.label)),rep(NA, length(tree$tip.label)),rep(NA, length(tree$tip.label)))
rownames(esbl_abx) = tree$tip.label  
colnames(esbl_abx) = c('Toy_a','Toy_b','Toy_c')

mlst_media_esbl_esblastabx_tree = gheatmap(mlst_media_esbl_esblast_tree2, esbl_abx, colnames_angle = 90, font.size = 6, colnames_offset_y = -2.5,
                                        offset = 0.0195, width = 0.12)+
    labs(fill = 'Patient antibiotic exposure in last 12 weeks')+
    scale_fill_manual(values = abx_colors, breaks = names(abx_colors))

mlst_media_esbl_esblastabx_tree2 = mlst_media_esbl_esblastabx_tree + new_scale_fill()

mlst_media_esbl_esblastabx_quin_tree = gheatmap(mlst_media_esbl_esblastabx_tree2, ordered_quinoloner, colnames_angle = 90, font.size = 6, colnames_offset_y = -1,
                                                offset = 0.0243, width = 0.18)+
  labs(fill = 'Quinolone resistance gene/mutation')+
  scale_fill_manual(values = quinolone_colors, breaks = names(quinolone_colors))+
  theme(legend.box = 'vertical')

mlst_media_esbl_esblastabx_quin_tree2 = mlst_media_esbl_esblastabx_quin_tree  + new_scale_fill()

mlst_media_esbl_esblastabx_quinast_tree = gheatmap(mlst_media_esbl_esblastabx_quin_tree2, quinolone_asts, colnames_angle = 90, font.size = 6, colnames_offset_y = -1.7,
                                                   offset = 0.0301, width = 0.03)+
  labs(fill = 'Susceptibility')+
  scale_fill_manual(values = ast_colors, breaks = names(ast_colors))

mlst_media_esbl_esblastabx_quinast_tree2 = mlst_media_esbl_esblastabx_quinast_tree + new_scale_fill()

# Patient antibiotic exposure data will not be published; generating empty dataframe here so the following code will run for others
quinolone_abx = data.frame(rep(NA, length(tree$tip.label)))
rownames(quinolone_abx) = tree$tip.label
colnames(quinolone_abx) = 'Toy_d'

mlst_media_esbl_esblastabx_quinastabx_tree = gheatmap(mlst_media_esbl_esblastabx_quinast_tree2, quinolone_abx, colnames_angle = 90, font.size = 6, colnames_offset_y = -1.5,
                                                      offset = 0.0313, width = 0.03)+
  labs(fill = 'Patient antibiotic exposure in last 12 weeks')+
  scale_fill_manual(values = abx_colors, breaks = names(abx_colors))+
  theme(legend.box = 'vertical')

#pdf(file = paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig3b.pdf'), width = 20, height = 20)
mlst_media_esbl_esblastabx_quinastabx_tree
#dev.off()

## 3c: Quinolone + ESBL bubble plot and histograms; histograms aligned in InkScape
# ESBL histogram
esbl_hist_data = esbl_data
esbl_hist_data = melt(esbl_hist_data, id = 'Name')
esbl_hist_data = esbl_hist_data[esbl_hist_data$value != 'not present',]
esbl_hist_data = esbl_hist_data[esbl_hist_data$value %in% c('blaCTX-M-3','blaCTX-M-14','blaCTX-M-15','blaCTX-M-27',
                                                            'blaCTX-M-55','blaEC-15','blaEC-19','blaSHV-12'),]

#pdf(paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig3c_ESBL_histogram.pdf'), width = 12, height = 4)
ggplot(esbl_hist_data, aes(x = factor(value, levels = c('blaCTX-M-3','blaCTX-M-14','blaCTX-M-15','blaCTX-M-27',
                                                        'blaCTX-M-55','blaEC-15','blaEC-19','blaSHV-12'))))+
  geom_histogram(stat='count', fill = 'violetred3', alpha = 0.5)+
  labs(x = 'ESBL resistance gene', y = 'Isolate count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('ESBL genes')
#dev.off() 

# quinolone histogram
quin_hist_data = ordered_quinoloner
quin_hist_data$Name = rownames(quin_hist_data)
quin_hist_data = melt(quin_hist_data, id = 'Name')
quin_hist_data = quin_hist_data[quin_hist_data$value != 'not present',]
quin_hist_data = quin_hist_data[quin_hist_data$value %in% c('gyrA (D87N, D87Y, S83L)','parC (A56T, E84G, E84V, S57T, S80I, S80R)',
                                                            'parE (E460K, I529L, L416F, L445H, S458A)','qnrB2',
                                                            'qnrS1',"aac(6')-Ib-cr5"),]

#pdf(paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig3c_quinolone_histogram.pdf'), width = 15, height = 8)
ggplot(quin_hist_data, aes(x = factor(value, levels = c('gyrA (D87N, D87Y, S83L)','parC (A56T, E84G, E84V, S57T, S80I, S80R)',
                                                        'parE (E460K, I529L, L416F, L445H, S458A)','qnrB2',
                                                        'qnrS1',"aac(6')-Ib-cr5"))))+
  geom_histogram(stat='count', fill = 'violetred3', alpha = 0.5)+
  labs(x = 'Quinolone resistance gene/mutation', y = 'Isolate count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('Quinolone resistance genes/mutations')
#dev.off()

# bubble plot
colnames(esbl_hist_data) = c('Name','element','ESBL_gene')
colnames(quin_hist_data) = c('Name','gene','QuinoloneR_element')
esbl_quin_count_data = full_join(esbl_hist_data, quin_hist_data, relationship = 'many-to-many')
esbl_quin_count_data$duplicate_check = paste0(esbl_quin_count_data$Name, esbl_quin_count_data$ESBL_gene, esbl_quin_count_data$QuinoloneR_element)
esbl_quin_count_data = esbl_quin_count_data[!duplicated(esbl_quin_count_data$duplicate_check),]
count_data = data.frame(table(esbl_quin_count_data[,c('ESBL_gene','QuinoloneR_element')]))



#pdf(paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig3b_bubbleplot.pdf'), width = 15, height = 12)
ggplot(count_data, aes(x = factor(ESBL_gene, levels = c('blaCTX-M-3','blaCTX-M-14','blaCTX-M-15','blaCTX-M-27',
                                                      'blaCTX-M-55','blaEC-15','blaEC-19','blaSHV-12')), 
                       y = factor(QuinoloneR_element, levels = c("aac(6')-Ib-cr5",'qnrS1','qnrB2','qnrB4','qnrB19',
                                                      'parE (E460K, I529L, L416F, L445H, S458A)','parC (A56T, E84G, E84V, S57T, S80I, S80R)',
                                                      'gyrA (D87N, D87Y, S83L)')),
                       size = Freq))+
  geom_point(color = 'violetred3', alpha = 0.85)+
  scale_size_area(max_size = 30)+
  geom_text(aes(label = Freq),  # Label by gene  # Place text in the middle of each segment
            color = "white", size = 4, fontface = 'bold')+  
  labs(y = 'Quinolone resistance gene/mutation', x = 'ESBL gene', size = 'Count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('ESBL genes and quinolone resistance genes/mutations')
#dev.off()


#### Figure 4 a - d ####
## 4a, ESBL genes in ST131 vs. non-ST131 isolates
table_esbl_data = data.frame(c('0','1','2'),c('24','15','0'),c('12','8','16'))
colnames(table_esbl_data) = c('Count of ESBL genes per isolate','ST131','Other STs')

kbl(table_esbl_data) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left")

## 4b, quinolone resistance elements in ST131 vs. non-ST131 isolates
table_quin_data = data.frame(c('0','1','2','3','4'),c('0','0','0','36','3'),c('4','5','4','18','5'))
colnames(table_quin_data) = c('Count of quinolone resistance\n genes/mutations per isolate','ST131','Other STs')

kbl(table_quin_data) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left")

## 4c, ST131 ESBL gene and quinolone resistance elements co-occurrence bubble plot; histograms aligned in InkScape
st131_list = rownames(mlst_df)[mlst_df$MLST == '131']
st131_count_data = esbl_quin_count_data[esbl_quin_count_data$Name %in% st131_list,]
st131_esbl_hist = esbl_hist_data[esbl_hist_data$Name %in% st131_list,]
st131_esbl_hist$duplicatecheck = paste0(st131_esbl_hist$Name, st131_esbl_hist$ESBL_gene)
st131_esbl_hist = st131_esbl_hist[!duplicated(st131_esbl_hist$duplicatecheck),]
st131_quin_hist = quin_hist_data[quin_hist_data$Name %in% st131_list,]
st131_quin_hist$duplicatecheck = paste0(st131_quin_hist$Name, st131_quin_hist$QuinoloneR_element)
st131_quin_hist = st131_quin_hist[!duplicated(st131_quin_hist$duplicatecheck),]

st131_count_data = data.frame(table(st131_count_data[,c('ESBL_gene','QuinoloneR_element')]))


#pdf(paste0('./figures/', Sys.Date(),'_ESBLCiproRe_esbl_vs_cipro_bubbleplot_st131.pdf'), width = 12)
ggplot(st131_count_data, aes(x = factor(ESBL_gene, levels = c('blaCTX-M-3','blaCTX-M-14','blaCTX-M-15','blaCTX-M-27',
                                                            'blaCTX-M-55','blaEC-15','blaEC-19','blaSHV-12')), 
                             y = factor(QuinoloneR_element, levels = c("aac(6')-Ib-cr5",'qnrS1','qnrB2','qnrB4','qnrB19',
                                                            'parE (E460K, I529L, L416F, L445H, S458A)','parC (A56T, E84G, E84V, S57T, S80I, S80R)',
                                                            'gyrA (D87N, D87Y, S83L)')),
                             size = Freq))+
  geom_point(color = '#6984B3', alpha = 0.85)+
  scale_size_area(max_size = 30)+
  geom_text(aes(label = Freq),  # Label by gene  # Place text in the middle of each segment
            color = "white", size = 4, fontface = 'bold')+  
  labs(y = 'Quinolone resistance gene/mutation', x = 'ESBL gene', size = 'Count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('ST131 ESBL genes and quinolone resistance genes/mutations')
#dev.off()

### 4c ST131 ESBL histogram
#pdf(paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig4c_esbl_histogram_st131.pdf'), width = 12, height = 4)
ggplot(st131_esbl_hist, aes(x = factor(ESBL_gene, levels = c('blaCTX-M-3','blaCTX-M-14','blaCTX-M-15','blaCTX-M-27',
                                                         'blaCTX-M-55','blaEC-15','blaEC-19','blaSHV-12'))))+
  geom_histogram(stat='count', fill = '#6984B3', alpha = 0.5)+
  labs(x = 'ESBL resistance gene', y = 'Isolate count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('ESBL genes')
#dev.off() 

### 4c ST131 quinolone histogram
#pdf(paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig4c_quinolone_histogram_st131.pdf'), width = 12, height = 8)
ggplot(st131_quin_hist, aes(x = factor(QuinoloneR_element, levels = c('gyrA (D87N, D87Y, S83L)','parC (A56T, E84G, E84V, S57T, S80I, S80R)',
                                                         'parE (E460K, I529L, L416F, L445H, S458A)','qnrB19','qnrB4','qnrB2',
                                                         'qnrS1',"aac(6')-Ib-cr5"))))+
  geom_histogram(stat='count', fill = '#6984B3', alpha = 0.5)+
  labs(x = 'Quinolone resistance gene/mutation', y = 'Isolate count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('Quinolone resistance genes/mutations')
#dev.off() 

## 4d, non-ST131 ESBL gene and quinolone resistance elements co-ocurrence bubble plot; histograms aligned in InkScape
nonst131_list = rownames(mlst_df)[mlst_df$MLST != '131']
nonst131_count_data = esbl_quin_count_data[esbl_quin_count_data$Name %in% nonst131_list,]
nonst131_esbl_hist = esbl_hist_data[esbl_hist_data$Name %in% nonst131_list,]
nonst131_esbl_hist$duplicatecheck = paste0(nonst131_esbl_hist$Name, nonst131_esbl_hist$ESBL_gene)
nonst131_esbl_hist = nonst131_esbl_hist[!duplicated(nonst131_esbl_hist$duplicatecheck),]
nonst131_quin_hist = quin_hist_data[quin_hist_data$Name %in% nonst131_list,]
nonst131_quin_hist$duplicatecheck = paste0(nonst131_quin_hist$Name, nonst131_quin_hist$QuinoloneR_element)
nonst131_quin_hist = nonst131_quin_hist[!duplicated(nonst131_quin_hist$duplicatecheck),]

nonst131_count_data = data.frame(table(nonst131_count_data[,c('ESBL_gene','QuinoloneR_element')]))


#pdf(paste0('./figures/', Sys.Date(),'_ESBLCiproRe_esbl_vs_cipro_bubbleplot_nonst131.pdf'), width = 12)
ggplot(nonst131_count_data, aes(x = factor(ESBL_gene, levels = c('blaCTX-M-3','blaCTX-M-14','blaCTX-M-15','blaCTX-M-27',
                                                              'blaCTX-M-55','blaEC-15','blaEC-19','blaSHV-12')), 
                             y = factor(QuinoloneR_element, levels = c("aac(6')-Ib-cr5",'qnrS1','qnrB2','qnrB4','qnrB19',
                                                                       'parE (E460K, I529L, L416F, L445H, S458A)','parC (A56T, E84G, E84V, S57T, S80I, S80R)',
                                                                       'gyrA (D87N, D87Y, S83L)')),
                             size = Freq))+
  geom_point(color = 'aquamarine4', alpha = 0.85)+
  scale_size_area(max_size = 30)+
  geom_text(aes(label = Freq),  # Label by gene  # Place text in the middle of each segment
            color = "white", size = 4, fontface = 'bold')+  
  labs(y = 'Quinolone resistance gene/mutation', x = 'ESBL gene', size = 'Count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('ST131 ESBL genes and quinolone resistance genes/mutations')
#dev.off()


### 4d non-ST131 ESBL histogram
#pdf(paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig4c_esbl_histogram_nonst131.pdf'), width = 12, height = 4)
ggplot(nonst131_esbl_hist, aes(x = factor(ESBL_gene, levels = c('blaCTX-M-3','blaCTX-M-14','blaCTX-M-15','blaCTX-M-27',
                                                                'blaCTX-M-55','blaEC-15','blaEC-19','blaSHV-12'))))+
  geom_histogram(stat='count', fill = 'aquamarine4', alpha = 0.5)+
  labs(x = 'ESBL resistance gene', y = 'Isolate count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('ESBL genes')
#dev.off() 

#### 4d non-ST131 quinolone histogram
#pdf(paste0('./figures/', Sys.Date(),'_BenedictAgee_CAEC_Fig4c_quinolone_histogram_st131.pdf'), width = 12, height = 8)
ggplot(nonst131_quin_hist, aes(x = factor(QuinoloneR_element, levels = c('gyrA (D87N, D87Y, S83L)','parC (A56T, E84G, E84V, S57T, S80I, S80R)',
                                                                         'parE (E460K, I529L, L416F, L445H, S458A)','qnrB19','qnrB4','qnrB2',
                                                                         'qnrS1',"aac(6')-Ib-cr5"))))+
  geom_histogram(stat='count', fill = 'aquamarine4', alpha = 0.5)+
  labs(x = 'Quinolone resistance gene/mutation', y = 'Isolate count')+
  theme_bw()+
  theme(axis.text.x = element_text(size = 14, angle = 90), axis.text.y = element_text(size = 14),
        axis.ticks = element_blank(),
        title = element_text(size = 16, face = 'bold'), legend.text = element_text(size = 12))+
  ggtitle('Quinolone resistance genes/mutations')
#dev.off() 



#### Supplemental Figure 1 ####
## core-genome SNPs histograms; histograms combined in InkScape
#pdf(paste0('./figures/',Sys.Date(),'_BenedictAgee_SupFig1_allSNPs_histogram.pdf'), width = 8, height = 8)
ggplot(snps, aes(x = core_snps))+
  geom_histogram(binwidth = 5, fill = natparks.pals('Arches',6)[3],color = natparks.pals('Arches',6)[4])+
  #geom_vline(xintercept = 22, lty = 2, color = 'navy', size = 1.5)+
  #xlim(c(0,250))+
  #ylim(c(0,6))+
  theme_bw()+
  labs(x = 'Pairwise core genome SNPs', y = 'Isolate count')+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'), legend.text = element_text(size = 12),
        title = element_text(size = 16, face = 'bold'))+
  ggtitle('Core genome SNPs')
#dev.off()

#pdf(paste0('./figures/',Sys.Date(),'_BenedictAgee_SupFig1_under250SNPs_histogram.pdf'), width = 8, height = 8)
ggplot(snps, aes(x = core_snps))+
  geom_histogram(binwidth = 1, fill = natparks.pals('Arches',6)[3],color = natparks.pals('Arches',6)[4])+
  geom_vline(xintercept = 22, lty = 2, color = 'navy', size = 1.5)+
  xlim(c(0,250))+
  ylim(c(0,30))+
  theme_bw()+
  labs(x = 'Pairwise core genome SNPs', y = 'Isolate count')+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'), legend.text = element_text(size = 12),
        title = element_text(size = 16, face = 'bold'))+
  ggtitle('Core genome SNPs')
#dev.off()

#### Supplemental Figure 2 ####
## ARG-annotated maximum likelihood tree; antibiotic class labels added in InkScape
thick_mlst_media_tree = gheatmap(tip_treeplot, mlst_df, colnames_angle = 90, colnames_offset_y = -2.5,
                                 offset = 0.003, width = 1, font.size = 4)+
  labs(fill = 'MLST')+
  scale_fill_manual(values = mlst_colors, breaks = names(mlst_colors))+ggtree::vexpand(0.2,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))

thick_mlst_media_tree2 = thick_mlst_media_tree + new_scale_fill()

allamr_class_mlst_media_tree = gheatmap(thick_mlst_media_tree2, ordered_allamr_class, colnames_angle = 90, font.size = 4, colnames_offset_y = -5,
                                        offset = 0.05, width = 50)+
  labs(fill = 'Resistance class')+
  scale_fill_manual(values = amrfinder_colors, breaks = names(amrfinder_colors))

# to combine with plasmid results
#pdf(file = paste0('./manuscript/figures/', Sys.Date(),'_BenedictAgee_CAEC_SupFig2.pdf'), width = 30, height = 20)
allamr_class_mlst_media_tree
#dev.off()

#### Supplemental Figure 3 ####
## Maximum likelihood tree of isolates from within-patient co-colonization cases; patient number added in InkScape
  ## Patient list will not be published, omitting the patient ring from the following code
# subsetting tree to those with high SNPs from same patient; the data will not be published, so assembling list manually
samept_difec = c('EC23','EC24','EC37','EC38','EC121','EC122','EC71','EC72','EC12','EC13','EC3','EC4','EC43','EC44','EC16','EC17',
                 'EC89','EC90','EC55','EC56','EC68','EC69')
highsnps_tree = drop.tip(tree, setdiff(tree$tip.label, samept_difec))
highsnps_tree = midpoint.root(highsnps_tree)
# adding mlst to circle treeplot
highsnps_circle_treeplot = ggtree(highsnps_tree, layout = 'fan', open.angle = 15) %<+% media_df
# generating mlst-tip treeplot
highsnps_tip_circle_treeplot = highsnps_circle_treeplot + geom_tiplab(align = TRUE, size = 0)+geom_tippoint(aes(color = Media), size = 3.7) +
  scale_color_manual(values = media_colors, breaks = names(media_colors))

highsnps_media_mlst_tree = gheatmap(highsnps_tip_circle_treeplot, mlst_df, colnames_angle = 90,
                                    offset = 0.001, width = 0.09, font.size = 5)+
  labs(fill = 'MLST')+
  scale_fill_manual(values = mlst_colors, breaks = names(mlst_colors))+ggtree::vexpand(0.2,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))

highsnps_media_mlst_tree2 = highsnps_media_mlst_tree + new_scale_fill()

#pdf(paste0('./manuscript/figures/',Sys.Date(),'_BenedictAgee_CAEC_SupFig3.pdf'), width = 20, height = 18)
gheatmap(highsnps_media_mlst_tree2, ordered_amr_classes[,2:ncol(ordered_amr_classes)], colnames_angle = 90,
         offset = 0.02, width = 1.8, font.size = 3)+
  labs(fill = 'Resistance class')+
  scale_fill_manual(values = amrfinder_colors, breaks = names(amrfinder_colors))+ggtree::vexpand(0.2,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))
#dev.off()

#### Supplemental Figure 4 ####
## Maximum likelihood tree of same-strain isolates from different patients; patient number added in InkScape for true figure
## Patient list will not be published, omitting the patient ring from the following code
lowsnps_tree = drop.tip(tree, highsnps_tree$tip.label)
lowsnps_tree = midpoint.root(lowsnps_tree)
# adding mlst to circle treeplot
lowsnps_circle_treeplot = ggtree(lowsnps_tree, layout = 'fan', open.angle = 15) %<+% media_df
# generating mlst-tip treeplot
lowsnps_tip_circle_treeplot = lowsnps_circle_treeplot + geom_tiplab(align = TRUE, size = 0)+geom_tippoint(aes(color = Media), size = 3.7) +
  scale_color_manual(values = media_colors, breaks = names(media_colors))

lowsnps_media_mlst_tree = gheatmap(lowsnps_tip_circle_treeplot, mlst_df, colnames_angle = 90,
                                    offset = 0.001, width = 0.09, font.size = 5)+
  labs(fill = 'MLST')+
  scale_fill_manual(values = mlst_colors, breaks = names(mlst_colors))+ggtree::vexpand(0.2,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))

lowsnps_media_mlst_tree2 = lowsnps_media_mlst_tree + new_scale_fill()

#pdf(paste0('./manuscript/figures/',Sys.Date(),'_BenedictAgee_CAEC_SupFig4.pdf'), width = 20, height = 18)
gheatmap(lowsnps_media_mlst_tree2, ordered_amr_classes[,2:ncol(ordered_amr_classes)], colnames_angle = 90,
         offset = 0.02, width = 1.8, font.size = 3)+
  labs(fill = 'Resistance class')+
  scale_fill_manual(values = amrfinder_colors, breaks = names(amrfinder_colors))+ggtree::vexpand(0.2,-1)+
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14, face = 'bold'))
#dev.off()
