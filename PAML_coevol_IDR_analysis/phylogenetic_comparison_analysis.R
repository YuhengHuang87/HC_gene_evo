required.libraries <- c("ape", "picante", "pez", "phytools",
                        "vegan", "adephylo", "phylobase",
                        "geiger", "mvMORPH", "OUwie", 
                        "hisse", "BAMMtools",
                        "phylosignal", "Biostrings",
                        "devtools",
                        "ggplot2",
                        "fulltext")

needed.libraries <- required.libraries[!(required.libraries %in% installed.packages()[,"Package"])]

if(length(needed.libraries)) install.packages(needed.libraries)

# Load all required libraries at once
lapply(required.libraries, require, character.only = TRUE)

library(dplyr)
library(treeio);library("tidytree")
#full<-read.nexus("/Users/yuhenghuang/Documents/Postdoc_UCI/selection_heterochromatin_pathway/drosophila_introgression_data/tree/schemeA.tre")

full<-read.iqtree("/Users/yuhenghuang/Documents/Postdoc_UCI/selection_heterochromatin_pathway/drosophila_introgression_data/tree/iqtree.tre")
full_tree<-as.phylo(full)
full_tree$tip.label
tree2 = tree_subset(full_tree, "D_ananassae", levels_back=5)  
tree2$tip.label
gTree <- ape::read.tree(text='(((((((D_melanogaster,(D_simulans,D_sechellia)),(D_erecta,D_yakuba)),D_eugracilis),((D_biarmipes, D_suzukii), D_takahashii)),(D_elegans,D_rhopaloa)),D_ficusphila),(D_serrata,D_kikkawai),(D_ananassae, D_bipectinata));')
gTree$tip.label
to_drop <- tree2$tip.label[!tree2$tip.label %in% gTree$tip.label]
focal_tree <- drop.tip(tree2, to_drop)
focal_tree$tip.label
focal_tree$edge.length

idr<-read.csv("/Users/yuhenghuang/Documents/Postdoc_UCI/selection_heterochromatin_pathway/IDR_Analysis_name_fixed.csv",header=T)

df<-subset(idr, select = -c(Gene.Class))

df<-read.table("/Users/yuhenghuang/Documents/Postdoc_UCI/selection_heterochromatin_pathway/fIDP_percentage_16_species.txt",header=T)
df<-df%>%filter(rowSums(across(where(is.numeric)))!=0)
df2 <- data.frame(t(df[-1]))

colnames(df2) <- df[, 1]

K_idr=c(); P_idr=c();
for (j in 1:ncol(df2)){
  
idr.K<- phylosig(focal_tree, 
                 setNames(as.numeric(df2[,j]), rownames(df2)),
                     method = "K", test = T)

K_idr=c(K_idr,idr.K$K);P_idr=c(P_idr,idr.K$P);
}
plot(K_idr,P_idr)
#outputdf <- data.frame(K_idr,P_idr)
outputdf <- data.frame(df,K_idr)
write.table(outputdf, file="/Users/yuhenghuang/Documents/Postdoc_UCI/selection_heterochromatin_pathway/IDR_K_test_random_genes.txt",sep="\t",eol="\n",row.names = F,col.names = F)
write.tree(focal_tree)


