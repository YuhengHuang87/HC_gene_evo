#install.packages("ape");library(ape)
tr <- ape::read.tree(text='((D_miranda,(D_persimilis,D_pseudoobscura)),(D_obscura,(D_guanche,D_subobscura)),((((D_ananassae,D_pandora),(D_bipectinata,D_pseudoananassae_pseudoananassae)),D_setifemur),(((D_jambulina,(D_birchii,(D_bunnanda,D_serrata))),D_kikkawai),(((((D_biarmipes,(D_subpulchrella,D_suzukii)),(D_pseudotakahashii,D_takahashii)),(((D_erecta,(D_teissieri_273_3,(D_santomea,D_yakuba))),(((D_mauritiana,D_simulans),D_sechellia),D_melanogaster)),D_eugracilis)),(D_rhopaloa,(D_gunungcola,D_elegans))),D_ficusphila))));')
tip.label<-c("D_miranda","D_persimilis","D_pseudoobscura","D_obscura","D_guanche","D_subobscura","D_ananassae","D_pandora","D_bipectinata","D_pseudoananassae_pseudoananassae","D_setifemur","D_jambulina","D_birchii","D_bunnanda","D_serrata","D_kikkawai","D_biarmipes","D_subpulchrella","D_suzukii","D_pseudotakahashii","D_takahashii","D_erecta","D_teissieri_273_3","D_santomea","D_yakuba","D_mauritiana","D_simulans","D_sechellia","D_melanogaster","D_eugracilis","D_rhopaloa","D_gunungcola","D_elegans","D_ficusphila")

files <- list.files(path = "/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/Clustal_input_23species", pattern = "_species\\.txt$", full.names = TRUE)
for (file in files) {
  data <- read.table(file)
  #data_list[[file]] <- data
  to_drop <- tip.label[!tip.label %in% data$V1]
  missing <- setdiff(data$V1, tip.label)
  if (length(missing) > 0) warning(file, ": species not in the tree: ", paste(missing, collapse = ", "))
  focal_tree <- ape::drop.tip(tr, to_drop)
  # PAML expects an unrooted tree: pruning can leave a bifurcating root, so collapse it
  if (ape::Ntip(focal_tree) > 2 && ape::is.rooted(focal_tree)) focal_tree <- ape::unroot(focal_tree)
  # gene name = file name without "_species.txt" (keeps underscores, e.g. Su_var_2-10)
  id = sub("_species\\.txt$", "", basename(file))
  filename = paste("/dfs7/grylee/yuhenh3/Heterochromatic_repeat_Novagenes_17species/MACSE_PRANK_23species/species_num/",id,"_tree.txt", sep = "")
  tree_num = paste(nrow(data),"1",sep=" ")
  write(tree_num, file=filename,append=FALSE)
  #write(c("\n\n", focal_tree), file=filename,append=TRUE)
  write.tree(focal_tree, file=filename,append=TRUE)
}