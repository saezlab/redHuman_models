library(readr)

bigg_Recon3D_chemicals <- 
  read_delim("Dropbox/redHuman_models/bigg_Recon3D.chemicals.tsv",
             "\t", 
             escape_double = FALSE, 
             trim_ws = TRUE)

bigg_Recon3D_chemicals <- bigg_Recon3D_chemicals[-1,]

xrefs_list <- apply(bigg_Recon3D_chemicals,1,function(x){
  metab <- x
  xrefs <-strsplit(metab[7],";")[[1]]
  hmdb <- xrefs[grepl("hmdb",xrefs)]
  kegg <- xrefs[grepl("kegg.compound",xrefs)]

  if(length(hmdb) > 0 & length(kegg) > 0)
  {
    xrefs_df <- merge(hmdb,kegg)
    names(xrefs_df) <- c("HMDB","KEGG")
    return(xrefs_df)
  }
  return(NA)
  
})

xrefs_df <- do.call(rbind,xrefs_list)
xrefs_df <- xrefs_df[complete.cases(xrefs_df),]
xrefs_df$HMDB <- gsub("hmdb:","",xrefs_df$HMDB)
xrefs_df$KEGG <- gsub("kegg.compound:","",xrefs_df$KEGG)

xrefs_df <- unique(xrefs_df)

library(ocean)

mapping_table <- mapping_table


all_pathways <- unique(recon2_redhuman$pathway)
sub_network <- model_to_pathway_sif(pathway_to_keep = all_pathways$X1)
metabs <- unique(c(sub_network$source, sub_network$target))
metabs <- metabs[grepl("cpd:",metabs)]
metabs <- unique(gsub("cpd:","",gsub("_[a-z]$","",metabs)))

xrefs_df_filtered <- xrefs_df[xrefs_df$KEGG %in% metabs,]
names(xrefs_df_filtered) <- names(mapping_table)

mapping_table <- mapping_table[!grepl("HMDB",mapping_table$metab),]

mapping_table <- rbind(mapping_table,xrefs_df_filtered)
save(mapping_table, file = "Dropbox/ocean/data/mapping_table.RData")
