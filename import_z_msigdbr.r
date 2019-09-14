MSigDBr_import<-function(gs_category){
  #funkcja pobierajaca z bazy MSigDB geneset z wybranej przez uzytkownika kategorii oraz zapisujaca
  #przetworzony geneset do pliku .gmt
  
  require(msigdbr)
  require(GSEABase)
  
  if (gs_category==1){
    Gene_Sets=msigdbr(species='Homo sapiens', category='H')
  }
  else if (gs_category==2){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C1')
  }
  else if (gs_category==3){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C2')
  }
  else if (gs_category==4){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C2', subcategory='CGP')
  }
  else if (gs_category==5){
    Gene_Sets=msigdbr(species='Homo sapiens',category='C2', subcategory='CP')
  }
  else if (gs_category==6){
    Gene_Sets=msigdbr(species='Homo sapiens',category='C2', subcategory='CP:BIOCARTA')
  }
  else if (gs_category==7){
    Gene_Sets=msigdbr(species='Homo sapiens',category='C2', subcategory='CP:KEGG')
  }
  else if (gs_category==8){
    Gene_Sets=msigdbr(species='Homo sapiens',category='C2', subcategory='CP:REACTOME')
  }
  else if (gs_category==9){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C3')
  }
  else if (gs_category==10){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C3', subcategory='MIR')
  }
  else if (gs_category==11){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C3', subcategory='TFT')
  }
  else if (gs_category==12){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C4')
  }
  else if (gs_category==13){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C4', subcategory='CGN')
  }
  else if (gs_category==14){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C4', subcategory='CM')
  }
  else if (gs_category==15){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C5')
  }
  else if (gs_category==16){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C5', subcategory='BP')
  }
  else if (gs_category==17){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C5', subcategory='CC')
  }
  else if (gs_category==18){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C5', subcategory='MF')
  }
  else if (gs_category==19){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C6')
  }
  else if (gs_category==20){
    Gene_Sets=msigdbr(species='Homo sapiens', category='C7')
  }
  
  geneset.names=unique(Gene_Sets$gs_name)
  # genesets<-vector("list",length(geneset.names))
  # for(i in 1:length(geneset.names)){
  #   nazwa<-geneset.names[i]
  #   geny<-Gene_Sets$gene_symbol[Gene_Sets$gs_name==nazwa]
  #   genesets[[i]] <- GeneSet(setName=nazwa, geneIds=as.character(geny), geneIdType=SymbolIdentifier())
  # }
  # 
  # Names_Col=GeneSetCollection(genesets)
  # toGmt(Names_Col, "GeneSetNames.gmt")
  
  genesets_e<-vector("list",length(geneset.names))
  for(i in 1:length(geneset.names)){
    nazwa<-geneset.names[i]
    geny<-Gene_Sets$entrez_gene[Gene_Sets$gs_name==nazwa]
    genesets_e[[i]] <- GeneSet(setName=nazwa, geneIds=as.character(geny), geneIdType=EntrezIdentifier())
  }
  
  Entrez_Col=GeneSetCollection(genesets_e)
  toGmt(Entrez_Col, "EntrezIDs.gmt")
}

