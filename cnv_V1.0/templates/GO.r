library(GO.db);library('biomaRt');library('affy')
ids=<insertpidshere>
unimart=useMart("ENSEMBL_MART_ENSEMBL",host='www.ensembl.org')
martdataset = useDataset("<insertspecieshere>_gene_ensembl",mart=unimart)
filters = listFilters(martdataset);attributes = listAttributes(martdataset)
ids<-unlist(sapply(as.character(rownames(a)), function(x) do.call(rbind,strsplit(x, ";"))))
t<-getBM(attributes=c("<insertfilterhere>", "wikigene_name","go_id") , filters = "<insertfilterhere>", values = ids, mart = martdataset)
for (i in 1:nrow(a)){
  id<-as.vector(unlist(sapply(as.character(rownames(a)[i]), function(x) do.call(rbind,strsplit(x, ";")))))
  a$short.name[i]<-unique(t[t[,"<insertfilterhere>"]%in%id, "wikigene_name"])[1]
  a$accession[i]<-unique(t[t[,"<intertfilterhere>"]%in%id, "<insertfilterhere>"])[1]
}
nucccgo = GOCCOFFSPRING[["GO:0005634"]];nucccgo=c(c('GO:0005634'),nucccgo)
cytccgo = GOCCOFFSPRING[["GO:0005737"]];cytccgo=c(c('GO:0005737'),cytccgo)
mitccgo = GOCCOFFSPRING[["GO:0005739"]];mitccgo=c(c('GO:0005739'),mitccgo)
extccgo1 = GOCCOFFSPRING[["GO:0005576"]];extccgo2 = GOCCOFFSPRING[["GO:0031012"]];extccgo3 = GOCCOFFSPRING[["GO:0044421"]];extccgo4 = GOCCOFFSPRING[["GO:0044420"]];extccgo=union(union(union(extccgo1,extccgo2),extccgo3),extccgo4);extccgo=c(c('GO:0005576','GO:0031012',"GO:0044421","GO:0044420"),extccgo)

write.table(ids[which(a$accession %in% subset(t,go_id %in% nucccgo)$uniprotswissprot),],'templates/nucids.txt',sep=' ',quote=FALSE)
write.table(ids[which(a$accession %in% subset(t,go_id %in% cytccgo)$uniprotswissprot),],'templates/cytids.txt',sep=' ',quote=FALSE)
write.table(ids[which(a$accession %in% subset(t,go_id %in% extccgo)$uniprotswissprot),],'templates/extids.txt',sep=' ',quote=FALSE)
write.table(ids[which(a$accession %in% subset(t,go_id %in% mitccgo)$uniprotswissprot),],'templates/mitids.txt',sep=' ',quote=FALSE)
