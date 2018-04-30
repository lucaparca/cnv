library("ggplot2")
library('GO.db')
library('biomaRt')
library('fdrtool')
library('affy')
library('preprocessCore')
library('limma')
a=read.table('auto_cnv_input.txt',sep='\t',quote='',header=TRUE)
rownames(a)=a[,1]

unimart=useMart("ENSEMBL_MART_ENSEMBL",host='www.ensembl.org')
martdataset = useDataset("<insertspecieshere>_gene_ensembl",mart=unimart)
filters = listFilters(martdataset);attributes = listAttributes(martdataset)
ids<-unlist(sapply(as.character(rownames(a)), function(x) do.call(rbind,strsplit(x, ";"))))
t<-getBM(attributes=c("<insertfilterhere>", "wikigene_name","go_id") , filters = "<insertfilterhere>", values = ids, mart = martdataset)
for (i in 1:nrow(a)){
  id<-as.vector(unlist(sapply(as.character(rownames(a)[i]), function(x) do.call(rbind,strsplit(x, ";")))))
  a$short.name[i]<-unique(t[t[,"<insertfilterhere>"]%in%id, "wikigene_name"])[1]
  a$accession[i]<-unique(t[t[,"<insertfilterhere>"]%in%id, "<insertfilterhere>"])[1]
}

nucccgo = GOCCOFFSPRING[["GO:0005634"]];nucccgo=c(c('GO:0005634'),nucccgo)
cytccgo = GOCCOFFSPRING[["GO:0005737"]];cytccgo=c(c('GO:0005737'),cytccgo)
extccgo1 = GOCCOFFSPRING[["GO:0005576"]];extccgo2 = GOCCOFFSPRING[["GO:0031012"]];extccgo3 = GOCCOFFSPRING[["GO:0044421"]];extccgo4 = GOCCOFFSPRING[["GO:0044420"]];extccgo=union(union(union(extccgo1,extccgo2),extccgo3),extccgo4);extccgo=c(c('GO:0005576','GO:0031012',"GO:0044421","GO:0044420"),extccgo)
mitccgo = GOCCOFFSPRING[["GO:0005739"]];mitccgo=c(c('GO:0005739'),mitccgo)
perccgo = GOCCOFFSPRING[["GO:0005777"]];perccgo=c("GO:0005777",perccgo)
lysccgo = GOCCOFFSPRING[["GO:0005764"]];lysccgo=c("GO:0005764",lysccgo)
relccgo= GOCCOFFSPRING[["GO:0005783"]];relccgo=c("GO:0005783",relccgo)
golccgo = GOCCOFFSPRING[["GO:0005794"]];golccgo = c(c("GO:0005794"),golccgo)
memccgo = GOCCOFFSPRING[["GO:0005886"]];memccgo=c("GO:0005886",memccgo)
nmeccgo = GOCCOFFSPRING[["GO:0031965"]];nmeccgo=c("GO:0031965",nmeccgo)

compartments=list(nucccgo,cytccgo,extccgo,mitccgo,perccgo,lysccgo,relccgo,golccgo,memccgo,nmeccgo)
compartment_names=c('Nucleus','Cytoplasm','Extracellular','Mitochondrion','Peroxisome','Lysosome','Endoplasmic reticulum','Golgi','Cell membrane','Nuclear membrane')
out_table=c()
for (i in 1:length(compartment_names)){
  selected_proteins=which(a$accession %in% subset(t,go_id %in% unlist(compartments[i]))$uniprotswissprot)
  if (length(selected_proteins)<20) {next}
  grep_data=a[selected_proteins,]
  grep_data$compartment=compartment_names[i]
  lin_regr=lm(grep_data$sample2_abundance~grep_data$sample1_abundance)
  grep_data$residuals=lin_regr$residuals
  grep_data$cnv_value=rstandard(lin_regr)
  grep_data$cnv_fdrtool.pval=fdrtool(grep_data$cnv_value,plot=FALSE)$pval
  grep_data$cnv_fdrtool.qval=fdrtool(grep_data$cnv_value,plot=FALSE)$qval
  out_table=rbind(out_table,grep_data)
}
write.table(out_table,'output_auto_cnv.txt',sep='\t',quote=FALSE,row.names=FALSE)
