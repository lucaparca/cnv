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

sel1=which(a$accession %in% subset(t,go_id %in% nucccgo)$<insertfilterhere>)
sel2=which(a$accession %in% subset(t,go_id %in% cytccgo)$<insertfilterhere>)
sel3=which(a$accession %in% subset(t,go_id %in% extccgo)$<insertfilterhere>)
selm=which(a$accession %in% subset(t,go_id %in% mitccgo)$<insertfilterhere>)
selp=which(a$accession %in% subset(t,go_id %in% perccgo)$<insertfilterhere>)
sell=which(a$accession %in% subset(t,go_id %in% lysccgo)$<insertfilterhere>)
selr=which(a$accession %in% subset(t,go_id %in% relccgo)$<insertfilterhere>)
selg=which(a$accession %in% subset(t,go_id %in% golccgo)$<insertfilterhere>)
selmem=which(a$accession %in% subset(t,go_id %in% memccgo)$<insertfilterhere>)
selnme=which(a$accession %in% subset(t,go_id %in% nmeccgo)$<insertfilterhere>)

nuc=a[sel1,]
cyt=a[sel2,]
ext=a[sel3,]
mit=a[selm,]
per=a[selp,]
lys=a[sell,]
ret=a[selr,]
gol=a[selg,]
mem=a[selmem,]
nme=a[selnme,]

nuc$compartment="Nucleus"
cyt$compartment="Cytoplasm"
ext$compartment="Extracellular"
mit$compartment="Mitochondrion"
per$compartment="Peroxisome"
lys$compartment="Lysosome"
ret$compartment="Endoplasmic reticulum"
gol$compartment="Golgi"
mem$compartment="Cell membrane"
nme$compartment="Nuclear membrane"

regrnuc=lm(nuc$sample2_abundance~nuc$sample1_abundance)
regrcyt=lm(cyt$sample2_abundance~cyt$sample1_abundance)
regrext=lm(ext$sample2_abundance~ext$sample1_abundance)
regrmit=lm(mit$sample2_abundance~mit$sample1_abundance)
regrper=lm(per$sample2_abundance~per$sample1_abundance)
regrlys=lm(lys$sample2_abundance~lys$sample1_abundance)
regrret=lm(ret$sample2_abundance~ret$sample1_abundance)
regrgol=lm(gol$sample2_abundance~gol$sample1_abundance)
regrmem=lm(mem$sample2_abundance~mem$sample1_abundance)
regrnme=lm(nme$sample2_abundance~nme$sample1_abundance)

nuc$residuals=regrnuc$residuals
cyt$residuals=regrcyt$residuals
ext$residuals=regrext$residuals
mit$residuals=regrmit$residuals
per$residuals=regrper$residuals
lys$residuals=regrlys$residuals
ret$residuals=regrret$residuals
gol$residuals=regrgol$residuals
mem$residuals=regrmem$residuals
nme$residuals=regrnme$residuals

nuc$cnv_value=rstandard(regrnuc)
cyt$cnv_value=rstandard(regrcyt)
ext$cnv_value=rstandard(regrext)
mit$cnv_value=rstandard(regrmit)
per$cnv_value=rstandard(regrper)
lys$cnv_value=rstandard(regrlys)
ret$cnv_value=rstandard(regrret)
gol$cnv_value=rstandard(regrgol)
mem$cnv_value=rstandard(regrmem)
nme$cnv_value=rstandard(regrnme)

nuc$std_fdrtool.pval=fdrtool(nuc$cnv_value, plot=FALSE)$pval
cyt$std_fdrtool.pval=fdrtool(cyt$cnv_value, plot=FALSE)$pval
ext$std_fdrtool.pval=fdrtool(ext$cnv_value, plot=FALSE)$pval
mit$std_fdrtool.pval=fdrtool(mit$cnv_value, plot=FALSE)$pval
per$std_fdrtool.pval=fdrtool(per$cnv_value, plot=FALSE)$pval
lys$std_fdrtool.pval=fdrtool(lys$cnv_value, plot=FALSE)$pval
ret$std_fdrtool.pval=fdrtool(ret$cnv_value, plot=FALSE)$pval
gol$std_fdrtool.pval=fdrtool(gol$cnv_value, plot=FALSE)$pval
mem$std_fdrtool.pval=fdrtool(mem$cnv_value, plot=FALSE)$pval
nme$std_fdrtool.pval=fdrtool(nme$cnv_value, plot=FALSE)$pval

write.table(rbind(nuc,cyt,ext,mit,per,lys,ret,gol,mem,nme),'output_auto_cnv.txt',sep='\t',quote=FALSE,row.names=FALSE)
