library(preprocessCore)
library(limma)
library(ggplot2)

data2 = normalize.quantiles(as.matrix(a[,c(<insertcolumns>)]))
rownames(data2) = rownames(a)
Control=c(<insertcontrol>)
Treated=c(<inserttreated>)
design=cbind(Control, Treated)
fit = lmFit(data2, design=design)
contrastsMatrix = makeContrasts("Control-Treated",levels = design)
fit2 = contrasts.fit(fit, contrasts = contrastsMatrix)
fit2 = eBayes(fit2)
hits_all<-topTable(fit2, coef = "Control-Treated",number=20000)

unimart=useMart("ENSEMBL_MART_ENSEMBL",host='www.ensembl.org')
martdataset = useDataset("<insertspecieshere>_gene_ensembl",mart=unimart)
filters = listFilters(martdataset);attributes = listAttributes(martdataset)
ids<-unlist(sapply(as.character(rownames(hits_all)), function(x) do.call(rbind,strsplit(x, ";"))))
t<-getBM(attributes=c("<insertfilterhere>", "wikigene_name","go_id") , filters = "<insertfilterhere>", values = ids, mart = martdataset)
for (i in 1:nrow(hits_all)){
  id<-as.vector(unlist(sapply(as.character(rownames(hits_all)[i]), function(x) do.call(rbind,strsplit(x, ";")))))
  hits_all$short.name[i]<-unique(t[t[,"<insertfilterhere>"]%in%id, "wikigene_name"])[1]
  hits_all$accession[i]<-unique(t[t[,"<insertfilterhere>"]%in%id, "<insertfilterhere>"])[1]
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

sel1=which(hits_all$accession %in% subset(t,go_id %in% nucccgo)$<insertfilterhere>)
sel2=which(hits_all$accession %in% subset(t,go_id %in% cytccgo)$<insertfilterhere>)
sel3=which(hits_all$accession %in% subset(t,go_id %in% extccgo)$<insertfilterhere>)
selm=which(hits_all$accession %in% subset(t,go_id %in% mitccgo)$<insertfilterhere>)
selp=which(hits_all$accession %in% subset(t,go_id %in% perccgo)$<insertfilterhere>)
sell=which(hits_all$accession %in% subset(t,go_id %in% lysccgo)$<insertfilterhere>)
selr=which(hits_all$accession %in% subset(t,go_id %in% relccgo)$<insertfilterhere>)
selg=which(hits_all$accession %in% subset(t,go_id %in% golccgo)$<insertfilterhere>)
selmem=which(hits_all$accession %in% subset(t,go_id %in% memccgo)$<insertfilterhere>)
selnme=which(hits_all$accession %in% subset(t,go_id %in% nmeccgo)$<insertfilterhere>)


z1<-hits_all[sel1,"logFC"]
z2<-hits_all[sel2,"logFC"]
z3<-hits_all[sel3,"logFC"]
z4<-hits_all[selm,"logFC"]
z5<-hits_all[selp,"logFC"]
z6<-hits_all[sell,"logFC"]
z7<-hits_all[selr,"logFC"]
z8<-hits_all[selg,"logFC"]
z9<-hits_all[selmem,"logFC"]
z10<-hits_all[selnme,"logFC"]

dataset1<-as.data.frame(z1);colnames(dataset1)[1]<-"logFC";dataset1$cond<-"Nucleus"
dataset2<-as.data.frame(z2);colnames(dataset2)[1]<-"logFC";dataset2$cond<-"Cytoplasm";
dataset3<-as.data.frame(z3);colnames(dataset3)[1]<-"logFC";dataset3$cond<-"Extracellular";
dataset4<-as.data.frame(z4);colnames(dataset4)[1]<-"logFC";dataset4$cond<-"Mitochondrion"
dataset5<-as.data.frame(z5);colnames(dataset5)[1]<-"logFC";dataset5$cond<-"Peroxisome"
dataset6<-as.data.frame(z6);colnames(dataset6)[1]<-"logFC";dataset6$cond<-"Lysosome"
dataset7<-as.data.frame(z7);colnames(dataset7)[1]<-"logFC";dataset7$cond<-"Endoplasmic reticulum"
dataset8<-as.data.frame(z8);colnames(dataset8)[1]<-"logFC";dataset8$cond<-"Golgi apparatus"
dataset9<-as.data.frame(z9);colnames(dataset9)[1]<-"logFC";dataset9$cond<-"Cell membrane"
dataset10<-as.data.frame(z10);colnames(dataset10)[1]<-"logFC";dataset10$cond<-"Nuclear membrane"

dataset_all<-rbind(dataset1,dataset2,dataset3,dataset4,dataset5,dataset6,dataset7,dataset8,dataset9,dataset10)
pdf('output_compartment_shifts_density.pdf')
ggplot(dataset_all, aes(x=logFC, fill=as.factor(cond))) + geom_density(alpha=.5) + scale_x_continuous() + theme_bw() + scale_fill_manual(values=c("#4900FF", "#FB7706", "#53FF05", "#15A9FE", "#DF04FF", "#22FFAB", "#E3FF09", "#22FF24", "#FB0007", "#0001FF"))  + theme(legend.justification=c(0,1), legend.position=c(0, 1), legend.background = element_rect()) + theme(panel.grid.minor = element_line(colour = NA))+
theme(axis.text=element_text(size=7),
      axis.title=element_text(face = "bold",size=7),
      plot.title=element_text(size=7),
      panel.background=element_rect(fill='white',color='black'),
      panel.grid=element_blank())
dev.off()
pdf('output_compartment_shifts_boxplot.pdf')
ggplot(dataset_all, aes(factor(cond),logFC, fill=as.factor(cond))) +
      geom_boxplot() +
      geom_hline(yintercept=0,color='red',linetype='dashed') +
      scale_fill_manual(values=c("#4900FF", "#FB7706", "#53FF05", "#15A9FE", "#DF04FF", "#22FFAB", "#E3FF09", "#22FF24", "#FB0007", "#0001FF"))+
      theme_bw()+
      theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
dev.off()

p1=t.test(z1,c(z2,z3,z4,z5,z6,z7,z8,z9,z10))$p.value
p2=t.test(z2,c(z1,z3,z4,z5,z6,z7,z8,z9,z10))$p.value
p3=t.test(z3,c(z2,z1,z4,z5,z6,z7,z8,z9,z10))$p.value
p4=t.test(z4,c(z2,z3,z1,z5,z6,z7,z8,z9,z10))$p.value
p5=t.test(z5,c(z2,z3,z4,z1,z6,z7,z8,z9,z10))$p.value
p6=t.test(z6,c(z2,z3,z4,z5,z1,z7,z8,z9,z10))$p.value
p7=t.test(z7,c(z2,z3,z4,z5,z6,z1,z8,z9,z10))$p.value
p8=t.test(z8,c(z2,z3,z4,z5,z6,z7,z1,z9,z10))$p.value
p9=t.test(z9,c(z2,z3,z4,z5,z6,z7,z8,z1,z10))$p.value
p10=t.test(z10,c(z2,z3,z4,z5,z6,z7,z8,z9,z1))$p.value

shift_pvalue=p.adjust(c(p1,p2,p3,p4,p5,p6,p7,p8,p9,10),method='fdr')
shift_mean=c(mean(z1),mean(z2),mean(z3),mean(z4),mean(z5),mean(z6),mean(z7),mean(z8),mean(z9),mean(z10))

outshift=cbind(cbind(c('Nucleus','Cytoplasm','Extracellular','Mitochondrion','Peroxisome','Lysosome','Endoplasmic reticulum','Golgi apparatus','Cell membrane','Nuclear Membrane'),shift_mean),shift_pvalue)
write.table(outshift,'output_compartment_shift_analysis.txt',sep='\t',quote=FALSE,row.names=FALSE)
