library('GO.db','biomaRt','affy','fdrtool')
a=read.table('custom_cnv_input.r',sep='\t',quote='',header=TRUE)
rownames(a)=a$V1

sel1=which(a$accession %in% subset(t,go_id %in% nucccgo)$<insertfilterhere>)
sel2=which(a$accession %in% subset(t,go_id %in% cytccgo)$<insertfilterhere>)
sel3=which(a$accession %in% subset(t,go_id %in% extccgo)$<insertfilterhere>)
selm=which(a$accession %in% subset(t,go_id %in% mitccgo)$<insertfilterhere>)


subset(aa,anno %in% grep('c',anno,value=TRUE))

nuc=a[sel1,]
cyt=a[sel2,]
ext=a[sel3,]
mit=a[selm,]

regrnuc=lm(nuc$sample2_abundance~nuc$sample1_abundance)
regrcyt=lm(cyt$sample2_abundance~cyt$sample1_abundance)
regrext=lm(ext$sample2_abundance~ext$sample1_abundance)
regrmit=lm(mit$sample2_abundance~mit$sample1_abundance)

nuc$residuals=regrnuc$residuals
cyt$residuals=regrcyt$residuals
ext$residuals=regrext$residuals
mit$residuals=regrmit$residuals

nuc2$cnv_value=rstandard(regrnuc)
cyt2$cnv_value=rstandard(regrcyt)
ext2$cnv_value=rstandard(regrext)
mit2$cnv_value=rstandard(regrmit)

nuc$std_fdrtool.pval=fdrtool(nuc$cnv_value, plot=FALSE)$pval
cyt$std_fdrtool.pval=fdrtool(cyt$cnv_value, plot=FALSE)$pval
ext$std_fdrtool.pval=fdrtool(ext$cnv_value, plot=FALSE)$pval
mit$std_fdrtool.pval=fdrtool(mit$cnv_value, plot=FALSE)$pval

write.table(rbind(nuc,cyt,ext,mit),'output_auto_cnv.txt',sep='\t',quote=FALSE)
