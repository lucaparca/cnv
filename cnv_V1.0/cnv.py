import os,sys,string,cPickle,time
from lib import load_parameters,load_data,ctc_replicates

#parameters section
replicates,s1col,s2col,toleratedNA,annotation,species,idtype,input,output,dolog2=load_parameters()

#check variables for shift analysis
if len(s1col)==1 or len(s2col)==1:
	print 'Limma analysis cannot be done (no multiple samples per condition).';time.sleep(3);
	shiftan=False
else:
	shiftan=True
	shiftcol=[]
	for col in s1col+s2col:shiftcol.append(str(col+1))
	shiftcol=string.join(shiftcol,',')
	control=string.join(['0']*len(s1col)+['1']*len(s2col),',')
	treated=string.join(['1']*len(s1col)+['0']*len(s2col),',')

#load data
content=load_data(input)

#check for replicates, tolerated missing data and collapse replicates for every sample
content2=ctc_replicates(replicates,content,toleratedNA,s1col,s2col,dolog2)

#automatic protein grouping
if annotation=='automatic':
	out=open('auto_cnv_input.txt','w')
	for line in content2:out.write(string.join(line,'\t')+'\n')
	out.close()
	temp=open('templates/auto_cnv.r').read()
	if shiftan==True:
		temp+=open('templates/shift_analysis.r').read()
		temp=temp.replace('<insertcolumns>',shiftcol);temp=temp.replace('<insertcontrol>',control);temp=temp.replace('<inserttreated>',treated)
	temp=temp.replace('<insertspecieshere>',species);temp=temp.replace('<insertfilterhere>',idtype);temp=temp.replace('output',output)
	outr=open('auto_cnv.r','w');outr.write(temp);outr.close()
	os.system('Rscript auto_cnv.r')
#provided protein grouping
else:
	compartments={}
	annotation=content2[0].index('annotation')
	for line in content2[1:]:
		compartment=line[annotation].split(';')
		for c in compartment:
			if c not in compartments:compartments[c]=[]
			compartments[c].append(line)
	out=open('custom_cnv_input.txt','w')
	for line in content2:out.write(string.join(line,'\t')+'\n')
	out.close()

	out=open('custom_cnv.r','w')
	out.write("library('GO.db');library('biomaRt');library('affy');library('fdrtool')\na=read.table('custom_cnv_input.txt',quote='',header=TRUE)\nrownames(a)=a[,1]\n")

	out.write('out_table=c()\n')
	out.write('compartment_names=as.vector(unique(a$annotation))\n')
	out.write('for (i in 1:%d){\n'%(len(compartments)))
	out.write("   grep_data=subset(a,annotation %in% grep(compartment_names[i],annotation,value=TRUE))\n")
	out.write('   if (length(grep_data$annotation)<20) {next}\n')
	out.write('   lin_regr=lm(grep_data$sample2_abundance~grep_data$sample1_abundance)\n')
	out.write('   grep_data$residuals=lin_regr$residuals\n')
	out.write('   grep_data$cnv_value=rstandard(lin_regr)\n')
	out.write('   grep_data$cnv_fdrtool.pval=fdrtool(grep_data$cnv_value,plot=FALSE)$pval\n')
	out.write('   grep_data$cnv_fdrtool.qval=fdrtool(grep_data$cnv_value,plot=FALSE)$qval\n')
	out.write('   out_table=rbind(out_table,grep_data)\n')
	out.write('}\n')
	out.write("write.table(out_table,'%s_custom_cnv.txt',sep='\\t',quote=FALSE)\n"%(output))

#	for compartment in compartments:
#		out.write("%s=subset(a,annotation %%in%% grep('%s',annotation,value=TRUE))\n"%(compartment,compartment))
#		out.write("regr%s=lm(%s$sample2_abundance~%s$sample1_abundance)\n"%(compartment,compartment,compartment))
#		out.write("%s$residuals=regr%s$residuals\n"%(compartment,compartment))
#		out.write("%s$cnv_value=rstandard(regr%s)\n"%(compartment,compartment))
#		out.write("%s$std_fdrtool.pval=fdrtool(%s$cnv_value)$pval\n"%(compartment,compartment))
#	out.write("write.table(rbind(%s),'%s_custom_cnv.txt',sep='\t',quote=FALSE)\n"%(string.join(compartments.keys(),','),output))
	if shiftan==True:
		temp=open('templates/shift_analysis.r').read()
		temp=temp.replace('<insertcolumns>',shiftcol);temp=temp.replace('<insertcontrol>',control);temp=temp.replace('<inserttreated>',treated)
		temp=temp.replace('<insertspecieshere>',species);temp=temp.replace('<insertfilterhere>',idtype);temp=temp.replace('output',output)
		out.write(temp)
	out.close()
	os.system('R -q -f custom_cnv.r')
