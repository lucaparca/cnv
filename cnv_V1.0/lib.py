import os,sys,string,cPickle,numpy,math

#parameters section

def load_parameters():
	paramfile=open('parameters.txt').read()
	if '\r' in paramfile:paramfile=paramfile.split('\r')
	else:paramfile=paramfile.split('\n')
	s1col=[]
	s2col=[]
	for line in paramfile:
		data=line.split('=')
		if data[0]=='replicates':
			if data[1]=='False':replicates=False
			else:replicates=True
		elif data[0] in ['sample1_columns','sample2_columns']:
			columndata=data[1]
			if ',' in columndata:
				for columnserie in columndata.split(','):
					if ':' in columnserie:
						if '1' in data[0]:s1col+=range(int(columnserie.split(':')[0]),int(columnserie.split(':')[1])+1)
						else:s2col+=range(int(columnserie.split(':')[0]),int(columnserie.split(':')[1])+1)
					elif '-' in columnserie:
						if '1' in data[0]:s1col+=range(int(columnserie.split('-')[0]),int(columnserie.split('-')[1])+1)
						else:s2col+=range(int(columnserie.split('-')[0]),int(columnserie.split('-')[1])+1)
					else:
						if '1' in data[0]:s1col.append(int(columnserie))
						else:s2col.append(int(columnserie))
			elif ':' in columndata:
				if '1' in data[0]:s1col=range(int(columndata.split(':')[0]),int(columndata.split(':')[1])+1)
				else:s2col=range(int(columndata.split(':')[0]),int(columndata.split(':')[1])+1)
			elif '-' in columndata:
				if '1' in data[0]:s1col=range(int(columndata.split('-')[0]),int(columndata.split('-')[1])+1)
				else:s2col=range(int(columndata.split('-')[0]),int(columndata.split('-')[1])+1)
			else:
				if '1' in data[0]:s1col=[int(columndata)]
				else:s2col=[int(columndata)]
		elif data[0]=='toleratedNA':toleratedNA=int(data[1])
		elif data[0]=='annotation':
			if data[1]=='automatic':annotation='automatic'
			else:annotation='custom'
		elif data[0]=='species':species=data[1]
		elif data[0]=='idtype':
			if data[1]=='uniprot':idtype='uniprotswissprot'
			elif data[1]=='ensembl':idtype='ensembl_peptide_id'
		elif data[0]=='inputfile':input=data[1]
		elif data[0]=='projectname':output=data[1]
		elif data[0]=='dolog2':
			if data[1]=='False':dolog2=False
			else:dolog2=True
	return replicates,s1col,s2col,toleratedNA,annotation,species,idtype,input,output,dolog2

def load_data(input):
	inputfile=open(input).read();inputfile.replace('"','')
	content=[]
	if '\r' in inputfile:contentlines=inputfile.split('\r')
	else:contentlines=inputfile.split('\n')
	for line in contentlines:
		if line.split()==[]:continue
		content.append(line.split())
	return content

def ctc_replicates(replicates,content,toleratedNA,s1col,s2col,dolog2):
	content2=[]
	outline=[]
	for x in range(len(content[0])):
		#if replicates==False and (x in s1col or x in s2col):continue
		outline.append(content[0][x])
	outline.append('sample1_abundance');outline.append('sample2_abundance')
	content2.append(outline)
	for line in content[1:]:
		s1values=[]
		s2values=[]
		for col in s1col:
			try:
				if dolog2==True:s1values.append(math.log(float(line[col]),2))
				else:s1values.append(float(line[col]))
			except:continue
		for col in s2col:
			try:
				if dolog2==True:s2values.append(math.log(float(line[col]),2))
				else:s2values.append(float(line[col]))
			except:continue
		if replicates==True and (len(s1values)<(len(s1col)-toleratedNA) or len(s2values)<(len(s2col)-toleratedNA)):continue
		s1value=numpy.mean(s1values)
		s2value=numpy.mean(s2values)
		outline=[]
		for x in range(len(line)):
		#	if replicates==False and (x in s1col or x in s2col):continue
			if x in s1col or x in s2col:
				try:
					if dolog2==True:outline.append(str(math.log(float(line[x]),2)))
					else:outline.append(line[x])
				except:outline.append(line[x])
			else:outline.append(line[x])
		outline.append(str(s1value));outline.append(str(s2value))
		content2.append(outline)
	return content2
