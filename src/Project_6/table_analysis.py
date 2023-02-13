import pandas as pd
import math

infile='tmp_GCF_fixed.csv'

df=pd.read_csv(infile, header=None)
masterlist=[]

GT_NCBI = []
GT_GeneMark = []
GT_Prodigal = []
GT_StORF = []


for i in range(len(df)):
#	print (df[0][i], df.iloc[i])
	rown = df[0][i]
	nlist = rown.split('|')
	tool = nlist[1].split(',')
	#masterlist.append(tool)
	if "GCF_000009045" in tool:
		series1 = df.iloc[i]
		for y in range(1, len(series1)):	
			element = series1.iloc[y]
			if isinstance(element, float):
				continue
			else:		
				element = element.split("|")
				if float(element[1]) >= 0.5:
					GT_NCBI.append(element[0])

	if "GeneMark_S_2" in tool:
		series1 = df.iloc[i]
		for y in range(1, len(series1)):
			element = series1.iloc[y]
			if isinstance(element, float):
				continue
			else:
				element = element.split("|")
				if float(element[1]) >= 0.5:
					GT_GeneMark.append(element[0])

	if "Prodigal" in tool:
		series1 = df.iloc[i]
		for y in range(1, len(series1)):
			element = series1.iloc[y]
			if isinstance(element, float):
				continue
			else:
				element = element.split("|")
				if float(element[1]) >= 0.5:
					GT_Prodigal.append(element[0])

	if "StORF-Reporter" in tool:
		series1 = df.iloc[i]
		for y in range(1, len(series1)):
			element = series1.iloc[y]
			if isinstance(element, float):
				continue
			else:
				element = element.split("|")
				if float(element[1]) >= 0.5:
					GT_StORF.append(element[0])


GT_NCBI_set = set(GT_NCBI)
GT_GeneMark_set = set(GT_GeneMark)
GT_Prodigal_set = set(GT_Prodigal)
GT_StORF_set = set(GT_StORF)

NCBI_total = len(GT_NCBI_set)
GeneMark_total = len(GT_GeneMark_set)
Prodigal_total = len(GT_Prodigal_set)
StORF_total = len(GT_StORF_set)

NCBI_diff_all = len(GT_NCBI_set.difference(GT_GeneMark_set, GT_Prodigal_set, GT_StORF_set))
Prodigal_diff_all = len(GT_Prodigal_set.difference(GT_GeneMark_set, GT_NCBI_set, GT_StORF_set))
GeneMark_diff_all = len(GT_GeneMark_set.difference(GT_NCBI_set, GT_Prodigal_set, GT_StORF_set))
StORF_diff_all = len(GT_StORF_set.difference(GT_GeneMark_set, GT_Prodigal_set, GT_NCBI_set))

NCBI_and_Prodigal = len(GT_NCBI_set.intersection(GT_Prodigal_set))
#NCBI_and_StORF = len(GT_NCBI_set.intersection(GT_StORF_set))
Prodigal_and_GeneMark = len(GT_Prodigal_set.intersection(GT_GeneMark_set))
#GeneMark_and_StORF = len(GT_GeneMark_set.intersection(GT_StORF_set))
NCBI_and_GeneMark = len(GT_NCBI_set.intersection(GT_GeneMark_set))

Inter_of_all=len(GT_NCBI_set.intersection(GT_GeneMark_set, GT_Prodigal_set))


print(f'NCBI len: {NCBI_total}\nGeneMark len: {GeneMark_total}\nProdigal len: {Prodigal_total}\nStORF len: {StORF_total}\n')
print(f'NCBI diff all: {NCBI_diff_all}\nGeneMark diff all: {GeneMark_diff_all}\nProdigal diff all: {Prodigal_diff_all}\nStORF diff all: {StORF_diff_all}\n')
print(f'NCBI and Prodigal: {NCBI_and_Prodigal}\nProdigal and GeneMark: {Prodigal_and_GeneMark}\nGeneMark and NCBI:{NCBI_and_GeneMark}\n')
print(f'Intersection of all: {Inter_of_all}\n')
#print(GT_GeneMark_set.difference(GT_NCBI_set, GT_Prodigal_set, GT_StORF_set))
#print(GT_NCBI_set.difference(GT_GeneMark_set, GT_Prodigal_set, GT_StORF_set))
