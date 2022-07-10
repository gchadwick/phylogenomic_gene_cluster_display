#this script just iterates through the HMM hits from step3.1 and parses them into more useful formats
import os
from simplehmmer.simplehmmer import HMMERParser

#manually curated list of all associated hmms to be used for this step
hmmlist = ['K00330','K00331','K00332','K00333','K00334','K00335','K00336','K00337','K00338','K00339','K00340','K00341','K00342','K00343','K00437','PF00374','TIGR00940','TIGR00941','TIGR00942','TIGR00943','TIGR00944','DUF2309','ChpXY','TIGR00073','TIGR00100','TIGR00074','TIGR00075','pfam03063','TIGR02951','pfam01656','pfam13244','pfam03334','pfam04066','TIGR00143','TIGR02124']
#make a list of all genomes in the "archaea" folder
archaea_genomes = []
for file in os.listdir('./archaea/'):
	if file.endswith(".fna"):
		archaea_genomes.append(file.split('_genomic.fna')[0])
#this bit of code relies on the HMMERParser function from simplehmmer to parse the output files
#generated in step3.1.py. This outputs more useful .csv format files into the hmm directories.
#each gene that was hit with the HMM search is recorded along with the HMM that is the best hit
#as well as the evalue associated with that hit
counter = 0
for archaea in archaea_genomes:
		counter += 1
		print counter
		besthmmhit = {}
		bestevalhit = {}
		outfile = open('./archaea/'+archaea+'/hmm.hits.csv','w')
		outfile2 = open('./archaea/'+archaea+'/hmm.hits.evals.csv','w')
		for hmm in hmmlist:
			with open('./archaea/'+archaea+'/'+hmm+'/'+hmm+'_out.txt','r') as fh:
				HP = HMMERParser(fh)
				while True:
					result = HP.next()
					if result:
						gene = str(result).split('\t')[0]
						if gene in besthmmhit.keys():
							if bestevalhit[gene] > float(str(result).split('\t')[6]):
								bestevalhit[gene] = float(str(result).split('\t')[6])
								besthmmhit[gene] = hmm
							else:
								pass
						else:
							if float(str(result).split('\t')[6]) < 1e-3:
								bestevalhit[gene] = float(str(result).split('\t')[6])
								besthmmhit[gene] = hmm
							else:
								pass
					else:
						break
		for gene in besthmmhit.keys():
			outfile.write(gene+','+besthmmhit[gene]+'\n')
			outfile2.write(gene+','+str(bestevalhit[gene])+'\n')
		outfile.close()
		outfile2.close()