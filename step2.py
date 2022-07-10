#before you start, get a hidden markov model from pfam, tigrfam, etc. or generate your own. 
#store these in a folder called "hmms" and insure they have the extension ".hmm"
#make sure all of your genomes are in a folder, in this case called "archaea"


import os
from simplehmmer.simplehmmer import HMMERRunner
HR = HMMERRunner(prefix='antiporter')

#make a list of all genomes in the "archaea" folder
archaea_genomes =[]
for file in os.listdir('./archaea/'):
	if file.endswith(".fna"):
		archaea_genomes.append(file.split('_genomic.fna')[0])

#iterate through the archaea genomes searching each one of their protein .faa files
#for hits to the 'Proton_antipo_M.hmm' HMM of pfam00361. Output of this search is saved
#as *archaea*.hmm.txt.
counter = 0
for archaea in archaea_genomes:
	counter += 1
	print counter,"of:",len(archaea_genomes),archaea
	HR.search('Proton_antipo_M.hmm','./archaea/'+archaea+'/'+archaea+'.faa','./archaea/'+archaea+'/'+archaea+'.hmm.txt')
	