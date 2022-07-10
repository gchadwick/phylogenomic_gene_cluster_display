#before you start, collect a set of hidden markov models from pfam, tigrfam, etc. or generate your own. 
#store these in a folder called "hmms" and insure they have the extension ".hmm"
#make sure all of your genomes are in a folder, in this case called "archaea"

import os
from simplehmmer.simplehmmer import HMMERRunner

#make a list of all genomes in the "archaea" folder
archaea_genomes =[]
for file in os.listdir('./archaea/'):
	if file.endswith(".fna"):
		archaea_genomes.append(file.split('_genomic.fna')[0])

#make a list of all hmms to be searched for
hmms = []
for file in os.listdir('./hmms/'):
	if file.endswith(".hmm"):
		hmms.append(file.split('.hmm')[0])

#iterate through each hmm searching against each genome. writes outputs to new file named after the hmm  
#within the genome folder (each genome should have a corresponding folder if you ran Step 1 with Prokka)
for hmm in hmms:
	print hmm
	HR = HMMERRunner(prefix=hmm)
	for archaea in archaea_genomes:
		HR.search('./hmms/'+hmm+'.hmm','./archaea/'+archaea+'/'+archaea+'.faa','./archaea/'+archaea+'/'+hmm)
