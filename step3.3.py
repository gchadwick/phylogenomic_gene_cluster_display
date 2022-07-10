#this script iterates through all of the genomes, counting how many hits there are for the antiporter
#genes and extracting all of the genes into a single combined fasta file for the entire dataset. In addition,
#for each genome the script extracts all the hits and stores them in their own directory within the genome
#directory.
import os
from simplehmmer.simplehmmer import HMMERParser
from Bio import SeqIO

outfile = open('./archaea/archaea.antiporter.counts.txt','w')
outfile_fasta = open('./archaea/archaea.all.antiporters.fasta','w')

#make a list of all genomes in the "archaea" folder
archaea_genomes = []
for file in os.listdir('./archaea/'):
    if file.endswith(".fna"):
        archaea_genomes.append(file.split('_genomic.fna')[0])

#iterate through the each genome, populating the large files created above
#as well as created genome-specific files and populating them with the gene
#hits
counter = 0
for archaea in archaea_genomes:
    counter +=1
    with open('./archaea/'+archaea+'/'+archaea+'.hmm.txt/'+'antiporter_out.txt','r') as fh:
        HP = HMMERParser(fh)
        antiporter_list = []
        while True:
            result = HP.next()
            if result:
                #print str(result).split('\t')[0], str(result).split('\t')[-1]
                antiporter_list.append(str(result).split('\t')[0])
            else:
                break
        antiporter_list_unique = []
        outfile2 = open('./archaea/'+archaea+'/'+archaea+'.hmm.txt/'+'antiporter_gene_list.txt','w')
        outfile3 = open('./archaea/'+archaea+'/'+archaea+'.hmm.txt/'+'antiporter_out.fasta','w')
        for x in antiporter_list:
            if x not in antiporter_list_unique:
                antiporter_list_unique.append(x)
        for x in antiporter_list_unique:
            outfile2.write(x+'\n')
            for seq_record in SeqIO.parse('./archaea/'+archaea+'/'+archaea+'.faa','fasta'):
                if seq_record.id == x:
                    outfile_fasta.write('>')
                    outfile_fasta.write(x)
                    outfile_fasta.write('\n')
                    outfile_fasta.write(str(seq_record.seq))
                    outfile_fasta.write('\n')
                    outfile3.write('>')
                    outfile3.write(x)
                    outfile3.write('\n')
                    outfile3.write(str(seq_record.seq))
                    outfile3.write('\n')
        outfile2.close()
        outfile3.close()
        print counter,len(archaea_genomes),archaea, len(antiporter_list_unique)
        outfile.write(archaea+','+str(len(antiporter_list_unique))+'\n')
outfile.close()
outfile_fasta.close()