#this step takes all of the genes identified in steps 2 and 3, extracts their surrounding genes from the 
#genbank files generated in step 1, and then uses the hard-coded color definitions to generate gene context
#images as .png files
import os
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv

#define colors that will be used for the various associated genes
nuoL_color = colors.Color(0.4,0.76,0.647,1)
nuoM_color = colors.Color(0.988,0.553,0.384,1)
nuoN_color = colors.Color(0.553,0.6275,0.796,1)
nuoAHJK_color = colors.Color(0.906,0.541,0.7647,1)
nuoBCDI_color = colors.Color(0.651,0.847,0.3294,1)
nuoEFG_color = colors.Color(1,0.851,0.1843,1)
mrpBCE_color = colors.Color(0.6,0.4,0.1,1)
mrpAD_color = colors.Color(0.6,0.6,0.1,1)
cupAB_color = colors.Color(1,0.6,0.6,1)
DUF2309_color = colors.Color(0.6,1,0.6,1)
lightgray = colors.Color(0.8,0.8,0.8,1)

# makes a dictionary that converts between hmm numbers (i.e. K00330) to gene (i.e. NuoA)
hmmconverter = {}
with open('./hmm.gene.coverter.csv') as csvfile:
	for line in csvfile:
		hmmconverter[line.split(',')[0]] = line.split(',')[1].rstrip()
		
# make a list of all archaea genome names
archaea_genomes = []
for file in os.listdir('./archaea/'):
	if file.endswith(".fna"):
		archaea_genomes.append(file.split('_genomic.fna')[0])
		
dummy_genbank = SeqIO.read('./break.gb','genbank')
fails = []
big_counter = 1
for archaea in archaea_genomes:
			print big_counter, 'of', len(archaea_genomes), archaea
			big_counter += 1
			# make a list of gene ids (locus tags) that have been hit by the hmm for antiporter subunits (pfam Proton_antipo_M.hmm)
			gene_ids = []
			with open('./archaea/'+archaea+'/'+archaea+'.hmm.txt/'+'antiporter_gene_list.txt','r') as fh:
				for line in fh:
					gene_ids.append(line.rstrip())
			max_length = 0
			combined_record2 = 'None'
			combine_counter = 0
			# make an empty GenomeDiagram object
			gd_diagram = GenomeDiagram.Diagram()
			#hmmhits is a dict of all genes that have hit with one of our supplemental
			#hmms (not the antiporter, but NuoA,J,K,MrpC,etc.)
			hmmhits = {}
			with open('./archaea/'+archaea+'/hmm.hits.csv') as csvfile:
				for line in csvfile:
					hmmhits[line.split(',')[0]] = line.split(',')[1].rstrip()
			#loop through all the stuff in the genbank file
			try:
				for records in SeqIO.parse('./archaea/'+archaea+'/'+archaea+'.gbf','genbank'):
					starts = []
					ends = []
					orders = []
					hits = []
					hit_counter = 0
					#for each feature in a record we will ask if it has a locus tag, if so
					#then we as if that locus tage is in our list of things that hit with
					#the pfam Proton_antipo_M.hmm.  if so, then we add its start to the starts
					#list. we make the starts and ends 10kb on either side of the feature
					#if another feature is within this area, we modify the start and end to
					#contain both
					for feature in records.features:
						if 'locus_tag' in feature.qualifiers.keys():
							if feature.qualifiers['locus_tag'][0] in gene_ids:
								if len(starts)>0:
									if abs(feature.location.start-(starts[-1]+10000))<10000:
										starts[-1]=min(starts[-1],(max(0, feature.location.start-10000)))
										ends[-1]=max(ends[-1],(min(len(records), feature.location.end+10000)))
										orders[-1]=feature.location.strand
										hit_counter += 1
										hits[-1]=hit_counter
									else:
										starts.append(max(0, feature.location.start-10000))
										ends.append(min(len(records),feature.location.end+10000))
										orders.append(feature.location.strand)
										hit_counter = 1
										hits.append(hit_counter)
								else:
									starts.append(max(0, feature.location.start-10000))
									ends.append(min(len(records),feature.location.end+10000))
									orders.append(feature.location.strand)
									hit_counter = 1
									hits.append(hit_counter)
					counter = 0
					#once we have the starts and ends of all the regions of the genbank file
					#that contain Proton_antipo_M.hmm hits, we will slice the genbank record
					#so that we add only that region to our genome diagram.
					for st in starts:	
						sub_record = records[st:ends[counter]]
						if orders[counter] == -1:
							sub_record = sub_record.reverse_complement(name=True)
						if combine_counter == 0:
							combined_record2 = sub_record
						else:
							combined_record2 = combined_record2+dummy_genbank+sub_record
						counter += 1
						combine_counter += 1
						max_length = max(max_length,len(sub_record))
				if type(combined_record2) != str:
					gd_track_for_features = gd_diagram.new_track(1,name='',greytrack=False,start=0,end=len(combined_record2))
					gd_feature_set = gd_track_for_features.new_set()
					for sub_feature in combined_record2.features:
						if 'gene' in sub_feature.qualifiers.keys():
							if sub_feature.qualifiers['locus_tag'][0] in hmmhits.keys():
								if hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'NuoL':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoL_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['product'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'NuoM':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoM_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'NuoN':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoN_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['NuoA','NuoH','NuoJ','NuoK']:
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoAHJK_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['NuoB','NuoC','NuoD','NuoI']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoBCDI_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['NuoE','NuoF','NuoG']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoEFG_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['MrpB','MrpC','MrpE','MrpG','MrpF','MrpDUF']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=mrpBCE_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['MrpA','MrpD']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=mrpAD_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'CupAB':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.green,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'DUF2309':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.black,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['HypA','HypB','HypC','HypD','HypE','HypF']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.purple,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['CooC','CooS','CooF']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.blue,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['NiFe1','NiFe2']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.red,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
								else:
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.gray,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_('+sub_feature.qualifiers['gene'][0].split('_')[0]+')'+'_'+sub_feature.qualifiers['locus_tag'][0])
							else:
								gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=lightgray,label=False,name=sub_feature.qualifiers['gene'][0])
						elif 'locus_tag' in sub_feature.qualifiers.keys():
							if sub_feature.qualifiers['locus_tag'][0] in hmmhits.keys():
								if hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'NuoL':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoL_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'NuoM':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoM_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'NuoN':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoN_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['NuoA','NuoH','NuoJ','NuoK']:
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoAHJK_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['NuoB','NuoC','NuoD','NuoI']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoBCDI_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['NuoE','NuoF','NuoG']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=nuoEFG_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['MrpB','MrpC','MrpE','MrpG','MrpF','MrpDUF']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=mrpBCE_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['MrpA','MrpD']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=mrpAD_color,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'CupAB':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.green,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] == 'DUF2309':	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.black,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['HypA','HypB','HypC','HypD','HypE','HypF']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.purple,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['CooC','CooS','CooF']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.blue,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								elif hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]] in ['NiFe1','NiFe2']:	
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.red,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
								else:
									gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=colors.gray,label=False,name=hmmconverter[hmmhits[sub_feature.qualifiers['locus_tag'][0]]]+'_'+sub_feature.qualifiers['locus_tag'][0])
							elif sub_feature.qualifiers['locus_tag'][0]=="BREAKBREAK":
								gd_feature_set.add_feature(sub_feature,sigil="BOX", color=colors.black,label=False,name='BREAK',strand=False)
							else:
								gd_feature_set.add_feature(sub_feature,sigil="BIGARROW", arrowshaft_height=1.0,color=lightgray,label=False)
					gd_diagram.draw(format="linear",orientation="landscape", pagesize=(5*cm,len(combined_record2)/500*cm),fragments=1,start=0,end=len(combined_record2))
					gd_diagram.write('./archaea/'+archaea+'/'+'allantiporter.png', 'png')
			except:
				print big_counter,archaea,"FAIL"
				fails.append(archaea)
print fails

