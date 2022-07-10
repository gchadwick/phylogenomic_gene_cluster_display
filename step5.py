#step 5 uses the r95_taxonomy.csv file from GTBD, the phylogenomic tree from GTDB,
#and the gene context images made in step 4 to generate pdf files with trees where every
#node is labelled with the full taxonomy string and the gene clusteres.
from ete3 import Tree, faces
import os
from shutil import copyfile
taxonomy_converter = {}
with open('./r95_taxonomy.csv') as csvfile:
	for line in csvfile:
		taxonomy_converter[line.split(',')[0]] = line.split(',')[1].rstrip()
t = Tree("ar122_r95_renamed.geneious.tre")
counter = 0
outfile = open('arc_leaves.csv','w')
for node in t:
	counter += 1
	leaf_name = str(node).lstrip().split('--')[1].replace("\'","")
	taxonomy_label = faces.TextFace(taxonomy_converter[leaf_name])
	print counter, leaf_name, taxonomy_converter[leaf_name]
	node.add_face(taxonomy_label,0,"aligned")
	if os.path.isfile('./archaea/'+leaf_name+'/allantiporter.pdf'):
		node.add_face(faces.ImgFace('./archaea/'+leaf_name+'/allantiporter.png',height=10),1,"aligned")
		copyfile('./archaea/'+leaf_name+'/allantiporter.png', './named_pngs/'+taxonomy_converter[leaf_name].split(' ')[0]+"_"+leaf_name+'.png')
		outfile.write(leaf_name+","+taxonomy_converter[leaf_name]+",pass\n")
	else:
		print counter, 'Fail'
		outfile.write(leaf_name+","+taxonomy_converter[leaf_name]+",fail\n")
t.render("arc_tree.pdf")