#sometimes a full domain tree is just too big, so you might want to tree based on some phylogenetic
#level. This optional step will iterate through all phylum-level designations within bacteria or archaea
#and make trees specifically for that group. For some phyla (like the proteobacteria) these trees will still
#be too big, so this can be extended to finer phylogenetic resolution groups.
from ete3 import Tree, faces
import os
from shutil import copyfile

#iterate through all the genomes in the archaea taxonomy file and assign phylogeny to each
#genome in the dictionary "taxonomy_converter"
taxonomy_converter = {}
with open('./ar122_taxonomy_r95.csv') as csvfile:
	for line in csvfile:
		taxonomy_converter[line.split(',')[0]] = line.split(',')[1].rstrip()

#make a list of all phyla in the domain
phyla = []
for taxon in taxonomy_converter.values():
	if taxon.split(";")[1] in phyla:
		pass
	else:
		phyla.append(taxon.split(";")[1])
#load tree and store all leaves in a list
leaves = []
t = Tree("ar122_r95_renamed.geneious.fig.tre")
for node in t:
	leaves.append(str(node).lstrip().split('--')[1])
counter = 0
#iterate through each phyla, store all of the genomes within that taxa, then use the 
#get_common_ancestor function to generate the subtree just containing that phyla. labelling proceeds
#similarly to step 5.
for phylum in phyla:
    taxa = []
    counter +=1
    #print phylum, counter, "of", len(phyla)
    for leaf in taxonomy_converter.keys():
        if phylum+';' in taxonomy_converter[leaf]:
            leafstring = "\'"+leaf+"\'"
            if leafstring in leaves:
                taxa.append(leafstring)
    print phylum, len(taxa)
    if len(taxa) < 2000:
        if len(taxa) > 1:
            t2 = t.get_common_ancestor(taxa)
            for node in t2:
                leaf_name = str(node).lstrip().split('--')[1].replace("\'","")
                taxonomy_label = faces.TextFace(taxonomy_converter[leaf_name],fsize=20)
                #print counter, leaf_name, taxonomy_converter[leaf_name]
                node.add_face(taxonomy_label,0,"aligned")
                if os.path.isfile('./archaea/'+leaf_name+'/allantiporter.png'):
                    node.add_face(faces.ImgFace('./archaea/'+leaf_name+'/allantiporter.png',height=30),1,"aligned")
                    #copyfile('./bacteria/'+leaf_name+'/allantiporter.png', './named_pngs/'+taxonomy_converter[leaf_name].split(' ')[0]+"_"+leaf_name+'.png')
                    #outfile.write(leaf_name+","+taxonomy_converter[leaf_name]+",pass\n")
                else:
                    pass
                    #print counter, 'Fail'
                    #outfile.write(leaf_name+","+taxonomy_converter[leaf_name]+",fail\n")
            t2.render("phyla_trees/archaea_"+str(phylum)+"_tree.pdf")
