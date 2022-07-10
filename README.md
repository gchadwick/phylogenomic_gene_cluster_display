# Phylogenomic Gene Cluster Display
This series of python scripts takes the genomes collected by the Genome Taxonomy Database (GTDB), extracts proteins from each genome based on hits to a Hidden Markov Model, and then displays these genes and their genomic neighborhood on a phylogenomic tree.  These scripts have been a work in progress, helping with the analysis done in the following publications, as well as several ongoing projects:

Chadwick, G. L., Hemp, J., Fischer, W. W., & Orphan, V. J. (2018). Convergent evolution of unusual complex I homologs with increased proton pumping capacity: energetic and ecological implications. The ISME journal, 12(11), 2668.

Yu, H., Chadwick, G. L., Lingappa, U. F., & Leadbetter, J. R. (2022). Comparative Genomics on Cultivated and Uncultivated Freshwater and Marine “Candidatus Manganitrophaceae” Species Implies Their Worldwide Reach in Manganese Chemolithoautotrophy. mBio, 13(2), e03421-21.

<img src="https://github.com/gchadwick/gchadwick.github.io/blob/master/images/fulls/17.png?raw=true" width="1000">

## Step 0: Get some things
This analysis relies on a variety of data products from the GTDB. You can find the latest release [here](https://data.gtdb.ecogenomic.org/releases/latest/). You will also need a gene calling/annotation tool, here I have used [prokka](https://github.com/tseemann/prokka). Other requirements are [Biopython](https://biopython.org/), [ete3](https://github.com/etetoolkit/ete), and [SimpleHMMER](https://github.com/minillinim/SimpleHMMER).
## Step 1: Gene Calls and Annotations
GTDB as of this work only supplies genomes as nucleotide fasta files, not genbank files with the associated gene calls and annotations. The first step of this pipeline is to iterate through all of the genomes and do this. Using the guidance of the prokka devs [here](https://github.com/tseemann/prokka/issues/187):

    for F in *_genomic.fna; do 
    N=$(basename $F _genomic.fna) ; 
    prokka --locustag $N --outdir $N --prefix $N  $F ;
    done

This iterates through all of the genomes in a given directory, assuming that GTDB is still using the "_genomic.fna" suffix on their genome files. This process might take an excessive ammount of time depending on your machine. You could try using the --noanno flag, to just have prokka do gene calls and not annotations which is a longer process. The HMM searches in the subsequent steps does not rely on the prokka annotations. This will make a directory for each individual genome with the prokka outputs.
## Step 2: Find all genes matching a certain Hidden Markov Model
For this example workflow I will search for all of the genes belonging to to broad family of proton transporters represented by the "Proton_antipo_M (PF00361)" pfam. For whatever gene you are interested in you can either make your own HMM from an alignment, or find existing HMMs from a variety of different databases.  step2.py simply iterates through all the genomes annotated in step 1 doing an HMM search for the input HMM and saves the output in the genome directory.
## Step 3: Identify the common associated genes for display purposes
This multi-part step uses a pre-defined set of HMMs to search against all of the genomes. These genes are thought of as accessory or associated genes for the main PF00361. Depending on the host of associated genes with the PF00361 homologs, we would predict a different possible function of the gene cluster. Each HMM can be colored differently in the final plots allowing for quick vizualization of operon changes through evolution
### Step 3.1:
This step does the initial itration through all of the selected HMMs and genomes and saves the hits.
### Step 3.2:
This step uses the HMMParser function from simplehmmer to parse the output of step 3.1 into a more useful format.
### Step 3.3:
This step uses the output of 3.2 to extract gene identifiers as well as sequences from protein files and store them in new fasta files.
## Step 4: Make image files for all of the identified gene clusters
The information produced in the previous steps about gene identifiers, HMM hits, etc. are all combined in this step, along with the gene size and orientation information from the genbank files to produce .png files of the gene contexts with functions encoded in colors. The syntax for this code come from GenomeDiagram and SeqFeature/FeatureLocation from Bio.Graphics and Bio.SeqFeature packages (more info can be found in the BioPython documentation [here](https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec340)). Lots of specifics to an application need to be hard-coded in here, like the color of the different HMM groups, size of the regions to extract around the specific hits, etc.
## Step 5: Make combined phylogenomic tree and gene cluster figures
The images from step 4 are now combined with the phylogenomic tree calculated with each GTDB release and the standardized taxonomy string to generate the final figure that is a tree where every leaf is labeled with the unique genome id, taxonomy string, and the image of all the gene clusters that contain proton pumping subunits.
## Step 6 (optional): Phylum or sub-phylum trees
Running step 5 on the entire bacterial tree from GTDB will work but produces and absurdly large pdf file that are functionally unusable. The optional step 6 allows you to iterate through various phylogenetic levels and make sub-trees that are more manageable files.