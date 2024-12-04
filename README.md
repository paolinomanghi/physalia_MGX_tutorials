## physalia_metagenomics_2024
#### lecture 1 - preprocessing

#### Download (and install) anaconda
```
wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh
path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/tools/anaconda3course/bin"
```

#### Raw data pre-processing (folder "1_pre-processing")
```
source ${path}/activate
conda create -n trimmomatic -c bioconda trimmomatic
conda create -n bowtie2 -c bioconda bowtie2
conda create -n samtools -c bioconda samtools
```

#### Getting fastq example files "seq_1.fastq.gz" and "seq_2.fastq.gz" from https://github.com/biobakery/biobakery/wiki/kneaddata
```
wget https://github.com/biobakery/kneaddata/files/4703820/input.zip
```

#### Define variable "s" with the sampleID
```
s="seq"
```

#### Run trimmomatic
```
source ${path}/activate trimmomatic

trimmomatic PE -threads 8 -phred33 -trimlog ${s}_trimmomatic.log ${s}_1.fastq.gz ${s}_2.fastq.gz \
${s}_filtered_1.fastq.gz ${s}_unpaired_1.fastq.gz ${s}_filtered_2.fastq.gz ${s}_unpaired_2.fastq.gz \
ILLUMINACLIP:${path}/../envs/trimmomatic/share/trimmomatic/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:75

for i in *.gz; do echo -ne "${i}\t"; zcat "$i" | wc -l; done
```

#### Getting human genome and generate bowtie2 index
#### Getting the file GCF_009914755.1_T2T-CHM13v2.0.fna from https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/
```
source ${path}/activate bowtie2

#bowtie2-build human_genome/GCF_009914755.1_T2T-CHM13v2.0.fna human_genome/GCF_009914755.1_T2T-CHM13v2.0 ### DON'T RUN IT! IT TAKES A FEW HOURS TO BE EXECUTED

bowtie2 -x human_genome/GCF_009914755.1_T2T-CHM13v2.0 -1 ${s}_filtered_1.fastq.gz -2 ${s}_filtered_2.fastq.gz -S ${s}.sam --very-sensitive-local -p 8 > ${s}_bowtie2.log 2>&1

source ${path}/activate samtools

samtools view -bS ${s}.sam > ${s}.bam
samtools view -b -f 12 -F 256 ${s}.bam > ${s}.bothunmapped.bam
samtools sort -n -m 5G -@ 2 ${s}.bothunmapped.bam -o ${s}.bothunmapped.sorted.bam
samtools fastq ${s}.bothunmapped.sorted.bam -1 >(gzip > ${s}_filtered.final_1.fastq.gz) -2 >(gzip > ${s}_filtered.final_2.fastq.gz) -0 /dev/null -s /dev/null -n
#rm ${s}.sam; rm ${s}.bam; rm ${s}.bothunmapped.bam; rm ${s}.bothunmapped.sorted.bam ### IF YOU WANT TO REMOVE THE INTERMEDIATE FILES

for i in *.gz; do echo -ne "${i}\t"; zcat "$i" | wc -l; done
```

##### -- end of lecture 1 - preprocessing
#### Lecture 2 - MetaPhlAn profiling
```
path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/tools/anaconda3course/bin"
```

#### MetaPhlAn for taxonomic profiling (folder "2_metaphlan")
```
source ${path}/activate
conda create -n mpa -c conda-forge -c bioconda python=3.7 metaphlan=4.1.0
source ${path}/activate mpa
```
 
#### Getting example files (6 fasta files) from https://github.com/biobakery/biobakery/wiki/metaphlan4
```
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014476-Supragingival_plaque.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014494-Posterior_fornix.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014459-Stool.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014464-Anterior_nares.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014470-Tongue_dorsum.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014472-Buccal_mucosa.fasta.gz

s="SRS014476-Supragingival_plaque"
metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --stat_q 0.1 --nproc 8 
# If this is your first time running MetaPhlAn, the first step involves downloading, decompressing, and indexing the MetaPhlAn marker gene database. This process may take ~30 minutes, but only needs to be performed once.

s="SRS014494-Posterior_fornix"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --stat_q 0.1 --nproc 8
s="SRS014459-Stool"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --stat_q 0.1 --nproc 8
s="SRS014464-Anterior_nares"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --stat_q 0.1 --nproc 8
s="SRS014470-Tongue_dorsum"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8
s="SRS014472-Buccal_mucosa"; metaphlan ${s}.fasta.gz --input_type fasta --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8

merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
```

#### Getting another example file (the fastq file SRS013951.fastq.bz2) from here: http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS013951.fastq.bz2

s="SRS013951";
metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --stat_q 0.1 --nproc 8
metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}_unclas.bowtie2.bz2 --samout ${s}_unclas.sam.bz2 -o ${s}_unclas_profile.txt --stat_q 0.1 --nproc 8 --unclassified_estimation
metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}_sub.bowtie2.bz2 --samout ${s}_sub.sam.bz2 -o ${s}_sub_profile.txt --stat_q 0.1 --nproc 8 --subsampling 1000
```

#### Generate heatmap with hclust2
```
conda create -n hclust2 -c bioconda hclust2 python=2.7
source ${path}/activate hclust2
grep -E "s__|SRS" merged_abundance_table.txt | grep -v "t__" | sed "s/^.*|//g" | sed "s/SRS[0-9]*-//g" > merged_abundance_table_species.txt

hclust2.py \
-i merged_abundance_table_species.txt \
-o metaphlan4_abundance_heatmap_species.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 0.5 \
--log_scale \
--flabel_size 10 --slabel_size 10 \
--max_flabel_len 100 --max_slabel_len 100 \
--minv 0.1 \
--dpi 300
```

##### End of Lecture - MetaPhlAn profiling
#### Lecture 3 - What is GraPhlAn ?
```
path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/tools/anaconda3course/bin"
```

#### GraPhlAn for circular representations of taxonomic and phylogenetic trees (folder "3_graphlan")
```
source ${path}/activate
conda create -n graphlan -c bioconda graphlan
source ${path}/activate graphlan
```

#### Getting example files from https://github.com/biobakery/graphlan/tree/master/examples/guide
```
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/guide.txt
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/annot_0.txt
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/annot_1.txt
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/annot_2.txt
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/guide/annot_3.txt

graphlan.py guide.txt step_0.png --dpi 300 --size 3.5
graphlan.py guide.txt step_0.svg --dpi 300 --size 3.5

graphlan_annotate.py --annot annot_0.txt guide.txt guide_1.xml
graphlan.py guide_1.xml step_1.png --dpi 300 --size 3.5
graphlan.py guide_1.xml step_1.svg --dpi 300 --size 3.5

graphlan_annotate.py --annot annot_1.txt guide_1.xml guide_2.xml
graphlan.py guide_2.xml step_2.png --dpi 300 --size 3.5
graphlan.py guide_2.xml step_2.svg --dpi 300 --size 3.5

graphlan_annotate.py --annot annot_2.txt guide_2.xml guide_3.xml
graphlan.py guide_3.xml step_3.png --dpi 300 --size 3.5
graphlan.py guide_3.xml step_3.svg --dpi 300 --size 3.5

graphlan_annotate.py --annot annot_3.txt guide_3.xml guide_4.xml
graphlan.py guide_4.xml step_4.png --dpi 300 --size 3.5 --pad 0.0
graphlan.py guide_4.xml step_4.svg --dpi 300 --size 3.5 --pad 0.0
```

#### Getting another example (PhyloPhlAn) from https://github.com/biobakery/graphlan/tree/master/examples/PhyloPhlAn
```
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/PhyloPhlAn/ppa_tol.xml 
wget -O annot_PhyloPhlAn.txt https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/PhyloPhlAn/annot.txt

graphlan_annotate.py --annot annot_PhyloPhlAn.txt ppa_tol.xml ppa_tol.annot.xml 
graphlan.py ppa_tol.annot.xml ppa_tol.png --dpi 200 --size 15 --pad 0.6
```

#### Getting another example (HMP_tree) from https://github.com/biobakery/graphlan/tree/master/examples/HMP_tree
```
wget https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/HMP_tree/hmptree.xml
wget -O annot_HMP_tree.txt https://raw.githubusercontent.com/biobakery/graphlan/refs/heads/master/examples/HMP_tree/annot.txt

graphlan_annotate.py --annot annot_HMP_tree.txt hmptree.xml hmptree.annot.xml 
graphlan.py hmptree.annot.xml hmptree.png --dpi 150 --size 14 
```

##### End of lecture 4 - What is GraPhlAn
#### Lecture 5 - What is StrainPhlAn 
```
path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/tools/anaconda3course/bin"
```

#### StrainPhlAn for strain-level profiling (folder "3_strainphlan")
```
source ${path}/activate mpa
```

#### Getting example files (6 fastq files) from https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4.1
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS013951.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS014613.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS019161.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS022137.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS055982.fastq.bz2
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/SRS064276.fastq.bz2

s="SRS013951"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8
s="SRS014613"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8
s="SRS019161"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8
s="SRS022137"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8
s="SRS055982"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8
s="SRS064276"; metaphlan ${s}.fastq.bz2 --input_type fastq --bowtie2out ${s}.bowtie2.bz2 --samout ${s}.sam.bz2 -o ${s}_profile.txt --nproc 8

sample2markers.py -i *.sam.bz2 -o ./ -n 8

mkdir -p db_markers
extract_markers.py -c t__SGB1877 -o db_markers/
```

#### Getting a reference genome ("GCF000273725")
```
mkdir -p reference_genomes
wget -P reference_genomes/ http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/reference_genomes/G000273725.fna.bz2
```

#### Build the multiple sequence alignment and the phylogenetic tree:
```
mkdir -p output
strainphlan -s *.json.bz2 -m db_markers/t__SGB1877.fna -r reference_genomes/G000273725.fna.bz2 -o output -c t__SGB1877 -n 8
```

#### Getting the metadata file ("metadata.txt")
```
wget http://cmprod1.cibio.unitn.it/biobakery4/github_strainphlan4/fastq/metadata.txt
add_metadata_tree.py -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre -f metadata.txt -m subjectID --string_to_remove .fastq.bz2

source ${path}/activate graphlan
${path}/../envs/mpa/bin/plot_tree_graphlan.py -t output/RAxML_bestTree.t__SGB1877.StrainPhlAn4.tre.metadata -m subjectID
```
##### End of Lecture 4 - what is StrainPhlAn 
#### Lecture 5 - What is PanPhlAn 
```
path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/tools/anaconda3course/bin"
```

#### PanPhlAn for Pangenome-based Phylogenomic Analysis (folder "5_panphlan")
```
source ${path}/activate

conda create -n panphlan -c bioconda panphlan
conda install -c conda-forge matplotlib

source ${path}/activate panphlan
```

#### Getting fastq example files from https://github.com/SegataLab/panphlan/wiki/Tutorial-3_0
```
wget https://www.dropbox.com/s/oi26jg0v7ktlavc/panphlan_tutorial_samples.tar.bz2
tar -xvjf panphlan_tutorial_samples.tar.bz2

panphlan_download_pangenome.py -i Eubacterium_rectale -o ./

mkdir -p map_results
s="CCMD34381688ST-21-0"
panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
```

#### The same script must be run for the other samples (fastq files)
```
s="G78505"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
s="G88884"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
s="G88970"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
s="G89027"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
s="H2M514903"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
s="HD-1"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
s="T2D-063"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8
s="T2D-105"; panphlan_map.py -i samples_fastq/${s}.fastq --indexes Eubacterium_rectale/Eubacterium_rectale -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv -o map_results/${s}_erectale.tsv --nproc 8

panphlan_profiling.py -i map_results/ --o_matrix ./result_profile_erectale.tsv -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv --o_covplot ./erectale_covplot
panphlan_profiling.py -i map_results/ --o_matrix ./result_profile_erectale_annotation.tsv -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv --func_annot Eubacterium_rectale/panphlan_Eubacterium_rectale_annot.tsv --field 2
panphlan_profiling.py -i map_results/ --o_matrix ./result_profile_erectale_annotation_withref.tsv -p Eubacterium_rectale/Eubacterium_rectale_pangenome.tsv --func_annot Eubacterium_rectale/panphlan_Eubacterium_rectale_annot.tsv --field 2 --add_ref
```

#### Generate the file "metadata_erectale.txt" by taking information from https://github.com/SegataLab/panphlan/wiki/Tutorial-3_0
#### Script in R to generate Heatmap and MDS
##### End of Lecture 5 - What is PanPhlAn 

#### Lecture 6 - HUMAnN 3
```
path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/tools/anaconda3course/bin"
```

#### HUMAnN for profiling the abundance of microbial metabolic pathways and other molecular functions (folder "6_humann")
```
source ${path}/activate

conda create -n humann -c bioconda python=3.9
source ${path}/activate humann
conda install -c biobakery humann
```

#### Test the local HUMAnN environment
```
humann_test
humann_config
```

#### Download databases
```
humann_databases --download chocophlan full humann_dbs --update-config yes
humann_databases --download uniref uniref90_diamond humann_dbs --update-config yes
humann_databases --download utility_mapping full humann_dbs --update-config yes
```

#### Getting example of fastq file from EBI
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/096/SRR15408396/SRR15408396.fastq.gz
```

#### Run humann to perform functional profiling
```
s="SRR15408396"
humann --input ${s}.fastq.gz --output ${s} --threads 8
```

#### Manipulating HUMAnN output tables
```
humann_renorm_table -i ${s}/${s}_genefamilies.tsv -o ${s}/${s}_genefamilies-relab.tsv -u relab
humann_renorm_table -i ${s}/${s}_pathabundance.tsv -o ${s}/${s}_pathabundance-relab.tsv -u relab
```

#### Regrouping genes to other functional categories
```
humann_regroup_table -i ${s}/${s}_genefamilies-relab.tsv -o ${s}/${s}_rxn-relab.tsv --groups uniref90_rxn
```

#### Getting another fastq file from EBI and run humann
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/098/SRR15408398/SRR15408398.fastq.gz
s="SRR15408398"
humann --input ${s}.fastq.gz --output ${s} --threads 8
humann_renorm_table -i ${s}/${s}_genefamilies.tsv -o ${s}/${s}_genefamilies-relab.tsv -u relab
humann_renorm_table -i ${s}/${s}_pathabundance.tsv -o ${s}/${s}_pathabundance-relab.tsv -u relab
```

#### Merge multiple outputs
```
mkdir -p merged
cp SRR15408396/SRR15408396_genefamilies-relab.tsv merged/
cp SRR15408398/SRR15408398_genefamilies-relab.tsv merged/
cp SRR15408396/SRR15408396_pathabundance-relab.tsv merged/
cp SRR15408398/SRR15408398_pathabundance-relab.tsv merged/
cp SRR15408396/SRR15408396_pathcoverage.tsv merged/
cp SRR15408398/SRR15408398_pathcoverage.tsv merged/

humann_join_tables -i merged -o merged/merged_genefamilies-relab.tsv --file_name genefamilies-relab
humann_join_tables -i merged -o merged/merged_pathabundance-relab.tsv --file_name pathabundance-relab
humann_join_tables -i merged -o merged/merged_pathcoverage.tsv --file_name pathcoverage

```

##### End of Lecture 6 - HUMAnN 3
#### Lecture 7 - Metagenomic assembly
```
path="/shares/CIBIO-Storage/CM/scratch/users/e.pasolli/tools/anaconda3course/bin"
```

#### Megahit for de novo metagenomic assembly (folder "7_assembly")
```
source ${path}/activate

conda create -n megahit -c bioconda megahit
source ${path}/activate megahit
```

#### Getting example of fastq file from https://github.com/voutcn/megahit/wiki/An-example-of-real-assembly
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/SRR341725/SRR341725_[12].fastq.gz
```

#### Run megahit to generate metagenomic assembly
```
s="SRR341725"
megahit -1 ${s}_1.fastq.gz -2 ${s}_2.fastq.gz -o ${s}.megahit_asm -t 8
```

#### Do some post-processing on the contigs file
```
conda install -c conda-forge biopython
conda install -c anaconda pandas

python megahit2spades.py SRR341725.megahit_asm/final.contigs.fa SRR341725.megahit_asm/contigs.fasta
python filter_contigs.py SRR341725.megahit_asm/contigs.fasta SRR341725.megahit_asm/contigs_filtered.fasta
python filter_contigs.py SRR341725.megahit_asm/contigs.fasta SRR341725.megahit_asm/contigs_filtered_50000.fasta -l 50000
```

#### Run prokka for rapid prokaryotic genome annotation
```
conda create -n prokka -c bioconda prokka
source ${path}/activate prokka
prokka --outdir ${s}_prokka --centre CDC --compliant --cpus 8 SRR341725.megahit_asm/contigs.fasta

```
