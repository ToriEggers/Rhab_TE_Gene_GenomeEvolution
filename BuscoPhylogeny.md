## Comparative Genomics and Phylogenomics of Rhabdiditae Nematodes ###

***Making a Phylogeny with BUSCO shared single copy orthologs***

Code/methods derived from [Manni et al. 2021](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.323). For more information please read their paper.

To download the data and set up the directories see: [Data.md]()

<details>
    
<summary><b>1. Run BUSCO on all genomes in parallel</b></summary>

```
vi busco.sh
```

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_busco.log
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 40

max_jobs=10
job_count=0

while read -r line; do
    busco -f -c 4 -m genome -i ./RhabditinaPhylogeny_NCBI/"${line}"/*.fna -o busco_"${line}" --offline --lineage_dataset /home/data/jfierst/veggers/nematoda_odb10 &

    job_count=$((job_count + 1))

        if [ "$job_count" -ge "$max_jobs" ]; then
                wait
                job_count=0
        fi

done < busco_list.txt
```
busco_list.txt is a list of all nematode species/strains I'm interested in. It looks like:
```
TWN1964
TWN1984
SX3368
SB194
QG2083
PX534
PX506
PX439
PX356
PS2068
PS1017
```

```
sbatch busco.sh
```
Takes less than a day to get through all 70 genomes.

Make sure that all the buscos worked. If they didn't why? probably because fna files are still gzipped or something. 

</details> 

<details>
    <summary><b>2. Concatenate BUSCO summaries for all genomes</b></summary>
    
```
vi busco_summary.sh
```

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_busco_summary.log

WORKING_DIR=/home/data/jfierst/veggers/RhabditinaPhylogeny/

cd ${WORKING_DIR}

#make header line 
echo -e "ID\tbusco\tsingle_copy_count" > busco_summary.txt

#loop through busco_list.txt 
while read -r line; do
    cd RhabditinaPhylogeny_Buscos/busco_${line}/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/

#list all single copy orthologs into a list (note: the paths might be messed up here. If you run it and don't see the output file, search around a bit)
    ls *.fna | sort > ./../../../${line}_single_copy_list.txt
    cd ${WORKING_DIR}

#count the number of single copy orthologs found
    single_copy_count=$(wc -l RhabditinaPhylogeny_Buscos/busco_${line}/${line}_single_copy_list.txt | awk '{print $1}')
    ID=${line}

#extract the summary line from short_summary.txt 
    busco=$(grep "C:" RhabditinaPhylogeny_Buscos/busco_${line}/short_summary*.txt)

#append all calculated values to busco_summary.txt
    echo -e "${ID}\t${busco}\t${single_copy_count}" >> busco_summary.txt

done < busco_list.txt
```
```
sbatch busco_summary.sh
```

The output will be busco_summary.txt and looks like:
```
ID      busco   single_copy_count
AF16            C:98.1%[S:97.8%,D:0.3%],F:1.1%,M:0.8%,n:3131            3063
AF72            C:64.5%[S:44.2%,D:20.3%],F:5.8%,M:29.7%,n:3131          1383
APS25           C:74.1%[S:73.8%,D:0.3%],F:5.8%,M:20.1%,n:3131           2310
APS4            C:76.3%[S:76.0%,D:0.3%],F:6.4%,M:17.3%,n:3131           2379
```
I used this to sure that the busco is good quality, although good quality is sort of subjective with organism of study. Here I define good quality as not too many duplicates. 3 species were removed from the analysis for having more than 5% duplicate BUSCOs, like you see in AF72 in the above output.

Another output from the above code is RhabditinaPhylogeny_Buscos/busco_${line}/${line}_single_copy_list.txt (replace ${line} with the species name. Thus, you would generate 70 of these files (one for each species). This output is important for the next step and should look like:
```
0at6231.fna
10010at6231.fna
10011at6231.fna
10018at6231.fna
10021at6231.fna
10025at6231.fna
10032at6231.fna
1003at6231.fna
10044at6231.fna
10045at6231.fna
10059at6231.fna
10063at6231.fna
10068at6231.fna
10076at6231.fna
10081at6231.fna
```

</details>

<details>
    <summary><b>3. Concatenate all shared single copy orthologous sequences</b></summary>

```
vi busco_parse.sh
```
Note: busco_parse.txt and busco_list.txt are the same thing (a list of species/strain names), but I've removed the three species with poor BUSCO scores, thus busco_parse.txt has 67 lines while busco_list.txt has 70.

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_busco_parse.log

WORKING_DIR=/home/data/jfierst/veggers/RhabditinaPhylogeny

cd ${WORKING_DIR}

#prep loop #1 (this code specifies the first species in the list as the start and copies it's single_copy_list.txt to com.txt)
startFile=$(head -1 busco_parse.txt)
cp ${WORKING_DIR}/RhabditinaPhylogeny_Buscos/busco_${startFile}/${startFile}_single_copy_list.txt com.txt

#loop through species to find buscos common between the current species and the file com.txt. All matches between the two files are saved to temp.txt. temp.txt is then moved to com.txt and the loop repeats with the next species.
while read -r file; do
    comm -12 com.txt ${WORKING_DIR}/RhabditinaPhylogeny_Buscos/busco_${file}/${file}_single_copy_list.txt > temp.txt
    mv temp.txt com.txt
done < busco_parse.txt

#make directory for organization 
mkdir -p busco_msa

#loop #2 makes a directory for each BUSCO_ID, copies the species BUSCO_ID gene sequence to that directory, and adds the species name to the fasta header
while read -r BUSCO_ID; do
    mkdir -p ./busco_msa/${BUSCO_ID}
    cd busco_msa/${BUSCO_ID}

    while read -r species; do
        cp ${WORKING_DIR}/RhabditinaPhylogeny_Buscos/busco_${species}/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/${BUSCO_ID} .
        sed -e "s/^>/>${species}_/" ${BUSCO_ID} | sed 's/-/_/g' | sed 's/\./_/g' | sed 's/\:/_/g' > ${species}_${BUSCO_ID}
        rm ${BUSCO_ID}
    done < ./../../busco_parse.txt

    cat * > ${BUSCO_ID}.fasta
    rm *.fna
    cd ${WORKING_DIR}

done < com.txt
```
Your output is ./busco_msa/${BUSCO_ID}/${BUSCO_ID}.fasta, generating a fasta file for each shared single copy ortholog present in all genomes included in this study (about 400 of the possible 3131 nematode BUSCOs). Each file should contain 67 sequences, one for each species. If you head one of the files, it should look like:
```
>AF16_NC_013486_2_6886414_6898962
ATGATACGCTGGAAGTACGGAATTCACTACCTCATATGGCTCCTTCTCGTGCTGCATTTG
TCGACGTGTCAATCAGATTCCTCTCTGACGACGTCGGCCGAGCAGCACGAGTTGTTTGCC
ATCAAGAAGGACTCGTTGTCTCCGTGGTCTCAAATTCTTGTTTCATTACCACGGAGGCAT
CCTTTGTATCAGAGTTTTGCTGCCAAGATTCAAGATGTAACTGAGAATATGtcggATGAT
GTTAGAGACGCGAACAAAACATTCGTTTCGTCCGATGATTCTCCATATAACATCAGGATA
CATGCTCTGAAACCAGGACACCGATActcgattgcTATCCACGGCCAAAAAGATGGATCG
ACCTCCTTGATAAAAGAAGAATCGGTTGTTATGGACCCTCGGGCTCCCGACTTCCGATCA
ATAGATTCTGATATCCAGGTGGCAGAGCACAACATCACAATGAGAACAATTAAGAACGAT
TCTTACTTACAAGACTCTTTCTCAATTGAATATCGTCAGATTAACCCGGATAAGAAGTTT
```

</details>

<details>
    <summary><b>4. Multiple sequence alignment and trimming</b></summary>

<br>

mafft.sh was created with the help of [Daniel Morales](https://github.com/dmora127) (arrays are still a bit over my head)
```
vi mafft.sh
```

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_mafft.log
#SBATCH --array=1-499%25

WORKING_DIR=/home/data/jfierst/veggers/RhabditinaPhylogeny

cd ${WORKING_DIR}

module load mafft-7.221-gcc-8.2.0-y6cgezm

# Define your input/output directories
INPUT_DIR="/home/data/jfierst/veggers/RhabditinaPhylogeny/busco_msa"
OUTPUT_DIR="/home/data/jfierst/veggers/RhabditinaPhylogeny/busco_msa"

# Define your taxonList
BUSCO_ID=""

# Read all BUSCO IDs into an array
readarray -t BUSCO_IDS < "com.txt"

# Use SLURM_ARRAY_TASK_ID to pick the corresponding BUSCO ID
BUSCO_ID=${BUSCO_IDS[$((SLURM_ARRAY_TASK_ID-1))]}

# The SLURM_ARRAY_TASK_ID variable helps in handling different parts in parallel
INPUT_FILE="${INPUT_DIR}/${BUSCO_ID}/${BUSCO_ID}.fasta"
OUTPUT_FILE="${OUTPUT_DIR}/${BUSCO_ID}/${BUSCO_ID}_mafftAligned.fasta"

mafft-linsi --thread 8 ${INPUT_FILE} > ${OUTPUT_FILE}
```
```
sbatch mafft.sh
```
The output is an alignment file for each BUSCO single copy ortholog, located at /home/data/jfierst/veggers/RhabditinaPhylogeny/busco_msa/${BUSCO_ID}/${BUSCO_ID}_mafftAligned.fasta

Now trim the alignments, removing any uninformative gaps (gaps in 90% of sequences). I'm using ClipKit here, which was conda installed:

```
module load mamba/23.1.0-4
```

```
conda create -n clipkit
```

```
source activate clipkit
```

```
mamba install bioconda::clipkit
```
Type `clipkit -h` and the manual should appear if you've installed it correctly.

```
vi clipkit.sh
```
```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_clipkit.log

WORKING_DIR=/home/data/jfierst/veggers/RhabditinaPhylogeny

cd ${WORKING_DIR}

# load software
module load mamba/23.1.0-4
source activate clipkit

while read -r BUSCO_ID; do
    cd busco_msa/${BUSCO_ID}
    clipkit -m smart-gap ${BUSCO_ID}_mafftAligned.fasta
    cd ${WORKING_DIR}
done < com.txt

#make directory 
mkdir -p busco_msa/total

#move all trimmed alignments to a directory called total
while read -r BUSCO_ID; do
    cd busco_msa/${BUCO_ID}
    cp ${BUSCO_ID}_mafftAligned.fasta.clipkit ./../total/.
    cd ${WORKING_DIR}
done < com.txt

# fix names . . . again, this changes >AF16_NC_013486_2_6886414_6898962 to >AF16 for each of the trimmed alignments in the total directory
while read -r BUSCO_ID; do
    cd busco_msa/total/
    cut -d "_" -f 1 ${BUSCO_ID}_mafftAligned.fasta.clipkit > temp
    mv temp ${BUSCO_ID}_mafftAligned.fasta.clipkit
    cd ${WORKING_DIR}
done < com.txt
```
```
sbatch clipkit.sh
```
The output is a directory called total with all trimmed alignments for each BUSCO_ID.

</details>

<details>
<summary><b>5. AMAS and IQTREE</b></summary>

```
# get partition file with AMAS
#module load mamba/23.1.0-4
#source activate AMAS

#python3 /home/data/jfierst/veggers/programs/AMAS/amas/AMAS.py concat -c 40 -f fasta -d dna --part-format raxml -i ${WORKING_DIR}/busco_msa/total/*

# iqtree
module load iqtree-2-gcc-8.2.0

iqtree2 -s concatenated.out -spp partitions.txt -m MFP+MERGE -bb 1000 -alrt 1000 -nt 40
```
</details>
