## Comparative Genomics and Phylogenomics of Rhabdiditae Nematodes ###

***Gene Annotation***

All preparation for RNA evidence for BRAKER3 was done in scratch directory.
This includes downloading RNA data from NCBI (aimed for total RNA seq illumina data), alignment with STAR, and indexing/sorting with samtools. Once the bam files were created successfully, they were moved to the directory containing the braker.sif file.

<details>

<summary>Download RNA data</summary>

The file: RNA_accessions.txt
contains 2 columns separated by a tab, ID and ACCESSION. It is the input to the following script which uses sratoolkit v.3.0.0 to (1)prefetch the accession of the RNA seq data and (2)fastq-dump the RNA reads into a directory named by the ID.

```
vi RNA_accessions.txt
```

```
PS1010  ERR13319081
PS1017  ERR10787778
PS2068  ERR10787775
PX356   SRR5837623
PX439   SRR5837881
BOV     ERR3610811
```

```
vi download.sh
```
 
```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=./logs/RNAseq_download_%j.log

module load proxy #needed to connect to internet
module load sratoolkit-3.0.0

INPUT_FILE=RNA_accessions.txt #input file: ID ACCESSION (2 columns seperated by a tab)

while read -r line; do
        species=$(echo $line | awk '{print $1}') #set variables
        accession=$(echo $line | awk '{print $2}')

        prefetch ${accession} #sratoolkit commands
        fasterq-dump ${accession} -O ${species}

        rm -r ${accession} #delete prefetch created directory
done < ${INPUT_FILE}
```
</details>

<details>
<summary>Align RNA reads to genome with STAR</summary>

STAR.txt has a list of IDs:
```
vi STAR.txt
```

```
QG555
NIC564
RS0144
RS5460
RS5133
```

STAR.sh is an array script to generate a genome index of each species and then map the rna reads to that index. --outSAMstrandField intronMotif and --outSAMtype BAM Unsorted are required for latter input into BRAKER3
```
vi STAR.sh
```
```
#!/bin/bash

#SBATCH --account acc_jfierst
#SBATCH --qos standard
#SBATCH --partition HighMem1
#SBATCH --array=1-5
#SBATCH --output=./logs/STAR_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH --mail-type=ALL

module load miniconda3/24.7.1-none-none-mjgmhio
source activate STAR

SPECIES=$(sed "${SLURM_ARRAY_TASK_ID}q;d" STAR.txt)

echo "$SPECIES"
mkdir ${SPECIES}_STAR

# Generate genome index
STAR \
    --runThreadN 12 --runMode genomeGenerate --genomeDir ${SPECIES}_STAR \
    --genomeSAindexNbases 12 --genomeFastaFiles /home/data/jfierst/veggers/RhabditinaPhylogeny/RhabditinaPhylogeny_repeatmasker/${SPECIES}/${SPECIES}.masked

# Map the reads
STAR \
    --runThreadN 12 --runMode alignReads --genomeDir ${SPECIES}_STAR --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --twopassMode Basic \
--readFilesIn /scratch/jfierst/tori/${SPECIES}/${SPECIES}_1.fastq /scratch/jfierst/tori/${SPECIES}/${SPECIES}_2.fastq --out
FileNamePrefix /scratch/jfierst/tori/${SPECIES}_STAR/${SPECIES}_
```
</details>
