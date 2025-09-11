## Comparative Genomics and Phylogenomics of Rhabdiditae Nematodes ###

***Gene Annotation***

All gene annotation work was done on the FIU HPC.

All preparation for RNA evidence for BRAKER3 was done in scratch directory, which deletes files last touched over 30 days.
Preparation of RNA evidence includes downloading RNA data from NCBI (aimed for total RNA seq illumina data), alignment with STAR, and indexing/sorting with samtools. Once the bam files were created successfully, they were moved to the directory containing the braker.sif file.

<details>

<summary><b>Download RNA data</b></summary>

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

Note: NIC394 does not have RNA seq available on NCBI

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
<summary><b>Align illumina RNA reads to genome with STAR</b></summary>

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

<details>

<summary><b>Align isoseq RNA reads to genome with minimap</b></summary>

Following directions from [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) github.

PX506 only had isoseq RNA reads available, which requires a bit of a different alignment process:

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=./logs/isoseq_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH --mail-type=ALL

module load samtools-1.15.1-gcc-8.2.0
module load minimap2-2.24

minimap2 -t 40 -ax splice:hq -uf /home/data/jfierst/veggers/RhabditinaPhylogeny/RhabditinaPhylogeny_repeatmasker/PX506/PX506.masked PX506_isoseq.fa > PX506_isoseq.sam     
samtools view -bS --threads 40 PX506_isoseq.sam -o PX506_isoseq.bam
```

move PX506_isoseq.bam to /home/data/jfierst/veggers/RhabditinaPhylogeny/.

</details>

<details>
<summary><b>Predict Genes with BRAKER3</b></summary>

```
#!/bin/bash

#SBATCH --job-name=braker3
#SBATCH --output=./logs/braker3.%j.out
#SBATCH --array=1
#SBATCH --account=iacc_jfierst
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --partition=highmem1
#SBATCH --qos=highmem1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vegge003@fiu.edu

#load modules
module load singularity-3.8.7
module load proxy

#export paths
export BRAKER_SIF=/home/data/jfierst/veggers/RhabditinaPhylogeny/braker3.sif

#set species variable
SPECIES=$(sed "${SLURM_ARRAY_TASK_ID}q;d" braker3.txt)

echo "$SPECIES"

cat ./RhabditinaPhylogeny_repeatmasker/${SPECIES}/${SPECIES}.masked | cut -f 1 -d " " > ${SPECIES}.masked

cp /scratch/jfierst/tori/${SPECIES}_STAR/${SPECIES}_Aligned.out.bam ./${SPECIES}_Aligned.out.bam

#sleep for a few seconds so braker doesn't try to name multiple species the same (causes write permission failures otherwise)
#sleep $((SLURM_ARRAY_TASK_ID * 25))

#organize and remove working directory if it already exists

wd=./RhabditinaPhylogeny_braker3/${SPECIES}_braker3

if [ -d $wd ]; then
    rm -r $wd
fi

#run braker
singularity exec  -B ${PWD}:${PWD}  ${BRAKER_SIF} braker.pl --genome=${SPECIES}.masked --prot_seq=refseq_db.faa --bam=${SPECIES}_Aligned.out.bam --workingdir=${wd} --GENEMARK_PATH=${ETP}/gmes --AUGUSTUS_CONFIG_PATH=/home/veggers/.augustus --threads 8 --softmasking --busco_lineage nematoda_odb10
```

BRAKER renames things the same thing and so if the jobs aren't spaced out enough you'll get an error about a species directory not existing or being writable. That's why I added the sleep command in the script above, however it does waste some computational resources. sloppy fix but it works sometimes. You will still get an error even with the sleep. Just modify the BRAKER3.txt file with the genomes that failed due to this issue and rerun, they'll all work eventually.

BOV, LJ9110, and LIT will all fail because there is not enough RNA seq evidence, likely because these are parasitic species and so the RNA available is from the infected organism rather than the nematodes themselves.

JU2585 and NIC534 fail because there is too little intron evidence.

</details>

<details>
<summary><b>Predict Genes with BRAKER2</b></summary>

</details>

<details>
<summary><b>Get longest isoform with AGAT</b></summary>

</details>

<details>
<summary><b>Predict Function with InterProScan</b></summary>

</details>
