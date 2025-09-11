## Comparative Genomics and Phylogenomics of Rhabdiditae Nematodes ###

***Gene Annotation***

All gene annotation work was done on the FIU HPC.

Before beginning, please complete up through repeatmasker on [repeatAnnotations.md](https://github.com/ToriEggers/Rhab_TE_Gene_GenomeEvolution/blob/main/repeatAnnotations.md)

All preparation for RNA evidence for BRAKER3 was done in scratch directory, which deletes files last touched over 30 days.
Preparation of RNA evidence includes downloading RNA data from NCBI (aimed for total RNA seq illumina data), alignment with STAR, and indexing/sorting with samtools. Once the bam files were created successfully, they were moved to the directory containing the braker.sif file.

<details>

<summary><b>Download RNA data</b></summary>

The file: RNA_accessions.txt
contains 2 columns separated by a tab, ID and ACCESSION. It is the input to the following script which uses sratoolkit v.3.0.0 to (1)prefetch the accession of the RNA seq data and (2)fastq-dump the RNA reads into a directory named by the ID.

```
vi RNA_accessions.txt
```

copy/paste the following:
```
AF72	ERR4264631
BOV	ERR3610811
DF5000	ERR10787778
DF5013	ERR13623740
DF5033	SRR18615345
EM464	SRR5837623
JU1373	SRR12623043
JU2817	ERR13319083
JU3284	ERR13319072
JU75	ERR13319099
LJ9110	SRR28370082
NKZ352	DRR252169
PB127	SRR18615296
QG2083	SRR25478168
SB194	ERR13319073
BRC20483	ERR13319078
QG555	ERR13623727
NIC564	ERR13623728
RS0144	ERR2235011
RS5460	ERR2019977
EX	ERR2652888
RS5133	ERR2235013
AF16	SRR30148319
APS25	ERR13319064
APS4	ERR3150274
APS7	ERR11177444
Aroian	SRR609884
BAKE	SRR20115893
CB4856	SRR29887579
CEW1	SRR19570657
CFB2252	ERR12997162
DF5070	ERR13319086
DF5081	SRR19570648
DF5083	ERR13319102
DF5112	SRR19570652
DF5120	ERR10787777
DF5173	ERR13623733
EG5942	ERR13319091
ISE	SRR27211025
JU1182	ERR13319057
JU1286	ERR1059227
JU1382	ERR10787776
JU1421	SRR3031146
JU1771	ERR1039279
JU1809	ERR13319052
JU1904	ERR13319080
JU1917	ERR13319067
JU1968	ERR1055674
JU2083	ERR1018630
JU2190	SRR8869244
JU2788	ERR13623725
JU2809	ERR13319082
JU800	ERR13319092
KANDY	ERR11518137
LIT	SRR30363359
MIMR	ERR2023646
MONO	ERR690851
N2	SRR32732686
NIC203	SRR12341282
NIC58	SRR12623044
NKZ35	DRR481199
OM	ERR13623726
PDL0010	ERR13319096
PF1305	ERR10787776
PLIC	ERR13623716
PS1010	ERR13319081
PS1017	ERR10787778
PS2068	ERR10787775
PX356	SRR5837623
PX439	SRR5837881
PX506	SRR10276661
PX534	SRR5831583
SX3368	SRR12868616
TWN1964	ERR13319056
TWN1984	ERR13319060
BRC20456	ERR13319076
CP168	SRR3031146
JU3283	ERR13319098
QG2077	ERR13362831
TT01	SRR29886265
VIVI	SRR1021578
TRI	SRR8636392
OST	SRR2567544
PS312	SRR23038884
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
...
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
vi isoseq_align.sh
```

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

cd to RhabditinaPhylogeny

Get container:

```
module load singularity-3.8.7
module load proxy
singularity build braker3.sif docker://teambraker/braker3:latest
```

Move the refseq_db.faa from nematoda_odb10 to the same directory as braker3.sif


```
vi fiu_array_singularity_braker3.sh
```

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
sleep $((SLURM_ARRAY_TASK_ID * 25))

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
<summary><b>Predict Genes with BRAKER3 using isoseq data</b></summary>

Following directions from [BRAKER3](https://github.com/Gaius-Augustus/BRAKER) github.

Get container:

```
module load singularity-3.8.7
module load proxy
singularity build braker3_lr.sif docker://teambraker/braker3:isoseq
```
Also make sure that you've moved PX506_isoseq.bam and refseq_db.faa to the same directory containing braker3_lr.sif

```
vi fiu_singularity_braker3_isoseq.sh
```

```
#!/bin/bash

#SBATCH --job-name=braker3_isoseq
#SBATCH --output=./logs/braker3_isoseq.%j.out
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
export BRAKER_SIF=/home/data/jfierst/veggers/RhabditinaPhylogeny/braker3_lr.sif

cat ./RhabditinaPhylogeny_repeatmasker/PX506/PX506.masked | cut -f 1 -d " " > PX506.masked

#organize and remove working directory if it already exists
wd=./RhabditinaPhylogeny_braker3/PX506_braker3

if [ -d $wd ]; then
    rm -r $wd
fi

#run braker
singularity exec  -B ${PWD}:${PWD}  ${BRAKER_SIF} braker.pl --genome=PX506.masked --prot_seq=refseq_db.faa --bam=PX506_isoseq.bam --workingdir=${wd} --GENEMARK_PATH=${ETP}/gmes --AUGUSTUS_CONFIG_PATH=/home/veggers/.augustus --threads 8 --softmasking --busco_lineage nematoda_odb10
```

</details>

<details>
<summary><b>Predict Genes with BRAKER2</b></summary>

BOV, LJ9110, LIT, JU2585, NIC534, and NIC394 were annotated using BRAKER2 due to difficulties with BRAKER3. LIT subsequently failed with BRAKER2 as well and was left out of all analyses.

The same container built and used for braker3 can also be used for braker2. They are the same program, braker2 just doesn't incorporate RNA evidence. Remember to copy refseq_db.faa from nematoda_odb10 to the same directory containing braker3.sif.

```
vi  fiu_array_singularity_braker2.sh
```

```
#!/bin/bash
#SBATCH --job-name=braker2
#SBATCH --output=./logs/braker2.%j.%N.out
#SBATCH --array=1
#SBATCH --account=iacc_jfierst
#SBATCH --cpus-per-task=1
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
SPECIES=$(sed "${SLURM_ARRAY_TASK_ID}q;d" braker2_list.txt)

echo "$SPECIES"

cp ./RhabditinaPhylogeny_repeatmasker/${SPECIES}/*.masked ./${SPECIES}.masked

#sleep for a few seconds so braker doesn't try to name multiple species the same (causes write permission failures otherwise)
sleep $((SLURM_ARRAY_TASK_ID * 25))

#organize and remove working directory if it already exists
wd=./RhabditinaPhylogeny_braker2/${SPECIES}_braker2

if [ -d $wd ]; then
    rm -r $wd
fi

#run braker
singularity exec  -B ${PWD}:${PWD}  ${BRAKER_SIF} braker.pl --genome=${SPECIES}.masked --prot_seq=refseq_db.faa --workingdir=${wd} --GENEMARK_PATH=${ETP}/gmes --AUGUSTUS_CONFIG_PATH=/home/veggers/.augustus --threads 8 --softmasking --busco_lineage nematoda_odb10
```

</details>

<details>
<summary><b>Get longest isoform with AGAT</b></summary>

</details>

<details>
<summary><b>Predict Function with InterProScan</b></summary>

</details>

<details>
<summary><b>Predict Orthologues with Orthofinder</b></summary>

</details>
