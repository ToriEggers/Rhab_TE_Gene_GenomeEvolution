#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%busco_msa.log

WORKING_DIR=/home/data/jfierst/veggers/RhabditinaPhylogeny

cd ${WORKING_DIR}

# msa for each gene
module load mafft-7.221-gcc-8.2.0-y6cgezm

while read -r gene; do
    cd busco_msa/${gene}
    linsi --thread 16 --maxiterate 1000 ${gene}.fasta > ${gene}_aligned.fasta
    cd ${WORKING_DIR}

done < com.txt


# trim alignments
module load mamba/23.1.0-4
source activate trimal

while read -r gene; do
    cd busco_msa/${gene}
    trimal -in ${gene}_aligned.fasta -keepheader -keepseqs -gappyout -out ${gene}_aligned_clean.fasta
    cd ${WORKING_DIR}
done < com.txt

# put all inputs into same directory
mkdir -p busco_msa/total

while read -r gene; do
    cd busco_msa/${gene}
    cp ${gene}_aligned_clean.fasta ./../total/.
    cd ${WORKING_DIR}
done < com.txt

# fix names . . . again
while read -r gene; do
    cd busco_msa/total/
    cut -d "_" -f 1 ${gene}_aligned_clean.fasta > temp
    mv temp ${gene}_aligned_clean.fasta
    cd ${WORKING_DIR}
done < com.txt

# get partition file with AMAS
module load mamba/23.1.0-4
source activate AMAS

python3 /home/data/jfierst/veggers/programs/AMAS/amas/AMAS.py concat -c 40 -f fasta -d dna --part-format raxml -i ${WORKING_DIR}/busco_msa/total/*

# iqtree
module load iqtree-2-gcc-8.2.0

/scratch/jfierst/tori/phylogeny/iqtree-1.6.12-Linux/bin/iqtree -s concatenated.out -spp partitions.txt -m MFP+MERGE -bb 1000 -alrt 1000 -nt 40
