#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%busco_parse_2.log

WORKING_DIR=/home/data/jfierst/veggers/RhabditinaPhylogeny

cd ${WORKING_DIR}

#prep loop
startFile=$(head -1 busco_parse.txt)
cp ${WORKING_DIR}/RhabditinaPhylogeny_Buscos/${startFile}/busco_${startFile}/${startFile}_single_copy_list.txt com.txt

#loop through species to find buscos found in all species (com.txt)
while read -r file; do
    comm -12 com.txt ${WORKING_DIR}/RhabditinaPhylogeny_Buscos/${file}/busco_${file}/${file}_single_copy_list.txt > temp.txt
    mv temp.txt com.txt
done < busco_parse.txt

#fix names and create file for each gene
mkdir -p busco_msa

while read -r gene; do
    mkdir -p ./busco_msa/${gene}
    cd busco_msa/${gene}

    while read -r species; do
        cp ${WORKING_DIR}/RhabditinaPhylogeny_Buscos/${species}/busco_${species}/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/${gene} .
        sed -e "s/^>/>${species}_/" ${gene} | sed 's/-/_/g' | sed 's/\./_/g' | sed 's/\:/_/g' > ${species}_${gene}
        rm ${gene}
    done < ./../../busco_parse.txt

    cat * > ${gene}.fasta
    rm *.fna
    cd ${WORKING_DIR}

done < com.txt

#msa for each gene with 12 jobs at a time
module load mafft-7.221-gcc-8.2.0-y6cgezm

        max_jobs=12
        job_count=0

while read -r gene; do
    linsi --thread 4 --maxiterate 1000 ${gene}.fasta > ${gene}_aligned.fasta &

        job_count=$((job_count + 1))
            if [ "$job_count" -ge "$max_jobs" ]; then
                wait
                job_count=0
            fi

done < com.txt
