#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%busco_parse.log

WORKING_DIR=/home/data/jfierst/veggers/RhabditinaPhylogeny/

cd ${WORKING_DIR}

echo -e "ID\tbusco\tsingle_copy_count" > busco_summary.txt

while read -r line; do
    cd RhabditinaPhylogeny_Buscos/${line}/busco_${line}/run_nematoda_odb10/busco_sequences/single_copy_busco_sequences/
    ls *.fna | sort > ./../../../${line}_single_copy_list.txt
    cd ${WORKING_DIR}

    single_copy_count=$(wc -l RhabditinaPhylogeny_Buscos/${line}/busco_${line}/${line}_single_copy_list.txt | awk '{print $1}')
    ID=${line}
    busco=$(grep "C:" RhabditinaPhylogeny_Buscos/${line}/busco_${line}/short_summary*.txt)
    echo -e "${ID}\t${busco}\t${single_copy_count}" >> busco_summary.txt

done < busco_list.txt
