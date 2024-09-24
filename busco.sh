#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%busco.log
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH --mail-type=ALL
#SBATCH -n 40

module load busco/5.4.7

export AUGUSTUS_CONFIG_PATH="/home/data/jfierst/veggers/programs/Augustus"
max_jobs=10
job_count=0

while read -r line; do
  busco -f -m genome -i ./RhabditinaPhylogeny_NCBI/"${line}"/*.fna -o busco_"${line}" --offline --lineage_dataset /home/data/jfierst/veggers/nematoda_odb10 &

    job_count=$((job_count + 1))

        if [ "$job_count" -ge "$max_jobs" ]; then
                wait
                job_count=0
        fi

done < busco_list.txt


#some of my files were gzipped and so they didn't work the first time. consequences of downloading genomes one by one instead of computationally. 
#I made the following loop to run through the ones that needed to be fixed after I unzipped the files

#while read -r line; do
# busco -f -c 4 -m genome -i ./RhabditinaPhylogeny_NCBI/${line}/*.fna -o busco_${line} --offline --lineage_dataset /home/data/jfierst/veggers/nematoda_odb10
#done < fix_list.txt
