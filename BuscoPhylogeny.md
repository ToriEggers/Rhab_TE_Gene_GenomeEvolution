## Comparative Genomics and Phylogenomics of Rhabdiditae Nematodes ###

***Making a Phylogeny with BUSCO shared single copy orthologs***

Code was derived from [Manni et al. 2021](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.323). For more information please read their paper.

<details>
    
<summary><b>Run BUSCO on all genomes in parallel</b></summary>

```
vi busco.sh
```

```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_%busco.log
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
```
sbatch busco.sh
```
Takes less than a day to get through all 70 genomes. 

</details> 
