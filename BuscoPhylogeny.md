## Comparative Genomics and Phylogenomics of Rhabdiditae Nematodes ###

***Making a Phylogeny with BUSCO shared single copy orthologs***

Code was derived from [Manni et al. 2021](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.323). For more information please read their paper.

To download the data and set up the directories see: [Data.md]()

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
    <summary><b>Concatenate BUSCO summaries for all genomes</b></summary>
    
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

</details>

<details>
    <summary><b>Concatenate all shared single copy orthologous sequences</b></summary>
</details>
