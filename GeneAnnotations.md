## Comparative Genomics and Phylogenomics of Rhabdiditae Nematodes ###

All preparation for RNA evidence for BRAKER3 was done in scratch directory.
This includes downloading RNA data from NCBI (aimed for total RNA seq illumina data), alignment with STAR, and indexing/sorting with samtools. Once the bam files were created successfully, they were moved to the directory containing the braker.sif file.

***Download RNA data***

[RNA_accessions.txt]()

download.sh 
```
#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=./logs/RNAseq_download_%N_%j.log

module load proxy #needed to connect to internet
module load sratoolkit-3.0.0

INPUT_FILE=refseq_RNA_accessions.txt #input file: ID ACCESSION (2 columns seperated by a tab)

while read -r line; do
        species=$(echo $line | cut -f 1) #set variables
        accession=$(echo $line | cut -f 2)

        prefetch $accession #sratoolkit commands
        fasterq-dump $accession --outdir $species

        rm -r $accession #delete prefetch created directory
done < $INPUT_FILE
```
