#!/bin/bash
#SBATCH --job-name=EDTA
#SBATCH --output=./logs/EDTA.%j.%N.out
#SBATCH --array=1-67%34
#SBATCH --mem=90G
#SBATCH --account=TG-BIO240239
#SBATCH --cpus-per-task=10
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=shared
#SBATCH -t 40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vegge003@fiu.edu

module load cpu
module load slurm
source activate EDTA

export wd=/expanse/lustre/projects/fiu112/veggers

SPECIES=$(sed "${SLURM_ARRAY_TASK_ID}q;d" list.txt)

echo "$SPECIES"

#unwrap the fasta
awk '/^>/ {print (NR==1?"":"\n") $0; next} {printf "%s", $0} END {print ""}' ./RhabditinaPhylogeny_repeatmasker/${SPECIES}/*.masked > unwrapped_${SPECIES}.masked

#remove spaces in header
sed -i 's/ /_/g' unwrapped_${SPECIES}.masked

/home/veggers/miniconda3/envs/EDTA/share/EDTA/EDTA.pl --genome unwrapped_${SPECIES}.masked --rmlib ./RhabditinaPhylogeny_repeatmodeler/${SPECIES}/${SPECIES}.repeats --sensitive 1 --force 1 --anno 1 --threads 10

mkdir -p ./RhabditinaPhylogeny_EDTA/${SPECIES}
mv unwrapped_${SPECIES}.masked.mod* ./RhabditinaPhylogeny_EDTA/${SPECIES}/. 
cp -r unwrapped_${SPECIES}.masked.mod* ./RhabditinaPhylogeny_EDTA/${SPECIES}/.
rm -r unwrapped_${SPECIES}.masked.mod.EDTA.raw
