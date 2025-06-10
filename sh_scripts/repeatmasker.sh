#!/bin/bash
#SBATCH --job-name=repeatmasker
#SBATCH --output=./logs/repeatmasker.%j.%N.out
#SBATCH --array=1-67%34
#SBATCH --mem=90G
#SBATCH --account=TG-BIO240239
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --partition=shared
#SBATCH -t 40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vegge003@fiu.edu

module load cpu
module load slurm
source activate repeatmodeler

export wd=/expanse/lustre/projects/fiu112/veggers

SPECIES=$(sed "${SLURM_ARRAY_TASK_ID}q;d" list.txt)

echo "$SPECIES"

#organize
cd RhabditinaPhylogeny_repeatmasker
mkdir -p ${SPECIES}
cd ${SPECIES}
	
GENOME=${wd}/RhabditinaPhylogeny_NCBI/${SPECIES}/*.fna
REPEATS=${wd}/RhabditinaPhylogeny_repeatmodeler/${SPECIES}/${SPECIES}.repeats

#Mask the genome of known repeats
RepeatMasker -lib ${REPEATS} -pa 8 -xsmall -nolow ${GENOME} 

cd ${wd}


