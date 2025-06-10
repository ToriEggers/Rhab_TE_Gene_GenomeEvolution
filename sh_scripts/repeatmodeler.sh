#!/bin/bash
#SBATCH --job-name=repeatmodeler
#SBATCH --output=./logs/repeatmodeler.%j.%N.out
#SBATCH --array=1-67%34
#SBATCH --mem=90G
#SBATCH --account=TG-BIO240239
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=3
#SBATCH --nodes=1
#SBATCH --partition=shared
#SBATCH -t 40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vegge003@fiu.edu

module load cpu
module load slurm
source actiavte repeatmodeler

export wd=/expanse/lustre/projects/fiu112/veggers


SPECIES=$(sed "${SLURM_ARRAY_TASK_ID}q;d" list.txt)

echo "$SPECIES"

#organize
cd RhabditinaPhylogeny_repeatmodeler
mkdir -p ${SPECIES}
cd ${SPECIES}
	
#Build the database
BuildDatabase -name ${SPECIES} ./../../RhabditinaPhylogeny_NCBI/${SPECIES}/*.fna

#Run RepeatModeler for de novo repeat identification and characterization. Takes long time.
RepeatModeler -pa 3 -database ${SPECIES}

#Combine the files to create a library of de novo and known repeats
cat RM*/consensi.fa.classified ./../../Rhab.repeats > ${SPECIES}.repeats

cd ${wd}


