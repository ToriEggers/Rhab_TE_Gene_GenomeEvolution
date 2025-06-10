#!/bin/bash
#SBATCH --job-name=quick_check
#SBATCH --account=TG-BIO240239
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=shared
#SBATCH -t 10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vegge003@fiu.edu

module load cpu
module load slurm

while read -r species; do
	mv ./RhabditinaPhylogeny_NCBI/${species}/*.fna.cat.gz ./RhabditinaPhylogeny_repeatmasker/${species}/.
	mv ./RhabditinaPhylogeny_NCBI/${species}/*.fna.masked ./RhabditinaPhylogeny_repeatmasker/${species}/.
	mv ./RhabditinaPhylogeny_NCBI/${species}/*.fna.out ./RhabditinaPhylogeny_repeatmasker/${species}/.
	mv ./RhabditinaPhylogeny_NCBI/${species}/*.fna.tbl ./RhabditinaPhylogeny_repeatmasker/${species}/. 
done < list.txt
