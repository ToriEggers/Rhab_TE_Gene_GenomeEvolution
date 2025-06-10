#!/bin/bash

#SBATCH --job-name=bed_intersect_annotation_quality_check
#SBATCH --account=TG-BIO240239
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=shared
#SBATCH -t 10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vegge003@fiu.edu

#run on SDSC expanse. remember to change the species and path variables
#runs in seconds

module load cpu/0.17.3b
module load gcc/10.2.0/npcyll4
module load bedtools2/2.30.0/5anbfkq

mkdir -p RhabditinaPhylogeny_annotation_comparison 
output_dir=./RhabditinaPhylogeny_annotation_comparison

#for N2 
#intersect v5 with v6 with EDTA

species=N2

edta=./RhabditinaPhylogeny_EDTA/${species}/${species}_edta_TEanno_sorted.bed
v5=./RhabditinaPhylogeny_EarlGrey_v5/${species}_earlgrey_v5_filteredRepeats_sorted.bed
v6=./RhabditinaPhylogeny_EarlGrey_v6/${species}_earlgrey_v6_filteredRepeats_sorted.bed

#overlapping bases for each file
bedtools intersect -wao -a ${v6} \
	-b ${v5} ${edta} \
	-names earlgrey_v5 edta \
	-sorted > ${output_dir}/${species}_repeat_annotations_bed_wao.txt

#number of overlaps in each b file for each feature in file a (detecting fragmentation?)
bedtools intersect -C -a ${v6} \
        -b ${v5} ${edta} \
        -names earlgrey_v5 edta \
        -sorted > ${output_dir}/${species}_repeat_annotations_bed_C.txt

#features that don't overlap (detecting FP or FN?)
bedtools intersect -v -a ${v5} \
	-b ${v6} -sorted > ${output_dir}/${species}_v5_not_in_v6.txt
bedtools intersect -v -a ${edta} \
	-b ${v6} -sorted > ${output_dir}/${species}_edta_not_in_v6.txt

#for NKZ352
#intersect refseq with braker2 with braker3 with tiberius

species=NKZ352

refseq=./Tiberius_denovo/gtfs/${species}_refseq_CDS_sorted.bed
braker2=./RhabditinaPhylogeny_braker2/${species}_braker2_CDS_sorted.bed
braker3=./RhabditinaPhylogeny_braker3/${species}_braker3_CDS_sorted.bed
tiberius=./Tiberius_mammalian/${species}_tiberius_mammalian_CDS_sorted.bed

#overlapping bases for each file
bedtools intersect -wao -a ${refseq} \
	-b ${braker2} ${braker3} ${tiberius} \
	-names braker2 braker3 tiberius \
	-sorted > ${output_dir}/${species}_gene_annotations_bed_wao.txt

#detecting fragmentation?
bedtools intersect -C -a ${refseq} \
	-b ${braker2} ${braker3} ${tiberius} \
	-names braker2 braker3 tiberius \
        -sorted > ${output_dir}/${species}_gene_annotations_bed_C.txt

#detecting FP and FN?
bedtools intersect -v -a ${braker2} \
        -b ${refseq} -sorted > ${output_dir}/${species}_braker2_not_in_refseq.txt
bedtools intersect -v -a ${braker3} \
        -b ${refseq} -sorted > ${output_dir}/${species}_braker3_not_in_refseq.txt
bedtools intersect -v -a ${tiberius} \
        -b ${refseq} -sorted > ${output_dir}/${species}_tiberius_not_in_refseq.txt

