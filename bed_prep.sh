#!/bin/bash

#SBATCH --job-name=bed_prep
#SBATCH --account=TG-BIO240239
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --partition=shared
#SBATCH -t 20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vegge003@fiu.edu

input_file=list.txt

#organize/create bed files
while read -r line; do

#paths
	export v5=./RhabditinaPhylogeny_EarlGrey_v5
	export v6=./RhabditinaPhylogeny_EarlGrey_v6
	export edta=./RhabditinaPhylogeny_EDTA/${line}
	export braker2=./RhabditinaPhylogeny_braker2
	export braker3=./RhabditinaPhylogeny_braker3
	export tiberius=./Tiberius_mammalian
	export refseq=./Tiberius_denovo/gtfs

#earlgrey_v5
	if [[ -f ${v5}/${line}.filteredRepeats.gff ]]; then
		cat ${v5}/${line}.filteredRepeats.gff | cut -f 1,4,5 > ${v5}/${line}_temp.bed
		awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${v5}/${line}_temp.bed > ${v5}/${line}_earlgrey_v5_filteredRepeats.bed
		sort -k1,1 -k2,2n ${v5}/${line}_earlgrey_v5_filteredRepeats.bed | uniq > ${v5}/${line}_earlgrey_v5_filteredRepeats_sorted.bed
		rm ${v5}/${line}_temp.bed
		echo "${line} earlgrey_v5 bed done"
	else
		echo "${v5}/${line}.filteredRepeats.gff not found"
	fi

#earlgrey_v6
	if [[ -f ${v6}/${line}.filteredRepeats.gff ]]; then
                cat ${v6}/${line}.filteredRepeats.gff | cut -f 1,4,5 > ${v6}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${v6}/${line}_temp.bed > ${v6}/${line}_earlgrey_v6_filteredRepeats.bed
                sort -k1,1 -k2,2n ${v6}/${line}_earlgrey_v6_filteredRepeats.bed | uniq > ${v6}/${line}_earlgrey_v6_filteredRepeats_sorted.bed
		rm ${v6}/${line}_temp.bed
		echo "${line} earlgrey_v6 bed done"
	else
                echo "${v6}/${line}.filteredRepeats.gff not found"
        fi

#EDTA
        if [[ -f ${edta}/unwrapped_${line}.masked.mod.EDTA.TEanno.gff3 ]]; then
                cat ${edta}/unwrapped_${line}.masked.mod.EDTA.TEanno.gff3 | grep -v "##" | cut -f 1,4,5 > ${edta}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${edta}/${line}_temp.bed > ${edta}/${line}_edta_TEanno.bed
                sort -k1,1 -k2,2n ${edta}/${line}_edta_TEanno.bed | uniq > ${edta}/${line}_edta_TEanno_sorted.bed
		rm ${edta}/${line}_temp.bed
		echo "${line} EDTA bed done"	
	else
                echo "${edta}/unwrapped_${line}.masked.mod.EDTA.TEanno.gff3 not found"
        fi

#braker2
        if [[ -f ${braker2}/${line}_braker2.gtf ]]; then
		#CDS
		cat ${braker2}/${line}_braker2.gtf | awk '{if ($3=="CDS") print}' | cut -f 1,4,5 > ${braker2}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${braker2}/${line}_temp.bed > ${braker2}/${line}_temp_longnames.bed
		awk -F '\t' -v OFS='\t' -v contig_names="./contig_IDs/${line}_contig_IDs.txt" 'BEGIN {while ((getline line < contig_names) > 0) {short[line] = 1}} {for (s in short) {pat = s "_[^ \t]*"; $1 = gensub(pat, s, "g", $1)} print }' ${braker2}/${line}_temp_longnames.bed > ${braker2}/${line}_braker2_CDS.bed
                sort -k1,1 -k2,2n ${braker2}/${line}_braker2_CDS.bed | uniq > ${braker2}/${line}_braker2_CDS_sorted.bed
		
		#exon
                cat ${braker2}/${line}_braker2.gtf | awk '{if ($3=="exon") print}' | cut -f 1,4,5 > ${braker2}/${line}_temp.bed
		awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${braker2}/${line}_temp.bed > ${braker2}/${line}_temp_longnames.bed
		awk -F '\t' -v OFS='\t' -v contig_names="./contig_IDs/${line}_contig_IDs.txt" 'BEGIN {while ((getline line < contig_names) > 0) {short[line] = 1}} {for (s in short) {pat = s "_[^ \t]*"; $1 = gensub(pat, s, "g", $1)} print }' ${braker2}/${line}_temp_longnames.bed > ${braker2}/${line}_braker2_exon.bed
                sort -k1,1 -k2,2n ${braker2}/${line}_braker2_exon.bed | uniq > ${braker2}/${line}_braker2_exon_sorted.bed

		#intron
                cat ${braker2}/${line}_braker2.gtf | awk '{if ($3=="intron") print}' | cut -f 1,4,5 > ${braker2}/${line}_temp.bed
		awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${braker2}/${line}_temp.bed > ${braker2}/${line}_temp_longnames.bed
		awk -F '\t' -v OFS='\t' -v contig_names="./contig_IDs/${line}_contig_IDs.txt" 'BEGIN {while ((getline line < contig_names) > 0) {short[line] = 1}} {for (s in short) {pat = s "_[^ \t]*"; $1 = gensub(pat, s, "g", $1)} print }' ${braker2}/${line}_temp_longnames.bed > ${braker2}/${line}_braker2_intron.bed
                sort -k1,1 -k2,2n ${braker2}/${line}_braker2_intron.bed | uniq > ${braker2}/${line}_braker2_intron_sorted.bed

		rm ${braker2}/${line}_temp.bed
		echo "${line} braker2 bed done"
	else
                echo "${braker2}/${line}_braker2.gtf not found"
        fi

#braker3
        if [[ -f ${braker3}/${line}_braker3.gtf ]]; then
		#CDS
		cat ${braker3}/${line}_braker3.gtf | awk '{if($3=="CDS") print}' | cut -f 1,4,5 > ${braker3}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${braker3}/${line}_temp.bed > ${braker3}/${line}_braker3_CDS.bed
                sort -k1,1 -k2,2n ${braker3}/${line}_braker3_CDS.bed | uniq > ${braker3}/${line}_braker3_CDS_sorted.bed
		#exon
                cat ${braker3}/${line}_braker3.gtf | awk '{if($3=="CDS") print}' | cut -f 1,4,5 > ${braker3}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${braker3}/${line}_temp.bed > ${braker3}/${line}_braker3_exon.bed
                sort -k1,1 -k2,2n ${braker3}/${line}_braker3_exon.bed | uniq > ${braker3}/${line}_braker3_exon_sorted.bed
		#intron
                cat ${braker3}/${line}_braker3.gtf | awk '{if($3=="intron") print}' | cut -f 1,4,5 > ${braker3}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${braker3}/${line}_temp.bed > ${braker3}/${line}_braker3_intron.bed
                sort -k1,1 -k2,2n ${braker3}/${line}_braker3_intron.bed | uniq > ${braker3}/${line}_braker3_intron_sorted.bed

		rm ${braker3}/${line}_temp.bed
		echo "${line} braker3 bed done"
	else
                echo "${braker3}/${line}_braker3.gtf not found"
        fi

#Tiberius
        if [[ -f ${tiberius}/${line}_tiberius_mammalian.gtf ]]; then
		cat ${tiberius}/${line}_tiberius_mammalian.gtf | awk '{if ($3=="CDS") print}' > ${tiberius}/${line}_tiberius_mammalian_CDS.gtf
		cat ${tiberius}/${line}_tiberius_mammalian.gtf | awk '{if ($3=="exon") print}' > ${tiberius}/${line}_tiberius_mammalian_exon.gtf
		cat ${tiberius}/${line}_tiberius_mammalian.gtf | awk '{if ($3=="intron") print}' > ${tiberius}/${line}_tiberius_mammalian_intron.gtf
		#CDS
		cut -f 1,4,5 ${tiberius}/${line}_tiberius_mammalian_CDS.gtf > ${tiberius}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${tiberius}/${line}_temp.bed > ${tiberius}/${line}_tiberius_mammalian_CDS.bed
                sort -k1,1 -k2,2n ${tiberius}/${line}_tiberius_mammalian_CDS.bed | uniq > ${tiberius}/${line}_tiberius_mammalian_CDS_sorted.bed
		#exon
		cut -f 1,4,5 ${tiberius}/${line}_tiberius_mammalian_exon.gtf > ${tiberius}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${tiberius}/${line}_temp.bed > ${tiberius}/${line}_tiberius_mammalian_exon.bed
                sort -k1,1 -k2,2n ${tiberius}/${line}_tiberius_mammalian_exon.bed | uniq > ${tiberius}/${line}_tiberius_mammalian_exon_sorted.bed
		#intron
		cut -f 1,4,5 ${tiberius}/${line}_tiberius_mammalian_intron.gtf > ${tiberius}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${tiberius}/${line}_temp.bed > ${tiberius}/${line}_tiberius_mammalian_intron.bed
                sort -k1,1 -k2,2n ${tiberius}/${line}_tiberius_mammalian_intron.bed | uniq > ${tiberius}/${line}_tiberius_mammalian_intron_sorted.bed
		
		rm ${tiberius}/${line}_temp.bed
		echo "${line} tiberius bed done"
	else
                echo "${tiberius}/${line}_tiberius_mammalian.gtf not found"
        fi

#RefSeq
        if [[ -f ${refseq}/${line}_refseq.gtf ]]; then
                cat ${refseq}/${line}_refseq.gtf | grep -v "#" | grep "gbkey \"CDS\"" | cut -f 1,4,5 > ${refseq}/${line}_temp.bed
                awk '{if($2>$3) {tmp=$2; $2=$3; $3=tmp} print}' ${refseq}/${line}_temp.bed > ${refseq}/${line}_refseq_CDS.bed
                sort -k1,1 -k2,2n ${refseq}/${line}_refseq_CDS.bed | uniq > ${refseq}/${line}_refseq_CDS_sorted.bed
		rm ${refseq}/${line}_temp.bed
		echo "${line} refseq bed done"
	else
                echo "${refseq}/${line}_refseq.gtf not found"
        fi

done < ${input_file}

