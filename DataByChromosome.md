## Comparative Genomics and Phylogenomics of Rhabdiditae Nematodes ###

***Data by Chromosome***

code to separate the data by chr (I am not in a coding state of mind. don't judge)

I took the time to go through each species and each chromosome. 
The difficulty is that naming isn't consistent; a chromosome could be named chrI, I, chr1, chromosome 1, chromosome:1, etc.
Tedious task but it'll speed up the code latter on. 


```
cd /home/data/jfierst/veggers/RhabditinaPhylogeny/RhabditinaPhylogeny_repeatmasker/${species}
grep -c ">" ${species}.masked
grep ">" ${species}.masked | head (or just grep without the head depending on how many contigs)
module load mamba/23.1.0-4
source activate seqkit
seqkit grep -n -p "NC_013486.2" ${species}.masked.prettyNames > ${species}_chr2_masked.fasta
	#repeat for all chromosomes
	#any not placed scaffolds larger than 1Mb were labelled chrU1, chrU2, chrU3 . . .
	#note: many unplaced scaffolds in the majority of assemblies, however most are less than 100Kb. Do the best with what we've got
mkdir by_chr
mv ${species}_chr* by_chr

#the output will be fastas containing 1 chr each in the dir /home/data/jfierst/veggers/RhabditinaPhylogeny/RhabditinaPhylogeny_repeatmasker/${species}/by_chr
#the name of the file will have chr# (1/2/3/4/5/X)
#the header line will be cleaned fasta header (ie. > NC_013486.2) which can be matched to annotation files
```

```
#!/bin/bash

while read -r species; do
        for file in ./RhabditinaPhylogeny_repeatmasker/${species}/by_chr/*; do
                ID=$(head -1 ${file} | sed 's/>//g')
                length=$(grep "${ID}" ./RhabditinaPhylogeny_repeatmasker/${species}/${species}_contig_lengths.txt | cut -f 3 -d " ")
                chr=$(basename ${file} | cut -f 2 -d "_")
                echo -e "${species}\t${chr}\t${ID}\t${length}" >> species_chr_ID_length.txt
        done
done < chr_lvl_species.txt
```

The output is species_chr_ID_length.txt, which looks like:

```
AF16    chr1    NC_013489.2     15455979
AF16    chr2    NC_013486.2     16627154
AF16    chr3    NC_013490.2     14578851
AF16    chr4    NC_013487.2     17485439
AF16    chr5    NC_013488.2     19495157
AF16    chrX    NC_013491.2     21540570
APS25   chr1    OZ181880.1      16121759
APS25   chr2    OZ181881.1      12655968
APS25   chr3    OZ181882.1      12260034
```

Number of chromosome level assemblies:
```
wc -l chr_lvl_species.txt
# 72
```

IDs of those I consider chromosome level:
```
AF16
APS25
APS4
APS7
Aroian
BAKE
BRC20456
BRC20483
CB4856
CEW1
CFB2252
CP168
DF5070
DF5081
DF5083
DF5112
DF5120
DF5173
EG5942
ISE
JU1182
JU1286
JU1382
JU1421
JU1667
JU1771
JU1809
JU1904
JU1917
JU1968
JU2083
JU2190
JU2585
JU2788
JU2809
JU3283
JU3284
JU75
JU800
KANDY
MIMR
N2
NIC203
NIC394
NIC534
NIC564
NIC58
NKZ35
OM
OST
PF1305
PS1010
PS1017
PS2068
PS312
PX356
PX439
PX506
PX534
QG2077
QG555
RS0144
RS5133
RS5460
RS5522B
SB355
SX3368
TRI
TT01
TWN1964
TWN1984
VIVI
```
