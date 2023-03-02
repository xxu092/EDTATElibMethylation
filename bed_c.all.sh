#!/usr/bin/bash -l
#SBATCH -p batch -n 1 --mem 64gb -N 1 --out logs/bed_EDTA%A.log

##this is to overlap all the methylation calls from megalodon to L1
module load bedtools

DNAMET=~/bigdata/M.cicadina/methylation/new.con.5mc.bed
#DNAMET=rs.5mc.bed.gz

OUT=methylsummary/TElib_con_5mc_c.tsv
#OUT=TE_rs_5mc_c.tsv
mkdir -p methylsummary

IN=M.cicadina.TElib.clean.gff3
bedtools intersect -a $IN -b $DNAMET -c > $OUT

