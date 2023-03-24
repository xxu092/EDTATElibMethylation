# EDTATElibMethylation
This repository is to analyze transposable element methylation

## in gff3 file, structrually intact TEs are represented as several records, just keep one record with whole length for these elements, filter out the ones with parent 

`grep -v "Parent" M.cicadina.TElib.gff3 > M.cicadina.TElib.clean.gff3`

## intersect TE library with all CpG sites from megalodon output
script: `bed_c.all.sh`

## filter methylation bed file, >20% confidence in methylation counts as positively methylated
`awk '$11 >= 20 {print}' new.con.5mc.bed > new.con.5mc.20.bed`

## intersect TE library with only positively methylated CpG sites 
script: `bed_c20.sh`

## head to R make the table for methylation and length for all TEs.
script: `TEmethylsummary.R`

## analyze LINE families length and methylation
script: `LINE/LINEsummary.R`
