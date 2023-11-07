#!/bin/bash

##Shell scripts
#Usage: ./clump.sh g1000_eur EUR LST 
./software/plink/plink --bfile ./work_dir/data/g1000/${1} \
--clump /media/yidie/Linuxdata/GWAS/${2}/${3}.tsv \
--clump-p1 5e-8 --clump-p2 1 --clump-r2 0.001 --clump-kb 10000 \
--out tmp

awk 'BEGIN{FS="\t"; OFS="\t"}FNR==NR{a[$1]=$0}FNR<NR&&($3 in a){print $0}' \
<(cat tmp.clumped | awk -v OFS="\t" '{print $3}') \
<(cat ./GWAS/${2}/${3}.tsv) \
> ./analysis0731_stroke/clumped_plink/${3}.clumped.tsv

rm -rf tmp*
