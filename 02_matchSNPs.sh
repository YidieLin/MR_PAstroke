#!/bin/bash

# Usage: ./02_matchV1.sh <exposure> <outcome> 
# Example: ./02_matchV1.sh LST AIS 
exposure=$1
outcome=$2

echo "BEGIN:matching $exposure & $outcome SNPs"
# extracting SNPs information
awk 'FNR==NR{a[$3];next} $3 in a' /media/yidie/Workspace/analysis0731_stroke/clumped_plink/${exposure}.clumped.tsv /media/yidie/Linuxdata/GWAS/EUR/${outcome}.tsv > /media/yidie/Workspace/analysis0731_stroke/match_V1/${exposure}_${outcome}.matchV1.tsv

line1=$(awk 'END {print NR}' /media/yidie/Workspace/analysis0731_stroke/match_V1/${exposure}_${outcome}.matchV1.tsv)-1
echo "step1 successfully matched: $line1 。"

# marked missing SNPs
awk 'FNR==NR{a[$3];next} !($3 in a){print $3}' /media/yidie/Workspace/analysis0731_stroke/match_V1/${exposure}_${outcome}.matchV1.tsv /media/yidie/Workspace/analysis0731_stroke/clumped_plink/${exposure}.clumped.tsv > /media/yidie/Workspace/analysis0731_stroke/missingSNPs/${outcome}_${exposure}.snpmiss

line2=$(awk 'END {print NR}' /media/yidie/Workspace/analysis0731_stroke/missingSNPs/${outcome}_${exposure}.snpmiss)
echo "step2 missing SNPs: $line2 。"
echo "END:---------------------------------------------------------------------------"
