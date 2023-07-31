#!/bin/bash

# Usage: ./match_snpmiss.sh <exposure> <outcome>
# Example: ./match_snpmiss.sh LST AIS
exposure=$1
outcome=$2
expodir=$3

# 使用awk命令从暴露数据集和结局数据集中提取数据，其中第三列是SNP
awk 'FNR==NR{a[$3];next} $3 in a' /media/yidie/Linuxdata/clumped/${exposure}.clumped.tsv /media/yidie/Linuxdata/GWAS/EUR/${outcome}.tsv > /media/yidie/Workspace/analysis0607_stroke/data/matched/${exposure}_${outcome}.clumped.tsv

# 使用awk命令找到结局数据集中缺失的SNP
awk 'FNR==NR{a[$3];next} !($3 in a){print $3}' /media/yidie/Workspace/analysis0607_stroke/data/matched/${exposure}_${outcome}.clumped.tsv /media/yidie/Linuxdata/clumped/${exposure}.clumped.tsv > /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/${outcome}_${exposure}.snpmiss
