#!/bin/bash

# Usage: ./fill_proxy.sh <exposure> <outcome> 
# Example: ./fill_proxy.sh LST AIS 
exposure=$1
outcome=$2

echo "BEGIN:开始获取 $exposure 与 $outcome 的SNP proxy"
###step1: 从LDlink获取缺失snp的proxy
while read snp; do
    curl -k -X GET "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=${snp}&pop=EUR&r2_d=r2&window=5000&genome_build=grch37&token=37c638ef14d5" > /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/SNP/${snp}.tmp
    awk -v OFS="\t" '$7 > 0.8&&NR>2{print "'"${snp}"'",$1,$2,$4,$5,$7,$8}' /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/SNP/${snp}.tmp > /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/SNP/${snp}.ld
done < /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/${outcome}_${exposure}.snpmiss

echo "step1: 已从LDlink完成proxy获取"

###step2: 为step3准备可匹配的proxy snp数据
# 将所有代理SNP的数据文件合并为一个文件
cat /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/SNP/*.ld > /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/${outcome}_${exposure}.tmp
# 创建获取到proxy的snp表格：去掉重复的snp,只保留获取到proxy的snp
awk '!seen[$1]++ {print $1}' /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/${outcome}_${exposure}.tmp > /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/${outcome}_${exposure}.psnpmiss

echo "step2: 已创建能够获取到proxy的snp信息文件，后缀为psnpmiss。"

##step3: 从outcome数据集中匹配候选proxy的信息（即筛选出outcome数据集中也包含的proxy snp）
awk -v OFS="\t" 'NR==FNR{a[$2]=$0;next}NR>FNR && ($3 in a){print a[$3],$0}' /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/${outcome}_${exposure}.tmp /media/yidie/Linuxdata/GWAS/EUR/${outcome}.tsv > /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}.proxy

echo "step3: 已筛选出 $outcome .tsv数据中也包含的snp proxy"

##step4: 筛选出LD关联度最高，且距离最近的Proxy snp

awk -v OFS="\t" '{$5 = ($5 < 0) ? -$5 : $5; print}' /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}.proxy | sort -k1,1 -k6,6nr -k5,5n > /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_sorted.proxy

awk -v OFS="\t" '!seen[$1]++ {print $1,$2,$5,$6,$7}' /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_sorted.proxy > /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_final.proxy

lines=$(wc -l < /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_final.proxy)
echo "step4: 确定最终的proxy,存于后缀为final.proxy的文件中，最终用于替代缺失SNP的proxy数量为: $line 。"

##step5: 从outcome数据集中匹配Proxy snp的数据
cat /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_final.proxy | awk -F '\t' 'BEGIN{OFS=FS}{split($5,a,"[=,=]"); print "TRUE",$2,$1,a[1],a[2],a[3],a[4]}' > /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_finaltrans.proxy

sed -i '1i proxy.outcome\ttarget_snp.outcome\tproxy_snp.outcome\tpa1\ta1\tpa2\ta2' /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_finaltrans.proxy

cat /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_final.proxy | awk -v OFS="\t" '{print $2}' | xargs -I{} awk -v OFS="\t" '$3=="{}"{print $0}' /media/yidie/Linuxdata/GWAS/EUR/${outcome}.tsv > /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_info.proxy

awk -v OFS="\t" 'NR==FNR{a[$3"_"$4"_"$5]=$0;next}{if($2"_"$5"_"$7 in a){print a[$2"_"$5"_"$7],$4,$6,$1,$2,$3,$4,$6,$5,$7}else if($2"_"$7"_"$5 in a){print a[$2"_"$7"_"$5],$6,$4,$1,$2,$3,$6,$4,$7,$5}}' /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_info.proxy /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_finaltrans.proxy > /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_allinfo.proxy

cat /media/yidie/Workspace/analysis0607_stroke/data/proxy/${outcome}_${exposure}_allinfo.proxy | awk -v OFS="\t" '{print $18,$15,$16,$6,$7,$8,$9,$13,$14,"TRUE","reported","'"${outcome}"'",$17,$18,$19,$20,$21,$22,$23}' > /media/yidie/Workspace/analysis0607_stroke/data/proxy/${exposure}_${outcome}_match.proxy

cat /media/yidie/Workspace/analysis0607_stroke/data/matched/${exposure}_${outcome}.matched.tsv /media/yidie/Workspace/analysis0607_stroke/data/proxy/${exposure}_${outcome}_match.proxy > /media/yidie/Workspace/analysis0607_stroke/data/matched/${exposure}_${outcome}.clumped2.tsv

echo "step5: $exposure 与 $outcome 的SNP proxy获取完毕"

rm -rf /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/*.tmp
rm -rf /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/SNP/*.ld
rm -rf /media/yidie/Workspace/analysis0607_stroke/data/missingSNPs/SNP/*.tmp
echo "END:----------------------------------分割线------------------------------------------"
