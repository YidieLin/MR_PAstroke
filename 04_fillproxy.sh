#!/bin/bash

# Usage: ./04_fillproxy.sh <exposure> <outcome> 
# Example: ./04_fillproxy.sh LST AIS 
exposure=$1
outcome=$2

echo "BEGIN: achieve $exposure vs $outcome SNP proxy"
###step1: finding proxy from LDlink
while read snp; do
    curl -k -X GET "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=${snp}&pop=EUR&r2_d=r2&window=500000&genome_build=grch37&token=37c638ef14d5" > /media/yidie/Workspace/analysis0731_stroke/missingSNPs/SNP/${snp}.tmp
    awk -v OFS="\t" '$7 > 0.8&&NR>2{print "'"${snp}"'",$1,$2,$4,$5,$7,$8}' /media/yidie/Workspace/analysis0731_stroke/missingSNPs/SNP/${snp}.tmp > /media/yidie/Workspace/analysis0731_stroke/missingSNPs/SNP/${snp}.ld
done < /media/yidie/Workspace/analysis0731_stroke/missingSNPs/${outcome}_${exposure}.snpmiss

echo "step1: proxy achieved"

###step2: finding proxy snps information

cat /media/yidie/Workspace/analysis0731_stroke/missingSNPs/SNP/*.ld > /media/yidie/Workspace/analysis0731_stroke/missingSNPs/${outcome}_${exposure}.tmp
# Keep available proxy
awk '!seen[$1]++ {print $1}' /media/yidie/Workspace/analysis0731_stroke/missingSNPs/${outcome}_${exposure}.tmp > /media/yidie/Workspace/analysis0731_stroke/missingSNPs/${outcome}_${exposure}.psnpmiss

echo "step2: proxy snp information .psnpmiss."

##step3: matching proxy snp
awk -v OFS="\t" 'NR==FNR{a[$2]=$0;next}NR>FNR && ($3 in a){print a[$3],$0}' /media/yidie/Workspace/analysis0731_stroke/missingSNPs/${outcome}_${exposure}.tmp /media/yidie/Linuxdata/GWAS/EUR/${outcome}.tsv > /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}.proxy

echo "step3:  $outcome snp proxy"

##step4: selecting the final proxy

awk -v OFS="\t" '{$5 = ($5 < 0) ? -$5 : $5; print}' /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}.proxy | sort -k1,1 -k6,6nr -k5,5n > /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_sorted.proxy

awk -v OFS="\t" '!seen[$1]++ {print $1,$2,$5,$6,$7}' /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_sorted.proxy > /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_final.proxy

lines=$(wc -l < /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_final.proxy)
echo "step4: save in final.proxy，SNP proxy No. of $lines 。"

##step5: tidy Proxy snp data
cat /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_final.proxy | awk -F '\t' 'BEGIN{OFS=FS}{split($5,a,"[=,=]"); print "TRUE",$2,$1,a[1],a[2],a[3],a[4]}' > /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_finaltrans.proxy

sed -i '1i proxy.outcome\tproxy_snp.outcome\ttarget_snp.outcome\ta1\tpa1\ta2\tpa2' /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_finaltrans.proxy

cat /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_final.proxy | awk -v OFS="\t" '{print $2}' | xargs -I{} awk -v OFS="\t" '$3=="{}"{print $0}' /media/yidie/Linuxdata/GWAS/EUR/${outcome}.tsv > /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_info.proxy

awk -v OFS="\t" 'NR==FNR{a[$3"_"$4"_"$5]=$0;next}{if($2"_"$5"_"$7 in a){print a[$2"_"$5"_"$7],$4,$6,$1,$2,$3,$4,$6,$5,$7}else if($2"_"$7"_"$5 in a){print a[$2"_"$7"_"$5],$6,$4,$1,$2,$3,$6,$4,$7,$5}}' /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_info.proxy /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_finaltrans.proxy > /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_allinfo.proxy

cat /media/yidie/Workspace/analysis0731_stroke/proxy/${outcome}_${exposure}_allinfo.proxy | awk -v OFS="\t" '{print $16,$12,$13,$6,$7,$8,$9,$10,$11,"TRUE","reported","'"${outcome}"'",$14,$16,$15,$17,$18,$19,$20}' > /media/yidie/Workspace/analysis0731_stroke/proxy/${exposure}_${outcome}_match.proxy

cat /media/yidie/Workspace/analysis0731_stroke/match_V2/${exposure}_${outcome}.matchV2.tsv /media/yidie/Workspace/analysis0731_stroke/proxy/${exposure}_${outcome}_match.proxy > /media/yidie/Workspace/analysis0731_stroke/match_V3/${exposure}_${outcome}.matchV3.tsv

echo "step5 $exposure vs $outcome SNP proxy: DONE"

rm -rf /media/yidie/Workspace/analysis0731_stroke/missingSNPs/*.tmp
rm -rf /media/yidie/Workspace/analysis0731_stroke/missingSNPs/SNP/*.ld
rm -rf /media/yidie/Workspace/analysis0731_stroke/missingSNPs/SNP/*.tmp
echo "END:--------------------------------------------------------------------------"
