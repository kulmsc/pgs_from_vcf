prefix=$1


zcat ${prefix}.lifted.vcf.gz | fgrep -v '#' | cut -f3,6 > qual
zcat ${prefix}.lifted.vcf.gz | fgrep -v '#' | cut -f8 | cut -f3 -d';' | cut -f2 -d'=' > dp
zcat ${prefix}.lifted.vcf.gz | fgrep -v '#' | cut -f8 | cut -f4 -d';' | cut -f2 -d'=' > mq
paste qual dp mq > qc_vcf

mean_dp=`awk '{ total += $2 } END { print total/NR }' qc_vcf`
dp_hi=`echo $mean_dp*3 | bc`

cat qc_vcf | fgrep rs | awk -v var="$dp_hi" '$3 < 5 || $2 < 15 || $3 > var || $4 < 30   {print $0}'

cat qc_vcf | fgrep rs | awk -v var="$dp_hi" '$3 < 5 || $2 < 15 || $3 > var || $4 < 30   {print $0}' | cut -f1 > ${prefix}.bad_qc_rsid

rm qual dp mq qc_vcf
