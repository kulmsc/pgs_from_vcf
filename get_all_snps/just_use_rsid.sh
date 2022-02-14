chr=$1
prefix=$2

rm all_rsids 
rm all_ids
ls ../mod_sets/ | fgrep .${chr}. | while read line;do
  cat ../mod_sets/$line | tail -n+2 | cut -f3 >> all_rsids
  cat ../mod_sets/$line | tail -n+2 | cut -f3-5  >> all_ids
done

ls ../pgs_sets | fgrep .${chr}. | while read line;do
  cat ../pgs_sets/$line | tail -n+2 | cut -f2 >> all_rsids
done

cat all_rsids | sort | uniq > temp; mv temp all_rsids
cat all_ids | sort | uniq > temp; mv temp all_ids


plink --vcf ../prep_vcf/${prefix}.lifted.vcf.gz --chr $chr --allow-extra-chr --extract all_rsids --make-bed --out split_files/${prefix}.${chr}.variable
cat split_files/${prefix}.${chr}.variable.bim | cut -f2 > variable_rsids

#need to only keep things where allele matches well

#I lose a very small number of rsids because there is chromosomal disagreement

cat all_rsids | fgrep -v -w -f variable_rsids > ref_rsids
zcat ~/athena/refs/Homo_sapiens_assembly19.dbsnp138.vcf.gz | head -10000 | fgrep '#' | head -n -1 > ref.vcf
echo "#CHROM    POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT FAKE_PERSON" | sed 's/  */ /g' | tr ' ' '\t' >> ref.vcf
zcat ~/athena/refs/Homo_sapiens_assembly19.dbsnp138.vcf.gz | fgrep -v '#' | awk -v var="$chr" '$1 == var {print $0}' | fgrep -w -f ref_rsids > start_ref
len=`cat start_ref | wc -l`
yes "GT:GQ:DP:MQ:PS:PQ:RFQUAL:FT" | head -$len > one_col
yes "0|0" | head -$len > two_col
paste start_ref one_col two_col >> ref.vcf
plink --vcf ref.vcf --make-bed --out split_files/${prefix}.${chr}.ref

rm ref.vcf one_col two_col start_ref ref_rsids variable_rsids all_rsids all_ids
rm split_files/*nosex
rm split_files/*log
