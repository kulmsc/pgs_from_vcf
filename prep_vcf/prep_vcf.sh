#filename=$1
ls ../raw_vcf/ | cut -f1 -d'.' > temp1
ls | fgrep vcf.gz | cut -f1 -d'.' > temp2
cat temp1 | fgrep -w -v -f temp2 > todo

ls ../raw_vcf/ | fgrep -w -f todo | while read filename;do

  prefix=`echo $filename | cut -f1 -d'.'`
  if [ ! -e ${prefix}.lifted.vcf.gz ];then

    #liftover to grch37
    zcat ../raw_vcf/$filename | fgrep -v "*" | gzip > ready_for_lift.vcf.gz
    java -jar ~/Programs/picard/build/libs/picard.jar LiftoverVcf I=ready_for_lift.vcf.gz O=temp.lifted.vcf  CHAIN=/home/kulmsc/athena/refs/for_picard/hg38ToHg19.over.chain REJECT=rejected_variants.vcf R=/home/kulmsc/athena/refs/for_picard/hg19.fa

    cat temp.lifted.vcf | awk '{gsub(/^chr/,""); print}' > new.vcf
    rm temp.lifted.vcf

    bgzip new.vcf

    bcftools index new.vcf.gz

    bcftools annotate -a /home/kulmsc/athena/refs/Homo_sapiens_assembly19.dbsnp138.vcf.gz -c ID -o ${prefix}.lifted.vcf new.vcf.gz

    gzip ${prefix}.lifted.vcf
    rm new.vcf.gz
    rm new.vcf.gz.csi
    rm temp.lifted.vcf.idx
    rm rejected_variants.vcf
    rm ready_for_lift.vcf.gz

    ./get_bad_rsids.sh ${prefix}
  fi
done

rm temp1 temp2 todo
