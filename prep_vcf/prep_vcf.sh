
#get the names of all the file prefixes that have not already been prepped
#those files need to be located in the raw_vcf directory
ls ../raw_vcf/ | cut -f1 -d'.' > temp1
ls | fgrep vcf.gz | cut -f1 -d'.' > temp2
cat temp1 | fgrep -w -v -f temp2 > todo

#iterating through all of the files that need to be prepped
ls ../raw_vcf/ | fgrep -w -f todo | while read filename;do

  #get just the prefix of the file
  prefix=`echo $filename | cut -f1 -d'.'`
  if [ ! -e ${prefix}.lifted.vcf.gz ];then  #double check that prepped file does not exist already

    #liftover to grch37
    zcat ../raw_vcf/$filename | fgrep -v "*" | gzip > ready_for_lift.vcf.gz
    java -jar ~/Programs/picard/build/libs/picard.jar LiftoverVcf I=ready_for_lift.vcf.gz O=temp.lifted.vcf  CHAIN=/home/kulmsc/athena/refs/for_picard/hg38ToHg19.over.chain REJECT=rejected_variants.vcf R=/home/kulmsc/athena/refs/for_picard/hg19.fa

    #change each variant's chromosome descript from "chr1" to just "1" (for example)
    cat temp.lifted.vcf | awk '{gsub(/^chr/,""); print}' > new.vcf
    rm temp.lifted.vcf

    #compress
    bgzip new.vcf

    #index
    bcftools index new.vcf.gz

    #annotate the file to add rsids
    bcftools annotate -a /home/kulmsc/athena/refs/Homo_sapiens_assembly19.dbsnp138.vcf.gz -c ID -o ${prefix}.lifted.vcf new.vcf.gz

    #compress
    gzip ${prefix}.lifted.vcf

    #remove temporary junk
    rm new.vcf.gz
    rm new.vcf.gz.csi
    rm temp.lifted.vcf.idx
    rm rejected_variants.vcf
    rm ready_for_lift.vcf.gz

    #get the rsids of variants that do not meet quality control specifications 
    #qual score > 15, read depth (dp) > 5 & read depth < mean_dp*3, mapping quality (mq) > 30 
    ./get_bad_rsids.sh ${prefix}
  fi
done

rm temp1 temp2 todo
