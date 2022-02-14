
ls output_files/ | cut -f1 -d'.' | sort | uniq > temp1
ls ../prep_vcf/ | fgrep vcf.gz | cut -f1 -d'.' | sort | uniq > temp2
cat temp2 | fgrep -w -v -f temp1 > prefix_list
rm temp1 temp2

cat prefix_list | while read pre;do
  if [ ! -e ${pre}.ref.bed ];then
    for chrom in {1..22};do

      ./just_use_rsid.sh $chrom $pre

    done

    ls split_files/${pre}*variable.bed | cut -f1-3 -d'.' > merge_list
    plink --merge-list merge_list --make-bed --out output_files/${pre}.variable

    ls split_files/${pre}*ref.bed | cut -f1-3 -d'.' > merge_list
    plink --merge-list merge_list --make-bed --out output_files/${pre}.ref

    rm split_files/*

  fi
done

./fix_fam.sh
