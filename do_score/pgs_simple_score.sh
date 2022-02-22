
ss_name=$1
prefix=$2
chr=$3
i=$4

echo $ss_name $prefix $chr $i

#get the length of the files
chr_col=`zcat ../raw_catalog_sets/${ss_name}.txt.gz | head -20 | fgrep chr_name | tr '\t' '\n' | fgrep -n chr_name | cut -f1 -d':'`
nsnps_raw=`zcat ../raw_catalog_sets/${ss_name}.txt.gz | fgrep -v '#' | fgrep -v "chr_name" | cut -f${chr_col} | fgrep -w $chr | wc -l`
nsnps_refine=`cat ../pgs_sets/${ss_name}.${chr}.ss | tail -n+2 | wc -l`


if [ ! -e full_small_score_files/score.${prefix}.${ss_name}.${chr}.profile.zst ];then


  #The VCF ###
  plink --memory 12000 --threads 12 --bfile use_done --keep-allele-order --score ../pgs_sets/${ss_name}.${chr}.ss 2 3 4 sum --allow-extra-chr --out small_score_files/score.${prefix}.full.${ss_name}.${chr}

  zstd --rm small_score_files/score.${prefix}.full.${ss_name}.${chr}.profile

  nsnps_scored=`cat small_score_files/score.${prefix}.full.${ss_name}.${chr}.log  | fgrep valid | cut -f2 -d' '`
  
  echo nsnps_raw $nsnps_raw > small_report_files/${prefix}.${ss_name}.${chr}.txt
  echo nsnps_refine $nsnps_refine >> small_report_files/${prefix}.${ss_name}.${chr}.txt
  echo nsnps_scored $nsnps_scored >> small_report_files/${prefix}.${ss_name}.${chr}.txt

else

  echo already exists

fi

