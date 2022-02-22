
ss_name=$1
prefix=$2
chr=$3
i=$4

echo $ss_name $prefix $chr $i



if [ ! -e full_small_score_files/score.${prefix}.${ss_name}.${chr}.profile.zst ];then

  nsnps_refine=`cat ../mod_sets/${ss_name}.${chr}.ss | wc -l`

  #The VCF ###
  plink --memory 12000 --threads 12 --bfile use_done --keep-allele-order --score ../mod_sets/${ss_name}.${chr}.ss 3 4 7 sum --allow-extra-chr --out small_score_files/score.${prefix}.full.${ss_name}.${chr}

  zstd --rm small_score_files/score.${prefix}.full.${ss_name}.${chr}.profile

  nsnps_scored=`cat small_score_files/score.${prefix}.full.${ss_name}.${chr}.log  | fgrep valid | cut -f2 -d' '`

  echo nsnps_refine $nsnps_refine >> small_report_files/${prefix}.${ss_name}.${chr}.txt
  echo nsnps_scored $nsnps_scored >> small_report_files/${prefix}.${ss_name}.${chr}.txt

else

  echo already exists

fi

