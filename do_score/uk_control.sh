rm comp_sizes bad_rsid
prefix=jc_full_run #of the small_score output

#watch out for the grep !!!
for chr in {1..12};do
#for chr in 12;do

  #--- get rsids in mod_sets
  rm all_rsids
  #ls ../mod_sets/ | fgrep .${chr}. | while read line;do
  cat ss_names | while read line;do
    cat ../mod_sets/${line}.${chr}.ss | tail -n+2 | cut -f3 >> all_rsids
  done
  cat pgs_names | while read line;do
    cat ../pgs_sets/${line}.${chr}.ss | tail -n+2 | cut -f2 >> all_rsids
  done

  #--- get ukbb data
  bgenix -g ~/athena/ukbiobank/imputed/ukbb.${chr}.bgen -incl-rsids all_rsids > temp.bgen
  plink2 --bgen temp.bgen ref-first --sample ~/athena/ukbiobank/imputed/ukbb.${chr}.sample --keep-fam ~/athena/doc_score/do_score/temp_files/brit_eid --make-bed --out use
  rm temp.bgen

  #--- make chomosome specific file
  rm ../get_all_snps/temp_for_scoring/*   
  ls ../get_all_snps/output_files/*bed | cut -f1-4 -d'.' | while read bigline;do
    x1=`echo $bigline | cut -f4 -d'/' | cut -f1 -d'.'`
    x2=`echo $bigline | cut -f4 -d'/' | cut -f2 -d'.'`
    plink --bfile $bigline --chr $chr --make-bed --out ../get_all_snps/temp_for_scoring/${x1}.${chr}.${x2}
  done

  #--- get list of all files to be merged
  echo use > to_merge
  ls ../get_all_snps/temp_for_scoring/*.${chr}.variable.bed | cut -f1-5 -d'.' >> to_merge
  ls ../get_all_snps/temp_for_scoring/*.${chr}.ref.bed | cut -f1-5 -d'.' >> to_merge

  #--- fix places pos or alleles  mismatches with same rsid between nygc and ukbb
  Rscript fix_pos.R
  cat bad_rsid | wc -l

  ls ../get_all_snps/temp_for_scoring/ | fgrep bim | cut -f1 -d'.' | sort | uniq | while read crp;do
    echo $crp
    if [ -e ../prep_vec/${crp}.bad_qc_rsid ];then
      cat ../${crp}.bad_qc_rsid >> bad_rsid
    fi
    blen=`cat bad_rsid | wc -l`
    echo bad qc snps: $blen

  #--- also need to remove snps that are not in nygc
    cat ../get_all_snps/temp_for_scoring/${crp}.${chr}.ref.bim | cut -f2 > ny_snps
    cat ../get_all_snps/temp_for_scoring/${crp}.${chr}.variable.bim | cut -f2 >> ny_snps
    cat all_rsids | fgrep -v -w -f ny_snps >> bad_rsid
    blen=`cat bad_rsid | wc -l`
    echo bad does not exist in vcf snps: $blen
    clen=`cat all_rsids | fgrep -w -f ny_snps | wc -l`
    echo not in: $clen
    echo
  done

  cat bad_rsid | sort | uniq  > tempy; mv tempy bad_rsid

  #--- remove the bad rsids - which come from fix_pos.R
  rm real_to_merge
  cat to_merge | while read line;do
    short_line=`echo $line | cut -f4 -d'/'`
    plink --bfile $line --exclude bad_rsid --make-bed --out temp_files/${short_line}
    echo temp_files/${short_line} >> real_to_merge
  done

  #--- try to merge
  plink --merge-list real_to_merge --make-bed --out use_done

  #--- try again to merge - this time without allele mismatches between nygc files
  cat use_done-merge.missnp >> bad_rsid 
  cat to_merge | while read line;do
    short_line=`echo $line | cut -f4 -d'/'`
    plink --bfile $line --exclude bad_rsid --make-bed --out temp_files/${short_line}
  done

  plink --merge-list real_to_merge --make-bed --out use_done

  rm use.bed use.bim use.fam temp_files/use.bed temp_files/use.bim temp_files/use.fam

  #./get_subset_rsids.sh $chr

  #Rscript align_sumstats.R $chr
  #mv big_mod_set big_mod_set.${chr}
  #num_cols=`head -1 big_mod_set | cut -f3-300 | tr '\t' '\n' | wc -l`
  
  #plink2 --bfile use_done --score big_mod_set 1 2 header-read cols=+scoresums --score-col-nums 3-${num_cols} --out big_small_scores/res.${chr}
  
  #gzip big_small_scores/res.${chr}.sscore

  cat ss_names | while read ss;do
      ./uk_simple_score.sh $ss $prefix $chr $chr
      rm full_small_score_files/*nosex
      sleep 30
  done

  cat pgs_names | while read ss;do
      ./pgs_simple_score.sh $ss $prefix $chr $chr
      rm full_small_score_files/*nosex
      sleep 30
  done



  #rm big_mod_set.${chr}
  rm use_done*
  #rm exp_done*  
  rm temp_files/*

done

