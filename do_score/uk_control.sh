rm comp_sizes bad_rsid
prefix=$1 #of the small_score output
#prefix=imp2

#watch out for the grep !!!
for chr in {1..22};do

  rm ../get_all_snps/temp_for_scoring/*

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

  #--- make chomosome specific file from the wcm data
  rm ../get_all_snps/temp_for_scoring/*   
  ls ../get_all_snps/output_files/*bed | cut -f1-4 -d'.' | while read bigline;do
    x1=`echo $bigline | cut -f4 -d'/' | cut -f1 -d'.'`
    x2=`echo $bigline | cut -f4 -d'/' | cut -f2 -d'.'`
    plink --bfile $bigline --chr $chr --extract all_rsids --make-bed --out ../get_all_snps/temp_for_scoring/${x1}.${chr}.${x2}
  done

  #--- get list of all files to be merged
  echo use > to_merge
  ls ../get_all_snps/temp_for_scoring/*.${chr}.variable.bed | cut -f1-5 -d'.' >> to_merge
  ls ../get_all_snps/temp_for_scoring/*.${chr}.ref.bed | cut -f1-5 -d'.' >> to_merge

  #--- fix places pos or alleles  mismatches with same rsid between nygc and ukbb
  Rscript fix_pos.R
  cat bad_rsid | wc -l > small_report_files/uk_mismatch_bad_rsid.${chr}.txt

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
  cat bad_rsid | wc -l > small_report_files/full_bad_rsid.${chr}.txt

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
  cat use_done-merge.missnp | wc -l > small_report_files/missnp_bad_rsid.${chr}.txt
  cat use_done-merge.missnp >> bad_rsid 
  cat to_merge | while read line;do
    short_line=`echo $line | cut -f4 -d'/'`
    plink --bfile $line --exclude bad_rsid --make-bed --out temp_files/${short_line}
  done

  plink --merge-list real_to_merge --make-bed --out use_done

  rm use.bed use.bim use.fam temp_files/use.bed temp_files/use.bim temp_files/use.fam


  # ----------------- Note this is an alternative way of calculating polygenic risk scores
  #It combines all of the variants sets (prepared from PGSCatalog or elsewhere) into a single file
  #And then scores that combined variant set file at one time
  #Thereby generating a file which contains multiple types of polygenic risk scores
  #The process is slightly more convoluted and has some round-off error (I think) so I stuck with
  #The direct approach, although if things get really slow and there are tons of scores
  #Then this could be a good thing to try
  #./get_subset_rsids.sh $chr

  #Rscript align_sumstats.R $chr
  #mv big_mod_set big_mod_set.${chr}
  #num_cols=`head -1 big_mod_set | cut -f3-300 | tr '\t' '\n' | wc -l`
  
  #plink2 --bfile use_done --score big_mod_set 1 2 header-read cols=+scoresums --score-col-nums 3-${num_cols} --out big_small_scores/res.${chr}
  
  #gzip big_small_scores/res.${chr}.sscore #Note this directory does not even exist right now
  #-----------------------------------------------------------------------------------------


  #--- iterate through polygenic risk scores that I derived
  cat ss_names | while read ss;do
      ./uk_simple_score.sh $ss $prefix $chr $chr
      rm full_small_score_files/*nosex
      sleep 30
  done

  #--- iterate through polygenic risk scores from pgs catalog
  cat pgs_names | while read ss;do
      ./pgs_simple_score.sh $ss $prefix $chr $chr
      rm full_small_score_files/*nosex
      sleep 30
  done



  rm big_mod_set.${chr}
  rm use_done*
  rm exp_done*  
  rm temp_files/*

done

rm ny_snps all_rsids bad_rsid to_merge use.log
rm real_to_merge
