rm small_score_files/*
rm small_report_files/*

file_starter=handoff
pheno_name="MS" #must be listed within ~/athena/doc_score/raw_ss/meta_stats

#--- create a polygenic risk score for each chromsome and variant set
./uk_control.sh $file_starter

#--- combine the single chromosome polygenic risk scores
Rscript uk_simple_assemble.R $file_starter

#--- add the phenotype
Rscript finish.R $file_starter $pheno_name
