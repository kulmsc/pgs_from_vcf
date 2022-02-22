
all_args <- commandArgs(trailingOnly=TRUE)
prefix_name <- all_args[1]
author_pheno <- all_args[2]

#convert the nice, interpreatble disease name
dict <- read.table("~/athena/doc_score/raw_ss/meta_stats", stringsAsFactors=F, header=T, sep=",")
dict$author <- tolower(dict$author)
if(author_pheno %in% dict$trait){
  author_pheno <- dict$author[dict$trait == author_pheno]
} else {
  print("your trait name is incorrect")
  print("your trait name must be listed within the file ~/athena/doc_score/raw_ss/meta_stats")
  stop()
}

#read in the scores
scores <- readRDS(paste0("final_scores/", prefix_name, "_score.RDS"))
eids <- read.table(paste0("final_scores/", prefix_name, ".eid_ids.txt"), stringsAsFactors=F)

#need to get the eids of the ukbb
construct_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
construct_eids <- construct_eids[order(construct_eids[,1]),]
construct_eids <- construct_eids[-length(construct_eids)]

scores <- scores[eids[,1] %in% construct_eids | grepl("new_", eids[,1]),,drop=F]
eids <- eids[eids[,1] %in% construct_eids | grepl("new_", eids[,1]),,drop=F]

#can in theory combine multiple phenotypes into the score file
for(author in c(author_pheno)){

  #read in the pheno, and only look at the columns corresponding to self-reports, icd, opcs
  pheno <- read.table(paste0("~/athena/doc_score/analyze_score/construct_defs/pheno_defs/diag.", tolower(author), ".txt.gz"), stringsAsFactors=F)
  pheno <- pheno[,1:4]
  pheno <- rowSums(pheno)
  pheno[pheno > 1] <- 1

  #need to get the eids of the phenotypes
  pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
  pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
  pheno_eids <- pheno_eids[-length(pheno_eids)]

  #line everything up
  subpheno <- pheno[pheno_eids %in% eids[,1]]
  subeids <- pheno_eids[pheno_eids %in% eids[,1]]
  subpheno <- subpheno[order(subeids)[rank(eids[1:length(subeids),1])]]
  subpheno <- c(subpheno, rep(NA, nrow(scores) - length(subpheno)))

  scores[paste0(author, "_pheno")] <- subpheno

}


saveRDS(eids, paste0("final_product/", prefix_name, ".eids_wpheno.RDS"))
saveRDS(scores, paste0("final_product/", prefix_name, ".score_wpheno.RDS"))

saveRDS(eids, paste0("../check_acc/score_files/", prefix_name, ".eids_wpheno.RDS"))
saveRDS(scores, paste0("../check_acc/score_files/", prefix_name, ".score_wpheno.RDS"))
