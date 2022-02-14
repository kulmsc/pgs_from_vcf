
prefix_name <- "ms_new_run"

scores <- readRDS(paste0("final_scores/", prefix_name, "_score.RDS"))

#somehow this bit got cut-off so need to make it from the components
#temp <- scores$shah.2822944P.variable_SUM + scores$shah.2822944P.ref_SUM
#scores$shah.5547789P.variable_SUM <- temp - scores$shah.5547789P.ref_SUM
#scores$shah.all.all_SUM <- temp

#extra_ids <- grep("GEN", scores$IID, value=T)

eids <- read.table(paste0("final_scores/", prefix_name, ".eid_ids.txt"), stringsAsFactors=F)
#eids <- read.table("final_scores/eid_ids.txt", stringsAsFactors=F)
#eids <- scores[,1,drop=F]
#scores <- scores[,grep("all", colnames(scores))]
#colnames(scores) <- unlist(lapply(strsplit(colnames(scores), ".", fixed = T), function(x) x[1]))

construct_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
construct_eids <- construct_eids[order(construct_eids[,1]),]
construct_eids <- construct_eids[-length(construct_eids)]

scores <- scores[eids[,1] %in% construct_eids | grepl("new_", eids[,1]),,drop=F]
eids <- eids[eids[,1] %in% construct_eids | grepl("new_", eids[,1]),,drop=F]

#for(author in colnames(scores)){
for(author in c("christophersen", "nikpay")){

  pheno <- read.table(paste0("~/athena/doc_score/analyze_score/construct_defs/pheno_defs/diag.", tolower(author), ".txt.gz"), stringsAsFactors=F)
  pheno <- pheno[,1:4]
  pheno <- rowSums(pheno)
  pheno[pheno > 1] <- 1

  pheno_eids <- read.table("~/athena/doc_score/analyze_score/construct_defs/eid.csv", header = T)
  pheno_eids <- pheno_eids[order(pheno_eids[,1]),]
  pheno_eids <- pheno_eids[-length(pheno_eids)]

  subpheno <- pheno[pheno_eids %in% eids[,1]]
  subeids <- pheno_eids[pheno_eids %in% eids[,1]]
  subpheno <- subpheno[order(subeids)[rank(eids[1:length(subeids),1])]]
  subpheno <- c(subpheno, rep(NA, nrow(scores) - length(subpheno)))

  scores[paste0(author, "_pheno")] <- subpheno

}


saveRDS(eids, paste0("final_scores/", prefix_name, ".eids_wpheno.RDS"))
saveRDS(scores, paste0("final_scores/", prefix_name, ".score_wpheno.RDS"))
