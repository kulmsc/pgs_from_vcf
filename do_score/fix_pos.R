
all_files <- read.table("to_merge", stringsAsFactors=F)
all_files <- all_files[all_files[,1] != "use",1]

ref_bim <- read.table("use.bim", stringsAsFactors = F)
bad_rsid <- list()
for(f in all_files){
  check <- read.table(paste0(f, ".bim"), stringsAsFactors = F)

  check <- check[order(check[,2]),]
  ref_bim <- ref_bim[order(ref_bim[,2]),]

  check[check[,2] %in% ref_bim[,2],4] <- ref_bim[ref_bim[,2] %in% check[,2],4]

  allele_bool <- (check[check[,2] %in% ref_bim[,2],5] == ref_bim[ref_bim[,2] %in% check[,2],5] &
                  check[check[,2] %in% ref_bim[,2],6] == ref_bim[ref_bim[,2] %in% check[,2],6]) |
                 (check[check[,2] %in% ref_bim[,2],5] == ref_bim[ref_bim[,2] %in% check[,2],6] &
                  check[check[,2] %in% ref_bim[,2],6] == ref_bim[ref_bim[,2] %in% check[,2],5])
  print(sum(!allele_bool))
  bad_rsid[[which(all_files == f)]] <- ref_bim[ref_bim[,2] %in% check[,2],][!allele_bool,2]

  to_comp <- read.table(paste0(f, ".bim"), stringsAsFactors = F)
  check <- check[order(check[,2])[rank(to_comp[,2])],]
  write.table(check, paste0(f, ".bim"), row.names = F, col.names = F, quote = F, sep = "\t")
}

write.table(unique(unlist(bad_rsid)), "bad_rsid", row.names=F, col.names=F, quote=F)
