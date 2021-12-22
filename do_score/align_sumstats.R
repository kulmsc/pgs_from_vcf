library(data.table)

args = commandArgs(trailingOnly=TRUE)
chr <- as.character(args[1])

ss_names <- read.table("ss_names", stringsAsFactors=F)
mod_sets <- list()

for(i in 1:nrow(ss_names)){
  mod_sets[[i]] <- as.data.frame(fread(paste0("../mod_sets/", ss_names[i,1], ".", chr, ".ss")))
}


sub_names <- list.files("subset_rsids")
subset_rsids <- list()

for(i in 1:length(sub_names)){
  subset_rsids[[i]] <- read.table(paste0("subset_rsids/", sub_names[i]), stringsAsFactors = F)
}

master_ms <- data.frame("rsid" = unlist(lapply(mod_sets, function(x) x[,3])), "A1" = unlist(lapply(mod_sets, function(x) x[,4])), stringsAsFactors=F)
master_ms <- master_ms[!duplicated(master_ms[,1]),]

k <- 1
betas <- list()
name_beta <- rep("", length(mod_sets) * (1+length(subset_rsids)))
use_sub_names <- c(sub_names, "all")
for(i in 1:length(mod_sets)){
  for(j in 1:(1+length(subset_rsids))){
    print(paste(i, "of", length(mod_sets)))

    betas[[k]] <- rep(0, nrow(master_ms))
    rele_rsids <- master_ms[master_ms[,1] %in% mod_sets[[i]][,3],1]
    curr_mod_set <- mod_sets[[i]][order(mod_sets[[i]][,3])[rank(rele_rsids)],]

    if(j > length(subset_rsids)){
      betas[[k]][master_ms[,1] %in% curr_mod_set[,3]] <- curr_mod_set[,7]
      name_beta[k] <- paste(ss_names[i,1], use_sub_names[j], use_sub_names[j], sep = "-")
    } else {
      betas[[k]][master_ms[,1] %in% curr_mod_set[,3] & master_ms[,1] %in% subset_rsids[[j]][,1]] <- curr_mod_set[curr_mod_set[,3] %in% subset_rsids[[j]][,1], 7]
      name_beta[k] <- paste(ss_names[i,1], strsplit(strsplit(use_sub_names[j], ".", fixed = T)[[1]][2], "-")[[1]][2], strsplit(use_sub_names[j], ".", fixed = T)[[1]][3], sep = "-")
    }
    
    k <- k + 1
  }
}

betas <- data.frame(do.call("cbind", betas))
colnames(betas) <- name_beta

master_ms <- data.frame(master_ms, betas)
write.table(master_ms, "big_mod_set", row.names = F, col.names = T, quote = F, sep = "\t")
