library(data.table)

prefix_name <- "full_fifty_run"

all_files <- list.files("small_score_files/", pattern = "profile")
all_files <- grep(prefix_name, all_files, value = T)

split_files <- strsplit(all_files, ".", fixed = T)
all_person <- unlist(lapply(split_files, function(x) x[2]))
all_type <- unlist(lapply(split_files, function(x) x[3]))
all_ss <- unlist(lapply(split_files, function(x) x[4]))
all_chr <- unlist(lapply(split_files, function(x) x[5]))

u_all_ss <- unique(all_ss)
u_all_person <- unique(all_person)

i <- 1
system(paste0("zstd -d small_score_files/", all_files[i]))
sub_score <- as.data.frame(fread(paste0("small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4))))
system(paste0("rm small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4)))
write.table(sub_score[,1], paste0("final_scores/", prefix_name, ".eid_ids.txt"), row.names = F, col.names = F, quote = F)


scores <- data.frame(matrix(0, nrow = nrow(sub_score), ncol = length(u_all_ss)))
comp_percs <- rep(0, length(all_files))
colnames(scores) <- u_all_ss

all_percs <- matrix(0, nrow = length(all_files), ncol = sum(grepl("new_", sub_score$FID)))

for(i in 1:length(all_files)){
  system(paste0("zstd -d small_score_files/", all_files[i]))
  sub_score <- as.data.frame(fread(paste0("small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4))))
  print(paste("sub_score", i, "size: ", dim(sub_score)))
  print(paste("new in sub_score", i, "size: ", sum(grepl("new_", sub_score$FID))))

  percentile <- ecdf(sub_score$SCORESUM[!grepl("new_", sub_score$FID)])
  all_percs[i,] <- percentile(sub_score$SCORESUM[grepl("new_", sub_score$FID)])

  system(paste0("rm small_score_files/", substr(all_files[i], 1, nchar(all_files[i])-4)))

  scores[,which(u_all_ss == all_ss[i])] <- scores[,which(u_all_ss == all_ss[i])] + sub_score$SCORESUM 
}

rownames(all_percs) <- all_files


write.table(scores, paste0("final_scores/", prefix_name, "_score.txt"), row.names = F, col.names = T, quote = F, sep = ' ')
saveRDS(scores, paste0("final_scores/", prefix_name, "_score.RDS"))


score_percs <- matrix(NA, nrow = ncol(scores), ncol = sum(grepl("new_", sub_score$FID)))
rownames(score_percs) <- u_all_ss
colnames(score_percs) <- grep("new_", sub_score[,1], value = T)
for(i in 1:ncol(scores)){
  percentile <- ecdf(scores[!grepl("new_", sub_score$FID),i])
  score_percs[i,] <- percentile(scores[grepl("new_", sub_score$FID),i])
}

saveRDS(list("total" = score_percs, "by_chromosome" = all_percs), paste0("final_scores/", prefix_name, "_percs.RDS"))

