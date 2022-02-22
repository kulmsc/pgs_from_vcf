library(reshape2)
library(ggplot2)
library(cowplot)
library(stringr)
library(pROC)
library(viridis)
library(epitools)
theme_set(theme_cowplot())

prefix_name <- "jc_full_run"
disease_identifier <- "new_5710|new_GEN21" #these are the prefixes of the individuals with disease
author_pheno <- "coronary artery disease" #again this must be in the meta_stats file

dict <- read.table("misc_files/meta_stats", stringsAsFactors=F, header=T, sep=",")
dict$author <- tolower(dict$author)
if(author_pheno %in% dict$trait){
  pick_name <- dict$author[dict$trait == author_pheno]
} else {
  print("your trait name is incorrect")
  print("your trait name must be listed within the file ~/athena/doc_score/raw_ss/meta_stats")
  stop()
}

save_name <- gsub(" ", "_", tolower(author_pheno))
save_name <- paste0(save_name, ".", prefix_name)


#########################

scores <- readRDS(paste0("score_files/", prefix_name, ".score_wpheno.RDS"))
eids <- readRDS(paste0("score_files/", prefix_name, ".eids_wpheno.RDS"))
train_eid <- read.table("misc_files/train_eid.0.6.txt", stringsAsFactors = F)

#must look at train_eid because my own polygenic risk scores were constructed from train_eid
#to prevent overfitting, although if only pgs catalog scores are used then we can get rid of this

scores <- scores[!(eids[,1] %in% train_eid[,1]),]
eids <- eids[!(eids[,1] %in% train_eid[,1]),]

#change the name of the score that I self-derived to a generic name such that it
#does not need to be specified differently
colnames(scores)[colnames(scores) == pick_name] <- "pick_name"
colnames(scores)[colnames(scores) == paste0(pick_name, "_pheno")] <- "pick_name_pheno"
scores <- scores[,grepl("pick", colnames(scores)) | grepl("PGS", colnames(scores))]


#look at only ukbb for the score I derived
ukbb <- scores[!is.na(scores$pick_name_pheno),]
ukbb$pick_name_pheno <- as.factor(ukbb$pick_name_pheno)
the_plot <- ggplot(ukbb, aes(pick_name, fill = pick_name_pheno)) + geom_density(alpha = 0.5) +
  labs(x = paste0("Self-Derived ", save_name, " PRS"), y = "Density", fill = paste0("Has ", save_name, "?")) +
  scale_fill_discrete(labels = c("Yes", "No"))
plot(the_plot)
ggsave(paste0("meta_plots/ukbb_dists_", save_name, ".png"), the_plot, width = 5, height = 4)


#evaluate each of the polygenic risk scores against the disease for individuals in the ukbb
#measuring accuracy through auc (could be different in the future)
all_auc <- as.data.frame(do.call("rbind", lapply(1:(ncol(ukbb)-1), function(x) as.numeric(ci.auc(roc(ukbb$pick_name_pheno ~ ukbb[,x]))))))
all_auc$score <- factor(colnames(ukbb)[1:(ncol(ukbb)-1)], levels = colnames(ukbb)[1:(ncol(ukbb)-1)][order(all_auc$V2)])
the_plot <- ggplot(all_auc, aes(V2, score)) + geom_point() +
  geom_errorbarh(aes(xmin = V1, xmax = V3, height = 0)) +
  labs(x = "AUC", y = "Score Name")
plot(the_plot)
ggsave(paste0("meta_plots/prs_auc.", save_name, ".png"), the_plot, width = 4.5, height = 4)
best_name <- as.character(all_auc$score[which.max(all_auc[,2])])
saveRDS(best_name, paste0("meta_plots/best_score.", save_name, ".RDS")) 

#change the name of the best score to just a genetic "best score"
scores$best_score <- scores[,colnames(scores) == best_name]


#establish the phenotypes of the non ukbb individuals
#3 means non-disease of wcm individuals
#4 means disease of wcm individuals
scores$all_pheno <- as.numeric(as.character(scores$pick_name_pheno))
scores$all_pheno[is.na(scores$all_pheno)] <- 3
for(curr_id in strsplit(disease_identifier, "|", fixed = T)[[1]]){
  scores$all_pheno[grep(curr_id, eids)] <- 4
}
scores$all_pheno <- as.factor(scores$all_pheno)

#if things went wrong previously in picking the effect allele then we may need to switch the sign
if(mean(scores$best_score[scores$all_pheno == 0]) > mean(scores$best_score[scores$all_pheno == 1])){
  print("FLIP")
  scores$best_score <- as.numeric(scale(-1*scores$best_score))  
}

#make density plots across all 4 phenotypes
the_plot <- ggplot(scores, aes(best_score, fill = all_pheno)) + geom_density(alpha = 0.5) +
  labs(x = "PRS", y = "Density", fill = "Has Disease?") +
  scale_fill_discrete(labels = c("UKBB No", "UKBB Yes", "WGS No", "WGS Yes"))
plot(the_plot)
ggsave(paste0("meta_plots/all_dists.", save_name, ".png"), the_plot, width = 5, height = 4)


#compare the mean and sd of each distribution plotted previously
#and calculate p-values, then plot
score_stats <- as.data.frame(do.call("rbind", lapply(c(0, 1, 3, 4), function(x) c(mean(scores$best_score[scores$all_pheno == x]),
                                 sd(scores$best_score[scores$all_pheno == x])/sqrt(sum(scores$all_pheno == x)) ))))
score_stats$group <- c("UKBB Controls", "UKBB Cases", "WGS Controls", "WGS Cases")
ukbb_pval <- t.test(scores$best_score[scores$all_pheno==1], scores$best_score[scores$all_pheno==0])$p.value
wgs_pval <- t.test(scores$best_score[scores$all_pheno==4], scores$best_score[scores$all_pheno==3])$p.value
the_plot <- ggplot(score_stats, aes(V1, group)) + geom_point() +
  labs(x = "Mean PRS", y = "") +
  geom_errorbarh(aes(xmin = V1-V2, xmax = V1+V2, height = 0)) +
  annotate(geom = "text", label = paste0("P = ", signif(wgs_pval, 2)),
           x = score_stats$V1[score_stats$group == "WGS Cases"], y = 3.2, size = 3) +
  annotate(geom = "text", label = paste0("P = ", as.character(signif(ukbb_pval,3)) ),
           x = score_stats$V1[score_stats$group == "UKBB Cases"], y = 1.2, size = 3)
plot(the_plot)
ggsave(paste0("meta_plots/all_prs_stats.", save_name, ".png"), the_plot, width = 4, height = 2.5)



#plot just the box of the wcm individuals
the_plot <- ggplot(scores[scores$all_pheno %in% 3:4,], aes(best_score, all_pheno)) + 
  geom_boxplot() + geom_point() + 
  labs(x = "Polygenic Risk Score", y = paste0("Has ", save_name, "?")) +
  scale_y_discrete(labels = c("No", "Yes"))
plot(the_plot)
ggsave(paste0("meta_plots/prs_spectrum.", save_name, ".png"), the_plot, width = 4, height = 2.5)


#######################################

#big function that calculates different accuracy statistics
#for each of the wcm individuals
#the accuracies are measured with respect to either ukbb or wcm individuals
get_or <- function(mid_quant, type = "or", compto = c(0,1), buffer = 0.01){

  #compto c(0,1) means ukbb is the comparative population
  #compto c(3,4) means wcm is the comparative population
  #we create an "healthy" population in lo
  #then the population that represents the actual individuals in hi 
  #  the size of the population is determined by buffer
  lo <- quantile(scores$best_score[scores$all_pheno %in% compto], max(c(mid_quant-buffer, 0)))
  hi <- quantile(scores$best_score[scores$all_pheno %in% compto], min(c(1, mid_quant+buffer)))
  comp <- quantile(scores$best_score[scores$all_pheno %in% compto], 0.5)
  phen1 <- scores$all_pheno[scores$all_pheno %in% compto & scores$best_score > lo & scores$best_score < hi]
  phen2 <- scores$all_pheno[scores$all_pheno %in% compto & scores$best_score < comp]
  if(max(phen1) > 1){
    phen1 <- phen1 - 3
    phen2 <- phen2 - 3
  }
  
  #from our two populations we calculate some simple accuracy statistics
  #we only return the statistic we are interested in
  orval <- (sum(phen1 == 1)/sum(phen1 == 0))/(sum(phen2 == 1)/sum(phen2 == 0))
  rr <- (sum(phen1 == 1)/length(phen1))/(sum(phen2 == 1)/length(phen2))
  prev <- (sum(phen1 == 1)/length(phen1))
  
  if(type == "or"){
    return(orval)
  } else if(type == "ci"){
    val <- c(sum(phen1 == 1), sum(phen1 == 0), sum(phen2 == 1), sum(phen2 == 0))
    return(oddsratio(matrix(val, nrow=2))$measure[2,c(2,1,3)])
  } else if(type == "rr"){
    val <- matrix(c(sum(phen1 == 1), sum(phen2 == 1), sum(phen1 == 0), sum(phen2 == 0)), nrow = 2)
    val <- val[c(2,1),c(2,1)]
    return(riskratio(matrix(val, nrow=2))$measure[2,c(2,1,3)])    
  } else if(type == "ci_prev"){
    val <- prop.test(sum(phen1==1), length(phen1))
    return(c(val$conf.int[1], prev, val$conf.int[2]))
  } else {
    return(prev)
  }
}

#get the reference prevalence of each percentile of the ukbb
scores$all_pheno <- as.numeric(as.character(scores$all_pheno))
temp <- scores[scores$all_pheno %in% c(0,1),]
hundo_prev <- rep(0, 100)
for(i in 1:100){
  sub_pheno <- temp$all_pheno[temp$best_score > quantile(temp$best_score, (i-1)/100) &
                              temp$best_score <= quantile(temp$best_score, i/100)]
  hundo_prev[i] <- sum(sub_pheno)/length(sub_pheno)
}


#iterate through each of the wcm individuals and make a whole bunch of accuracy statistics
#here the comp pop is ukbb
uk_risk_perc <- ecdf(scores$best_score[scores$all_pheno <= 1])(scores$best_score[scores$all_pheno == 4])
uk_odds_ratio <- unlist(lapply(uk_risk_perc, get_or))
uk_prev <- unlist(lapply(uk_risk_perc, get_or, "prev"))
uk_ciprev <- do.call("rbind", lapply(uk_risk_perc, get_or, "ci_prev"))
uk_rr <- do.call("rbind", lapply(uk_risk_perc, get_or, "rr"))

#here the comp pop is wcm
wgs_risk_perc <- ecdf(scores$best_score[scores$all_pheno > 1])(scores$best_score[scores$all_pheno == 4])
wgs_odds_ratio <- unlist(lapply(wgs_risk_perc, get_or, "or", c(3,4), 0.05))
wgs_prev <- unlist(lapply(wgs_risk_perc, get_or, "prev", c(3,4), 0.05))


######################################################################

#plot the percentile of all ukbb individual against their prevalence
plot_df <- data.frame("risk_perc" = c(uk_risk_perc, wgs_risk_perc), 
                      "or" = c(uk_odds_ratio, wgs_odds_ratio), 
                      "prev" = c(uk_prev, wgs_prev),
                      "id" = rep(c("UKBB", "WGS"), each = length(uk_risk_perc)), stringsAsFactors = F)

the_plot <- ggplot(plot_df, aes(risk_perc, or, color = id )) + geom_point() + xlim(c(0,1)) +
  labs(y = "Odds Ratio", x = "PRS Percentile", color = "Comp.\nPop")
plot(the_plot)
ggsave(paste0("meta_plots/wgs_or.", save_name, ".png"), the_plot, width = 5, height = 3.5)

the_plot <- ggplot(plot_df[plot_df$id == "UKBB",], aes(risk_perc, prev )) + geom_point() +
  xlim(c(0,1)) + scale_color_viridis() +
  labs(y = paste0(save_name, " Prevalence"), x = "PRS Percentile", caption = "Compared to UKBB")
plot(the_plot)
ggsave(paste0("meta_plots/wgs_prev.", save_name, ".png"), the_plot, width = 4.5, height = 3.5)


######################################################################

comp <- quantile(scores$best_score[scores$all_pheno %in% c(0,1)], 0.5)
val <- scores$all_pheno[scores$all_pheno %in% c(0,1) & scores$best_score < comp]
temp <- prop.test(sum(val==1), length(val))
val <- c(temp$conf.int[1], temp$est, temp$conf.int[2])

#now actually going through and making the plots and saving the statistics
#that will be apart of the prs reports
for(ind in 1:sum(scores$all_pheno == 4)){

  #create bar plot that compares disease prevalences
  prev_df <- data.frame(rbind(val, uk_ciprev[ind,]))
  prev_df$group <- c("baseline", "indiv")
  bval <- floor((uk_risk_perc[ind]-0.005) * 100)
  tval <- bval+1
  the_plot <- ggplot(prev_df, aes(group, p)) + geom_col() +
    labs(x = "", y = "Prevalence") +
    scale_x_discrete(labels = c("Healthy\nPopulation", "")) +
    geom_errorbar(aes(ymin = V1, ymax = V3, width = 0.1)) +
    annotate("text", x = 1, y = prev_df$p[1]/2, label = "PRS Percentile:\n1-50", color = "white") +
    annotate("text", x = 2, y = prev_df$p[1]/2, label = paste0("PRS Percentile:\n", bval, "-", tval), color = "white")
  plot(the_plot)
  ggsave(paste0("report_plots/bar_prev.", ind, ".", save_name, ".png"), the_plot, width = 4, height = 4)
  
  #create distribution plot with individual as vertical line
  the_plot <- ggplot(scores[scores$all_pheno == 3,], aes(best_score)) + geom_density(size=2) +
    geom_vline(aes(xintercept = scores$best_score[scores$all_pheno == 4][ind]), color = "red", size = 1.5) +
    labs(x = "Polygenic Risk Score", y = "Density")
  plot(the_plot)
  ggsave(paste0("report_plots/bell.", ind, ".", save_name, ".png"), the_plot, width = 3, height = 3)
  
  #create a prevalence plot including all percentiles of ukbb with indiv as big red dot
  prev_df <- data.frame("prev" = hundo_prev, "x" = 1:100)
  the_plot <- ggplot(prev_df, aes(x, prev)) + geom_point() +
    geom_point(aes(x = uk_risk_perc[ind]*100, y = uk_prev[ind]), color = "red", size = 5) +
    labs(x = "Polygenic Risk Score Percentile", y = "Prevalence of A. Fib.")
  plot(the_plot)
  ggsave(paste0("report_plots/prev.", ind, ".", save_name, ".png"), the_plot, width = 4, height = 3.5)
  
  #save all the stats
  saveRDS(uk_risk_perc[ind], paste0("report_plots/uk_risk_perc.", ind, ".", save_name, ".RDS"))
  saveRDS(uk_rr[ind,], paste0("report_plots/uk_rr.", ind, ".", save_name, ".RDS"))
  saveRDS(uk_prev[ind], paste0("report_plots/uk_prev.", ind, ".", save_name, ".RDS"))
  saveRDS(uk_ciprev[ind,], paste0("report_plots/uk_ciprev.", ind, ".", save_name, ".RDS"))
  saveRDS(as.numeric(get_or(uk_risk_perc[ind], "ci")), paste0("report_plots/uk_or.", ind, ".", save_name, ".RDS"))
  
  saveRDS(wgs_risk_perc[ind], paste0("report_plots/wgs_risk_perc.", ind, ".", save_name, ".RDS"))
  
  temp <- as.numeric(ci.auc(roc(scores$all_pheno[scores$all_pheno %in% c(0,1)] ~ scores$best_score[scores$all_pheno %in% c(0,1)])))
  saveRDS(temp, paste0("report_plots/score_auc.", ind, ".", save_name, ".RDS"))
  
  temp <- c(sum(scores$all_pheno==0), sum(scores$all_pheno==1), sum(scores$all_pheno==3), sum(scores$all_pheno==4))
  saveRDS(temp, paste0("report_plots/sample_size.", ind, ".", save_name, ".RDS"))
  
  saveRDS(best_name, paste0("report_plots/score_name.", ind, ".", save_name, ".RDS"))


}


