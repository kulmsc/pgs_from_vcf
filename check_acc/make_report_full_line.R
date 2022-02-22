library(magick)

report_ind <- 20
#report_ind <- commandArgs(trailingOnly=TRUE)

prefix_name <- "jc_full_run"
author_pheno <- "coronary artery disease"

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


##################################################################


eids <- read.table("misc_files/eid_order.txt", stringsAsFactors = F)

indiv_name <- strsplit(eids[report_ind,1], "_")[[1]][2]
indiv_sex <- "NA"
indiv_dob <- "NA/NA/NA"

test_order_by <- "Dr. Cheung"
order_date <- "NA/NA/NA"
disease_ordered <- "CAD"
author_pheno <- "coronary artery disease"

test_complete_by <- "Scott Kulm"
complete_day <- format(Sys.Date(), "%m/%d/%Y")
variants_employed <- readRDS(paste0("meta_plots/best_score.", save_name, ".RDS"))
if(variants_employed == "pick_name"){
  variants_employed <- "Custom"
  dict <- read.table("misc_files/meta_stats", stringsAsFactors=F, header=T, sep=",")
  dict$author <- tolower(dict$author)
  if(author_pheno %in% dict$trait){
    pick_name <- dict$author[dict$trait == author_pheno]
  } 

  num_snps <- 0
  for(f in list.files("../mod_sets/", pick_name)){
    num_snps <- num_snps + nrow(read.table(paste0("../mod_sets/", f), header=T))
  }
} else {
  num_snps <- 0
  for(f in list.files("../pgs_sets/", variants_employed)){
    num_snps <- num_snps + nrow(read.table(paste0("../pgs_sets/", f), header=T))
  }
}

risk_group <- "NA"
rank_group <- "NA"
rec1 <- "NA"
rec2 <- "NA"

or_val <- signif(readRDS(paste0("report_plots/uk_or.", report_ind, ".", save_name, ".RDS")),3)
pop <- readRDS(paste0("report_plots/sample_size.", report_ind, ".", save_name, ".RDS"))
uk_perc <- signif(readRDS(paste0("report_plots/uk_risk_perc.", report_ind, ".", save_name, ".RDS")),2)*100
uk_prev <-  signif(readRDS(paste0("report_plots/uk_prev.", report_ind, ".", save_name, ".RDS"))*100,3)

wgs_perc <- signif(readRDS(paste0("report_plots/wgs_risk_perc.", report_ind, ".", save_name, ".RDS"))*100, 2)

score_auc <- signif(readRDS(paste0("report_plots/score_auc.", report_ind, ".", save_name, ".RDS")), 3)


line_space <- 45
group_space <- 70

#####################################################


tiger <- image_read_pdf("~/Downloads/Template_PRS_Report.pdf")

#header field
tiger <- image_annotate(tiger, paste0("Name: ", indiv_name), size = 10, location = "+146+250", color = "black")
tiger <- image_annotate(tiger, paste0("Sex: ", indiv_sex), size = 10, location = "+146+310", color = "black")
tiger <- image_annotate(tiger, paste0("Date of Birth: ", indiv_dob), size = 10, location = "+146+370", color = "black")
tiger <- image_annotate(tiger, paste0("Sample Type: ", "Blood"), size = 10, location = "+146+430", color = "black")

tiger <- image_annotate(tiger, paste0("Name: ", test_order_by), size = 10, location = "+810+250", color = "black")
tiger <- image_annotate(tiger, paste0("Date: ", order_date), size = 10, location = "+810+310", color = "black")
tiger <- image_annotate(tiger, paste0("Disease to be Screened:"), size = 10, location = "+810+370", color = "black")
tiger <- image_annotate(tiger, paste0(disease_ordered), size = 10, location = "+810+430", color = "black")

tiger <- image_annotate(tiger, paste0("Name: ", test_complete_by), size = 10, location = "+1474+250", color = "black")
tiger <- image_annotate(tiger, paste0("Date: ", complete_day), size = 10, location = "+1474+310", color = "black")
tiger <- image_annotate(tiger, paste0("Variant Set Employed:"), size = 10, location = "+1474+370", color = "black")
tiger <- image_annotate(tiger, paste0(variants_employed), size = 10, location = "+1474+430", color = "black")

#box field
tiger <- image_annotate(tiger, paste0(indiv_name, "'s odds of ", disease_ordered),
               size = 15, location = "+320+850", color = "black")
tiger <- image_annotate(tiger, paste0("is ", or_val[2], " times that of a healthy population"),
               size = 15, location = "+320+970", color = "black")

#uk details
tiger <- image_annotate(tiger, "Population: ", size = 10, location = "+283+1140", color = "black", weight = 700)
tiger <- image_annotate(tiger, paste0(pop[1], " healthy and ", pop[2], " individuals with ", disease_ordered),
                        size = 10, location = "+530+1140", color = "black")

tiger <- image_annotate(tiger, "Score: ", size = 10, location = paste0("+283+", 1140+group_space),
                        color = "black", weight = 700)
tiger <- image_annotate(tiger, paste0("The polygenic risk score (PRS) is in the ", uk_perc, "th percentile"),
                        size = 10, location = paste0("+430+", 1140+group_space), color = "black")

tiger <- image_annotate(tiger, "Odds Ratio: ", size = 10, location = paste0("+283+", 1140+group_space*2),
                        color = "black", weight = 700)
tiger <- image_annotate(tiger, paste0("The odds of a UK Biobank individual in the ", uk_perc, "th PRS"),
                        size = 10, location = paste0("+520+", 1140+(group_space*2)), color = "black")
tiger <- image_annotate(tiger, paste0("percentile experiencing ", tolower(disease_ordered), " is ",
                                      or_val[2], "(95% CI ", or_val[1], "-", or_val[3], ")"),
                        size = 10, location = paste0("+283+", 1140+(group_space*2)+line_space), color = "black")
tiger <- image_annotate(tiger, paste0("greater than that of individuals below the 50th percentile"),
               size = 10, location = paste0("+283+", 1140+(group_space*2)+(line_space*2)), color = "black")

tiger <- image_annotate(tiger, "Prevalence: ", size = 10, location = paste0("+283+", 1140+(group_space*3)+(line_space*2)), color = "black", weight = 700)
tiger <- image_annotate(tiger, paste0("The prevalence of ", tolower(disease_ordered), " in the ", uk_perc, "th PRS"),
                        size = 10, location = paste0("+520+", 1140+(group_space*3)+(line_space*2)), color = "black")
tiger <- image_annotate(tiger, paste0("is ", uk_prev, "% whereas the prevalence in the 50th PRS percentile is 4.29%"),
               size = 10, location = paste0("+283+", 1140+(group_space*3)+(line_space*3)), color = "black")

#wcm results
tiger <- image_annotate(tiger, "Population: ", size = 10, location = "+283+1840", color = "black", weight = 700)
tiger <- image_annotate(tiger, paste0(pop[3], " healthy and ", pop[4], " individuals with ", disease_ordered, ","),
                        size = 10, location = "+530+1840", color = "black")
tiger <- image_annotate(tiger, paste0("the same size is too small to generate accuracy statistics"),
               size = 10, location = paste0("+283+", 1840+line_space), color = "black")

tiger <- image_annotate(tiger, "Score: ", size = 10, location = paste0("+283+", 1840+group_space+line_space), color = "black", weight = 700)
tiger <- image_annotate(tiger, paste0("The polygenic risk score is in the ", wgs_perc, "th percentile,"),
                        size = 10, location = paste0("+425+", 1840+group_space+line_space), color = "black")
tiger <- image_annotate(tiger, paste0("shown as the vertical bar in the plot on the left"),
               size = 10, location = paste0("+283+", 1840+group_space+line_space*2), color = "black")


#what this really means
tiger <- image_annotate(tiger, paste0(indiv_name, " is in the ", risk_group, " risk group"),
               size = 13, location = "+100+2520", color = "black")
tiger <- image_annotate(tiger, paste0("The ", risk_group, " is the ", rank_group,
                                      "st highest group. General recommendations for this risk group include"),
               size = 10, location = "+100+2585", color = "black")
tiger <- image_annotate(tiger, paste0(rec1, " and ", rec2, ". The care required in this case may vary."),
               size = 10, location = "+100+2630", color = "black")


#about the test
tiger <- image_annotate(tiger, paste0("Variants Scored:"),
               size = 10, location = "+283+2950", color = "black", weight = 700)
tiger <- image_annotate(tiger, paste0("The ", variants_employed, " set was used to construct the PRS, which is available"),
               size = 10, location = "+615+2950", color = "black")
tiger <- image_annotate(tiger, paste0("from the PGS Catalog. The PRS contains ", num_snps, " variants."),
               size = 10, location = paste0("+283+", 2950+line_space), color = "black")

tiger <- image_annotate(tiger, paste0("Performance: "),
                        size = 10, location = paste0("+283+", 2950+group_space+line_space), color = "black", weight = 700)
tiger <- image_annotate(tiger, paste0("The ", variants_employed, " PRS generates an unadjusted AUC of ",
                                      score_auc[2], " (95% CI ", score_auc[1], "-", score_auc[2], ")"),
                        size = 10, location = paste0("+570+", 2950+group_space+line_space), color = "black")


new_image <- image_read(paste0("report_plots/bell.", report_ind, ".", save_name, ".png"))
tiger <- image_composite(tiger, image_scale(new_image, 530), operator = "atop", offset = "+1740+1770")

new_image <- image_read(paste0("report_plots/prev.", report_ind, ".", save_name, ".png"))
tiger <- image_composite(tiger, image_scale(new_image, 780), operator = "atop", offset = "+1695+785")


image_write(tiger, paste0("actual_reports/report.", report_ind, ".pdf"), format = "pdf")

#rm(list=ls())

