
#downloaded from:
#https://github.com/privefl/paper-ancestry-matching/blob/5b23bb2a6d0f1f8e1b44b1626e3d9fea7b851c3e/code/prepare-UKBB.R

##### Prepare genotypes ####

# file.symlink("~/NCRR-PRS/faststorage/UKBB/", ".")

write(sapply(1:22, function(chr) {
  paste0("/home/kulmsc/athena/ukbiobank/calls/",
         c(paste0("ukbb.", chr, ".bed"),
           paste0("ukbb.", chr, ".bim"),
           paste0("ukbb.", chr, ".fam")))
}), tmp <- tempfile(), ncolumns = 3)



library(bigsnpr)
snp_plinkQC(
  plink.path = "/home/kulmsc/bin/plink",
  prefix.in = tmp,
  file.type = "--merge-list",
  prefix.out = "data/ukbb",
  geno = 0.01,
  autosome.only = TRUE,
  extra.options = "--memory 100000"
)

#### Prepare other data (self-reported ancestry and PCs) ####

library(bigreadr)
library(dplyr)

## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("ukbb_coding/coding1001.tsv")
unknown <- c("Prefer not to answer", "Do not know", "Mixed",
             "Other ethnic group", "Any other mixed background")
code_country <- filter(fread2("ukbb_coding/coding89.tsv"), selectable == "Y")
code_continent <- filter(fread2("ukbb_coding/coding89.tsv"), selectable == "N")

phen1 <- as.data.frame(fread("~/athena/ukbiobank/phenotypes/ukb26867.csv.gz"))
phen1 <- phen1[,colnames(phen1) %in% c("eid", "21000-0.0", "20115-0.0", paste0("22009-0.", 1:16))]
phen2 <- as.data.frame(fread("~/athena/ukbiobank/phenotypes/ukb41972.csv.gz"))
phen2 <- phen2[,colnames(phen2) %in% c("eid", "21000-0.0", "20115-0.0", paste0("22009-0.", 1:16))]
write.table(merge(phen1, phen2, by = "eid"), "ukbb_coding/needed_pheno.csv", row.names = F, col.names = T, quote = F, sep = ",")

df0 <- fread2(
  "ukbb_coding/needed_pheno.csv",
  select = c("eid", "21000-0.0", "20115-0.0", paste0("22009-0.", 1:16)),
  col.names = c("eid", "pop", "country", paste0("PC", 1:16))
) %>%
  mutate(
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning) %>%
      forcats::fct_recode("Other White" = "Any other white background",
                          "Other Asian" = "Any other Asian background",
                          "Other Black" = "Any other Black background") %>%
      forcats::fct_other(drop = unknown, other_level = NA) %>%
      droplevels(pop, exclude = NA) %>%
      forcats::fct_relevel(c(
        "British", "Irish", "White", "Other White",
        "Indian", "Pakistani", "Bangladeshi", "Chinese", "Other Asian",
        "Caribbean", "African", "Other Black")),

    continent = factor(country, levels = code_country$coding,
                       labels = code_country$parent_id) %>%
      factor(labels = code_continent$meaning),

    country = factor(country, levels = code_country$coding,
                     labels = code_country$meaning)
  )

str(df0)
df0$country[with(df0, is.na(country) & pop == "British")] <- "United Kingdom"
df0$country[with(df0, is.na(country) & pop == "Irish")]   <- "Ireland"
df0$continent[df0$country %in% c("United Kingdom", "Ireland")] <- "Europe"

bed_eid <- fread2("data/ukbb.fam")[[2]]
df <- df0[match(bed_eid, df0$eid), ]
mean(is.na(df$pop))  # 1.6%

saveRDS(df, "data/info_UKBB.rds")
