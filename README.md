# pgs_from_vcf

This repository contains all of the scripts necessary to convert a VCF (generated from whole genome sequencing) into a specific polygenic risk score report.  The entire process has been created to be highly understandable and will require minimal interaction.

The process starts by linking or moving all VCFs into a directory (not provided here) called raw_vcfs.
Then scripts in prep_vcf are run, which changes the genetic coordinates and does some simple quality control.
Next the scripts in get_all_snps are run, which pulls in all of the reference variants that are listed in any of the variant sets that are used to compute polygenic risk scores.
The variants sets that define/constitute the polygenic risk scores are initially listed within a raw_catalog_sets directory, not provided here.  Each of these variant sets should/must be downloaded from the PGS Catalog.
The raw variant sets must undergo processing that can be completed by the scripts of the pgs_sets directory.  This processing splits the variant sets by chromosomes and matches all the alleles properly.
The variant sets and VCFs are combined together to actually form polygenic risk scores within the do_score directory.  The scoring is applied to specified variant sets and all VCFs and all UK Biobank individuals are scored together.
The polygenic risk scores are assessed with the scripts listed within the check_acc directory.  These scripts can eventually create a polygenic risk score report, which depicts the risk of an individual according to a specific polygenic risk score in a highly interpretable manner.
A final piece to be integrated within this process is the ancestry of each individual, as defined in the ancestry_match directory.

In theory, all of these scripts can be freshly applied to any whole genome sequenced VCF, although the UK Biobank data must be carefully integrated along with disease labels of those UK Biobank individuals.

