
chr=$1

 ls ../get_all_snps/output_files/*bed | cut -f4 -d'/' | cut -f1 -d'.' | sort | uniq | while read newfix;do
    for subset in variable ref;do

      cat ../get_all_snps/output_files/${newfix}.${chr}.${subset}.bim | cut -f2 > subset_rsids/exrsid.${newfix}.${subset}

    done
 done
