Rscript plot_it.R

ls report_plots/prev.* | cut -f2 -d'/' | cut -f2 -d'.' | while read i;do
  Rscript make_report_full_line.R $i
done
