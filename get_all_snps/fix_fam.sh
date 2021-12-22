cd output_files

ls *fam | while read line;do

  id=`echo $line | cut -f1 -d'.'`
  echo "new_${id} new_${id} 0 0 0 -9" > $line

done

cd ..
