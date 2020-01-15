totalseq=$1 
for filename in *.raw 
do 
FFour=${filename:0:2} 
for spline in $(seq $totalseq)
do
  ./GoBaby $filename $spline 
 echo $filename $FFour cor=$(head -n 49 ${FFour}*OUT | tail -n 10 | cut -c 70-77 | head -n 5 | tail -n 1) spline=$spline >> $filename.txt
 # echo $filename $FFour spline=$spline cor=$(head -n 49 $FFour*.OUT | tail -n 10 | cut -c 70-77 | head -n 5 | tail -n 1) >> $filename.txt
 # head -n 49 $FFourCOF.OUT | tail -n 10 >> $filename.txt #can use to check actually pulling correct value 
done
  sort -nr $filename.txt > Sorted$filename.txt 
done
