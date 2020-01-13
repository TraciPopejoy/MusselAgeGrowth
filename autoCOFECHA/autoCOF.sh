
totalseq=$1
for filename in *.raw
do
for spline in $(seq $totalseq)
do
  ./GoBaby $filename $spline
  echo $filename spline=$spline cor=$(head -n 49 *COF.OUT | tail -n 10 | cut -c 70-77 | head -n 5 | tail -n 1) >> $filename.txt
 # head -n 49 *COF.OUT | tail -n 10 >> $filename.txt #can use to check actually pulling correct value 
done 
done
