for I in {0..80..4}
do
   nn=`echo $I/10 | bc -l`
   ./spin_pd $nn
done
