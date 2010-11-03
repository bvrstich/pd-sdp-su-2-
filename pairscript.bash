for I in {4..400..4}
do
   ./spin_pd `echo $I/100 | bc -l`
done
