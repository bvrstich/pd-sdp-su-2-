for I in {8..4000..8}
do
   ./spin_pd `echo $I/100 | bc -l`
done
