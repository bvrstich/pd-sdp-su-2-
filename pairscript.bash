for I in {4..200..4}
do
   ./spin_pd `echo $I/10 | bc -l`
done
