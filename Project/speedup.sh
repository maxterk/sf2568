#!/bin/bash -l
OMPI_MCA_=0
export OMPI_MCA_btl_openib_warn_nonexistent_if=0
export OMPI_MCA_btl_openib_warn_nonexistent_if

#module add i-compilers intelmpi
for j in {1..3}
do
for i in {1..24}
do
	mpirun -n $i ./program >> output_seq_iter_$j.txt
done
done
#mpirun -n 4 ./program
#echo Hello world

