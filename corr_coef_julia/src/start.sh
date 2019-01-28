#!/bin/bash


START=$1
STOP=$2
NPROCS=$3
FILE_OUT=$4

JULIA=`whereis julia | awk {'print $2;'}`

DELTA=$((($STOP-$START) / $NPROCS))

echo $DELTA

for i in $(seq 1 $3 ); do
		$JULIA -O3 Corr_Coef.jl $(($START+$DELTA*($i-1)))  $(($START+$DELTA*$i-1)) ../results/$4_$i.dat &
done







