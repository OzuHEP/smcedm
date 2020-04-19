#!/bin/bash
start=`date +%s`
for i in `seq 0 38665 2783939`; 
do  
	root -l "SMCEDM_Analysis.C($i, $((i+38665)))" &
	#echo "SMCEDM_Analysis.C($i, $((i+38665)))" &
done
wait
end=`date +%s`
runtime=$((end-start))

echo $runtime
