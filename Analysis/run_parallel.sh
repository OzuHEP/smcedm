#!/bin/bash
rm output.txt
root -l "get_number_of_events.C(\"$1\", \"Delphes\")"

n_events=$(<output.txt)
interval=$((n_events/72))

echo "$n_events, $interval" 

for (( i=0; i<=$n_events; i+=$interval ));
do  
	#root -l "SMCEDM_Analysis.C($i, $((i+38665)))" &
	if [ $((i+interval)) -lt $n_events ]; then
		root -l "SMCEDM_Analysis.C(\"$1\", $i, $((i+interval)))" &
		#echo "SMCEDM_Analysis.C($i, $((i+interval)))" &
	else
		root -l "SMCEDM_Analysis.C(\"$1\",$i, $((i+n_events)))" &
		#echo "SMCEDM_Analysis.C($i, $((i+n_events)))" &
	fi
done
wait
cd /mnt/harddisk4/scratch/
hadd -f $2_delphes.root output*.root
rm output*.root
cd -
