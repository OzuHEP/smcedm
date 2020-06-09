#!/bin/bash

if [ -f output.txt ]; then
	rm output.txt
fi

root -l "get_number_of_events.C(\"$1\", \"Delphes\")"

n_events=$(<output.txt)
interval=$((n_events/72))

echo "$n_events, $interval" 

for (( i=0; i<=$n_events; i+=$interval ));
do  
	#root -l "SMCEDM_Analysis.C($i, $((i+38665)))" &
	if [ $((i+interval)) -lt $n_events ]; then
		./SMCEDM_Analysis $1 $i $((i+interval)) &
	else
		./SMCEDM_Analysis $1 $i $((i+n_events)) &
	fi
done
wait
cd /mnt/harddisk4/scratch/
hadd -f $2_delphes.root output*.root
rm output*.root
#cd -
root -l "get_xsec.C(\"$1\", \"Delphes\")"