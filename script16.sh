#!/bin/bash
# use the first two arguements for the first configuration and the last configuration the third argument defines the number of runs
# obviously, $1 < $2!
# $1 FIRST = START OF CONFIGURATIONS
# $2 SECTION = END OF CONFIGURATIONS
# $3 THRID = NUMBER OF RUNS
maxConc=16 # maximum of 15 concurrent processes
configs="$2" # number of configurations to be processed
RUNS="$3" # number of runs per configuration
currentCount="$1" #running numebr of for-loop
totalCount="$1"
sleepdur=12

echo "BASH: $1 configurations perform $2 runs"


## start 'maxConc' processes, for each increase count

## while count < configs and wait, if  maxConc --> wait, 
while [ $totalCount -lt $configs ]; do
	currentCount=$(jobs -rp | wc -l) ## get current count of jobs
	if [ $currentCount -lt $maxConc ]; then
		echo "BASH: totalcount: $totalCount runs: $RUNS currentCout $currentCount maxConc $maxConc"
		./distr_5 "$totalCount" "$RUNS" &
		((totalCount+=1))
	else
		sleep 3 ## WAIT for certain time
	fi
done

