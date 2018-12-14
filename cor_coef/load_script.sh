#!/bin/bash

read TIME


while true; do

		LOAD=`uptime | awk {'print $8'} | awk -F "," {'print $1'} | awk -F "." {'print $1'}`
		if (( "$LOAD" < "17" )); then
				echo $LOAD
				nohup ./corr_coef.x 10000000 $TIME >data/time_$TIME.dat &
				TIME=$[$TIME + 1]
		fi

		sleep 20

done


