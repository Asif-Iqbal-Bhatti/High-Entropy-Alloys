#!/bin/bash

mkdir -p .test
i=0
hh=bestcorr
jj=out

for d in sqs_x_0.20_*; do
	f=${d%/}; 
	ext=${d##*-}
	
	for k in $f/$hh*.$jj; do
	
		echo -ne "$k! \033[0K\r"
		
		cc=`printf "%04d" $i`
		cp $k .test/$hh-$cc.$jj
		
		ov=`tail -1 .test/$hh-$cc.$jj | awk -F' ' '{printf $2}'`
		echo .test/$hh-$cc.$jj, $ov >> objective_Ti0.20.csv
		
		let i=i+1
	
	
	done
done	
