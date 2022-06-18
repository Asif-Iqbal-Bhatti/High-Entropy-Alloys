#!/bin/bash

mkdir ./test
i=0

for d in */; 
do
	f=${d%/}; ext=${d##*-}
	for k in $f/mcsqs*.log; 
	do
  
	echo $k
  cc=`printf "%03d" $i`
	cp $k .test/mcsqs-$cc.log
  let i=i+1
	
	done
done	
