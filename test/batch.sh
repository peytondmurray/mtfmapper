#!/bin/bash

if [ $# -ne 1 ]; then
    echo "you must specify the desired blur sigma"
    exit
fi

:>results.txt
:>diff_results.txt
for i in 4 10 30; do
    echo $i
    for j in {1..30}; do
    	./generate_rectangle -b $1 -s $j -a $i > gen.txt
    	target=`awk '$1 ~ /MTF50/ {printf "%s ", $3}' gen.txt`
    	awk '$1 ~ /MTF50/ {printf "%s ", $3}' gen.txt >> results.txt
    	./mtf_mapper rect.png  2> res.txt
    	echo -n $i " " >> results.txt
    	cat res.txt >> results.txt
    	awk '$1 ~ /MTF50/ {printf "%s ", $3}' gen.txt >> diff_results.txt
    	echo -n  $i " " >> diff_results.txt
    	awk -vt=$target '{ sum=0; for (i=1; i <= NF; i++) sum += sqrt(($i - t)*($i - t)); printf "%lf\n", sum/NF}' res.txt >> diff_results.txt
    done
done
