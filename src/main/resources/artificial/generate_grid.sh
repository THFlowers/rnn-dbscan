#!/bin/bash

nodes=25
for (( i=1; i<=$nodes; i++ ))
do
	for (( j=1; j<=$nodes; j++ ))
	do
		echo -n $(( 5 * $i))'.0'
		printf '\t'
		echo -n $(( 5 * $j))'.0'
		printf '\t1\n'
	done
done

for (( i=1; i<=$nodes; i++ ))
do
	for (( j=1; j<=$nodes; j++ ))
	do
		echo -n $(( 127 + $i))'.5'
		printf '\t'
		echo -n $j'.0'
		printf '\t2\n'
	done
done
