#!/bin/bash

filename=${1%.*}
mv ${filename}.txt ${filename}.orig.txt

echo $(wc -l ${filename}.orig.txt | awk '{print $1}') $(head -1 ${filename}.orig.txt | awk '{print NF}') > ${filename}.txt
#sed  -e 's/..$//' ${filename}.orig.txt >> ${filename}.txt

while read -r line; do
	echo ${line%	*} >> ${filename}.txt
done < ${filename}.orig.txt
