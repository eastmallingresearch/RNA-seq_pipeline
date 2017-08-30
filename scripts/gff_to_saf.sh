#!/bin/bash

GFF=$1

grep exon $GFF|awk -F"\t" '{split($9,a,".");s=substr(a[1],4);print s,$1,$4,$5,"-" }' OFS="\t"