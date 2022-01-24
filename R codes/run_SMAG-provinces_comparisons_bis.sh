#!/bin/bash
jobify -c 1 -q normal -o SMAG-provinces_comparisons_bis_$1_$2.o -e SMAG-provinces_comparisons_bis_$1_$2.e  --time=06:00:00 "Rscript SMAG-provinces_comparisons_bis.R $1 $2"&
