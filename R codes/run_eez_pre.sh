#!/bin/bash
models="model-mean"
for mod in $models;
do
	echo ${mod}
	jobify  -c 1 -o eez_pre.o -e eez_pre.e "Rscript eez_pre_treatment.R ${mod}"&
done
