#!/bin/bash
jobify -c 5 -o dimred_clust_${1}_${2}_${3}.o -e dimred_clust_${1}_${2}_${3}.e --time=72:00:00 -q xlarge "Rscript dimred_clustering.R $1 $2 $3"&
