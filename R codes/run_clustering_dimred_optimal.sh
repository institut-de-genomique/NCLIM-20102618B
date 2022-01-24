#!/bin/bash
jobify -c 1 -o clustering_dimred_optimal_$1_$2_$4.o -e clustering_dimred_optimal_$1_$2_$4.e --time=06:00:00 -q normal "Rscript clustering_dimred_optimal.R $1 $2 $3 $4"&
