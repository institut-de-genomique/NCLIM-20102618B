#!/bin/bash
jobify -c 1 -o phate_$1.o -e phate_$1.e --time=48:00:00 -q xlarge "Rscript phate_fit.R $1"&
