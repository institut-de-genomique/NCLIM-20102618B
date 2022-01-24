#!/bin/bash
jobify -q normal -c 1 -o run_arules_carbonexport_$1_$2.o -e run_arules_carbonexport_$1_$2.e --time=24:00:00 "Rscript arules_carbonexport.R $1 $2 $3"&
