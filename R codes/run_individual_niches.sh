#!/bin/bash
jobify -c 1 -q normal -o indiv_niches.o -e indiv_niches.e  --time=06:00:00 "Rscript individual_niches.R"&
