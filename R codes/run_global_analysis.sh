#!/bin/bash
jobify -c 1 -o glob_an.o -e glob_an.e --time=06:00:00 -q normal "Rscript global_analysis.R"&
