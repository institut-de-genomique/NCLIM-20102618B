#!/bin/bash
jobify -c 1 -o post_analysis.o -e post_analysis.e --time=06:0:00 "Rscript post_analysis.R"&
