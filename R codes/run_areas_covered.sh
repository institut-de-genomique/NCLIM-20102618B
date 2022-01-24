#!/bin/bash
jobify -q normal -c 1 -o areas_covered.o -e areas_covered.e --time=00:10:00 "Rscript areas_covered.R"&
