#!/bin/bash
jobify -c 1 -o bio_preproc.o -e bio_preproc.e  -q small "Rscript biogeographies_preprocess.R"&
