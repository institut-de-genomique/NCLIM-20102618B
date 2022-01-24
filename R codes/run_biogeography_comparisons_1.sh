#!/bin/bash
jobify --time=12:00:00 -c 1 -o bio_comp_1.o -e bio_com_1.e  -q normal "Rscript biogeography_comparisons_1.R"&
