#!/bin/bash
jobify -c 1 -o phate_fit_woa_${1}.o -e phate_fit_woa_${1}.e --time=06:00:00 "Rscript phate_fit_woa.R $1"&
