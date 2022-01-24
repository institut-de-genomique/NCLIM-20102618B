#!/bin/bash
jobify -c 1  -e dimred_fig_${1}_${2}.e --time=00:30:00 -q normal "Rscript dimred_woa_figs.R $1 $2 $3"&
