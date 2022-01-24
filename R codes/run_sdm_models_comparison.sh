#!/bin/bash
jobify  -c 1 --time=00:10:00 -q normal -o sdm_comp.o -e sdm_comp.e "Rscript sdm_models_comparison.R"&
