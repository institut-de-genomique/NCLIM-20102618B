#!/bin/bash
jobify  -c 1 --time=06:00:00 -q normal -o proj_data.o -e proj_data.e "Rscript generate_model_mean_predictions.R"&
