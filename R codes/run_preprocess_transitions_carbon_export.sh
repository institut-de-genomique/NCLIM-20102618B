#!/bin/bash
jobify -c 1 -e run_preprocess_transitions_carbon_export.e -q normal --time=06:00:00 "Rscript preprocess_transitions_carbon_export.R"&
