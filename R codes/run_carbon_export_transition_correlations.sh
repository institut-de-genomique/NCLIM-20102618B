#!/bin/bash
jobify -c 1 -e run_carbon_export_transition_correlations.e -q normal --time=06:00:00 "Rscript carbon_export_transition_correlations.R"&
