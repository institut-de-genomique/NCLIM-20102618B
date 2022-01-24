#!/bin/bash
jobify -c 1 -e run_carbon_export_assemblages_bilan_maps.e -o run_carbon_export_assemblages_bilan_maps.o -q normal --time=06:00:00 "Rscript carbon_export_assemblages_bilan_maps.R"&
