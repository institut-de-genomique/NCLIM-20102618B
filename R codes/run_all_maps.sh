#!/bin/bash
./run_bray-curtis.sh
./run_individual_niches.sh
./run_global_analysis.sh
./run_drivers_analysis.sh
./run_network_analysis.sh
./run_dimred_woa_figs.sh prob phate 4
./run_dimred_woa_figs.sh prob phate 7
./run_dimred_figs_0690.sh prob phate 4
./run_pres_abs_map.sh
./run_post_analysis.sh
./run_areas_covered.sh
