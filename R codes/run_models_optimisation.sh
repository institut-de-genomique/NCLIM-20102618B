#!/bin/bash
#models="gam nn bt rf"
#for mod in $models;
#do
#	jobify  -c 31 -p normal --time=24:00:00 -o models_opt_${mod}.o -e models_opt_${mod}.e "Rscript models_optimisation_${mod}.R models_opt_${mod}.out"&
#done
jobify  -c 16 -p normal --time=24:00:00 -o models_opt_nn.o  -e models_opt_nn.e  "Rscript models_optimisation_nn.R models_opt_nn.out"&
jobify  -c 31 -p normal --time=24:00:00 -o models_opt_bt.o  -e models_opt_bt.e  "Rscript models_optimisation_bt.R models_opt_bt.out"&
jobify  -c 6  -p normal --time=24:00:00 -o models_opt_gam.o -e models_opt_gam.e "Rscript models_optimisation_gam.R models_opt_gam.out"&
jobify  -c 16 -p normal --time=24:00:00 -o models_opt_rf.o  -e models_opt_rf.e  "Rscript models_optimisation_rf.R models_opt_rf.out"&
