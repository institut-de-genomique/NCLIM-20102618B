#!/bin/bash
jobify -q normal -c 1 -o presence_absence_map.o -e presence_absence_maps.R.e --time=00:10:00 "Rscript presence_absence_maps.R"&
