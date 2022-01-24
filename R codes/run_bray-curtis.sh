#!/bin/bash
jobify -q normal -c 1 -o bray_curtis_map.o -e bray_curtis_map.e --time=00:10:00 "Rscript bray_curtis_map.R"&
