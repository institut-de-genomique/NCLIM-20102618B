#!/bin/bash

rm -f relative_influences_all.txt
cat *relative_influences*.txt >> relative_influences_all.txt
mv *relative_influences*.txt Projections/

