#!/bin/bash

mod="model-mean"
fractions="0.8-5 180-2000 20-180 5-20 0.22-3 0-0.2"
fractions_bis="0.8-5 180-2000 20-180 43952 0.22-3 0-0.2"

for frac in $fractions;
do
	mv *${mod}*${frac}* Projections/${frac}/
done

for frac in $fractions;
do
	mv ${frac}*drivers* Drivers/${frac}/
done

mv *fish*pdf Fish_EEZ/
mv *eez*pdf Fish_EEZ/
mv *${mod}*drivers*pdf Drivers/
mv *${mod}*metacommunities_* Projections/
mv *sdm*pdf Correlations_sdm_models/	
mv *${mod}*_2090.pdf corrections/
mv *${mod}*_2090_neg.pdf corrections/
mv *${mod}*_2006_neg.pdf corrections/
mv *${mod}*_cor.pdf corrections/
mv *${mod}*_woa.pdf corrections/
mv *${mod}*_2090_2006.pdf corrections/
mv taylor*.pdf Taylor_diagram/
mv *shifts*.pdf Shifts/
mv Relative_influences*pdf Relative-influences/
mv relative_influence*pdf Relative-influences/
mv areas*pdf Areas/
mv *network*pdf Network_analysis/
mv number_of_changes.pdf Network_analysis/
mv *tsne*pdf tsne/
mv shannon*pdf shannon/
mv *alluvial*pdf Areas/
mv *phate*pdf phate/
mv *carbon*pdf Carbon_export/

for frac in $fractions_bis;
do
        if [ "$frac" == "43952" ];
        then
                mv relative_influences_${frac}* Climate_niches/5-20/
                mv niche_response_${frac}* Climate_niches/5-20/
        else
                mv relative_influences_${frac}.* Climate_niches/${frac}/
                mv niche_response_${frac}* Climate_niches/${frac}/
        fi
done

mv *pdf diverse/
