This file describes the pipeline of the ecological niches model study "Restructuring of genomic provinces of surface ocean 1 plankton under climate change" based on Biogeographies from Richter et al.
Written by Paul Frémont
Before running all commands (./*sh) ; do module load r/3.3.1 (terminal etna0.genoscope.cns.fr); some scripts are run on rstudio: module load rstudio/0.99.115.3
R packages required:
gbm; randomForest; mdcv; nnet; dismo; FactoMineR; factoextra; readxl; ggplot2;matlab; reshape2; gplots; plotly; stringr; caret; mapproj; mapplots; SDMTools; RColoBrewer; ncdf4; CDFt; plotrix; png; grid; DALEX; ggalluvial; stringr; isofor; parallel; scales; Rtsne; sm; scatterplot3d; imager; ingredients; VennDiagramm; tidygraph; ggraph; igraph; animation; VennDiagram

The pipeline below allows to perform all analysis from the study. All data needed is in the directory so any of the program can be launched in any order to generate the figures. Some side analysis are also available in the folder other_analysis/ (to be moved in the main directory if performed)

Start directly at point 2, the climate data is pre-generated in the data folder.
1. Extract climatic data: WOA and ESM models
	1.1 WOA (https://www.nodc.noaa.gov/OC5/woa13/woa13data.html)
	done on ciclad.ipsl.jussieu.fr in the directory /data/pfremont/data/woa_month
	for each variable in *variable*: aou, no3, o2, o2s, po4, Si, so, T
	 1.1.0 for file in *.nc
		do ncks --mk-rec-dmn time $file $file
	       done
	 1.1.1 concatenate files
		ncrcat -h *.nc output.nc
	 1.1.2 Rearrange longitude
		cdo -sellonlatbox 0,360,-90,90 i.nc o.nc
	 1.1.3 select surface variable
		ncea -d depth,0,0 output.nc output1.nc
	 1.1.4 select variable of interest 
		cdo -selname,variable output1.nc output2.nc
	 1.1.5 change variable name dimension
		ncrename -h -O -d (dimension) -a (attribute) -v (variable) old-variable, new-var file
	 (see extract.sh and launch)
	1.2 ESM
	 do ./new_model-mean.sh
 
2. Niche model generation 
  Compute and optimize niche models using 4 machine learning techniques: RandomForest (rf), Neural Networks (nn), Boosted regression Trees (bt) and Generalized Additive Models (gam)
	2.1. Clean the biogeography dataset: clean_dataset_genocenoses.R
		input: Genocenoses_env_parameters_all_woa.xlsx
		output: Genocenoses_env_parameters_woa.rds
		command: execute the R program manually in rstudio (module load rstudio/0.99.115.3 then rstudio)
		It removes Genocenoses with presences < 4 and 3 lagoon niches: 180-2000 4 and 6; 20-180 4
	2.2. Optimize parameters of machine learning models
		input:Genocenoses_env_parameters_woa.rds
		output: best_models_*model*.txt *model*: rf, nn, bt, gam
		command: ./run_models_optimization.sh
		It finds the best combination of parameters for each of the four models (highest mean AUC on 30 cross validation on 75/25 partitionings of the dataset)
	2.3 Models validation 
		input:best_models_*model*.txt
		output: models_cross_validation.xlsx, best_selected_models.xlsx
		command: manual using excel
		A model is considered valid if 3 out of four models have a mean AUC > 0.65
	2.4 Save models
		input:  best_selected_models.xlsx; Genocenoses_env_parameters_woa.rds
		output: models_*model*.rds; good_models.rds 
		command: ./run_save_models.sh
		It saves the valid models; This time models are trained on the whole dataset
	2.5 Model visualization 
		input: models_*model*.rds; 
		output: relative_influences_*fraction*.txt/.pdf *fraction*: 180-2000; 20-180; 5-20 (or 43952); 0.8-5; 0.22-3; 0-0.2
			niche_response_*fraction*.pdf 
	2.6 Model comparison
		input: models_*model*.rds; Genocenoses_env_parameters_woa.rds
		output: correlation_sdm_models.pdf
		Compare outputs of machine learning models (correlation)
  2.7 Relative influences using DALEX R package
    !!! needs to be with r/3.6.0 (DALEX pacjage not available for r/3.3.1)
    input: *model*: bt, nn, rf, gam
            models_*model*.rds
    output: Relative_influences_violin_dalex_*fraction*.pdf
            relative_influence_violin_sd_dalex.pdf
            relative_influence_violin_sdfrac_dalex.pdf
            data_glob_relinf.txt
            data_relinf_dalex.txt
    command: ./run_dalex.sh
   2.8 Analysis of the cross validation
     input:models_cross_validation.xlsx
     output: bar_plot_auc_cross_val.pdf
     command: run manually cross_val_analysis.R

3. Bias correction using CDFt
  Use a bias correction method (Michelangeli et al 2009) to correct physico-chemical variables values from ESM models with respect to WOA data (based on Cumulative Distributions of the variables)  
	3.1 Perform the bias correction and save the data
		input: t_so_si_no3_po4_woa13_annual_mean_2005_2012.nc; no3_mon_woa-13-an_2006_201512_clim.nc; FED_PISCES2_1y_FEMIP.nc
		       all_7_mon_model-mean_rcp85_200601_201512_clim.nc; all_7_mon_model-mean_rcp85_209001_209912_clim.nc
		output: data_w.rds: world ocean atlas datas + dfe from PISCES2
			model-mean_2006-15.rds: model mean data in rds format present (no correction)
			model-mean_2090-99.rds: model mean data in rds format future (no correction)
			data_p.rds: corrected model mean data in rds format present
			data_f.rds: corrected model mean data in rds format future
			data_p_n.rds: corrected model mean data in rds format present with negative data
                        data_f_n.rds: corrected model mean data in rds format future with negative data
		command: ./run_save_climatic_data.sh

4. Make global scale projections
	4.1 Perform the projections on WOA, present, future (corrected and not corrected)
		input: climate data:data_w.rds, data_p.rds, data_f.rds, model-mean_2006-15.rds, model-mean_2090-99.rds
		       machine learning models: models_*model*.rds
		output:predictions_new_stations_w_model-mean_1deg.rds %% projections on WOA
		       predictions_new_stations_model-mean_1deg.rds  %% projections with bias correction
		       new_stations_1deg.rds %% latitude and longitude in the bias corrected model
		       predictions_new_stations_model-mean_1deg_nc.rds %% projections without bias correction
		       new_stations_1deg_nc.rds %% latitude and longitude in the non corrected model
		       predictions_new_stations_noT_model-mean_1deg.rds %% projections without temperature change
		       *time*: woa, 2006, 2090
		       *err*: delta, sd
		       *model*: rf, nn, bt, gam
		       model-mean_bray_curtis.rds %% bray-curtis index (2090-2006) using all 30 niches (not only dominants)
		       model-mean_pred_*time*_list.rds %% projections on WOA, 2006 and 2090
		       model-mean_pred_2090_list_noT.rds %% projections without temperature change
		       model-mean_pred_*time*_list_*err*_sdm.rds %% measure of uncertainty/disagreement between machine learning techniques (delta=max-min or sd)
		       model-mean_pred_*time*_list_*model*.rds %% individual projection of each machine learning techniques
		       model-mean_drivers_al.rds %% file that estimate the importance of each variable in community changes
		       model-mean_weights_cos.rds %% file of cos(laitude) to weight the area of each grid points
		command:./run_generate_projections_data.sh

5. Projection analysis
	5.1 Analyze each niches separately
		input: all predictions/projections files and models files generated precedently
      		output: 180-2000_5.pdf; 180-2000_8.pdf (examples of presence/absence data
			*fraction*: 180-2000; 20-180; 5-20 (or 43952); 0.8-5; 0.22-3; 0-0.2
			*num*: number of the genocenose (see best_models.xlsx)
			*time0*: 2006-15; 2006-15_nc; 2090-99; 2090-99_nc; woa; delta_2090-2006; shifts_delta_2090-2006
			*err*: sd; delta 
			*time* 2006;2090;woa
		        projection_model-mean_*fraction*_metacommunity_*num*_*time0*.pdf
			projection_model-mean_*fraction*_metacommunity_*err*_sdm_model_*time0*.pdf
			model-mean_*fraction*_metacommunity_*num*_delta_sdm_distirbution_*time*.pdf 
			model-mean_delta_sdm_model_distribution_*time*.pdf 
			Distances_shifts_model-mean.txt
			correlation_sdm_models_model-mean_*time*.pdf
		command:./run_individual_niches.sh
	
	5.2 Analyze all niches together
     Take all communities (=niches=) together and compute bray curtis dissimilarity index between present day and end of the century projection (as described in Barton et al. 2015); Also dominant communities map are generated for each size fraction
		input: all predictions/projections files and models files generated precedently; Distances_shifts_model-mean.txt; data_fisheries.rds; model-mean_eez_points1.rds; model-mean_eez_points.rds (fisheries and eez files genrated using ./run_fisheries_pre.sh and ./run_eez_pre.sh)
		output: PDFs:
			model-mean_fish_vs_bc.pdf; model-mean_fish_vs_bc_violin.pdf
			projection_model-mean_metacommunities_bray-curtis_2006-2090.pdf
			projection_model-mean_metacommunities_bray-curtis_high_fish_2006-2090.pdf
			projection_model-mean_metacommunities_bray-curtis_eez_2006-2090.pdf
			*fraction*: 180-2000; 20-180; 5-20 (or 43952); 0.8-5; 0.22-3; 0-0.2
			*num*: 1, 2, contour
			projection_model-mean_drivers_max_alone_*num*.pdf
			projection_model-mean_drivers_max_alone.pdf
			projection_model-mean_drivers_second_max_alone_*num*.pdf
                        projection_model-mean_drivers_second_max_alone.pdf
			*num0*: ''; _contour; _contour_alpha
			*time*: woa, 2006-15; 2090-99; 2090-99_noT
			projection_model-mean_*fraction*_metacommunities_*time**num0*.pdf
			delta_dom_*fraction*.pdf
			delta_dom_mean*fraction*.pdf
			correlation_drivers_max.pdf
			shifts_map.pdf
			number_of_changes.pdf
			R files (RDS):
			model-mean_bray_curtis_cond_fish_eez.rds
			number_of_changes.rds
			*time0*: 2006; 2090; woa; 2090_noT
			dominant_communities_*time0*.rds
		command:./run_global_analysis.sh

	5.3 Drivers analysis
    This analysis calculate the relative impact of each physico chemical drivers on the changes of niche between 2006-15 and end of the century 2090-99 projected environmental variables
		input: all predictions/projections files and models files generated precedently
		output: PDF
      *fraction*: 180-2000; 20-180; 5-20 (or 43952); 0.8-5; 0.22-3; 0-0.2
      *num*: number of the genocenose (see best_models.xlsx)
      *version*: alpha_dom, dom_contour, alpha_dom_contour, alpha_dom_noT, dom_contour_noT, alpha_dom_contour_noT
      projection_model-mean_metacommunities_bray-curtis_high_fish_dom_2006-2090.pdf
      *fraction*_*num*_drivers_*version*.pdf
      *version0*: '', '_alpha', '_contour', '_contour_alpha', '_contour_bcfilter', 'contour_bcfilter_noT'
      projection_model-mean_drivers_max_alone_dom*version*.pdf
      *var*: T, Si, NO3, Phos, SI_NO3, Fe, Sal
      taylor_diagram*var*.pdf
      taylor_diagram_all.pdf
      *version1*: 2006_woa, 2090, 2006, 2006_cor, 2090_cor (cor= bias corrected), 2090_2006 (delta between 2090 and 2006), 2006_neg, 2090_neg (neg= when bias correction gives negative values)
      model-mean_*var*_*version1*.pdf
      R files (RDS)
      mapping_lo_lt.rds, rev_mapping_lo_lt.rds'
      pred_dom06.rds; pred_dom90.rds; pred_dom90_noT.rds (dom= projection of only dominant community in a size fraction)
      model-mean_bray-curtis_dom.rds; model-mean_bray-curtis_dom_noT.rds (bray curtis on solely dominant communities and without temperature change in 2090 (noT)
      drivers_glob.rds; drivers_solo_niche.rds; drivers_solo_niche_glob.rds; drivers_glob_noT.rds; drivers_solo_niche_noT.rds; drivers_solo_niche_glob_noT.rds;
      selections.rds, selections_noT.rds
    command: ./run_drivers_analysis.sh
    
  5.4 Post analysis (just before run ./rel_inf_files.sh and ./covered_areas_files.sh)
    Analysis of covered oceanic areas by the different ecological niches and their displacement with climate change: centroid migrations
    input: shifts_dominant.rds; Projections/covered_areas_all.txt, Projections/relative_influences_all.txt, Genocenoses_env_parameters_woa.rds, 
    model-mean_full_area.txt, dominant_communities_2006.rds, dominant_communities_2090.rds, model-mean_weights_cos.rds
    outputs: PDF
      shifts_distibutions.pdf
      Relative_influences.pdf
      Relative_influences_*fraction*.pdf
      distribution_shifts.pdf
      relative_influence_violin.pdf
      Relative_influences_violin_*fraction*.pdf
      areas_*fraction*.pdf
      *loc*: North_A, South_A, North_P, South_P, Indian
      shifts_*loc*_fraction.pdf
      TXT files
      changes_trop_com.txt
      changes_temp_com.txt
      shifts_stats_fractions.txt
      shifts_stats_pole_mean_med.txt
      global_stats_area.txt
    command: ./run_post_analysis.sh
             ./run_alluvial_plots.sh
    
  5.5 PHATE analysis
    5.5.1. Use PHATE algorithm to visualize projection data in a reduced 3 dimension space.
    3 cases:
     - all projected probabilities of all niches are taken: *type* = prob
     - solely dominant probabilities (one per size fraction at a given oceanic point) of all niches are taken: *type* = dom
     - solely dominant probabilities are set to one: *type*=discrete 
    input: pred_dom06.rds; pred_dom90.rds; new_stations_1deg.rds
    *type0*: prob dom discrete; prob is used in the article
    *type*: prob dom discrete 06_prob 90_prob all; prob is used in the article
    output: R files (RDS)
      WOA:phate_fit_woa_*type0*.rds
      2006-90:phate_fit_0690_*type*.rds if *type* in(prob dom discrete)
              phate_fit_0690woa_*type*.rds  if *type*=all
              phate_fit_*type*.rds if *type* in (06_prob 90_prob)
      PDF
       rgb_phate_1deg_*type*.pdf
       3D_phate_*type*.pdf 
       rgb_2006_phate_1deg_*type*.pdf if *type* in (prob dom discrete all)
       rgb_2090_phate_1deg_*type*.pdf if *type* in (prob dom discrete all)
       rgb_woa_phate_1deg_*type*.pdf if *type*=all
       rgb_phate_1deg_*type*.pdf if *type* in (06_prob 90_prob)
    command:WOA: ./run_phate_woa.sh *type* 
            2006-90:./run_phate.sh *type*
    5.5.2 Clustering PHATE using k-medoids
      *period*: 06 90 woa 0690 0690woa
      *type*:prob dom discrete 06_prob 90_prob all
      *N*: optimal number of cluster found in the first step cluster_silhouette_*type*_phate_*period*.pdf
      *year*: 06 90
      *biog*: McKingley Reygondeau_BGCP19 Reygondeau_BGCP Reygondeau_BIOME
      input: phate_fit_*period*_*type*.rds
      output: cluster_silhouette_*type*_phate_*period*.pdf
              phate_cluster_medoids_*type*_*N*_*period*.rds
              rgb_2006_phate_1deg_*type*_*N*_clusts.pdf 
              rgb_2090_phate_1deg_*type*_*N*_clusts.pdf
              rgb_2090-2006_phate_1deg_*type*_*N*_clusts.pdf
              rgb_2006_phate_1deg_*type*_*N*_*year*_clusts.pdf
              rgb_2090_phate_1deg_*type*_*N*_*year*_clusts.pdf
              rgb_woa_phate_1deg_*type*_phate_*N*.pdf
              rgb_woa_phate_1deg_*type*_phate_*N*_clusts.pdf
              rgb_woa_phate_1deg_*type*_phate_*N*_*biog*.pdf
      command: 1 ./run_dimred_cluster.sh *type* phate *period*
               2 ./run_clustering_dimred_optimal.sh *type* phate *N* *period*
               3 ./run_dimred_figs_0690.sh *type* phate *N*
               4 ./run_dimred_figs_0690_sep.sh *type* phate *N* *year*
               5 ./run_dimred_woa_figs.sh *type* phate *N*

  5.6 Network analysis
    Consider each grid point of the ocean as an assemblage of 6 dominant communities of plankton; one of each size classes. Changes of networks and novelty as well as disparition of networks are analysed betweend 2006 and 2090 period
    input: dominant_communities_2006.rds; dominant_communities_2090.rds; model-mean_bray_curtis.rds; number_of_changes.rds; rev_mapping_lo_lt.rds; model-mean_bray-curtis_dom.rds; model-mean_weights_cos.rds
    output: PDF
          changes_network.pdf
          barplot_network.pdf
          venn_diagram_network.pdf
          common_change_fractions_network.pdf
          *co*: 75 25 50 
          major_changes_all_frac_pie_*co*.pdf
          TXT files
          niches_annotation.txt
          stats_area_network.txt
          stats_2006_vs_2090_network.txt
          glob_stats_changes_frac_network.txt
    command: ./run_network_analysis.sh
  
  5.7 Diverse analysis
	output: PDF
	    *frac*: '', _180-2000, _20-180, _5-20, _0.8-5, _0.22-3, _0-0.2
	    shannon06*frac*.pdf; shannon90*frac.pdf, shannon_distribution*frac*.pdf
	    *set*: set1, set2, set3
	    projection_model-mean_metacommunities_bray-curtis_dom_2006-2090_*set*.pdf
            drivers_solo_niches_area_upwelling.txt 
            areas_reorganization_strengths.pdf
	command:./run_shannon.sh, ./run_bray_curtis.sh drivers_solo_niche_statistics.R diagram_reorganization_strengths.R
      
   5.8 Biogeography comparison (before run ./run_biogeography_preprocess.sh)
      *biog*: McKingley Reygondeau_BGCP19 Reygondeau_BGCP Reygondeau_BIOME
      *type*: prob dom discrete 
      input: biogeographies
      output: randIndex_biogeographies.pdf
             randIndex_biogeographies_all.pdf
             biogeography_comparisons_*biog*_*type*.txt
             
      command: ./run_biogeography_comparisons.sh *biog* *type*
               ./run_biogeography_comparisons_1.sh
   5.9 Compositional shifts analysis
      output: PDF
              *frac0*: 180-2000, 20-180, 5-20, 0.8-5, 0.22-3, 0-0.2
              dominant_community_changes_trans_*frac0*.pdf
              *frac*: 180-2000, 20-180, 5-20, 0.8-5, 0.22-3
              *frac1*: 180-2000, 20-180, 5-20
              circos_transitions_diaz_*frac*.pdf 
              circos_transitions_hexanauplia_*frac1*.pdf
              circos_transitions_phot_*frac*.pdf
              circos_transitions_carb_*frac*.pdf
      command:in this order: ./run_areas_covered.sh
              in rstudio run save_SMAGs_stat.R
              ./run_SMAG-provinces_comparisons.sh SMAGs new
	      ./run_SMAG-provinces_comparisons.sh MAGprok new
              ./run_SMAG-provinces_comparisons_bis.sh SMAGs new
              ./run_SMAG-provinces_comparisons_bis.sh MAGprok new
              in rstudio run circos_plots.R
  6.0 Carbon export analysis
     6.1 in rstudio: run carbon_export_preprocess.R
     6.2 output: PDF
             *name*: diaz(_bc), hex(_bc), photo_algae(_bc), photo_cyano(_bc), photo_diatoms(_bc)
	     transitions_*name*.pdf 
         command: run_preprocess_transitions_carbon_export.sh
     6.3 Flux and delta flux maps and bilan from assemblages mean/median export based on different datasets 
     output: PDF
       *flux*: flux, flux_present, delta_flux
       *stat*: mean(_bc), median(_bc)
       *organism*: cyanobacteria(_big), copepods,algae(_big), diatoms(_big), cyano-diazotrophs (for flux==delta_flux)
       *version*: all, no_5-20
       *type*: '', guidi, henson, extrapolated_henson, extrapolated_eppley, extrapolated_laws, extrapolated_schlitzer
       carbon_*ty*_assemblages_*stat*_*organism*_*version*_*type*.pdf
       txt
       *type*: '', guidi, henson, extrapolated_henson, extrapolated_eppley, extrapolated_laws, extrapolated_schlitzer
       carbon_export_statistics_*type*.txt
     command:./run_carbon_export_assemblages_bilan_maps.sh
     6.4 Association rules between carbon export and organisms' changes using Apriori algorithm
     output: PDF
     *region*:'equatorial', 'subtropical_north','subtropical_south', 'subpolar'
     *maxn*: 6, 7 (maximum length of the rule)
     *support*: 0.05, 0.02 (minimum support of the rule) 
     associations_carbon_export_*maxn*_*region*_*support*.pdf 
     command:./run_carbon_export_transition_correlations.sh (output necessary data table for arules)
             ./run_arules_carbonexport.sh $1 $2 $1:*maxn* $2:*support*

  7.0 tidy the directory
    7.1 move figures of the study in a specific folder: ./command_paper_figs_v3.sh
    7.2 move pdfs in folders: ./transfer_files.sh
