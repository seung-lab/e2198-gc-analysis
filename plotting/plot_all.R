source("plot.R")

#NOTE: REMEMBER TO READ plot_type_boxplot signature before plotting new things

#FIG1d&e
#d
fig1d <- plot_ca_trials(20233, roi_sums, ca_ids, 1)
#ggsave('fig1/fig1d_20233.pdf',plot=fig1d,width=8.5,height=11,units="in")
#talk version
fig1d <- plot_ca_trials(20233, roi_sums, ca_ids, 1) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank())
ggsave('fig1/fig1d_20233_talk.pdf',plot=fig1d,width=8.5,height=11,units="in", bg="transparent")

#e
fig1e <- plot_ca_trials_by_stim(20233, roi_sums, ca_ids, 1) #note temporary code within this fn
#as far as I know, this needs to be saved through the GUI
multiplot(fig1e, cols=8)
fig1e2 <- plot_avg_type_traces('37v',ca_ids, ca) #note temporary code within this fn too (lazy)
ggsave('fig1/fig1e2_20233.pdf',plot=fig1e2,width=8.5,height=11,units="in")

#EXTFIG5
#e5a
soma_plot <- plot_type_boxplot( all_gc_types, soma_sizes_for_ids, "Soma Size (um^3)" )
ggsave('extfig5/extfig5_soma.pdf',plot=soma_plot, width=7.89, height=8.53,units="in")
#e5b
max_diam_plot <- plot_type_boxplot( all_gc_types, diameters_for_ids, "Max Diameter (um)", 1/1000, NULL, ci, ci_ids )
ggsave('extfig5/extfig5_max_diam.pdf',plot=max_diam_plot, width=7.89, height=8.53,units="in")
#e5c
#mean_bw_plot <- plot_type_boxplot( all_gc_types, mean_branch_width_for_ids, "Mean Branch Width", NULL, "mean_bw.csv")
mean_bw_plot <- plot_type_file_boxplot( "mean_bw.csv", "Mean Branch Width" )
ggsave('extfig5/extfig5_mean_branch_width.pdf',plot=mean_bw_plot, width=7.89, height=8.53,units="in")
#e5d
pl_corr   <- plot_pl_corr()
ggsave('extfig5/extfig5_pl_corr.pdf',plot=pl_corr, width=7.89, height=4.25,units="in")
#e5e
mean_bw_corr  <- plot_soma_file_scatter( "mean_bw.csv", "Mean Branch Width" )
ggsave('extfig5/extfig5_mbw_corr.pdf',plot=mean_bw_corr, width=7.89, height=4.26,units="in")

#EXTFIG8
#SDI (a&b)
sac_di <- plot_type_boxplot( all_gc_types, dssi_for_ids, "SAC Directionality Index", NULL, NULL, gcsac, num_ids, denom_ids) + scale_y_continuous(limits=c(0,1), expand=c(0,0))
sdi_pvals <- type_dssi_significance_test(all_gc_types, gcsac, num_ids, denom_ids)
#63 and 73 use only ON SACs
on_sdi_pvals <- type_dssi_significance_test(all_gc_types, gcsac, num_ids, denom_ids, "on")
sdi_pvals['63'] <- on_sdi_pvals['63']
sdi_pval_plot <- plot_pvals(sdi_pvals)
ggsave('extfig8/extfig8_sac_di.pdf',plot=sac_di, width=7.89, height=8.53,units="in")
ggsave('extfig8/extfig8_sac_di_p.pdf',plot=sdi_pval_plot, width=7.89, height=8.53,units="in")
#DSI (c&d)
#on_dsi <- plot_type_boxplot( all_gc_types, on_dsi_for_ids, "ON DSI", NULL, dsos, dsos_ids )
#off_dsi <- plot_type_boxplot( all_gc_types, off_dsi_for_ids, "OFF DSI", NULL, dsos, dsos_ids )
overall_dsi <- plot_type_boxplot( all_gc_types, overall_dsi_for_ids, "Direction Selectivity Index", NULL, NULL, dsos, dsos_ids )
dsi_pvals <- type_dsi_significance_test( all_gc_types, dsos, dsos_ids )
dsi_pval_plot <- plot_pvals(dsi_pvals)
ggsave('extfig8/extfig8_dsi.pdf',plot=overall_dsi, width=7.89, height=8.53,units="in")
ggsave('extfig8/extfig8_dsi_p.pdf',plot=dsi_pval_plot, width=7.89, height=8.53,units="in")

#EXTFIG9
#Arbor Density (a)
npa_plot <- plot_type_file_boxplot( "npha.csv", "Arbor Density (mm-2)", 1e6 )
ggsave('extfig9/extfig9_arbor_density.pdf',plot=npa_plot, width=7.89, height=8.53,units="in")
#Arbor Complexity (b)
bpppl_plot <- plot_type_file_boxplot( "bpppl.csv", "Arbor Complexity (mm-1)", 1e6 )
ggsave('extfig9/extfig9_arbor_complexity.pdf',plot=bpppl_plot, width=7.89, height=8.53,units="in")
#Arbor Asymmetry (c&d)
asym_plot <- plot_type_boxplot( all_gc_types, asym_index_for_ids, "Arbor Asymmetry Index", NULL, NULL, asym, asym_ids, F ) + scale_y_continuous(limits=c(-1,1),expand=c(0,0))
asym_ps <- type_asymmetry_significance_test(all_gc_types, asym, asym_ids)
asym_p_plot <- plot_pvals(asym_ps)
ggsave('extfig9/extfig9_arbor_asym.pdf',plot=asym_plot, width=7.89, height=8.53,units="in")
ggsave('extfig9/extfig9_arbor_asym_p.pdf',plot=asym_p_plot, width=7.89, height=8.53,units="in")

#Supplementary
#Normalized SAC contact
on_sac_plot <- plot_type_boxplot( all_gc_types, on_contact_for_ids, "ON SAC Contact", NULL, NULL, gcsac, num_ids, denom_ids )
off_sac_plot <- plot_type_boxplot( all_gc_types, off_contact_for_ids, "OFF SAC Contact", NULL, NULL, gcsac, num_ids, denom_ids )
both_sac_plot <- plot_type_boxplot( all_gc_types, both_contact_for_ids, "Overall SAC Contact", NULL, NULL, gcsac, num_ids, denom_ids )
ggsave('supp/supp_overall_sac.pdf',plot=both_sac_plot, width=7.89, height=8.53,units="in")
ggsave('supp/supp_on_sac.pdf',plot=on_sac_plot, width=7.89, height=8.53,units="in")
ggsave('supp/supp_off_sac.pdf',plot=off_sac_plot, width=7.89, height=8.53,units="in")

path_length_plot <- plot_type_file_boxplot( "pl.csv", "Path Length (mm)", 1/1e6 )
path_length_plot <- plot_type_boxplot( c(all_gc_types,"cutoffs","weirdos","Path Length (mm)", 1/1e6 )
hull_area_plot <- plot_type_boxplot( all_gc_types, hull_areas_for_ids, "Hull Area (mm^2)", 1/1e12, NULL, ci, id_lookup)
ggsave('supp/supp_path_length.pdf',plot=path_length_plot, width=7.89, height=8.53,units="in")
ggsave('supp/supp_hull_area.pdf',plot=hull_area_plot, width=7.89, height=8.53,units="in")

on_sac_di <- plot_type_boxplot( all_gc_types, dssi_for_ids, "SAC Directionality Index", NULL, NULL, gcsac, num_ids, denom_ids, "on" ) + scale_y_continuous(limits=c(0,1), expand=c(0,0))
on_sdi_pval_plot <- plot_pvals(on_sdi_pvals)
ggsave('supp/supp_on_sac_di.pdf',plot=sac_di, width=7.89, height=8.53,units="in")
ggsave('supp/supp_on_sac_di_p.pdf',plot=sdi_pval_plot, width=7.89, height=8.53,units="in")

off_sac_di <- plot_type_boxplot( all_gc_types, dssi_for_ids, "SAC Directionality Index", NULL, NULL, gcsac, num_ids, denom_ids, "off" ) + scale_y_continuous(limits=c(0,1), expand=c(0,0))
off_sdi_pvals <- type_dssi_significance_test(all_gc_types, gcsac, num_ids, denom_ids, "off")
off_sdi_pval_plot <- plot_pvals(on_sdi_pvals)
ggsave('supp/supp_off_sac_di.pdf',plot=sac_di, width=7.89, height=8.53,units="in")
ggsave('supp/supp_off_sac_di_p.pdf',plot=sdi_pval_plot, width=7.89, height=8.53,units="in")



###############################
#POTENTIALLY OLD BEYOND THIS POINT
###############################

#dendritic thickness (currently flawed)
mbw_plot <- plot_type_file_boxplot( "mbw.csv", "Median Branch Width (units eventually)" )

#bp / path length
#bp / hull_area
bpphl_plot <- plot_type_boxplot( all_gc_types, bp_per_hull_area_for_ids, "Branch Points per Hull Area", "1/(nm^2)", ci, id_lookup )
#num nodes / hull_area
primary_dendrite <- plot_type_file_boxplot( "pdw.csv", "Primary Dendrite Width" )
pdw20 <- plot_type_boxplot( all_gc_types, primary_dendrite_width_for_ids, "Primary Dendrite Width", NULL, "pdw20", ci, id_lookup, 0.2)


#Asymmetry
asym_plot <- plot_type_boxplot( all_gc_types, asym_index_for_ids, "Arbor Asymmetry Index", NULL, NULL, asym, asym_ids, F )
asym_ps <- type_asymmetry_significance_test(all_gc_types, asym, asym_ids)
asym_p_plot <- plot_pvals(asym_ps)
#dv_asym_plot <- plot_type_boxplot( all_gc_types, dv_asym_for_ids, "Dorsoventral Asymmetry", NULL, ci, id_lookup )


#OSI
on_osi <- plot_type_boxplot( all_gc_types, on_osi_for_ids, "ON OSI", NULL, NULL, dsos, dsos_ids )
#off_osi <- plot_type_boxplot( all_gc_types, off_osi_for_ids, "OFF OSI", NULL, dsos, dsos_ids )
overall_osi <- plot_type_boxplot( all_gc_types, overall_osi_for_ids, "Overall OSI", NULL, dsos, dsos_ids )

#Polar Plots
#norm_on_ds <-    plot_polar_selectivity( all_gc_types, dsos, dsos_ids, T,T,T)
#norm_off_ds <-   plot_polar_selectivity( all_gc_types, dsos, dsos_ids, T,F,T)
#unnorm_on_ds <-  plot_polar_selectivity( all_gc_types, dsos, dsos_ids, T,T,F)
#unnorm_off_ds <- plot_polar_selectivity( all_gc_types, dsos, dsos_ids, T,F,F)

#norm_on_ds <-    plot_polar_selectivity( all_gc_types, dsos, dsos_ids, F,T,T)
#norm_off_ds <-   plot_polar_selectivity( all_gc_types, dsos, dsos_ids, F,F,T)
#unnorm_on_ds <-  plot_polar_selectivity( all_gc_types, dsos, dsos_ids, F,T,F)
#unnorm_off_ds <- plot_polar_selectivity( all_gc_types, dsos, dsos_ids, F,F,F)

#covariate_plots
#path_length_by_soma  <- plot_type_median_scatter( all_gc_types, soma_sizes_for_ids, path_lengths_for_ids, "Soma Size", "Path Length" )
#branch_width_by_soma <- plot_type_median_scatter( all_gc_types, soma_sizes_for_ids, median_branch_width_for_ids, "Median Soma Size", "Median Median Branch Width" )
#hull_area_by_soma    <- plot_type_median_scatter( all_gc_types, soma_sizes_for_ids, hull_areas_for_ids, "Soma Size", "Median Hull Area", ci, ci_ids )
#max_diam_by_soma     <- plot_type_median_scatter( all_gc_types, soma_sizes_for_ids, diameters_for_ids, "Soma Size", "Median Max Diameter", ci, ci_ids )



#Scatter Plots
mbw_corr  <- plot_soma_file_scatter( "mbw.csv", "Median Branch Width " )
pdw_corr <- plot_soma_file_scatter( "pdw.csv", "Primary Dendrite Width")
pdw20_corr <- plot_soma_file_scatter( "pdw20.csv", "Primary Dendrite Width(20)")
md_corr <- plot_soma_fn_scatter( diameters_for_ids, "Max Diameter (um^2)", 1/1000, ci, id_lookup )
hull_corr <- plot_soma_fn_scatter( hull_areas_for_ids, "Hull Area (um^2)", 1/1e6, ci, id_lookup )
