# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 16:35:07 2024

@author: sf02
"""
import sys
import glob
import numpy as np
from matplotlib import pyplot as plt
from nilearn import surface, plotting
from scipy import stats

analyses = ['log-LASSO','SOSLASSO','linear-LASSO_dimension-1','linear-LASSO_dimension-2','linear-LASSO_dimension-3','grOWL_dimension-1','grOWL_dimension-2','grOWL_dimension-3']

# set path
dirp='/imaging/projects/cbu/wbic-p00591-DAISY/main/'
sys.path.append(dirp)

# get pial surface mesh and curvature data
sys.path.append(dirp+'/scripts/fsaverage')
pial_left = dirp+'/scripts/fsaverage/pial_left.gii.gz'
curv_left = dirp+'/scripts/fsaverage/curv_left.gii.gz'
pial_right = dirp+'/scripts/fsaverage/pial_right.gii.gz'
curv_right = dirp+'/scripts/fsaverage/curv_right.gii.gz'

subcode = ['sub-001','sub-002','sub-003','sub-004','sub-007','sub-009','sub-010','sub-011','sub-012','sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-022','sub-023','sub-024','sub-026','sub-028','sub-029','sub-030','sub-031','sub-032']
   
# for each analysis
for current_analysis in analyses:
    # initialise
    Ltotals = np.zeros(163842)
    Lanimatecoefs = np.zeros(163842)
    Rtotals = np.zeros(163842)
    Ranimatecoefs = np.zeros(163842)
    Lpermtotals = np.zeros(163842)
    Lpermanimatecoefs = np.zeros(163842)
    Rpermtotals = np.zeros(163842)
    Rpermanimatecoefs = np.zeros(163842)
    # for each participant
    for s in subcode:
        # load final coefficients
        Lfc = surface.load_surf_data(dirp+'/work/'+s+'/coefficients/surface/L_'+s+'_'+current_analysis+'_final_coefficients.func.gii')
        Rfc = surface.load_surf_data(dirp+'/work/'+s+'/coefficients/surface/R_'+s+'_'+current_analysis+'_final_coefficients.func.gii')
        # update counts. Create a mesh with 1 where there is a coefficient and 0 otherwise
        tmp = Lfc != 0
        tmp = tmp.astype(int)
        Ltotals += tmp
        tmp = Rfc != 0
        tmp = tmp.astype(int)
        Rtotals += tmp
        # update sum of NEGATIVE coefficients (N.B. Chris and Tim code animal as 1 and therefore count positive coefficients)
        tmp = Lfc < 0
        tmp = tmp.astype(int)
        Lanimatecoefs += tmp
        tmp = Rfc < 0
        tmp = tmp.astype(int)
        Ranimatecoefs += tmp
        # get lists of perms (they are not in order, but it doesn't matter)
        Lperms = glob.glob(dirp+'/work/'+s+'/coefficients/surface/L_'+s+'_'+current_analysis+'_randomseed-*_perm_coefficients.func.gii')
        Rperms = glob.glob(dirp+'/work/'+s+'/coefficients/surface/R_'+s+'_'+current_analysis+'_randomseed-*_perm_coefficients.func.gii')
        # for each perm
        for p in range(len(Lperms)):
            # load coefficients
            Lpc = surface.load_surf_data(Lperms[p])
            Rpc = surface.load_surf_data(Rperms[p])
            # update counts
            tmp = Lpc != 0
            tmp = tmp.astype(int)
            Lpermtotals += tmp
            tmp = Rpc != 0
            tmp = tmp.astype(int)
            Rpermtotals += tmp
            # update sum of negative coefficients
            tmp = Lpc < 0
            tmp = tmp.astype(int)
            Lpermanimatecoefs += tmp
            tmp = Rpc < 0
            tmp = tmp.astype(int)
            Rpermanimatecoefs += tmp
    # set unselected vertices to nan
    Ltotals[Ltotals==0] = np.nan
    Rtotals[Rtotals==0] = np.nan
    Lpermtotals[Lpermtotals==0] = np.nan
    Rpermtotals[Rpermtotals==0] = np.nan
    # convert final totals to proportion of participants in which that vertex is selected
    Lfinalselection = Ltotals/27
    Rfinalselection = Rtotals/27
    # plot selection
    # Left hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_left, Lfinalselection, hemi='left',
    colorbar=False, cmap='jet', alpha=1, vmin = 0, vmax = 1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_finaldist.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_right, Rfinalselection, hemi='right',
    colorbar=False, cmap='jet', alpha=1, vmin = 0, vmax = 1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_finaldist.jpeg',dpi=300)
    plt.close(fig)
    # Left hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_left, Lfinalselection, hemi='left',view=[270,270],
    colorbar=False, cmap='jet', alpha=1, vmin = 0, vmax = 1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_finaldist.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_right, Rfinalselection, hemi='right',view=[270,270],
    colorbar=False, cmap='jet', alpha=1, vmin = 0, vmax = 1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_finaldist.jpeg',dpi=300)
    plt.close(fig)
    
    # calculate coefficient directions
    Ldirections = Lanimatecoefs/Ltotals
    Rdirections = Ranimatecoefs/Rtotals
    # plot unthresholded
    # Left hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_left, Ldirections, hemi='left',
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_finalunthresholded.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_right, Rdirections, hemi='right',
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_finalunthresholded.jpeg',dpi=300)
    plt.close(fig)
    # Left hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_left, Ldirections, hemi='left',view=[270,270],
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_finalunthresholded.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_right, Rdirections, hemi='right',view=[270,270],
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_finalunthresholded.jpeg',dpi=300)
    plt.close(fig)
    
    # convert perm totals to proportion of permutations in which that vertex is selected
    Lpermselection = Lpermtotals/2700
    Rpermselection = Rpermtotals/2700
    # plot permutation distribution
    # Left hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_left, Lpermselection, hemi='left',
    colorbar=False, cmap='jet', threshold=0,alpha=1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_permdist.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_right, Rpermselection, hemi='right',
    colorbar=False, cmap='jet', threshold=0.000001, alpha=1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_permdist.jpeg',dpi=300)
    plt.close(fig)
    # Left hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_left, Lpermselection, hemi='left',view=[270,270],
    colorbar=False, cmap='jet', threshold=0.000001, alpha=1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_permdist.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_right, Rpermselection, hemi='right',view=[270,270],
    colorbar=False, cmap='jet', threshold=0.000001, alpha=1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_permdist.jpeg',dpi=300)
    plt.close(fig)
    
    # calculate perm directions
    Lpermdirections = Lpermanimatecoefs/Lpermtotals
    Rpermdirections = Rpermanimatecoefs/Rpermtotals
    # plot perm directions
    # Left hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_left, Lpermdirections, hemi='left',
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_permunthresholded.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_right, Rpermdirections, hemi='right',
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_permunthresholded.jpeg',dpi=300)
    plt.close(fig)
    # Left hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_left, Lpermdirections, hemi='left',view=[270,270],
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_permunthresholded.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_right, Rpermdirections, hemi='right',view=[270,270],
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_permunthresholded.jpeg',dpi=300)
    plt.close(fig)
    
    # # threshold attempt 1 - if the proportion of final participants in which that vertex is selected
    # # is LESS than the proportion of permutations in which that vertex is selected, set it to 0
    # Ldirections[Lfinalselection < Lpermselection] = 0
    # Rdirections[Rfinalselection < Rpermselection] = 0
    # # plot thresholded
    # # Left hemisphere viewed from the side
    # fig = plotting.plot_surf_stat_map(
    # pial_left, Ldirections, hemi='left',
    # colorbar=False, cmap='jet', threshold=0.000001,alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_dimension-'+d+'_thresholded1.jpeg',dpi=300)
    # plt.close(fig)
    
    # # Right hemisphere viewed from the side
    # fig = plotting.plot_surf_stat_map(
    # pial_right, Rdirections, hemi='right',
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_dimension-'+d+'_thresholded1.jpeg',dpi=300)
    # plt.close(fig)
    
    # # Left hemisphere viewed from underneath
    # fig = plotting.plot_surf_stat_map(
    # pial_left, Ldirections, hemi='left',view=[270,270],
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_dimension-'+d+'_thresholded1.jpeg',dpi=300)
    # plt.close(fig)
    # # Right hemisphere viewed from underneath
    # fig = plotting.plot_surf_stat_map(
    # pial_right, Rdirections, hemi='right',view=[270,270],
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_dimension-'+d+'_thresholded1.jpeg',dpi=300)
    # plt.close(fig)
    # # threshold attempt 2 - binomial probability of selection
    # Lpval = np.ones(163842)
    # Rpval = np.ones(163842)
    # for i in range(len(Ltotals)):
    #     tmp = stats.binomtest(int(Ltotals[i]),27,Lpermselection[i],alternative="greater")
    #     Lpval[i] = tmp.pvalue
    #     tmp = stats.binomtest(int(Rtotals[i]),27,Rpermselection[i],alternative="greater")
    #     Rpval[i] = tmp.pvalue
    # Ldirections = Lanimatecoefs/Ltotals
    # Rdirections = Ranimatecoefs/Rtotals
    # Ldirections[Lpval > 0.05] = 0
    # Rdirections[Rpval > 0.05] = 0
    # # plot thresholded
    # # Left hemisphere viewed from the side
    # fig = plotting.plot_surf_stat_map(
    # pial_left, Ldirections, hemi='left',
    # colorbar=False, cmap='jet', threshold=0.000001,alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_dimension-'+d+'_thresholded2.jpeg',dpi=300)
    # plt.close(fig)
    # # Right hemisphere viewed from the side
    # fig = plotting.plot_surf_stat_map(
    # pial_right, Rdirections, hemi='right',
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_dimension-'+d+'_thresholded2.jpeg',dpi=300)
    # plt.close(fig)
    # # Left hemisphere viewed from underneath
    # fig = plotting.plot_surf_stat_map(
    # pial_left, Ldirections, hemi='left',view=[270,270],
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_dimension-'+d+'_thresholded2.jpeg',dpi=300)
    # plt.close(fig)
    # # Right hemisphere viewed from underneath
    # fig = plotting.plot_surf_stat_map(
    # pial_right, Rdirections, hemi='right',view=[270,270],
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_dimension-'+d+'_thresholded2.jpeg',dpi=300)
    # plt.close(fig)
    
    # threshold attempt 3 - as above but FDR-corrected
    Ltotals[np.isnan(Ltotals)] = 0
    Rtotals[np.isnan(Rtotals)] = 0
    Lpermselection[np.isnan(Lpermselection)] = 0
    Rpermselection[np.isnan(Rpermselection)] = 0
    Lpval = np.ones(163842)
    Rpval = np.ones(163842)
    for i in range(len(Ltotals)):
        tmp = stats.binomtest(int(Ltotals[i]),27,Lpermselection[i],alternative="greater")
        Lpval[i] = tmp.pvalue
        tmp = stats.binomtest(int(Rtotals[i]),27,Rpermselection[i],alternative="greater")
        Rpval[i] = tmp.pvalue
    Lpval = stats.false_discovery_control(Lpval)
    Rpval = stats.false_discovery_control(Rpval)
    Ltotals[Ltotals==0] = np.nan
    Rtotals[Rtotals==0] = np.nan
    Ldirections = Lanimatecoefs/Ltotals
    Rdirections = Ranimatecoefs/Rtotals
    Ldirections[Lpval > 0.05] = np.nan
    Rdirections[Rpval > 0.05] = np.nan
    # plot thresholded
    # Left hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_left, Ldirections, hemi='left',
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_thresholded.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from the side
    fig = plotting.plot_surf_stat_map(
    pial_right, Rdirections, hemi='right',
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_thresholded.jpeg',dpi=300)
    plt.close(fig)
    # Left hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_left, Ldirections, hemi='left',view=[270,270],
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_thresholded.jpeg',dpi=300)
    plt.close(fig)
    # Right hemisphere viewed from underneath
    fig = plotting.plot_surf_stat_map(
    pial_right, Rdirections, hemi='right',view=[270,270],
    colorbar=False, cmap='jet', alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    bg_on_data=True,
    )
    fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_thresholded.jpeg',dpi=300)
    plt.close(fig)
    
    del(Lanimatecoefs,Ldirections,Lfc,Lfinalselection,Lpc,Lpermanimatecoefs,Lpermdirections,Lperms,Lpermselection,Lpermtotals,Lpval,Ltotals,Ranimatecoefs,Rdirections,Rfc,Rfinalselection,Rpc,Rpermanimatecoefs,Rpermdirections,Rperms,Rpermselection,Rpermtotals,Rpval,Rtotals)
    
    # # threshold attempt 4 - binomial probability of consistent direction
    # # For each voxel, some participants will have a positive coefficient and others will have
    # # a negative coefficient. Find the larger of these totals
    # Lagreement=np.maximum(Lanimatecoefs,Ltotals-Lanimatecoefs)
    # Ragreement=np.maximum(Ranimatecoefs,Rtotals-Ranimatecoefs)
    # # calculate the same agreement for the perms, but as a probability
    # Lpermagreement=np.maximum(Lpermdirections,1-Lpermdirections)
    # Rpermagreement=np.maximum(Rpermdirections,1-Rpermdirections)
    
    # # Find the p-value for each of these
    # Lpval = np.ones(163842)
    # Rpval = np.ones(163842)
    # for i in range(len(Lagreement)):
    #     if Ltotals[i] > 0:
    #         tmp = stats.binomtest(int(Lagreement[i]),int(Ltotals[i]),Lpermagreement[i],alternative="greater")
    #         Lpval[i] = tmp.pvalue
    #     if Rtotals[i] > 0:
    #         tmp = stats.binomtest(int(Ragreement[i]),int(Rtotals[i]),Rpermagreement[i],alternative="greater")
    #         Rpval[i] = tmp.pvalue
    # Ldirections = Lanimatecoefs/Ltotals
    # Rdirections = Ranimatecoefs/Rtotals
    # Ldirections[Lpval > 0.05] = 0
    # Rdirections[Rpval > 0.05] = 0
    # # plot thresholded
    # # Left hemisphere viewed from the side
    # fig = plotting.plot_surf_stat_map(
    # pial_left, Ldirections, hemi='left',
    # colorbar=False, cmap='jet', threshold=0.000001,alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_dimension-'+d+'_thresholded4.jpeg',dpi=300)
    # plt.close(fig)
    # # Right hemisphere viewed from the side
    # fig = plotting.plot_surf_stat_map(
    # pial_right, Rdirections, hemi='right',
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_dimension-'+d+'_thresholded4.jpeg',dpi=300)
    # plt.close(fig)
    # # Left hemisphere viewed from underneath
    # fig = plotting.plot_surf_stat_map(
    # pial_left, Ldirections, hemi='left',view=[270,270],
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_dimension-'+d+'_thresholded4.jpeg',dpi=300)
    # plt.close(fig)
    # # Right hemisphere viewed from underneath
    # fig = plotting.plot_surf_stat_map(
    # pial_right, Rdirections, hemi='right',view=[270,270],
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_dimension-'+d+'_thresholded4.jpeg',dpi=300)
    # plt.close(fig)
    
    # # threshold attempt 5 - as above but FDR-corrected
    # Ldirections = Lanimatecoefs/Ltotals
    # Rdirections = Ranimatecoefs/Rtotals
    # Lpval = stats.false_discovery_control(Lpval)
    # Rpval = stats.false_discovery_control(Rpval)
    # Ldirections[Lpval > 0.05] = 0
    # Rdirections[Rpval > 0.05] = 0
    # # plot thresholded
    # # Left hemisphere viewed from the side
    # fig = plotting.plot_surf_stat_map(
    # pial_left, Ldirections, hemi='left',
    # colorbar=False, cmap='jet', threshold=0.000001,alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_'+current_analysis+'_dimension-'+d+'_thresholded5.jpeg',dpi=300)
    # plt.close(fig)
    # # Right hemisphere viewed from the side
    # fig = plotting.plot_surf_stat_map(
    # pial_right, Rdirections, hemi='right',
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_'+current_analysis+'_dimension-'+d+'_thresholded5.jpeg',dpi=300)
    # plt.close(fig)
    # # Left hemisphere viewed from underneath
    # fig = plotting.plot_surf_stat_map(
    # pial_left, Ldirections, hemi='left',view=[270,270],
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_left,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_'+current_analysis+'_dimension-'+d+'_thresholded5.jpeg',dpi=300)
    # plt.close(fig)
    # # Right hemisphere viewed from underneath
    # fig = plotting.plot_surf_stat_map(
    # pial_right, Rdirections, hemi='right',view=[270,270],
    # colorbar=False, cmap='jet', threshold=0.000001, alpha=1, vmin=0, vmax = 1, bg_map= curv_right,
    # bg_on_data=True,
    # )
    # fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_'+current_analysis+'_dimension-'+d+'_thresholded5.jpeg',dpi=300)
    # plt.close(fig)
    # # clear up
    #del(Lagreement,Lanimatecoefs,Ldirections,Lfc,Lfinalselection,Lpc,Lpermagreement,Lpermanimatecoefs,Lpermdirections,Lperms,Lpermselection,Lpermtotals,Lpval,Ltotals,Ragreement,Ranimatecoefs,Rdirections,Rfc,Rfinalselection,Rpc,Rpermagreement,Rpermanimatecoefs,Rpermdirections,Rperms,Rpermselection,Rpermtotals,Rpval,Rtotals)
  