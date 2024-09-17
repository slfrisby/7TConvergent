# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 18:35:23 2024

@author: sf02
"""

import sys
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from nilearn import surface, plotting

# set path
dirp='/imaging/projects/cbu/wbic-p00591-DAISY/main/'
sys.path.append(dirp)

# get pial surface mesh and curvature data
sys.path.append(dirp+'/scripts/fsaverage')
pial_left = dirp+'/scripts/fsaverage/pial_left.gii.gz'
curv_left = dirp+'/scripts/fsaverage/curv_left.gii.gz'
pial_right = dirp+'/scripts/fsaverage/pial_right.gii.gz'
curv_right = dirp+'/scripts/fsaverage/curv_right.gii.gz'

# PLOT UNIVARIATE

# load data
L_uni_animate = surface.load_surf_data(dirp+"/work/ROIsurf/L_univariate_animate.func.gii")
L_uni_inanimate = surface.load_surf_data(dirp+"/work/ROIsurf/L_univariate_inanimate.func.gii")
R_uni_animate = surface.load_surf_data(dirp+"/work/ROIsurf/R_univariate_animate.func.gii")
R_uni_inanimate = surface.load_surf_data(dirp+"/work/ROIsurf/R_univariate_inanimate.func.gii")

# combine into a single image
L_univariate = np.zeros(163842,)
L_univariate[:] = np.nan
L_univariate[np.where(~np.isnan(L_uni_animate))] = 2
L_univariate[np.where(~np.isnan(L_uni_inanimate))] = 1
R_univariate = np.zeros(163842,)
R_univariate[:] = np.nan
R_univariate[np.where(~np.isnan(R_uni_animate))] = 2
R_univariate[np.where(~np.isnan(R_uni_inanimate))] = 1

# define colourmaps
univariatemap = ListedColormap(["darkblue","darkred"])

# plot
# Left hemisphere viewed from the side
fig = plotting.plot_surf_roi(
pial_left, L_univariate, hemi='left',
colorbar=False, cmap=univariatemap, alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_univariate.jpeg',dpi=300)
plt.close(fig)
# Right hemisphere viewed from the side
fig = plotting.plot_surf_roi(
pial_right, R_univariate, hemi='right',
colorbar=False, cmap=univariatemap, alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_univariate.jpeg',dpi=300)
plt.close(fig)
# Left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_left, L_univariate, hemi='left', view=[270,270],
colorbar=False, cmap=univariatemap, alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_univariate.jpeg',dpi=300)
plt.close(fig)
# Right hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_right, R_univariate, hemi='right', view=[270,270],
colorbar=False, cmap=univariatemap, alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_univariate.jpeg',dpi=300)
plt.close(fig)

# PLOT tSNR

# load data
L_tsnr = surface.load_surf_data(dirp+"/work/ROIsurf/L_tSNR.func.gii")
R_tsnr = surface.load_surf_data(dirp+"/work/ROIsurf/R_tSNR.func.gii")

# Left hemisphere viewed from the side
fig = plotting.plot_surf_stat_map(
pial_left, L_tsnr, hemi='left',
colorbar=False, cmap='jet', alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_lateral_tSNR.jpeg',dpi=300)
plt.close(fig)
# Right hemisphere viewed from the side
fig = plotting.plot_surf_stat_map(
pial_right, R_tsnr, hemi='right',
colorbar=False, cmap='jet', alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_lateral_tSNR.jpeg',dpi=300)
plt.close(fig)
# Left hemisphere viewed from underneath
fig = plotting.plot_surf_stat_map(
pial_left, L_tsnr, hemi='left', view=[270,270],
colorbar=False, cmap='jet', alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_tSNR.jpeg',dpi=300)
plt.close(fig)
# Right hemisphere viewed from underneath
fig = plotting.plot_surf_stat_map(
pial_right, R_tsnr, hemi='right', view=[270,270],
colorbar=False, cmap='jet', alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_tSNR.jpeg',dpi=300)
plt.close(fig)


# PLOT ROIs 

L_ATL = surface.load_surf_data(dirp+'/work/ROIsurf/L_vATL_ROI.func.gii')
R_ATL = surface.load_surf_data(dirp+'/work/ROIsurf/R_vATL_ROI.func.gii')
L_ATL_ant = surface.load_surf_data(dirp+'/work/ROIsurf/L_vATL_ant_ROI.func.gii')
L_ATL_pos = surface.load_surf_data(dirp+'/work/ROIsurf/L_vATL_pos_ROI.func.gii')
R_ATL_ant = surface.load_surf_data(dirp+'/work/ROIsurf/R_vATL_ant_ROI.func.gii')
R_ATL_pos = surface.load_surf_data(dirp+'/work/ROIsurf/R_vATL_pos_ROI.func.gii')

# binarise
L_ATL[L_ATL!=0] = 1
R_ATL[R_ATL!=0] = 1
L_ATL_ant[L_ATL_ant!=0] = 1
L_ATL_pos[L_ATL_pos!=0] = 1
R_ATL_ant[R_ATL_ant!=0] = 1
R_ATL_pos[R_ATL_pos!=0] = 1

# whole left hemisphere

# define colourmap for LASSO
tmpcmap = ListedColormap(["mediumblue"])

# whole left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_left, L_ATL, hemi='left', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_LASSO_ROI.jpeg',dpi=300)
plt.close(fig)

# define colourmap for SOSLASSO
tmpcmap = ListedColormap(['mediumvioletred'])

# whole left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_left, L_ATL, hemi='left', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_SOSLASSO_ROI.jpeg',dpi=300)
plt.close(fig)

# define colourmap for grOWL
tmpcmap = ListedColormap(["darkred"])

# whole left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_left, L_ATL, hemi='left', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_grOWL_ROI.jpeg',dpi=300)
plt.close(fig)

# halves of left hemisphere

# define halves of left hemisphere
tmp = np.zeros(163842,)
tmp[:] = np.nan
tmp[L_ATL_ant == 1] = 1
tmp[L_ATL_pos == 1] = 2

# define colourmap for LASSO
tmpcmap = ListedColormap(["dodgerblue","lightskyblue"])

# halves of left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_left, tmp, hemi='left', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_LASSO_ROI_halves.jpeg',dpi=300)
plt.close(fig)

# define colourmap for SOSLASSO
tmpcmap = ListedColormap(["palevioletred","pink"])

# halves of left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_left, tmp, hemi='left', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_SOSLASSO_ROI_halves.jpeg',dpi=300)
plt.close(fig)

# define colourmap for grOWL
tmpcmap = ListedColormap(["crimson","lightcoral"])

# halves of left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_left, tmp, hemi='left', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_left,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/L_ventral_grOWL_ROI_halves.jpeg',dpi=300)
plt.close(fig)

# whole right hemisphere

# define colourmap for LASSO
tmpcmap = ListedColormap(["green"])

# whole right hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_right, R_ATL, hemi='right', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_LASSO_ROI.jpeg',dpi=300)
plt.close(fig)

# define colourmap for SOSLASSO
tmpcmap = ListedColormap(["rebeccapurple"])

# whole right hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_right, R_ATL, hemi='right', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_SOSLASSO_ROI.jpeg',dpi=300)
plt.close(fig)

# define colourmap for grOWL
tmpcmap = ListedColormap(["gold"])

# whole right hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_right, R_ATL, hemi='right', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_grOWL_ROI.jpeg',dpi=300)
plt.close(fig)

# halves of right hemisphere

# define halves of right hemisphere
tmp = np.zeros(163842,)
tmp[:] = np.nan
tmp[R_ATL_ant == 1] = 1
tmp[R_ATL_pos == 1] = 2

# define colourmap for LASSO
tmpcmap = ListedColormap(["limegreen","palegreen"])

# halves of left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_right, tmp, hemi='right', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_LASSO_ROI_halves.jpeg',dpi=300)
plt.close(fig)

# define colourmap for SLASSO
tmpcmap = ListedColormap(["mediumorchid","plum"])

# halves of left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_right, tmp, hemi='right', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_SOSLASSO_ROI_halves.jpeg',dpi=300)
plt.close(fig)

# define colourmap for grOWL
tmpcmap = ListedColormap(["yellow","lemonchiffon"])

# halves of left hemisphere viewed from underneath
fig = plotting.plot_surf_roi(
pial_right, tmp, hemi='right', view=[270,270],
colorbar=False, cmap=tmpcmap, alpha=1, bg_map= curv_right,
bg_on_data=True,
)
fig.savefig('/group/mlr-lab/Saskia/7T_decoding/surface_plots/R_ventral_grOWL_ROI_halves.jpeg',dpi=300)
plt.close(fig)



