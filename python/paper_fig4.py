import os
import numpy as np
import pandas as pd
from scipy.io import loadmat
import matplotlib.pyplot as plt
import pyvista as pv
from visualtools import visualtools

# form the dataset
windows = [0, 100, 200, 400]
dataset = {}
subj = "NY742"
path = f"./sample_data/sample_subjects/{subj}/AudRep"
for w in windows:
    dl = loadmat(os.path.join(path, f"{subj}_AudRep_{w}msec.mat"))
    dataset[w] = np.abs(dl['sig'][:128])
    dataset[w][dataset[w]<=0.2] = 0.0 # noise-level

# load coordinates
fn_csv = f"./sample_data/sample_subjects/{subj}/coordinates.csv"
csv = pd.read_csv(fn_csv)
se = np.arange(128) # first 128 electrodes are grid
x = np.array(csv['T1_x'][se])[:,None]
y = np.array(csv['T1_y'][se])[:,None]
z = np.array(csv['T1_z'][se])[:,None]
elec_coords = np.hstack((x,y,z))
x = np.array(csv['MNI_x'][se])[:,None]
y = np.array(csv['MNI_y'][se])[:,None]
z = np.array(csv['MNI_z'][se])[:,None]
elec_coords_mni = np.hstack((x,y,z))
elec_regions = csv['T1_AnatomicalRegion'][se] # str list of 128x1

fn_csv = f"./sample_data/sample_subjects/{subj}/coordinates_fs.csv"
csv = pd.read_csv(fn_csv)
x = np.array(csv['fs_coords_1'][se])[:,None]
y = np.array(csv['fs_coords_2'][se])[:,None]
z = np.array(csv['fs_coords_3'][se])[:,None]
elec_coords_fs = np.hstack((x,y,z))

vt = visualtools(subj=subj, hemi='lh',
                 flag_use_annot=True,
                 brain_file=f'./sample_data/sample_subjects/{subj}/{subj}_pial_surf.mat',
                 annot_file=f'./sample_data/sample_subjects/{subj}/{subj}_lh_aparc.annot')

vt_mni = visualtools(subj='MNI_FSL', hemi='lh',
                 flag_use_annot=True,
                 brain_file='./sample_data/MNI_FSL152/FSL_MNI152_lh_pial.mat',
                 annot_file='./sample_data/MNI_FSL152/FSL_MNI152.lh.aparc.split_STG_MTG.annot')
elec_coords_mni, _ = vt_mni.project_elec_verts(elec_coords_mni,
                                               regions = elec_regions)

vt_fs = visualtools(subj='fsaverge', hemi='lh',
                 flag_use_annot=True,
                 brain_file='./sample_data/fsaverage/lh.pial',
                 annot_file='./sample_data/fsaverage/lh.aparc.split_STG_MTG.annot')


pl = pv.Plotter(shape=(3,4))

for i, w in enumerate(windows):
    pl.subplot(0,i)
    vt.flag_use_annot = True
    vt.flag_merge_STG = True
    plb = vt.plot_weighted_surface(elec_coords,
                                   dataset[w],
                                   regions=elec_regions,
                                   sigma = 3.0,
                                   plotter=pl,
                                   flag_scalar_bar = False,
                                   mode="normalize",
                                   cmap='hot_r',
                                   clim=(0,2),
                                   show=False)
    vt.flag_use_annot = False
    vt.plot_elec_on_pial(elec_coords,
                         elec_colors=dataset[w].squeeze(),
                         radius=.75,
                         plotter=pl,
                         brain_plotter=plb,
                         flag_scalar_bar = False,
                         cmap='hot_r',
                         clim=(0,1.2),
                         show=False)
    pl.subplot(1,i)
    vt_mni.flag_use_annot = True
    vt_mni.flag_merge_STG = True
    plb = vt_mni.plot_weighted_surface(elec_coords_mni,
                                       dataset[w],
                                       regions=elec_regions,
                                       sigma = 3.0,
                                       plotter=pl,
                                       mode="normalize",
                                       cmap='hot_r',
                                       flag_scalar_bar = False,
                                       clim=(0,2),
                                       show=False)
    vt_mni.flag_use_annot = False
    vt_mni.plot_elec_on_pial(elec_coords_mni,
                             elec_colors=dataset[w].squeeze(),
                             radius=.75,
                             plotter=pl,
                             brain_plotter=plb,
                             flag_scalar_bar = False,
                             cmap='hot_r',
                             clim=(0,1.2),
                             show=False)
    pl.subplot(2,i)
    vt_fs.flag_use_annot = True
    vt_fs.flag_merge_STG = True
    plb = vt_fs.plot_weighted_surface(elec_coords_fs,
                                       dataset[w],
                                       regions=elec_regions,
                                       sigma = 3.0,
                                       plotter=pl,
                                       mode="simple",
                                       cmap='hot_r',
                                       flag_scalar_bar = False,
                                       clim=(0,2),
                                       show=False)
    vt_mni.flag_use_annot = False
    vt_fs.plot_elec_on_pial(elec_coords_fs,
                             elec_colors=dataset[w].squeeze(),
                             radius=.75,
                             plotter=pl,
                             brain_plotter=plb,
                             flag_scalar_bar = False,
                             cmap='hot_r',
                             clim=(0,1.2),
                             show=False)

pl.link_views()
pl.show()


