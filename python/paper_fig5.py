import os
import numpy as np
import pandas as pd
from scipy.io import loadmat
import matplotlib.pyplot as plt
import pyvista as pv
from visualtools import visualtools

subjs = ['NY717', 'NY742', 'NY749', 'NY857', 'NY869']
grid =  [128, 128, 128, 64]
windows = [0, 100, 200, 400]

for w in windows:
    dataset = []
    mni = []
    fsa = []
    regions = []
    for si, subj in enumerate(subjs):
        path = f"./sample_data/sample_subjects/{subj}/AudRep"
        dl = loadmat(os.path.join(path, f"{subj}_AudRep_{w}msec.mat"))
        data = np.abs(dl['sig'][:grid[si]])
        data[data<=0.2] = 0.0
        dataset.append(data)

        fn_csv = f"./sample_data/sample_subjects/{subj}/coordinates.csv"
        csv = pd.read_csv(fn_csv)
        se = np.arange(grid[si])
        x = np.array(csv['MNI_x'][se])[:,None]
        y = np.array(csv['MNI_y'][se])[:,None]
        z = np.array(csv['MNI_z'][se])[:,None]
        elec_coords_mni = np.hstack((x,y,z))
        elec_regions = csv['T1_AnatomicalRegion'][se]
        mni.append(elec_coords_mni)
        regions.append(elec_regions)

        fn_csv = f"./sample_data/sample_subjects/{subj}/coordinates_fs.csv"
        csv = pd.read_csv(fn_csv)
        x = np.array(csv['fs_coords_1'][se])[:,None]
        y = np.array(csv['fs_coords_2'][se])[:,None]
        z = np.array(csv['fs_coords_3'][se])[:,None]
        elec_coords_fs = np.hstack((x,y,z))
        fsa.append(elec_coords_fs)


    dataset = np.concatenate(dataset)
    elec_coords_mni = np.concatenate(mni)
    elec_coords_fs = np.concatenate(fsa)
    elec_regions = np.concatenate(regions)

    print(dataset.shape)
    print(elec_coords_mni.shape)

    vt_mni = visualtools(subj='MNI_FSL', hemi='lh',
                     flag_use_annot=True,
                     brain_file='./sample_data/MNI_FSL152/FSL_MNI152_lh_pial.mat',
                     annot_file='./sample_data/MNI_FSL152/FSL_MNI152.lh.aparc.split_STG_MTG.annot')
    elec_coords_mni, _ = vt_mni.project_elec_verts(elec_coords_mni,
                                                   regions = elec_regions)
    print(elec_coords_mni.shape)

    vt_fs = visualtools(subj='fsaverge', hemi='lh',
                     flag_use_annot=True,
                     brain_file='./sample_data/fsaverage/lh.pial',
                     annot_file='./sample_data/fsaverage/lh.aparc.split_STG_MTG.annot')


    save_name = "Results"

    pl = pv.Plotter(shape=(1,1))
    vt_mni.flag_use_annot = True
    vt_mni.flag_merge_STG = True
    plb = vt_mni.plot_weighted_surface(elec_coords_mni,
                                       dataset,
                                       regions=elec_regions,
                                       sigma = 3.0,
                                       plotter=pl,
                                       mode="normalize",
                                       cmap='hot_r',
                                       flag_scalar_bar = False,
                                       clim=(0,1),
                                       show=False)
    vt_mni.flag_use_annot = False
    vt_mni.plot_elec_on_pial(elec_coords_mni,
                             elec_colors=dataset.squeeze(),
                             radius=.5,
                             plotter=pl,
                             brain_plotter=plb,
                             flag_scalar_bar = False,
                             cmap='hot_r',
                             clim=(0,1.5),
                             show=False)
    pl.save_graphic(os.path.join(save_name, f"MNI_{w}.svg"))
    pl.clear()
    pl = pv.Plotter(shape=(1,1))
    vt_fs.flag_use_annot = True
    vt_fs.flag_merge_STG = True
    plb = vt_fs.plot_weighted_surface(elec_coords_fs,
                                       dataset,
                                       regions=elec_regions,
                                       sigma = 3.0,
                                       plotter=pl,
                                       mode="normalize",
                                       cmap='hot_r',
                                       flag_scalar_bar = False,
                                       clim=(0,1),
                                       show=False)
    vt_mni.flag_use_annot = False
    vt_fs.plot_elec_on_pial(elec_coords_fs,
                             elec_colors=dataset.squeeze(),
                             radius=.5,
                             plotter=pl,
                             brain_plotter=plb,
                             flag_scalar_bar = False,
                             cmap='hot_r',
                             clim=(0,1.5),
                             show=False)
    pl.save_graphic(os.path.join(save_name, f"FS_{w}.svg"))
    pl.clear()


