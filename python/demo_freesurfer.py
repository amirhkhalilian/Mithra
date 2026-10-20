import os, glob, argparse

import numpy as np
import pandas as pd
import pyvista as pv
from matplotlib import colormaps

from visualtools import visualtools

FS_COLS = ['fs_x', 'fs_y', 'fs_z']


def load_mapped_elecs(data_dir, subjects=None,
                      fn_name='coordinates_fs_atlas.csv'):
    '''
    loads the hcp_mapping output of each subject in data_dir and
    keeps the electrodes that are mapped to fsaverage (fs_x/y/z not NaN).
    returns a single dataframe with a 'subj' column.
    '''
    if subjects:
        fns = [os.path.join(data_dir, s, fn_name) for s in subjects]
    else:
        fns = sorted(glob.glob(os.path.join(data_dir, '*', fn_name)))
    dfs = []
    for fn in fns:
        df = pd.read_csv(fn)
        df['subj'] = os.path.basename(os.path.dirname(fn))
        dfs.append(df[df[FS_COLS].notna().all(axis=1)])
    return pd.concat(dfs, ignore_index=True)


def demo_mapped_elecs_on_fsaverage(data_dir, fs_dir, subjects=None,
                                   radius=1.5, fn_save=None):
    '''
    This example plots the electrodes mapped to fsaverage by
    hcp_mapping.py on the fsaverage pial surface. Each subject
    gets its own color. Left and right hemispheres are shown
    in separate linked subplots.
    '''
    df = load_mapped_elecs(data_dir, subjects=subjects)
    subjs = sorted(df['subj'].unique())
    # one color per subject from the tab10/tab20 colormap
    cmap = colormaps['tab10' if len(subjs) <= 10 else 'tab20']
    colors = {s: tuple(cmap(i)) for i, s in enumerate(subjs)}
    # create the plotting grid, one subplot per hemisphere
    pl = pv.Plotter(shape=(1, 2), off_screen=fn_save is not None)
    for col, hemi in enumerate(['lh', 'rh']):
        pl.subplot(0, col)
        # plot the fsaverage pial surface of this hemisphere
        vt = visualtools(subj='fsaverage', hemi=hemi,
                         brain_file=os.path.join(fs_dir, f'{hemi}.pial'))
        vt.plot_pial_surface(plotter=pl)
        # plot the electrodes of each subject on top
        legend = []
        for s in subjs:
            m = (df['subj'] == s) & (df['hemi'] == hemi)
            if not m.any():
                continue
            vt.plot_elec_on_pial(df.loc[m, FS_COLS].to_numpy(dtype=float),
                                 elec_colors= colors[s],
                                 radius=radius,
                                 plotter=pl,
                                 brain_plotter=pl,
                                 show=False)
            legend.append([f'{s} (n={m.sum()})', colors[s]])
        if legend:
            pl.add_legend(legend, bcolor='white', face='circle')
        pl.add_text(f'fsaverage {hemi}', font_size=10,
                    font='times', position='upper_edge')
    pl.set_background(color='white')
    if fn_save is not None:
        pl.screenshot(fn_save, window_size=(1600, 700))
        print(f'saved {fn_save}')
    else:
        pl.show()
    return pl


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='plot electrodes mapped to fsaverage')
    parser.add_argument('--data_dir', default='./Data')
    parser.add_argument('--fs_dir', default='./sample_data/fsaverage')
    parser.add_argument('--subjects', nargs='*', default=None)
    parser.add_argument('--radius', type=float, default=1.5)
    parser.add_argument('--save', default=None,
                        help='save a screenshot to this file instead of showing')
    args = parser.parse_args()
    demo_mapped_elecs_on_fsaverage(args.data_dir, args.fs_dir,
                                   subjects=args.subjects,
                                   radius=args.radius,
                                   fn_save=args.save)
