import os

import numpy as np
import pyvista as pv

from visualtools import visualtools
from nibabel.freesurfer.io import read_geometry

from mat73 import loadmat
from misc import read_coords_txt

def rgb_to_ansi(r, g, b):
    """Convert RGB in [0,1] to ANSI 24-bit color escape code."""
    R, G, B = int(r * 255), int(g * 255), int(b * 255)
    return f"\033[38;2;{R};{G};{B}m"

# 10 color-blind–friendly colors (as RGB tuples)
colors = [
    (0.0, 0.45, 0.70),  # blue
    (0.80, 0.47, 0.65), # reddish purple
    (0.90, 0.60, 0.0),  # orange
    (0.0, 0.60, 0.50),  # bluish green
    (0.95, 0.90, 0.25), # yellow
    (0.56, 0., 0.),      # dark red
    (0.0, 0.85, 0.90),  # cyan
    (0.5, 0.5, 0.5),    # gray
    (0.8, 0.4, 0.0),    # vermilion
    (0.9, 0.6, 0.6),    # pinkish
    (0.35, 0.70, 0.90), # sky blue
]
colors = [(*c, 1.0) for c in colors]

coords = read_coords_txt('./sample_data/sample_subjects/NY912/NY912_T1.txt')

elec_coords = []
elec_colors = []
seen_shaft = {}
color_ind = -1
for lbl, coord, region in coords:
    if region.startswith('ctx-'):
        shaft_name = ''.join([c for c in lbl if c.isalpha()])
        if shaft_name=="SPI" or shaft_name=="SAI":
            continue
        elec_coords.append(coord)
        if shaft_name in seen_shaft:
            elec_colors.append(seen_shaft[shaft_name])
        else:
            color_ind += 1
            seen_shaft[shaft_name] = colors[color_ind]
            elec_colors.append(seen_shaft[shaft_name])
            r, g, b, _ = colors[color_ind]
            color_code = rgb_to_ansi(r, g, b)
            print(f"{color_code}Added new shaft: {shaft_name}\033[0m")
elec_coords = np.array(elec_coords)
elec_colors = np.array(elec_colors)
print(type(elec_colors))

# creatre the visualization object
vt = visualtools(subj='fs_avg', hemi='lh',
                 flag_use_annot = False,
                 brain_file='./sample_data/sample_subjects/NY912/NY912_lh_pial.mat',
                 annot_file='./sample_data/sample_subjects/NY912/lh.aparc.annot')
verts, faces = vt._get_trisurf()
print(verts.shape)

# define a plotter object
pl = pv.Plotter()
# plot and show the pial surface
pl = vt.plot_pial_surface(plotter=pl, brain_color=(0.9,0.9,0.9,0.95))

pl = vt.plot_elec_on_pial(elec_coords,
                          elec_colors=elec_colors,
                          radius=0.75,
                          brain_plotter=pl,
                          flag_scalar_bar = False,
                          show = False)
# pl.camera.Azimuth(-45)
# pl.camera.Roll(30)
# pl.show()
# quit()
vt._save_all_views(pl, 'NY912')

fn_fs_pial = "./sample_data/fsaverage/lh.pial"
fn_fs_sphr = "./sample_data/fsaverage/lh.sphere.reg"
fn_subj_pial = './sample_data/sample_subjects/NY912/NY912_lh_pial.mat'
fn_subj_sphr ='./sample_data/sample_subjects/NY912/lh.sphere.reg'

import numpy as np
import pyvista as pv
from nibabel.freesurfer import read_geometry
from visualtools import visualtools  # assuming your class is in visualtools.py

# create the visualization object
vt = visualtools(subj='fs_avg', hemi='lh',
                 flag_use_annot=False,
                 brain_file='./sample_data/sample_subjects/NY912/NY912_lh_pial.mat',
                 annot_file='./sample_data/sample_subjects/NY912/lh.aparc.annot')

verts_subj, faces_subj = vt._get_trisurf()

# -------------------------------------------------------------------------
# Conversion from subject → fsaverage using surface registration
# -------------------------------------------------------------------------
fn_fs_pial = "./sample_data/fsaverage/lh.pial"
fn_fs_inf = "./sample_data/fsaverage/lh.inflated"
fn_fs_sphr = "./sample_data/fsaverage/lh.sphere.reg"
fn_subj_pial = './sample_data/sample_subjects/NY912/NY912_lh_pial.mat'
fn_subj_sphr = './sample_data/sample_subjects/NY912/lh.sphere.reg'

verts_fs_pial, faces_fs_pial = read_geometry(fn_fs_pial)
verts_fs_inf, faces_fs_inf = read_geometry(fn_fs_inf)
verts_fs_sphr, faces_fs_sphr = read_geometry(fn_fs_sphr)
verts_subj_pial = verts_subj.copy()
verts_subj_sphr, faces_subj_sphr = read_geometry(fn_subj_sphr)

# assuming `elec_coords` are electrode coordinates in subject space
fs_elec_coords = np.zeros_like(elec_coords)
inf_elec_coords = np.zeros_like(elec_coords)
for i in range(elec_coords.shape[0]):
    subj_ind = np.argmin(np.sum((verts_subj_pial - elec_coords[i, :]) ** 2, axis=1))
    fs_ind = np.argmin(np.sum((verts_fs_sphr - verts_subj_sphr[subj_ind, :]) ** 2, axis=1))
    fs_elec_coords[i, :] = verts_fs_pial[fs_ind, :]
    inf_elec_coords[i, :] = verts_fs_inf[fs_ind, :]


vt = visualtools(subj='fs_avg', hemi='lh',
                 flag_use_annot = False,
                 brain_file = fn_fs_pial)
verts, faces = vt._get_trisurf()

# define a plotter object
pl = pv.Plotter()
# plot and show the pial surface
pl = vt.plot_pial_surface(plotter=pl, brain_color=(0.9,0.9,0.9,0.95))

pl = vt.plot_elec_on_pial(fs_elec_coords,
                          elec_colors=elec_colors,
                          radius=1.0,
                          brain_plotter=pl,
                          flag_scalar_bar = False,
                          show = False)
# pl.camera.Azimuth(-45)
# pl.camera.Elevation(-10)
# pl.show()
vt._save_all_views(pl, 'fspial')

vt = visualtools(subj='fs_avg', hemi='lh',
                 flag_use_annot = False,
                 brain_file = fn_fs_inf)

# define a plotter object
pl = pv.Plotter()
# plot and show the pial surface
pl = vt.plot_inflated_surf(morph_file="./sample_data/fsaverage/lh.curv",
                           plotter=pl)

pl = vt.plot_elec_on_pial(inf_elec_coords,
                          elec_colors=elec_colors,
                          radius=1.0,
                          brain_plotter=pl,
                          flag_scalar_bar = False,
                          show = False)
# pl.camera.Azimuth(-45)
# pl.camera.Elevation(-10)
# pl.show()
vt._save_all_views(pl, 'fsinf')
