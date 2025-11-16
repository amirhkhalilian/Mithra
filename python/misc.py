import sys, os, glob

import numpy as np
import scipy.io
import pandas as pd
import pyvista as pv
from nibabel.freesurfer.io import read_annot
from nibabel.freesurfer.io import read_geometry
from matplotlib import pyplot as plt

from visualtools import visualtools

import re

def read_coords_txt(filepath):
    result = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('%'):
                continue

            parts = line.split()
            if len(parts) < 5:
                continue

            name = parts[0]
            coords = list(map(float, parts[1:4]))  # three numeric columns

            # Find all "number% label" pairs
            matches = re.findall(r'([\d.]+)%\s+(\S+)', line)
            if not matches:
                continue

            # Sort by descending percentage
            percentages = sorted(
                [(float(p), lbl) for p, lbl in matches],
                key=lambda x: -x[0]
            )

            # Choose first valid label
            chosen_label = None
            for p, lbl in percentages:
                if p>=30 and lbl != "Unknown" and not lbl.endswith("White-Matter"):
                    chosen_label = lbl
                    break

            if chosen_label is None:
                chosen_label = percentages[0][1]

            result.append((name, coords, chosen_label))
    return result

def project_elecs_pial(elec_locs, elec_regions, BrainFile, AnnotFile):
    '''
    function to project the electrodes to the closest vertex (Euclidean
    distance) with the same annotation.
    '''
    # laod the BrainFile
    ext = os.path.splitext(BrainFile)[1]
    if ext=='.mat':
        # load using scipy
        Brain = scipy.io.loadmat(BrainFile)
        # convert brain faces and verts to PolyData in pyvista
        verts = np.array(Brain['coords'])
        faces = np.hstack((3*np.ones((Brain['faces'].shape[0],1),dtype=np.int32),\
                           np.array(Brain['faces']-1)))
    elif ext=='.pial':
        # load using read_geometry
        dl = read_geometry(BrainFile)
        verts = np.array(dl[0])
        faces = np.hstack((3*np.ones((dl[1].shape[0],1),dtype=np.int32),\
                           np.array(dl[1])))
    else:
        raise ValueError(f'file extension {ext} not valid!')

    # load the AnnotFile
    # Annot[0] -- list of region index per vertex
    # Annot[2] -- region name
    Annot = read_annot(AnnotFile)
    elec_locs_proj = np.zeros_like(elec_locs)
    elec_rgb = np.ones((elec_locs.shape[0],4))

    # Not efficient code...
    for ei in range(elec_locs.shape[0]):
        # find region index in the list
        ri = Annot[2].index(elec_regions[ei].encode('utf-8'))
        elec_rgb[ei,:] = Annot[1][ri,0:4]
        # find vertices in region
        vi = np.argwhere(Annot[0]==ri)
        if len(vi)==0:
            verts_r = verts
        else:
            verts_r = verts[Annot[0]==ri,:]
        # compute the shortest dist
        dist = verts_r - elec_locs[ei,:]
        dist = np.sum(dist**2, axis=1)
        elec_locs_proj[ei,:] = verts_r[np.argmin(dist),:]
    elec_rgb /= 255
    elec_rgb[:,3] = 1.0
    return elec_locs_proj, elec_rgb

if __name__=="__main__":
    # first we load the "coordinates.csv" file that has electrode
    # locations in subject T1 and MNI spaces.
    fn_csv = "./SampleData/NY717/coordinates.csv"
    csv = pd.read_csv(fn_csv)
    # we select the grid electrodes for now
    se = np.arange(128) # first 128 electrodes are grid
    # we get the x,y,z locations for the elctrodes in T1 space
    x = np.array(csv['MNI_x'][se])[:,None]
    y = np.array(csv['MNI_y'][se])[:,None]
    z = np.array(csv['MNI_z'][se])[:,None]
    elec_locs = np.hstack((x,y,z)) # a 128 x 3 matrix of coordinates
    # get the electrode locations in subject space
    elec_regions = csv['T1_AnatomicalRegion'][se] # str list of 128x1
    # set MNI brain and annot files
    AnnotFile='./SampleData/MNI/ch2_template.lh.aparc.split_STG_MTG.annot'
    BrainFile='./SampleData/MNI/ch2_template_lh_pial_120519.mat'
    # call projection function
    elec_locs_proj, elec_rgb = project_elecs_pial(elec_locs,
                                                  elec_regions,
                                                  BrainFile,
                                                  AnnotFile)
    # visualize the results
    VT = visualtools(Subj='MNI', HS='lh',
                     flag_UseAnnots = False,
                     BrainFile=BrainFile)
    pl = pv.Plotter(shape=(1,2)) # a 2x2 grid to plot different visuals
    ## plot the electrodes on the brain -- elecs are region color-coded
    pl.subplot(0,0) # select first subplot
    # use PlotElecOnBrain function from visualtools to plot the elecs
    VT.PlotElecOnBrain(elec_locs, ElecColor=elec_rgb,\
                       radius=1.0, BrainPlotter=pl)
    pl.subplot(0,1) # select first subplot
    # use PlotElecOnBrain function from visualtools to plot the elecs
    VT.PlotElecOnBrain(elec_locs_proj, ElecColor=elec_rgb,\
                       radius=1.0, BrainPlotter=pl)
    pl.link_views()
    pl.show()




