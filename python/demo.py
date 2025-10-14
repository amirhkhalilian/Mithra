import numpy as np
import pandas as pd
import pyvista as pv

from visualtools import visualtools

def demo_plot_brain_pial():
    '''
    This example plots the brain surface with annotation
    '''
    # creatre the visualization object
    vt = visualtools(subj='NY717', hemi='lh',
                     flag_use_annot = True,
                     brain_file='./sample_data/sample_subjects/NY717/NY717_pial_surf.mat',
                     annot_file='./sample_data/sample_subjects/NY717/NY717_lh_aparc.annot')
    # define a plotter object
    pv.set_plot_theme('dark')
    pl = pv.Plotter(notebook=False)
    # plot and show the pial surface
    pl = vt.plot_pial_surface(plotter=pl)
    pl.show()

def demo_elces_on_subject_brain():
    '''
    This example plots the grid electrodes from NY717 on
    the NY717 brain surface. Each electrode is shown as
    a red sphere.
    '''
    # first we load the "coordinates.csv" file that has electrode
    # locations in subject T1 and MNI spaces.
    fn_csv = "./sample_data/sample_subjects/NY717/coordinates.csv"
    csv = pd.read_csv(fn_csv)
    # we select the grid electrodes for now
    se = np.arange(128) # first 128 electrodes are grid
    # we get the x,y,z locations for the elctrodes in T1 space
    x = np.array(csv['T1_x'][se])[:,None]
    y = np.array(csv['T1_y'][se])[:,None]
    z = np.array(csv['T1_z'][se])[:,None]
    elec_coords = np.hstack((x,y,z)) # a 128 x 3 matrix of coordinates

    # creatre the visualization object
    vt = visualtools(subj='NY717', hemi='lh',
                     brain_file='./sample_data/sample_subjects/NY717/NY717_pial_surf.mat')
    # plot the electrodes on the brain
    pl = vt.plot_elec_on_pial(elec_coords,
                              elec_colors=(1.0,0.0,0.0,1.0),
                              radius=1.0,
                              show=True)

def demo_elec_project_mni():
    '''
    this example plots the electrodes on subject pial surface
    and MNI surface before and after projection
    '''
    # first we load the "coordinates.csv" file that has electrode
    # locations in subject T1 and MNI spaces.
    fn_csv = "./sample_data/sample_subjects/NY863/coordinates.csv"
    csv = pd.read_csv(fn_csv)
    # we select the grid electrodes for now
    se = np.arange(88) # first 128 electrodes are grid
    # we get the x,y,z locations for the elctrodes in T1 space
    x = np.array(csv['T1_x'][se])[:,None]
    y = np.array(csv['T1_y'][se])[:,None]
    z = np.array(csv['T1_z'][se])[:,None]
    elec_coords_T1 = np.hstack((x,y,z)) # a 128 x 3 matrix of coordinates
    # we get the x,y,z locations for the elctrodes in MNI space
    x = np.array(csv['MNI_x'][se])[:,None]
    y = np.array(csv['MNI_y'][se])[:,None]
    z = np.array(csv['MNI_z'][se])[:,None]
    elec_coords_mni = np.hstack((x,y,z)) # a 128 x 3 matrix of coordinates
    elec_regions = csv['T1_AnatomicalRegion'][se] # str list of 128x1

    # creatre the visualization object
    vt = visualtools(subj='NY863', hemi='lh',
                     flag_use_annot=True,
                     brain_file='./sample_data/sample_subjects/NY863/NY863_pial_surf.mat',
                     annot_file='./sample_data/sample_subjects/NY863/NY863_lh_aparc.annot')
    verts, faces = vt._get_trisurf()
    elec_coords_T1, elec_rgb = vt.project_elec_verts(elec_coords_T1,
                                                     regions = elec_regions)
    # plot the electrodes on the subject pial brain
    vt.flag_use_annot = False
    pl = vt.plot_elec_on_pial(elec_coords_T1,
                              elec_colors=(1.0,0.0,0.0,1.0),
                              radius=1.5,
                              show=False)


    # plot the electrodes on the brain
    vt.flag_use_annot = False
    pl = vt.plot_elec_on_pial(elec_coords_T1,
                              elec_colors=elec_rgb,
                              radius=0.9,
                              flag_scalar_bar=False,
                              show=False)

    # plot the electrodes on MNI brain
    vt = visualtools(subj='MNI_FSL', hemi='lh',
                     flag_use_annot=False,
                     brain_file='./sample_data/MNI_FSL152/FSL_MNI152_lh_pial.mat',
                     annot_file='./sample_data/MNI_FSL152/FSL_MNI152.lh.aparc.split_STG_MTG.annot')
    pl = vt.plot_elec_on_pial(elec_coords_mni,
                              elec_colors=elec_rgb,
                              radius=0.9,
                              flag_scalar_bar=False,
                              show=False)
    vt.flag_use_annot = False
    elec_coords_mni, _ = vt.project_elec_verts(elec_coords_mni,
                                               regions = elec_regions)
    vt.flag_use_annot = False
    pl = vt.plot_elec_on_pial(elec_coords_mni,
                              elec_colors=elec_rgb,
                              radius=0.9,
                              flag_scalar_bar=False,
                              show=False)
    pl.set_background(color='white')
    visualtools._save_all_views(pl, 'test')
    quit()


def demo_weighted_elecs_projection(make_gif=False):
    '''
    This example plots the grid electrodes from NY717 on the
    subject brain surface. Each electrode has an associated
    weight (random numbers here) that will be shown as a colored
    sphere or projected on the brain surface with a Gaussian
    kernel. To show the effect of region boundary limitations
    we assign random weights to only mSTG electrodes.
    '''
    # first we load the "coordinates.csv" file that has electrode
    # locations in subject T1 and MNI spaces.
    fn_csv = "./sample_data/sample_subjects/NY717/coordinates.csv"
    csv = pd.read_csv(fn_csv)
    # we select the grid electrodes for now
    se = np.arange(128) # first 128 electrodes are grid
    # we get the x,y,z locations for the elctrodes in T1 space
    x = np.array(csv['T1_x'][se])[:,None]
    y = np.array(csv['T1_y'][se])[:,None]
    z = np.array(csv['T1_z'][se])[:,None]
    elec_coords = np.hstack((x,y,z)) # a 128 x 3 matrix of coordinates
    # get the electrode locations in subject space
    # for each electrode
    elec_regions = csv['T1_AnatomicalRegion'][se] # str list of 128x1
    # for simulation we assign some random weights for
    # electrodes in mSTG
    # set seed for same result between different runs
    np.random.seed(65777382) # why this seed :-)
    # find index of mSTG electrodes and set rendom weights
    mSTG_inds = np.where(elec_regions=="mSTG")[0]
    elec_weights = np.zeros(elec_coords.shape[0]) # array of 128x1
    elec_weights[mSTG_inds] = np.random.rand(mSTG_inds.shape[0])
    # creatre the visualization object
    # we provide subject pial surface and annotation
    vt = visualtools(subj='NY717', hemi='lh',
                     flag_use_annot=True,
                     brain_file='./sample_data/sample_subjects/NY717/NY717_pial_surf.mat',
                     annot_file='./sample_data/sample_subjects/NY717/NY717_lh_aparc.annot')
    # create the plotting grid
    pl = pv.Plotter(shape=(2,2)) # a 2x2 grid to plot different visuals
    ## plot the electrodes on the brain -- Just electrodes with no weight
    pl.subplot(0,0) # select first subplot
    # use PlotElecOnBrain function from visualtools to plot the elecs
    # all electrodes are colored as red (1,0,0)
    # provide BrainPlotter argument to plot in the subplot
    vt.plot_elec_on_pial(elec_coords,
                         elec_colors=(0.0,0.0,0.0,1.0),
                         radius=1.0,
                         plotter=pl,
                         show = False)
    # set title for subplot
    pl.add_text("elecs on subj pial", font_size=10,
                font='times', position='upper_edge')
    ## plot electrodes weigthed by their colors
    pl.subplot(0,1)
    # set flag_use_annot to False to avoid region annotation
    # coloring when plotting the surface
    vt.flag_use_annot = False
    # set arguments for plottting the colorbar
    scalar_bar_args = dict(font_family='times', fmt='%.1f',
                           label_font_size=20)
    # use plot_elec_on_pial function from visualtools to plot the elecs
    # ElecColor input can be set to elec_weights to color each electrode
    # proportional to the weight value
    vt.plot_elec_on_pial(elec_coords,
                         elec_colors=elec_weights,
                         radius=1.0,
                         plotter=pl,
                         cmap='hot_r',
                         clim=(0,1),
                         scalar_bar_args=scalar_bar_args,
                         show=False)
    # set title for subplot
    pl.add_text("mSTG elecs weighted",
                font_size=10, font='times', position='upper_edge')
    ## plot projection of the weigths on the brain without
    # region boundary
    pl.subplot(1,0)
    # use PlotWeightedBrain from visualtools to project the electrode
    # weights on the brain.
    # setting flag_use_annot to False removes the dependancy on the
    # annotations and the boundary of Gaussian kernel is not limited
    vt.flag_use_annot = False
    vt.plot_weighted_surface(elec_coords,
                             elec_weights,
                             regions=elec_regions,
                             sigma = 2.5,
                             plotter=pl,
                             mode="simple",
                             cmap='hot_r',
                             show=False)
    pl.add_text("weights projected on brain w/o boundary constraint",
                font_size=10, font='times', position='upper_edge')
    # plot projection of the weigths on the brain with
    # region boundary
    pl.subplot(1,1)
    # use PlotWeightedBrain from visualtools to project the electrode
    # weights on the brain.
    # setting flag_use_annot to True makes the kernel depend on the
    # annotations and the boundary of Gaussian kernel is limited to
    # regions
    vt.flag_use_annot = True
    vt.plot_weighted_surface(elec_coords,
                             elec_weights,
                             regions=elec_regions,
                             sigma = 2.5,
                             plotter=pl,
                             mode="simple",
                             cmap='hot_r',
                             show=False)
    pl.add_text("weights projected on brain with boundary constraint",
                font_size=10, font='times', position='upper_edge')
    # link the views so the brains move together
    pl.link_views()
    # show the results
    cpos = pl.show(return_cpos=True, auto_close=False)
    if make_gif:
        ## create a gif of the results
        pl.camera_position = cpos
        pl.open_gif("linked.gif", subrectangles=True)
        nframe = 15
        for i in range(nframe):
            pl.camera.azimuth -= 3
            pl.camera.roll += 1
            pl.write_frame()
        for i in range(nframe):
            pl.camera.azimuth += 3
            pl.camera.roll -= 1
            pl.write_frame()
        pl.close()

if __name__ == "__main__":
    demo_weighted_elecs_projection(make_gif=False)
    quit()
    demo_elec_project_mni()
    demo_plot_brain_pial()
    demo_elces_on_subject_brain()
