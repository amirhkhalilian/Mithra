import sys, os, glob

import numpy as np
import scipy.io
import scipy.spatial as scp
from nibabel.freesurfer.io import read_annot
from nibabel.freesurfer.io import read_geometry
from matplotlib import pyplot as plt

import pyvista as pv
pv.set_plot_theme('dark')

class visualtools(object):
    """
    visualtools for plotting ECoG and sEEG electrodes on a pial surface.
    These tools can load triangulated pial surfaces from .pial files
    (Freesurfer format) or .mat files. To plot the brain surface and
    annotate regions, use .annot files from the Freesurfer output.
    --------------------
    Methods:
    --------------------
        - PlotBrainSurface
        - PlotElecOnBrain
        - PlotWeightedBrain
    """

    def __init__(self, subj, hemi,
                 brain_file=None,
                 annot_file=None,
                 flag_use_annot=False,
                 flag_merge_STG=False):
        super().__init__()
        self.subj = subj
        self.hemi = hemi
        if brain_file is None:
            raise ValueError('please provide path to pial surface path brain_file')
        self.brain_file = brain_file
        self.flag_use_annot = flag_use_annot
        self.flag_merge_STG = flag_merge_STG
        if self.flag_use_annot:
            if annot_file is None:
                raise ValueError('please provide annot_file')
        self.annot_file = annot_file

    @staticmethod
    def _get_trisurf_from_file(brain_file):
        '''
        loads triangulated pial surfaces from .pial or .mat file
        '''
        ext = os.path.splitext(brain_file)[1]
        if ext=='.mat':
            # load the brain file
            brain = scipy.io.loadmat(brain_file)
            verts = np.array(brain['coords'])
            faces = np.hstack((3*np.ones((brain['faces'].shape[0],1),dtype=np.int32),\
                               np.array(brain['faces']-1)))
        elif ext=='.pial':
            dl = read_geometry(brain_file)
            verts = np.array(dl[0])
            faces = np.hstack((3*np.ones((dl[1].shape[0],1),dtype=np.int32),\
                               np.array(dl[1])))
        else:
            raise ValueError(f'file extension {ext} not valid!')
        return verts, faces

    def _get_trisurf(self):
        '''
        loads triangulated pial surfaces from .pial or .mat file
        '''
        ext = os.path.splitext(self.brain_file)[1]
        if ext=='.mat':
            # load the brain file
            brain = scipy.io.loadmat(self.brain_file)
            verts = np.array(brain['coords'])
            faces = np.hstack((3*np.ones((brain['faces'].shape[0],1),dtype=np.int32),\
                               np.array(brain['faces']-1)))
        elif ext=='.pial':
            dl = read_geometry(self.brain_file)
            verts = np.array(dl[0])
            faces = np.hstack((3*np.ones((dl[1].shape[0],1),dtype=np.int32),\
                               np.array(dl[1])))
        else:
            raise ValueError(f'file extension {ext} not valid!')
        return verts, faces

    @staticmethod
    def _merge_annot_STG(annot):
        '''
        merges the three portions of STG and MTG annotations into one
        '''
        # change the verts assigned to mSTG and rSTG to cSTG
        cSTG_in_list = annot[2].index('cSTG'.encode('utf-8'))
        mSTG_in_list = annot[2].index('mSTG'.encode('utf-8'))
        mSTG_verts = np.where(annot[0]==mSTG_in_list)[0]
        annot[0][mSTG_verts] = cSTG_in_list
        rSTG_in_list = annot[2].index('rSTG'.encode('utf-8'))
        rSTG_verts = np.where(annot[0]==rSTG_in_list)[0]
        annot[0][rSTG_verts] = cSTG_in_list
        # chnage the verts assigned to mMTG and rMTG to cMTG
        cMTG_in_list = annot[2].index('cMTG'.encode('utf-8'))
        mMTG_in_list = annot[2].index('mMTG'.encode('utf-8'))
        mMTG_verts = np.where(annot[0]==mMTG_in_list)[0]
        annot[0][mMTG_verts] = cMTG_in_list
        rMTG_in_list = annot[2].index('rMTG'.encode('utf-8'))
        rMTG_verts = np.where(annot[0]==rMTG_in_list)[0]
        annot[0][rMTG_verts] = cMTG_in_list
        return annot

    @staticmethod
    def _project_elec_verts(elec_locs, verts,
                           flag_subject_anatomy = False,
                           elec_regions = None,
                           annot = None):
        '''
        function to project the elec on the pial surface
        elec_locs: array of electrode locations
        verts: array of pial surface verts
        annot[0] -- list of region index per vertex
        annot[2] -- region name
        '''
        elec_locs_proj = np.zeros_like(elec_locs)
        elec_rgb = 255*np.ones((elec_locs.shape[0],4))

        # Not efficient code...
        for ei in range(elec_locs.shape[0]):
            if flag_subject_anatomy:
                # use subject anatomy to find closest vert
                # within region
                # find region index in the list
                ri = annot[2].index(elec_regions[ei].encode('utf-8'))
                elec_rgb[ei,:] = annot[1][ri,0:4]
                # find vertices in region
                vi = np.argwhere(annot[0]==ri).squeeze()
                if len(vi)==0:
                    verts_r = verts
                else:
                    verts_r = verts[vi,:]
            else:
                verts_r = verts
            # compute the shortest dist
            dist = verts_r - elec_locs[ei,:]
            dist = np.sum(dist**2, axis=1)
            elec_locs_proj[ei,:] = verts_r[np.argmin(dist),:]
        elec_rgb /= 255
        elec_rgb[:,3] = 1.0
        return elec_locs_proj, elec_rgb

    def project_elec_verts(self,
                           elec_locs,
                           regions = None):
        verts, faces = self._get_trisurf()
        if self.flag_use_annot:
            annot = read_annot(self.annot_file)
        else:
            annot = None
        elec_locs_proj, elec_rgb = visualtools._project_elec_verts(
            elec_locs,
            verts,
            flag_subject_anatomy = self.flag_use_annot,
            elec_regions = regions,
            annot = annot)
        return elec_locs_proj, elec_rgb

    @staticmethod
    def _save_all_views(plotter, save_name):
        '''
        function to save the views
        '''
        plotter.camera_position = 'yz'
        plotter.camera.azimuth = 180
        plotter.save_graphic(save_name+'_sagittal.svg')
        plotter.camera_position = 'xz'
        plotter.camera.azimuth = 180
        plotter.save_graphic(save_name+'_coronal.svg')
        plotter.camera_position = 'xy'
        plotter.camera.azimuth = 0
        plotter.camera.elevation = 180
        plotter.save_graphic(save_name+'_horizontal.svg')


    def plot_pial_surface(self,
                          brain_color=(0.9,0.9,0.9,1.0),
                          flag_smooth = False,
                          smooth_iter = 10,
                          plotter = None,
                          flag_return_verts = False):
        '''
        function to plot a Brain pial surface
        '''
        verts, faces = self._get_trisurf()
        # convert brain faces and verts to PolyData in pyvista
        surf = pv.PolyData(verts, faces)
        if self.flag_use_annot:
            annot = read_annot(self.annot_file)
            brain_color = np.ones((verts.shape[0],4))
            if self.flag_merge_STG:
                annot = visualtools._merge_annot_STG(annot)
            brain_color[:,0:3] = annot[1][annot[0],0:3]/255.0
        if flag_smooth:
            surf = surf.smooth(n_iter=smooth_iter)
        # plot the poly data using
        if plotter is None:
            pl = pv.Plotter()
        else:
            pl = plotter
        if self.flag_use_annot:
            pl.add_mesh(surf, scalars=brain_color,
                        rgb = True,
                        use_transparency=True,
                        smooth_shading=True,
                        opacity=0.0)
        else:
            pl.add_mesh(surf, color=brain_color,
                        use_transparency=True,
                        smooth_shading=True,
                        opacity=1.0-brain_color[3])
        if self.hemi=='lh':
            pl.camera_position = 'yz'
            pl.camera.azimuth = 180
        elif self.hemi=='rh':
            pl.camera_position = 'yz'
        if flag_return_verts:
            return pl, verts
        else:
            return pl

    def plot_elec_on_pial(self,
                          elec_locs,
                          elec_colors=(1.0,0.0,0.0,1.0),
                          radius=1,
                          phi_resolution=10,
                          theta_resolution=10,
                          scale = False,
                          smooth_shading = True,
                          flag_scalar_bar = True,
                          scalar_bar_args = None,
                          cmap = 'coolwarm',
                          clim = None,
                          opacity = 1.0,
                          plotter = None,
                          brain_plotter = None,
                          brain_plotter_args = None,
                          show = True):

        if plotter is None:
            plotter = pv.Plotter()
        if brain_plotter is None:
            # plot the brain surface
            if brain_plotter_args is None:
                brain_plotter_args = dict(
                                      brain_color=(0.9,0.9,0.9,1.0),
                                      flag_smooth = False,
                                      smooth_iter = 10,
                                      plotter = plotter,
                                      flag_return_verts = False
                                      )
            pl = self.plot_pial_surface(**brain_plotter_args)
        else:
            pl = brain_plotter

        # locations to pv
        pdata = pv.PolyData(elec_locs)

        # single color vs colormap for elecs
        if type(elec_colors) is tuple:
            sphere = pv.Sphere(radius=radius,
                               phi_resolution=phi_resolution,
                               theta_resolution=theta_resolution)
            pc = pdata.glyph(scale=scale, geom=sphere, orient=False)
            pl.add_mesh(pc, color=elec_colors, smooth_shading=smooth_shading)
        elif type(elec_colors) is np.ndarray:
            if scalar_bar_args is None:
                scalar_bar_args = dict(title='electrode weights',
                                       font_family='times',
                                       fmt='%.1f',
                                       label_font_size=20,
                                       title_font_size=20)
            if elec_colors.ndim==1:
                pl.add_mesh(pdata, scalars=elec_colors,
                            point_size=15.0*radius,
                            render_points_as_spheres=True,
                            cmap = cmap,
                            clim = clim,
                            opacity = opacity,
                            smooth_shading = True,
                            show_scalar_bar=False)
                if flag_scalar_bar:
                    pl.add_scalar_bar(**scalar_bar_args)
            elif elec_colors.ndim==2 and elec_colors.shape[1]==4:
                pl.add_mesh(pdata, scalars=elec_colors,
                            point_size=15.0*radius,
                            render_points_as_spheres=True,
                            rgb=True,
                            smooth_shading=True,
                            show_scalar_bar=False)
                if flag_scalar_bar:
                    pl.add_scalar_bar(**scalar_bar_args)
            else:
                raise ValueError(f"Expected ElecColor to be a 4-tuple or ndarray")
        else:
            raise ValueError(f"Expected ElecColor to be a 4-tuple or ndarray")
        if show:
            pl.show()
        return pl

    def plot_weighted_surface(self,
                              elec_locs,
                              weights,
                              regions = None,
                              annot = None,
                              mode = "simple",
                              sigma = 2.0,
                              tol = 1e-3,
                              cmap = 'coolwarm',
                              clim = None,
                              plotter = None,
                              scalar_bar_args=None,
                              show=True):
        # load the brain
        verts, faces = self._get_trisurf()
        surf = pv.PolyData(verts, faces)
        color_values = np.zeros(surf.n_points)
        if mode == "normalize":
            norm_values = np.zeros(surf.n_points)

        # build a tree, query verts within 3*sigma of elec_locs
        tree = scp.cKDTree(verts)
        indices = tree.query_ball_point(elec_locs, 3*sigma)

        # loop to apply Gaussian filter
        for el, elec in enumerate(elec_locs):

            if not indices[el]:  # Skip if no points found
                continue

            if self.flag_use_annot: # Filter by region
                if annot is None:
                    annot = read_annot(self.annot_file)
                if regions is None:
                    raise ValueError('provide regions for elec_locs')
                region_id = annot[2].index(regions[el].encode('utf-8'))
                valid_indices = np.asarray(indices[el])[annot[0][indices[el]] == region_id]
            else:
                valid_indices = indices[el]

            if np.asarray(valid_indices).size == 0:  # Skip if no valid points
                continue

            # Apply Gaussian weights
            distances = np.linalg.norm(verts[valid_indices] - elec, axis=1)
            _gauss = np.exp(-0.5 * (distances / sigma)**2)

            # Accumulate weights into the surface values
            color_values[valid_indices] += weights[el] * _gauss
            if mode=="normalize":
                norm_values[valid_indices] += _gauss

        if mode=="normalize": # normalize the color_values
            color_values /= norm_values
            color_values[norm_values<=tol] = 0.0

        if plotter is None:
            pl = pv.Plotter()
        else:
            pl = plotter
        if scalar_bar_args is None:
            scalar_bar_args = dict(title='weighted surface',
                                   font_family='times',
                                   fmt='%.1f',
                                   label_font_size = 20,
                                   title_font_size = 20)

        pl.add_mesh(surf, scalars=color_values,\
                    cmap=cmap, clim=clim, smooth_shading=True,\
                    show_scalar_bar=False)

        pl.add_scalar_bar(**scalar_bar_args)

        if self.hemi=='lh':
            pl.camera_position = 'yz'
            pl.camera.azimuth = 180
        elif self.hemi=='rh':
            pl.camera_position = 'yz'

        if show:
            pl.show()
        return pl

if __name__=="__main__":
    pass

