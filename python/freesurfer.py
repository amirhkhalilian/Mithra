import sys, os, glob

import numpy as np
import scipy.io
from nibabel.freesurfer.io import read_annot
from nibabel.freesurfer.io import read_geometry
from nibabel.freesurfer.io import read_morph_data


def _convert_subject_T1_from_loc_to_fsaverage(subj_T1coords, HS,
                                              root_path=None,
                                              fsaverage_path=None,
                                              subj_path=None,
                                              flag_write=False,
                                              flag_visualize=False):
    '''
    function to load subject T1 locations from loc and covert
    to fs-average.
    This function needs the 'surf' directory in freesurfer output
    and uses the {HS}.pial and {HS}.sphere.reg
    '''
    pass

def get_mesh(file_name):
    dl = read_geometry(file_name)
    verts = np.array(dl[0])
    faces = np.hstack((3*np.ones((dl[1].shape[0],1),dtype=np.int32),\
                       np.array(dl[1])))
    return verts, faces

def morph_to_color(morph):
    color = np.zeros((len(morph), 3))
    color[morph >= 0, :] = [0.5, 0.5, 0.5]
    color[morph < 0, :] = [0.65, 0.65, 0.65]
    return color

def read_inflated(fn_inflated, fn_curv=None):
    verts, faces = get_mesh(fn_inflated)
    if fn_curv is not None:
        curv = read_morph_data(fn_curv)
        color = morph_to_color(curv)
        return verts, faces, color
    else:
        return verts, faces

# Handle freesurfer read curvature --- needed for .curv and .sulc files
def _fread3(fobj):
    """Read 3 bytes and adjust.
    Modified from mne:
    https://raw.githubusercontent.com/mne-tools/mne-python/maint/1.3/mne/surface.py
    """
    b1, b2, b3 = np.fromfile(fobj, ">u1", 3)
    return (b1 << 16) + (b2 << 8) + b3


def _fread3_many(fobj, n):
    """Read 3-byte ints from an open binary file object.
    Modified from mne:
    https://raw.githubusercontent.com/mne-tools/mne-python/maint/1.3/mne/surface.py
    """
    b1, b2, b3 = np.fromfile(fobj, ">u1",
                             3 * n).reshape(-1, 3).astype(np.int64).T
    return (b1 << 16) + (b2 << 8) + b3


def read_curvature(filepath, binary=True):
    """Load in curvature values from the ?h.curv file.
    and ?h.sulc files
    Modified from mne:
    https://raw.githubusercontent.com/mne-tools/mne-python/maint/1.3/mne/surface.py

    Parameters
    ----------
    filepath : str
        Input path to the .curv file.
    binary : bool
        Specify if the output array is to hold binary values. Defaults to True.

    Returns
    -------
    curv : array, shape=(n_vertices,)
        The curvature values loaded from the user given file.
    """
    with open(filepath, "rb") as fobj:
        magic = _fread3(fobj)
        if magic == 16777215:
            vnum = np.fromfile(fobj, ">i4", 3)[0]
            curv = np.fromfile(fobj, ">f4", vnum)
        else:
            vnum = magic
            _fread3(fobj)
            curv = np.fromfile(fobj, ">i2", vnum) / 100
    if binary:
        return 1 - np.array(curv != 0, np.int64)
    else:
        return curv


class fsaveragetools(object):
    """
    funcitonal tools for converting electrode locations from subject T1
    to fs-average space. The general approach is to find the closest vertex
    in sphere space.
    -------------------------
    Methods
    -------------------------
    """
    def __init__(self, Subj, HS, *args, **kwargs):
        self.Subj = Subj
        self.HS = HS
        if 'fsaverage_pial' in kwargs:
            self.fn_fs_pial = kwargs.get('fsaverage_pial')
        else:
            self.fn_fs_pial = f'./SampleData/fsaverage/{HS}.pial'
        if 'fsaverage_sphere' in kwargs:
            self.fn_fs_sphr = kwargs.get('fsaverage_sphere')
        else:
            self.fn_fs_sphr = f'./SampleData/fsaverage/{HS}.sphere.reg'
        if 'subj_pial' in kwargs:
            self.fn_subj_pial = kwargs.get('subj_pial')
        else:
            self.fn_subj_pial = f'./SampleData/{Subj}/{Subj}_pial_surf.mat'
        if 'subj_sphere' in kwargs:
            self.fn_subj_sphr = kwargs.get('subj_sphere')
        else:
            self.fn_subj_sphr = f'./SampleData/{Subj}/{Subj}_sphere.reg'
        self.fs_pial_verts = None
        self.fs_pial_faces = None
        self.fs_sphr_verts = None
        self.fs_sphr_faces = None
        self.subj_pial_verts = None
        self.subj_pial_faces = None
        self.subj_sphr_verts = None
        self.subj_sphr_faces = None
        return

    def load_fs_pial(self):
        ''' loads the fsaverage pial surface
        given the file given in self.fn_fs_pial
        '''
        ext = os.path.splitext(self.fn_fs_pial)[1]
        if ext=='.pial':
            dl = read_geometry(self.fn_fs_pial)
            self.fs_pial_verts = np.array(dl[0])
            self.fs_pial_faces = np.hstack((3*np.ones((dl[1].shape[0],1),dtype=np.int32),\
                                           np.array(dl[1])))
        elif ext=='.mat':
            # load the brain file
            Brain = scipy.io.loadmat(self.fn_fs_pial)
            # convert brain faces and verts to PolyData in pyvista
            self.fs_pial_verts = np.array(Brain['coords'])
            self.fs_pial_faces = np.hstack((3*np.ones((Brain['faces'].shape[0],1),dtype=np.int32),\
                                           np.array(Brain['faces']-1)))
        else:
            raise ValueError(f'file extension {ext} not supported!');
        return

    def load_subj_pial(self):
        ''' loads the subject pial surface
        given the file given in self.fn_subj_pial
        '''
        ext = os.path.splitext(self.fn_subj_pial)[1]
        if ext=='.pial':
            dl = read_geometry(self.fn_subj_pial)
            self.subj_pial_verts = np.array(dl[0])
            self.subj_pial_faces = np.hstack((3*np.ones((dl[1].shape[0],1),dtype=np.int32),\
                                           np.array(dl[1])))
        elif ext=='.mat':
            # load the brain file
            Brain = scipy.io.loadmat(self.fn_subj_pial)
            # convert brain faces and verts to PolyData in pyvista
            self.subj_pial_verts = np.array(Brain['coords'])
            self.subj_pial_faces = np.hstack((3*np.ones((Brain['faces'].shape[0],1),dtype=np.int32),\
                                           np.array(Brain['faces']-1)))
        else:
            raise ValueError(f'file extension {ext} not supported!');
        return

    def load_fs_sphere(self):
        ''' loads the fsaverage sphere surface
        given the file given in self.fn_fs_sphere
        '''
        ext = os.path.splitext(self.fn_fs_sphr)[1]
        if ext=='.reg':
            dl = read_geometry(self.fn_fs_sphr)
            self.fs_sphr_verts = np.array(dl[0])
            self.fs_sphr_faces = np.hstack((3*np.ones((dl[1].shape[0],1),dtype=np.int32),\
                                           np.array(dl[1])))
        elif ext=='.mat':
            # load the brain file
            Brain = scipy.io.loadmat(self.fn_fs_sphr)
            # convert brain faces and verts to PolyData in pyvista
            self.fs_sphr_verts = np.array(Brain['coords'])
            self.fs_sphr_faces = np.hstack((3*np.ones((Brain['faces'].shape[0],1),dtype=np.int32),\
                                           np.array(Brain['faces']-1)))
        else:
            raise ValueError(f'file extension {ext} not supported!');
        return

    def load_subj_shpere(self):
        ''' loads the subject shpere surface
        given the file given in self.fn_subj_sphr
        '''
        ext = os.path.splitext(self.fn_subj_sphr)[1]
        if ext=='.reg':
            dl = read_geometry(self.fn_subj_sphr)
            self.subj_sphr_verts = np.array(dl[0])
            self.subj_sphr_faces = np.hstack((3*np.ones((dl[1].shape[0],1),dtype=np.int32),\
                                           np.array(dl[1])))
        elif ext=='.mat':
            # load the brain file
            Brain = scipy.io.loadmat(self.fn_subj_sphr)
            # convert brain faces and verts to PolyData in pyvista
            self.subj_sphr_verts = np.array(Brain['coords'])
            self.subj_sphr_faces = np.hstack((3*np.ones((Brain['faces'].shape[0],1),dtype=np.int32),\
                                           np.array(Brain['faces']-1)))
        else:
            raise ValueError(f'file extension {ext} not supported!');
        return

    def convert_T1_to_fsaverage(self, ElecLocs, *args, **kwargs):
        '''
        converts the given coordinates in subject T1 to fsaverage
        by first finding the closest vertex on subjuct pial
        then converting the subject pial to subject sphere
        then find the closest fsaverage vert on sphere
        then convert fsaverage sphere to fsaverage pail
        '''
        # check if the pial and sphere verts are loaded properly
        # if not loaded we call the load function
        # set these object properties before call to avoid load
        # or provide in kwargs
        if 'fs_pial_verts' in kwargs:
            self.fs_pial_verts = kwargs.get('fs_pial_verts')
        else:
            if self.fs_pial_verts is None:
                self.load_fs_pial()
        if 'fs_sphr_verts' in kwargs:
            self.fs_sphr_verts = kwargs.get('fs_sphr_verts')
        else:
            if self.fs_sphr_verts is None:
                self.load_fs_sphere()
        if 'subj_pial_verts' in kwargs:
            self.subj_pial_verts = kwargs.get('subj_pial_verts')
        else:
            if self.subj_pial_verts is None:
                self.load_subj_pial()
        if 'subj_sphr_verts' in kwargs:
            self.subj_sphr_verts = kwargs.get('subj_sphr_verts')
        else:
            if self.subj_sphr_verts is None:
                self.load_subj_shpere()
        # init the fs locations
        fs_ElecLocs = np.zeros_like(ElecLocs)
        for i in range(ElecLocs.shape[0]):
            # get closest subj pial vert to ElecLoc
            subj_ind = np.argmin(np.sum((self.subj_pial_verts-ElecLocs[i,:])**2,axis=1))
            # get closest fs sphere to subj sphere in subj_ind
            fs_ind = np.argmin(np.sum((self.fs_sphr_verts-self.subj_sphr_verts[subj_ind,:])**2,axis=1))
            # set the fs location based on this index
            fs_ElecLocs[i,:] = self.fs_pial_verts[fs_ind,:]
        return fs_ElecLocs
