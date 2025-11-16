# Matlab implementation for electrode visualization

## Installation and running demos

You can add `matlab` directory to your Matlab path by running
```
addpath(genpath(./matlab))
```

## demos
Please review the [`demo.m`](./demo.m) file to see some of the functionality. The demo file shows plotting electrodes on subject or MNI pial surfaces as well as projection of electrode on the surface using a Gaussian kernel.

# Function Docs

## visualtools class

The class that implements visualization tools in Matlab. The functions here mostly use `tirsurf` for drawing the brain surfaces and electrodes. 

Properties:
- `obj.Subj`: subject ID. You can set to 'MNI-FS' to load the MNI brains automatically.
- `obj.HS`: indicator for which brain hemesphere is being plotted. can set to 'lh' or 'rh'. Default 'lh'.
- `obj.BrainFile`: name and address of the pial surface file in `*.mat` format (support for *.pial* will be added with fsaveragetools).
- `obj.AnnotFile`: name and address of the annotation file in `*.annot` format.
- `obj.flag_UseAnnots`: flag to indicate if we use the Annotations. Default False.
- `obj.flag_MergeSTG`: (NOT Implemented in Matlab yet) if this flag is set we merge rSTG, cSTG, and mSTG in the Annot file and Electorde regions into one region for plotting. Same with MTG. Default False. 

Methods:
- `__init__`: initialization
- `PlotBrainSurface`: function to plot a Brain surface
- `PlotElecOnBrain`: function to plot electrodes on a brain surface. each electrode shown as a small sphere.
- `PlotWeightedBrain`: function to plot a brain surface with electrode weights projected with a Gaussian kernel onto the surface.
- `read_annot`: function to read the freesurfer annotation files (from ntools). 

### PlotBrainSurface
The function used for plotting only the brain surface. This function loads the surface vertices and faces from `obj.BrainFile`. If `obj.flag_UseAnnots` is set to `True`, the function reads the RoI annotations from `obj.AnnotFile` using `read_annot` function.

Inputs:
- `obj`: visualtools object.
- `BrainColor`: 3-vector indicating the rgb values for brain surface color. Default `[0.7,0.7,0.7]`.
- `flag_AddFigure`: option to create a new figure. Default `True`.
- `FaceAlpha`: transparency value for brain surface. Default 1.

Returns:
- plotted figure.

### PlotElecOnBrain
The function used for plotting the electrodes on the brain surface. This function allows for plotting the electrodes as spheres with colors set to specific RGBA values or weighted based on the inputed weights. 

Inpute:
- `obj`: visualtools object.
- `ElecLoc`: required input of shape Nx3 for (x,y,z) electrode coordinates.
- `ElecColor`: if type is vector with (R,G,B) values for all electrodes. if  size is Nx1, interpreted as weights and colors are computed based on `cmap`. if type is `ndims==2`of size Nx3, interrupted as (RGB) values for each electrode. Default is (1,0,0) for red electrode spheres.
- `radius`: size of spheres. Default is 1.0
- `cmap`: kwarg passed to pyvista for colormap setting. Default `coolwarm`.
- `clim`: kwarg passed to pyvista for color limits. Default `None`.
- `flag_AddFigure`: option to create a new figure. Default `True`.
- `BrainColor`: 3-vector indicating the rgb values for brain surface color. Default `[0.7,0.7,0.7]`.
- `FaceAlpha`: transparency value for brain surface. Default 1.

Returns:
- plotted figure.

### PlotWeightedBrain
The function for plotting the brain surface such that the value associated with each electrode is projected on the surface by a Gaussian kernel.

Inputs:
- `obj`: visualtools object.
- `ElecLoc`: required input of shape Nx3 for (x,y,z) electrode coordinates.
- `Weights`: weighting values for each electrode vector of size Nx1.
- `ElecRegions`: the T1 regions assigned to each electrode cell of str. Default empty.
- `Sigma`: the Gaussian kernel spread. Default 2.0
- `Mode`: projection mode valid: `['SimpleProject', 'SimpleProjectNorm', 'NNProject', 'NNProjectNorm']`. Default `'SimpleProject'`.
- `FaceAlpha`: transparency value for brain surface. Default 1.

Returns:
- `color_weights`: weights per vertex on the brain surface

#### plotting modes
Let's assume we have $N$ electrodes with coordinates $(x_n,y_n,z_n)$ for $n\in\{1,\cdots,N\}$. Each electrode has an associated value $v_n$ that we want to project on the brain surface. In the cases when `self.flag_UseAnnots` is set to `True`, we also have the associated region for each electrode in Subject T1 space. Depending on the selected `Mode` and `self.flag_UseAnnots`, we compute the coloring weight for any location $(x,y,z)$ as follows:

- `Mode=="SimpleProject"` and `self.flag_UseAnnots==False`: In this case the color weight $c(x,y,z)$ can be written as $$c(x,y,z)=\sum_{n=1}^{N} v_n \times \exp\left(-\dfrac{\|(x,y,z)-(x_n,y_n,z_n)\|_F^2}{2\sigma^2}\right).$$
- `Mode=="SimpleProject"` and `self.flag_UseAnnots==True`: In this case the color weight $c(x,y,z)$ can be written as $$c(x,y,z)=\sum_{n=1}^{N} v_n \times \exp\left(-\dfrac{\|(x,y,z)-(x_n,y_n,z_n)\|_F^2}{2\sigma^2}\right) \times I(r(x,y,z)\in r_n)$$ where $I(\mathrm{condition})$ is the indicator function and is equal to 1 if condition is True and 0 otherwise.
- `Mode=="SimpleProject"` and `self.flag_UseAnnots==False`: In this case the color weight $c(x,y,z)$ can be written as $$c(x,y,z)=\dfrac{1}{Z} \sum_{n=1}^{N} v_n \times \exp\left(-\dfrac{\|(x,y,z)-(x_n,y_n,z_n)\|_F^2}{2\sigma^2}\right)$$ 

where $$Z = \sum_{n=1}^{N} \exp\left(-\dfrac{\|(x,y,z)-(x_n,y_n,z_n)\|_F^2}{2\sigma^2}\right).$$.
- `Mode=="SimpleProjectNorm"` and `self.flag_UseAnnots==True`: In this case the color weight $c(x,y,z)$ can be written as $$c(x,y,z)=\dfrac{1}{Z}\sum_{n=1}^{N} v_n \times \exp\left(-\dfrac{\|(x,y,z)-(x_n,y_n,z_n)\|_F^2}{2\sigma^2}\right) \times I(r(x,y,z)\in r_n).$$
