clc; clear; close all;
% demo file for matlab visualization tools

%% Plot NY717 grid electrodes on the T1 brain surface
% first we load the "coordinates.csv" file that has electrode
% location in subject T1 space.
fn_csv = './sample_data/sample_subjects/NY717/coordinates.csv';
csv = readtable(fn_csv);
% we select the grid electrodes 
se = 1:128; % first 128 electrodes are grid
% we get T1 x,y,z coordinates
elec_coords = table2array(csv(se,6:8)); % 6:8 is for T1... matrix of size 128x3
% create a visualization object
VT = visualtools('Subj', 'NY717', 'HS', 'lh',...
                 'BrainFile', './sample_data/sample_subjects/NY717/NY717_pial_surf.mat',...
                 'AnnotFile', []);
% plot electrodes on the brain
VT.PlotElecOnBrain(elec_coords, 'ElecColor',[1.0,0.0,0.0],...
                   'radius',1.0);

%% Plot NY717 grid electrodes on the MNI brain surface
% first we load the "coordinates.csv" file that has electrode
% location in subject T1 space.
fn_csv = './sample_data/sample_subjects/NY717/coordinates.csv';
csv = readtable(fn_csv);
% we select the grid electrodes 
se = 1:128; % first 128 electrodes are grid
% we get T1 x,y,z coordinates
elec_coords = table2array(csv(se,2:4)); % 2:4 is for MNI... matrix of size 128x3
VT = visualtools('Subj', 'MNI-FS', 'HS', 'lh',...
                 'flag_UseAnnots', true,...
                 'BrainFile','./sample_data/MNI_FSL152/FSL_MNI152_lh_pial.mat',...
                 'AnnotFile','./sample_data/MNI_FSL152/FSL_MNI152.lh.aparc.split_STG_MTG.annot');
% plot electrodes on the brain
VT.PlotElecOnBrain(elec_coords, 'ElecColor', [0.0,0.0,0.0],...
                   'radius',1.0);

%% Plot Weighted Electrode Projections
% This example plots the grid electrodes from NY717 on the
% subject brain surface. Each electrode has an associated
% weight (random numbers here) that will be shown as a colored
% sphere or projected on the brain surface with a Gaussian
% kernel. To show the effect of region boundary limitations
% we assign random weights to only mSTG electrodes.

% load T1 coordinates
fn_csv = './sample_data/sample_subjects/NY717/coordinates.csv';
csv = readtable(fn_csv);
% we select the grid electrodes 
se = [1:128]'; % first 128 electrodes are grid
% we get T1 x,y,z coordinates
elec_coords = table2array(csv(se,6:8)); % 6:8 is for T1... matrix of size 128x3
elec_regions = csv{se,11}; % T1_AnatomicalRegion for electrodes ... cell of 128x1

% for simulation we assign some random weights for
% electrodes in mSTG
% set seed
rng(65777382) % why this seed? :-)
mSTG_inds = find(contains(elec_regions, 'mSTG'));
elec_weights = zeros(size(se));
elec_weights(mSTG_inds) = rand(size(mSTG_inds));

% create a visualization object
VT = visualtools('Subj', 'NY717', 'HS', 'lh',...
                 'flag_UseAnnots', true,...
                 'flag_MergeSTG', false,...
                 'BrainFile', './sample_data/sample_subjects/NY717/NY717_pial_surf.mat',...
                 'AnnotFile', './sample_data/sample_subjects/NY717/NY717_lh_aparc.annot');
% new figure
SCR_SZ = get(0, 'Screensize');
fig = figure('Position',SCR_SZ);

% plot the electrodes on the brain -- Just electrodes with no weigh
ax1 = subplot(2,2,1);
VT.PlotElecOnBrain(elec_coords, 'ElecColor', [0.0,0.0,0.0],...
                   'flag_AddFigure', false, 'radius',1.0);
title('Elecs on Subj Brain');

% plot the electrodes colored by the elec_weights
ax2 = subplot(2,2,2);
VT.flag_UseAnnots = false;
VT.PlotElecOnBrain(elec_coords, 'ElecColor', elec_weights,...
                   'flag_AddFigure', false, 'radius',1.0,...
                   'cmap', 'hot_r', 'clim', [0,1]);
title('mSTG Elecs Weighted');

% plot projection of the weigths on the brain without
% region boundary
ax3 = subplot(2,2,3);
VT.flag_UseAnnots = false;
VT.PlotWeightedBrain(elec_coords, elec_weights,...
                     'Mode', 'SimpleProjectNorm',...
                     'Sigma', 3,...
                     'cmap', 'hot_r', 'clim', [0,1]);
title('Weights projected on brain w/o boundary constraint');

% plot projection of the weigths on the brain without
% region boundary
ax4 = subplot(2,2,4);
VT.flag_UseAnnots = true;
VT.PlotWeightedBrain(elec_coords, elec_weights,...
                     'ElecRegions', elec_regions,...
                     'Mode', 'SimpleProjectNorm',...
                     'Sigma', 3,...
                     'cmap', 'hot_r', 'clim', [0,1]);
title('Weights projected on brain with boundary constraint');

% link the axis
linkprop([ax1,ax2,ax3,ax4], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});

%% Plot NY717 grid electrodes on subject pial surface and fs-average pial surface
% this demo shows how to convert the T1 coordinates to fs-average
% and plot on the fs-average brain

% load T1 coordinates
fn_csv = './sample_data/sample_subjects/NY717/coordinates.csv';
csv = readtable(fn_csv);
% we select the grid electrodes 
se = [1:128]'; % first 128 electrodes are grid
% we get T1 x,y,z coordinates
elec_coords = table2array(csv(se,6:8)); % 6:8 is for T1... matrix of size 128x3

% get an fs-average tools object to convert T1 coords
% to fs-average
FST = fsaveragetools('HS', 'lh',...
                     'fn_fsaverage_pial', './sample_data/fsaverage/lh.pial',...
                     'fn_fsaverage_sphere', './sample_data/fsaverage/lh.sphere.reg',...
                     'fn_subj_pial', './sample_data/sample_subjects/NY717/NY717.lh.pial',...
                     'fn_subj_sphere', './sample_data/sample_subjects/NY717/NY717.lh.sphere.reg');
fs_elec_coords = FST.convert_T1_to_fsaverage(elec_coords);

% new figure
SCR_SZ = get(0, 'Screensize');
fig = figure('Position',SCR_SZ);

% create a visualization object for subject surface
VT = visualtools('Subj', 'NY717', 'HS', 'lh',...
                 'flag_UseAnnots', true,...
                 'BrainFile', './sample_data/sample_subjects/NY717/NY717_pial_surf.mat',...
                 'AnnotFile', './sample_data/sample_subjects/NY717/NY717_lh_aparc.annot');
ax1 = subplot(1,2,1);
VT.PlotElecOnBrain(elec_coords, 'ElecColor', [1.0,1.0,1.0],...
                   'flag_AddFigure', false, 'radius',1.0);
title('Elecs on Subj Brain');

% create a visualization object for fs-average surface
VT2 = visualtools('Subj', 'fs-average', 'HS', 'lh',...
                 'flag_UseAnnots', true,...
                 'BrainFile', './sample_data/fsaverage/lh.pial',...
                 'AnnotFile', './sample_data/fsaverage/lh.aparc.annot');
ax2 = subplot(1,2,2);
VT2.PlotElecOnBrain(fs_elec_coords, 'ElecColor', [1.0,1.0,1.0],...
                   'flag_AddFigure', false, 'radius',1.0);
title('Elecs on Subj Brain');

% link all axis
linkprop([ax1,ax2], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});

%% Plot NY717 grid electrodes on the T1 brain surface with random connectivity 
% first we load the "coordinates.csv" file that has electrode
% location in subject T1 space.
fn_csv = './sample_data/sample_subjects/NY717/coordinates.csv';
csv = readtable(fn_csv);
% we select the grid electrodes 
se = 1:128; % first 128 electrodes are grid
% we get T1 x,y,z coordinates
elec_coords = table2array(csv(se,6:8)); % 6:8 is for T1... matrix of size 128x3

elec_coords = [elec_coords; elec_coords(128,:)];
elec_coords(129,1) = elec_coords(129,1)+20;

% create an adjacency matrix
Adj = zeros(length(se)+1, length(se)+1);
Adj(1,8) = 0.8;
Adj(8,1) = 0.8;
Adj(1,4) = 0.9;
Adj(4,1) = 0.9;
Adj(1,48) = 1;
Adj(48,1) = 1;
Adj(3,32) = 0.1;
Adj(32,3) = 0.1;

Adj(3,129) = 1;
Adj(129,3) = 1;

Adj(3,128) = 1;
Adj(128,3) = 1;


% create a visualization object
VT = visualtools('Subj', 'NY717', 'HS', 'lh',...
                 'BrainFile', './sample_data/sample_subjects/NY717/NY717_pial_surf.mat',...
                 'AnnotFile', []);
% plot connectivity on the brain
VT.PlotElecConnectivityOnBrain(elec_coords, Adj, 'ElecColor',[0.0,0.0,0.0],...
                   'radius',1.0);

%% Projecting the electrodes on the pial surface with closest vert on same area
fn_csv = './sample_data/sample_subjects/NY717/coordinates.csv';
csv = readtable(fn_csv);
% we select the grid electrodes
se = 1:128; % first 128 electrodes are grid
% we get MNI x,y,z coordinates
elec_locs = table2array(csv(se,2:4)); % 2:4 is for MNI... matrix of size 128x3
elec_regions = csv{se,11};

AnnotFile='./sample_data/MNI_canonical/ch2_template.lh.aparc.split_STG_MTG.annot';
BrainFile='./sample_data/MNI_canonical/ch2_template_lh_pial_120519.mat';

[elec_locs_proj, elec_rgb] = project_elecs_pial(elec_locs, elec_regions, BrainFile, AnnotFile);

VT = visualtools('Subj', 'MNI', 'HS', 'lh',...
                 'flag_UseAnnots', false,...
                 'BrainFile', BrainFile,...
                 'AnnotFile', AnnotFile);
% plot electrodes on the brain
VT.PlotElecOnBrain(elec_locs, 'ElecColor', elec_rgb,...
                   'radius',2.0);
% plot electrodes on the brain
VT.PlotElecOnBrain(elec_locs_proj, 'ElecColor', elec_rgb,...
                   'radius',2.0);


%% Coloring regions on brian surface

% define visualization tools
VT = visualtools('Subj', 'MNI-FS', 'HS', 'lh',...
                 'flag_UseAnnots', true,...
                 'BrainFile','./sample_data/MNI_FSL152/FSL_MNI152_lh_pial.mat',...
                 'AnnotFile','./sample_data/MNI_FSL152/FSL_MNI152.lh.aparc.split_STG_MTG.annot');
% select name of areas you want
area_names = {'cSTG', 'precentral', 'parsopercularis', 'rostralmiddlefrontal'};
area_colors = [1 0 0;
               0 1 0;
               0 0 1;
               1 0 1];
VT.PlotRoIonSurface(area_names, area_colors);


%% Plotting freesurfer inflated brain

[sulc, ~] = fsaveragetools.read_curv( './sample_data/fsaverage/rh.curv');
% create a visualization object for fs-average surface
VT2 = visualtools('Subj', 'fs-average', 'HS', 'rh',...
                 'flag_UseAnnots', false,...
                 'BrainFile', './sample_data/fsaverage/rh.inflated',...
                 'AnnotFile', './sample_data/fsaverage/rh.aparc.annot');
%set color
color = zeros(length(sulc),3);
color(sulc>=0,:) = repmat([0.4,0.4,0.4], sum(sulc>=0),1);
color(sulc<0,:) = repmat([0.55,0.55,0.55], sum(sulc<0),1);
VT2.PlotBrainSurface('BrainColor', color);

%% plotting the subcortical regions in MNI-FS

names = {'subcort_lHipp', 'subcort_lAmgd', 'subcort_lThal',...
         'subcort_rHipp', 'subcort_rAmgd', 'subcort_rThal'};
colors = {[220,216,20]/256; [103,255,255]/256; [0, 118, 14]/256};
colors = [colors; colors];
figure;
VT = visualtools('Subj', 'MNI-FS', 'HS', 'lh',...
                 'flag_UseAnnots', false,...
                 'BrainFile',[]);
for i = 1:length(names)
    VT.BrainFile = fullfile('./sample_data/MNI_FSL152/', [names{i},'.mat']);
    VT.PlotBrainSurface('flag_AddFigure', false,...
                        'BrainColor', colors{i},...
                        'FaceAlpha', 1.0);
    hold on;
end


