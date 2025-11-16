function [fs_coords, fs_inds, hemi] = T1_to_fsaverage(subj, coords, varargin)
    % function to convert subject coordinates to fs-average coordinates
    % Inputs:
    %   - subj: subject "NY***" number
    %   - coords: struct containing
    %       - coords.MNI: mni coordinates
    %       - coords.T1:  T1 coordinates
    %       - coords.labels: electrode labels to remove depth

    % parse inputs 
    p = inputParser;
    p.KeepUnmatched = false;
    path_fsaverage_DEFAULT = './visualization-tools/SampleData/fsaverage';
    addParameter(p, 'path_fsaverage', path_fsaverage_DEFAULT, @(s) isdir(s));
    path_fs_recon_DEFAULT = './';
    addParameter(p, 'path_fs_recon', path_fs_recon_DEFAULT, @(s) isdir(s));
    path_subj_DEFAULT = [];
    addParameter(p, 'path_subj', path_subj_DEFAULT);
    addParameter(p, 'flag_Debug_mode', true);
    addParameter(p, 'flag_auto_rm_depth', true);
    % parse
    parse(p, varargin{:});
    % set values
    path_fsaverage  = p.Results.path_fsaverage;
    path_fs_recon   = p.Results.path_fs_recon;
    path_subj       = p.Results.path_subj;
    flag_Debug_mode = p.Results.flag_Debug_mode;
    flag_auto_rm_depth = p.Results.flag_auto_rm_depth;

    % check for subject auto recon dir
    if isempty(path_subj) % set path automatically
        path_subj = fullfile(path_fs_recon, sprintf('*%s*',subj));
        dir_with_subj_name = dir(path_subj);
        dirFlags = [dir_with_subj_name.isdir];
        dir_with_subj_name = dir_with_subj_name(dirFlags);
        if length(dir_with_subj_name)>1
            error('failed auto fill for subject recon path, multiple directories exist!... consider inputing the correct path in path_subj');
        end
        path_subj = fullfile(dir_with_subj_name(1).folder,...
                             dir_with_subj_name(1).name,...
                             'surf');
        if ~isdir(path_subj)
            error('failed auto fill for subject recon path, surf directory not found!');
        end
    end


    % set fsaverage files
    if ~isdir(path_fsaverage)
        error('failed to loacte fsavaverage recon path');
    end
    fs_lh_pial = fullfile(path_fsaverage, 'lh.pial');
    fs_rh_pial = fullfile(path_fsaverage, 'rh.pial');

    fs_lh_sphere = fullfile(path_fsaverage, 'lh.sphere.reg');
    fs_rh_sphere = fullfile(path_fsaverage, 'rh.sphere.reg');

    % set data path for subj
    if ~isdir(path_subj)
        error('failed to loacte subject recon path');
    end
    subj_lh_pial = fullfile(path_subj, 'lh.pial');
    subj_rh_pial = fullfile(path_subj, 'rh.pial');

    subj_lh_sphere = fullfile(path_subj, 'lh.sphere.reg');
    subj_rh_sphere = fullfile(path_subj, 'rh.sphere.reg');

    % find index of elecs on each hemisphere
    lh_ind = find(coords.MNI(:,1)<=0);
    rh_ind = find(coords.MNI(:,1)>0);

    if flag_auto_rm_depth
        % find depth electrodes (assuming labels starts with D)
        depth_ind = find(startsWith(coords.labels, 'd', 'IgnoreCase', true));
        surf_ind = setdiff([1:size(coords.MNI,1)], depth_ind);
    else
        surf_ind = [1:size(coords.MNI,1)];
    end

    % fs-average coords place holder
    fs_coords = nan(size(coords.T1));

    if nargout>=2
        fs_inds = nan(size(coords.T1,1),1);
    end

    % find indices on each hemisphere that are surface
    surf_lh_ind = intersect(surf_ind, lh_ind);
    surf_rh_ind = intersect(surf_ind, rh_ind);
    if nargout>=3
        hemi = nan(size(coords.T1,1),1);
        hemi(surf_lh_ind) = 1; % 1 for lh
        hemi(surf_rh_ind) = 2; % 2 for rh
    end

    % create fsaveragetools object for each hemisphere
    fstools_lh = fsaveragetools('HS', 'lh',...
                                'fn_fsaverage_pial', fs_lh_pial,...
                                'fn_fsaverage_sphere', fs_lh_sphere,...
                                'fn_subj_pial', subj_lh_pial,...
                                'fn_subj_sphere', subj_lh_sphere);

    fstools_rh = fsaveragetools('HS', 'rh',...
                                'fn_fsaverage_pial', fs_rh_pial,...
                                'fn_fsaverage_sphere', fs_rh_sphere,...
                                'fn_subj_pial', subj_rh_pial,...
                                'fn_subj_sphere', subj_rh_sphere);

    % convert coords using convert_T1_to_fsaverage
    [temp_coords_lh, fs_lh_vert_ind] = fstools_lh.convert_T1_to_fsaverage(coords.T1(surf_lh_ind,:));
    fs_coords(surf_lh_ind,:) = temp_coords_lh;
    if nargout>=2 & ~isempty(surf_lh_ind)
        fs_inds(surf_lh_ind) = fs_lh_vert_ind;
    end

    [temp_coords_rh, fs_rh_vert_ind] = fstools_rh.convert_T1_to_fsaverage(coords.T1(surf_rh_ind,:));
    fs_coords(surf_rh_ind,:) = temp_coords_rh;
    if nargout>=2 & ~isempty(surf_rh_ind)
        fs_inds(surf_rh_ind) = fs_rh_vert_ind;
    end

    % for debug mode plot the results
    if flag_Debug_mode
        % make a new figure
        SCR_SZ = get(0, 'Screensize');
        fig = figure('Position',SCR_SZ);

        % create a visualization object for subject surface
        VT = visualtools('Subj', subj, 'HS', 'lh',...
                         'flag_UseAnnots', false,...
                         'BrainFile', subj_lh_pial);
        ax1 = subplot(2,2,1);
        ElecColor_lh = tab10(length(lh_ind));
        VT.PlotElecOnBrain(coords.T1(lh_ind,:), 'ElecColor', ElecColor_lh, 'BrainColor', [0.9,0.9,0.9],...
                           'flag_AddFigure', false, 'radius',1.5, 'FaceAlpha', 0.9);
        title('Elecs on Subj Brain lh');

        % create a visualization object for subject surface
        VT = visualtools('Subj', subj, 'HS', 'rh',...
                         'flag_UseAnnots', false,...
                         'BrainFile', subj_rh_pial);
        ax3 = subplot(2,2,3);
        ElecColor_rh = tab10(length(rh_ind));
        VT.PlotElecOnBrain(coords.T1(rh_ind,:), 'ElecColor', ElecColor_rh, 'BrainColor', [0.9,0.9,0.9],...
                           'flag_AddFigure', false, 'radius',1.5, 'FaceAlpha', 0.9);
        title('Elecs on Subj Brain rh');

        % create a visualization object for fs-average surface
        VT2 = visualtools('Subj', 'fs-average', 'HS', 'lh',...
                         'flag_UseAnnots', false,...
                         'BrainFile', fs_lh_pial);
        ax2 = subplot(2,2,2);
        VT2.PlotElecOnBrain(fs_coords(lh_ind,:), 'ElecColor', ElecColor_lh, 'BrainColor', [0.9,0.9,0.9],...
                           'flag_AddFigure', false, 'radius',1.5, 'FaceAlpha', 0.9);
        title('Elecs on fsaverage Brain lh');

        % create a visualization object for fs-average surface
        VT2 = visualtools('Subj', 'fs-average', 'HS', 'rh',...
                         'flag_UseAnnots', false,...
                         'BrainFile', fs_rh_pial);
        ax4 = subplot(2,2,4);
        VT2.PlotElecOnBrain(fs_coords(rh_ind,:), 'ElecColor', ElecColor_rh, 'BrainColor', [0.9,0.9,0.9],...
                           'flag_AddFigure', false, 'radius',1.5, 'FaceAlpha', 0.9);
        title('Elecs on fsaverage Brain lh');

        % link all axis
        linkprop([ax3,ax4], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
        linkprop([ax1,ax2], {'CameraUpVector', 'CameraPosition', 'CameraTarget'});
    end

end


