classdef visualtools
    % Visualization tools for Ecog data in matlab
    properties
        Subj = 'MNI-FS';
        HS = 'lh';
        flag_UseAnnots = false;
        flag_MergeSTG = false;
        RootDir = '.';
        BrainFile = ['.', filesep, 'SampleData', filesep, 'MNI-FS', ...
                     filesep, 'FSL_MNI152_lh_pial.mat'];
        AnnotFile = ['.', filesep, 'SampleData', filesep, 'MNI-FS', ...
                     filesep, 'FSL_MNI152.lh.aparc.split_STG_MTG.annot'];
    end

    methods
        function [obj] = visualtools(varargin)
            % init construction function
            % parser for the input parameters
            % default to properties
            p = inputParser;
            addParameter(p, 'Subj', obj.Subj);
            addParameter(p, 'HS', obj.HS);
            addParameter(p, 'flag_UseAnnots', obj.flag_UseAnnots);
            addParameter(p, 'flag_MergeSTG', obj.flag_MergeSTG);
            addParameter(p, 'BrainFile', obj.BrainFile);
            addParameter(p, 'RootDir', obj.RootDir);
            addParameter(p, 'AnnotFile', obj.AnnotFile);
            % parse the input
            parse(p, varargin{:});
            % set parsed to obj
            obj.Subj = p.Results.Subj;
            obj.HS = p.Results.HS;
            obj.RootDir = p.Results.RootDir;
            obj.BrainFile = p.Results.BrainFile;
            obj.AnnotFile = p.Results.AnnotFile;
            obj.flag_UseAnnots = p.Results.flag_UseAnnots;
            obj.flag_MergeSTG = p.Results.flag_MergeSTG;
            if strcmp(obj.HS, 'rh') && strcmp(obj.Subj, 'MNI-FS')
                % correcting the default file to rh
                obj.BrainFile = ['.', filesep, 'SampleData', filesep, 'MNI-FS', ...
                                 filesep, 'FSL_MNI152_rh_pial.mat'];
                obj.AnnotFile = ['.', filesep, 'SampleData', filesep, 'MNI-FS', ...
                                 filesep, 'FSL_MNI152.rh.aparc.split_STG_MTG.annot'];
            end
        end

        % Plot a brain
        handle = PlotBrainSurface(obj, varargin);

        handle = PlotRoIonSurface(obj, area_names, area_colors, varargin);

        % Plot Electrodes on a Brain
        handle = PlotElecOnBrain(obj, ElecLocs, varargin);

        % Plot Weighted brian surface 
        [handle, color_weights] = PlotWeightedBrain(obj, ElecLocs, Weights, varargin);

        % Plot elecs connectivty on the brain
        handle = PlotElecConnectivityOnBrain(obj, ElecLocs, Adjacency, varargin);
    end

    methods (Static)
        % read annotation file
        [vertices, label, colortable] = read_annot(filename);

        % custom colormaps from 
        % https://www.mathworks.com/matlabcentral/fileexchange/69470-custom-colormap
        J = customcolormap(positions, colors, m)

        % function to convert value to color
        function [col, rgb_inds, cmap] = value2color(value, c_map, c_range)
            if ~exist('c_map','var') || isempty(c_map)
                c_map = 'hot';
            end
            if ~exist('c_range','var') || isempty(c_range)
                c_range = [min(value,[],'all'), max(value,[],'all')];
            end
            if isstring(c_map) || ischar(c_map)
                if strcmp(c_map(end), 'r')
                    cmap = flipud(eval([c_map(1:end-2),'(256)']));
                else
                    cmap = eval([c_map,'(256)']);
                end
            else
                cmap = c_map;
            end
            % clip the out of range values
            value(value<c_range(1)) = c_range(1);
            value(value>c_range(2)) = c_range(2);
            % normalize the values
            rgb_inds = round(((value-c_range(1))/(c_range(2)-c_range(1)))*255)+1;
            rgb_inds(isnan(rgb_inds)) = 1;
            % set colors
            col = cmap(rgb_inds,:);
        end

    end

end
