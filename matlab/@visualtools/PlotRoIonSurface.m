function handle = PlotRoIonSurface(obj, area_names, area_colors, varargin)
    % function to plot a brain surface
    p = inputParser;
    addParameter(p, 'flag_AddFigure', true);
    addParameter(p, 'FaceAlpha', 1);
    addParameter(p, 'BrainColor', [0.7, 0.7, 0.7]);
    % parse inputs
    parse(p, varargin{:});
    flag_AddFigure = p.Results.flag_AddFigure;
    FaceAlpha = p.Results.FaceAlpha;
    BrainColor = p.Results.BrainColor;

    % load the brian
    if exist(obj.BrainFile, 'file');
        % load the brain file using the load_surface function
        % supports both freesurfer files and mat file.
        [verts, faces] = fsaveragetools.load_surface(obj.BrainFile);
    else
        error(['Brain file: ', obj.BrainFile, ' does not exist!']);
    end
    
    % set default brain_color
    brain_color = repmat(BrainColor(:)', [size(verts,1),1]);

    % load annotations
    if obj.flag_UseAnnots
        [~, annt_lbls, annt_tbls] = visualtools.read_annot(obj.AnnotFile);
        % loop over brain colors
        for i = 1:length(area_names)
            [~,area_in_list] = ismember(area_names{i}, annt_tbls.struct_names);
            verts_ind = find(annt_lbls==annt_tbls.table(area_in_list,5));
            brain_color(verts_ind,:) = repmat(area_colors(i,:),[length(verts_ind),1]);
        end
    end

    % plot the brain
    if flag_AddFigure; figure; end
    % plot the brain surface
    handle = trisurf(faces,...
                     verts(:,1),...
                     verts(:,2),...
                     verts(:,3),...
                     'FaceVertexCData', brain_color,...
                     'FaceColor', 'interp',...
                     'FaceAlpha', FaceAlpha);

    % make it pretty with lighting and rotate based on HS
    shading interp; lighting gouraud; material dull;
    if strcmp(obj.HS, 'lh')
        view([-1 0 0]);
        set(light,  'position', [-1 0 0], ...
            'color', 0.8 * [1 1 1]);
    elseif strcmp(obj.HS, 'rh')
        view([1 0 0]);
        set(light,  'position', [1 0 0], ...
            'color', 0.8 * [1 1 1]);
    else
        error('is it right or left?')
    end
    alpha(handle, FaceAlpha);
    axis image; axis off; hold on;

    % stretch to axis outer edge
    if flag_AddFigure
        outerpos = get(gca, 'OuterPosition');
        left = outerpos(1); bottom = outerpos(2);
        ax_width = outerpos(3); ax_height = outerpos(4);
        set(gca, 'Position', [left bottom ax_width ax_height]);
    end
end

