function handle = PlotBrainSurface(obj, varargin)
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
    % load annotations
    if obj.flag_UseAnnots
        [~, annt_lbls, annt_tbls] = visualtools.read_annot(obj.AnnotFile);
        if obj.flag_MergeSTG
            % change verts assigned to mSTG and rSTG to cSTG
            [~,cSTG_in_list] = ismember('cSTG', annt_tbls.struct_names);
            [~,mSTG_in_list] = ismember('mSTG', annt_tbls.struct_names);
            annt_lbls(annt_lbls==annt_tbls.table(mSTG_in_list,5)) = annt_tbls.table(cSTG_in_list,5);
            annt_tbls.table(mSTG_in_list,:) = annt_tbls.table(cSTG_in_list,:);
            annt_tbls.struct_names{mSTG_in_list} = 'cSTG';
            [~,rSTG_in_list] = ismember('rSTG', annt_tbls.struct_names);
            annt_lbls(annt_lbls==annt_tbls.table(rSTG_in_list,5)) = annt_tbls.table(cSTG_in_list,5);
            annt_tbls.table(rSTG_in_list,:) = annt_tbls.table(cSTG_in_list,:);
            annt_tbls.struct_names{rSTG_in_list} = 'cSTG';
            % change verts assigned to mMTG and rMTG to cMTG
            [~,cMTG_in_list] = ismember('cMTG', annt_tbls.struct_names);
            [~,mMTG_in_list] = ismember('mMTG', annt_tbls.struct_names);
            annt_lbls(annt_lbls==annt_tbls.table(mMTG_in_list,5)) = annt_tbls.table(cMTG_in_list,5);
            annt_tbls.table(mMTG_in_list,:) = annt_tbls.table(cMTG_in_list,:);
            annt_tbls.struct_names{mMTG_in_list} = 'cMTG';
            [~,rMTG_in_list] = ismember('rMTG', annt_tbls.struct_names);
            annt_lbls(annt_lbls==annt_tbls.table(rMTG_in_list,5)) = annt_tbls.table(cMTG_in_list,5);
            annt_tbls.table(rMTG_in_list,:) = annt_tbls.table(cMTG_in_list,:);
            annt_tbls.struct_names{rMTG_in_list} = 'cMTG';
        end
        [~, annot] = ismember(annt_lbls, annt_tbls.table(:,5));
        annot(annot==0) = 1;
        brain_color = annt_tbls.table(annot, 1:3)./255;
    else
        if all(size(BrainColor)==[1,3]) | all(size(BrainColor)==[3,1])
            brain_color = repmat(BrainColor(:)', [size(verts,1),1]);
        else
            brain_color = BrainColor;
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

