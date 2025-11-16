function [handle, color_weights] = PlotWeightedBrain(obj, ElecLocs, Weights, varargin)
    % function to visualize a weighted brain surface based on a Gaussian projection
    % of associated weights to each electrode

    % variables to change to varargin
    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'ElecRegions', []);
    addParameter(p, 'Sigma', 2.0);
    addParameter(p, 'Mode', 'SimpleProject');
    addParameter(p, 'tol_norm', 1e-3);
    addParameter(p, 'FaceAlpha', 1);
    addParameter(p, 'cmap', 'hot');
    addParameter(p, 'clim', []);
    addParameter(p, 'flag_AddFigure', false);
    addParameter(p, 'flag_rm_sparse', false);
    addParameter(p, 'sparse_thresh', 20);
    addParameter(p, 'sparse_ratio', 0.1);
    % parse inputs
    parse(p, varargin{:});
    ElecRegions = p.Results.ElecRegions;
    Sigma = p.Results.Sigma;
    Mode = p.Results.Mode;
    tol_norm = p.Results.tol_norm;
    FaceAlpha = p.Results.FaceAlpha;
    cmap = p.Results.cmap;
    clim = p.Results.clim;
    flag_AddFigure = p.Results.flag_AddFigure;
    flag_rm_sparse = p.Results.flag_rm_sparse;
    sparse_thresh = p.Results.sparse_thresh;
    sparse_ratio = p.Results.sparse_ratio;

    % get brain verts and faces
    if exist(obj.BrainFile, 'file');
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
        if isempty(ElecRegions)
            error('flag_UseAnnots is set True but ElecRegions is not provided!');
        end
        if obj.flag_MergeSTG
            % change elec areas from mSTG, rSTG to cSTG to match the annotation
            mSTG_in_elec = find(contains(ElecRegions, 'mSTG'));
            ElecRegions(mSTG_in_elec) = {'cSTG'};
            rSTG_in_elec = find(contains(ElecRegions, 'rSTG'));
            ElecRegions(rSTG_in_elec) = {'cSTG'};
            % change elec areas from mMTG, rMTG to cMTG to match the annotation
            mMTG_in_elec = find(contains(ElecRegions, 'mMTG'));
            ElecRegions(mMTG_in_elec) = {'cMTG'};
            rMTG_in_elec = find(contains(ElecRegions, 'rMTG'));
            ElecRegions(rMTG_in_elec) = {'cMTG'};
        end
    end


    % create color based on weights
    if strcmp(Mode, 'SimpleProject')
        % simple projection of weights on surface with sum
        color_weights = zeros(size(verts,1),1);
        for i = 1:size(ElecLocs,1)
            if any(isnan(ElecLocs(i,:)))
                continue;
            end
            bx = abs(verts(:,1)-ElecLocs(i,1));
            by = abs(verts(:,2)-ElecLocs(i,2));
            bz = abs(verts(:,3)-ElecLocs(i,3));
            if obj.flag_UseAnnots
                [~, inds2color] = ismember(ElecRegions{i},...
                                           annt_tbls.struct_names);
                color_weights = color_weights + ...
                                double(inds2color==annot) .* ...
                                Weights(i) .* ...
                                exp(-(bx.^2 + by.^2 + bz.^2)./(2*Sigma^2));
            else
                color_weights = color_weights + ...
                                Weights(i) .* ...
                                exp(-(bx.^2 + by.^2 + bz.^2)./(2*Sigma^2));
            end
        end
    elseif strcmp(Mode, 'SimpleProjectNorm')
        % simple projection of weights on surface with avg normalization
        color_weights = zeros(size(verts,1),1);
        verts_weights = zeros(size(verts,1),1); 
        if flag_rm_sparse
            elec_dist = pdist(ElecLocs);
            elec_dist = squareform(elec_dist);
            elec_dist = sum(elec_dist<=(3*Sigma),1);
        end
        for i = 1:length(ElecLocs)
            if any(isnan(ElecLocs(i,:)))
                continue;
            end
            bx = abs(verts(:,1)-ElecLocs(i,1));
            by = abs(verts(:,2)-ElecLocs(i,2));
            bz = abs(verts(:,3)-ElecLocs(i,3));
            if obj.flag_UseAnnots
                [~, inds2color] = ismember(ElecRegions{i},...
                                           annt_tbls.struct_names);
                if flag_rm_sparse
                    if elec_dist(i)<=sparse_thresh;
                        norm_fact = sparse_ratio;
                    else
                        norm_fact = 1;
                    end
                else
                    norm_fact = 1;
                end
                color_weights = color_weights + ...
                                norm_fact.*double(inds2color==annot) .* ...
                                Weights(i) .* ...
                                exp(-(bx.^2 + by.^2 + bz.^2)./(2*Sigma^2));
                verts_weights = verts_weights + ...
                                exp(-(bx.^2 + by.^2 + bz.^2)./(2*Sigma^2));
            else
                color_weights = color_weights + ...
                                Weights(i) .* ...
                                exp(-(bx.^2 + by.^2 + bz.^2)./(2*Sigma^2));
                verts_weights = verts_weights + ...
                                exp(-(bx.^2 + by.^2 + bz.^2)./(2*Sigma^2));
            end
        end
        % remove devsion of small vals
        color_weights(abs(color_weights)<=tol_norm) = 0;
        verts_weights(abs(color_weights)<=tol_norm) = 1;
        color_weights = color_weights ./ verts_weights;
    elseif strcmp(Mode, 'NNProject')
        msg = 'NNProject simply moves the electrodes to the closest vertex, ';
        msg = [msg, 'consider using ''project_elecs_pial.m'' as an alternative that '];
        msg = [msg, 'uses within subject anatomy for projection.'];
        warning(msg);
        % first move elecs to closest vertex
        % then simple projection of weights on surface with sum
        color_weights = zeros(size(verts,1),1);
        for i = 1:length(ElecLocs)
            if any(isnan(ElecLocs(i,:)))
                continue;
            end
            bx = abs(verts(:,1)-ElecLocs(i,1));
            by = abs(verts(:,2)-ElecLocs(i,2));
            bz = abs(verts(:,3)-ElecLocs(i,3));
            dist = bx.^2+by.^2+bz.^2;
            [~, NN_ind] = min(dist);
            bnx = abs(verts(:,1)-verts(NN_ind,1));
            bny = abs(verts(:,2)-verts(NN_ind,2));
            bnz = abs(verts(:,3)-verts(NN_ind,3));
            if obj.flag_UseAnnots
                [~, inds2color] = ismember(ElecRegions{i},...
                                           annt_tbls.struct_names);
                color_weights = color_weights + ...
                                double(inds2color==annot) .* ...
                                Weights(i) .* ...
                                exp(-(bnx.^2 + bny.^2 + bnz.^2)./(2*Sigma^2));
            else
                color_weights = color_weights + ...
                                Weights(i) .* ...
                                exp(-(bnx.^2 + bny.^2 + bnz.^2)./(2*Sigma^2));
            end
        end
    elseif strcmp(Mode, 'NNProjectNorm')
        msg = 'NNProject simply moves the electrodes to the closest vertex, ';
        msg = [msg, 'consider using ''project_elecs_pial.m'' as an alternative that '];
        msg = [msg, 'uses within subject anatomy for projection.'];
        warning(msg);
        % first move elecs to closest vertex
        % then simple projection of weights on surface with avg normalization
        color_weights = zeros(size(verts,1),1);
        verts_weights = zeros(size(verts,1),1); 
        for i = 1:length(ElecLocs)
            if any(isnan(ElecLocs(i,:)))
                continue;
            end
            bx = abs(verts(:,1)-ElecLocs(i,1));
            by = abs(verts(:,2)-ElecLocs(i,2));
            bz = abs(verts(:,3)-ElecLocs(i,3));
            dist = bx.^2+by.^2+bz.^2;
            [~, NN_ind] = min(dist);
            bnx = abs(verts(:,1)-verts(NN_ind,1));
            bny = abs(verts(:,2)-verts(NN_ind,2));
            bnz = abs(verts(:,3)-verts(NN_ind,3));
            if obj.flag_UseAnnots
                [~, inds2color] = ismember(ElecRegions{i},...
                                           annt_tbls.struct_names);
                color_weights = color_weights + ...
                                double(inds2color==annot) .* ...
                                Weights(i) .* ...
                                exp(-(bnx.^2 + bny.^2 + bnz.^2)./(2*Sigma^2));
                verts_weights = verts_weights + ...
                                exp(-(bnx.^2 + bny.^2 + bnz.^2)./(2*Sigma^2));
            else
                color_weights = color_weights + ...
                                Weights(i) .* ...
                                exp(-(bnx.^2 + bny.^2 + bnz.^2)./(2*Sigma^2));
                verts_weights = verts_weights + ...
                                exp(-(bnx.^2 + bny.^2 + bnz.^2)./(2*Sigma^2));
            end
        end
        % remove devsion of small vals
        color_weights(abs(color_weights)<=tol_norm) = 0;
        verts_weights(abs(color_weights)<=tol_norm) = 1;
        color_weights = color_weights ./ verts_weights;
    else
        error(['Mode: ', Mode, 'Not Implemented!']);
    end

    % Plot the brain
    if flag_AddFigure; figure; end
    handle = trisurf(faces,...
                     verts(:,1),...
                     verts(:,2),...
                     verts(:,3),...
                     'FaceVertexCData', color_weights,...
                     'FaceAlpha', FaceAlpha,...
                     'FaceColor', 'interp');

    % set colormap
    if isstring(cmap) || ischar(cmap)
        if strcmp(cmap(end), 'r')
            cmap = flipud(eval([cmap(1:end-2),'(256)']));
        else
            cmap = eval([cmap,'(256)']);
        end
    end
    colormap(cmap);
    % set colorlim
    if isempty(clim)
        cval_min = min(color_weights, [], 'all');
        cval_max = max(color_weights, [], 'all');
        clim = [cval_min, cval_max];
    end
    caxis(clim);
    colorbar();


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
