function [elec_locs_proj, elec_rgb, vt_inds] = project_elecs_pial(elec_locs, elec_regions, BrainFile, AnnotFile);
    % function to project the electrodes to the 
    % closest vertex to within the same anatomical region
    % ----
    % load the brian
    if exist(BrainFile, 'file');
        % load the brain file using the load_surface function
        % supports both freesurfer files and mat file.
        [verts, faces] = fsaveragetools.load_surface(BrainFile);
    else
        error(['Brain file: ', BrainFile, ' does not exist!']);
    end
    % load annotation file
    [~, annt_lbls, annt_tbls] = visualtools.read_annot(AnnotFile);
    % quick check of sizes
    if size(annt_lbls,1)~=size(verts,1);
        error('mismatch between annotation and verts sizes!');
    end

    % pre-allocate vecs
    elec_locs_proj = zeros(size(elec_locs));
    elec_rgb = zeros(size(elec_locs,1),3);
    if nargout>2
        vt_inds = zeros(size(elec_locs,1),1);
    end

    % not efficient code... TODO: make the search by region
    for ei = 1:length(elec_regions)
        % find region in table
        [~,re_in_list] = ismember(elec_regions{ei},...
                                  annt_tbls.struct_names);
        elec_rgb(ei,:) = annt_tbls.table(re_in_list,1:3);
        ri = find(annt_lbls==annt_tbls.table(re_in_list,5));
        verts_r = verts(annt_lbls==annt_tbls.table(re_in_list,5),:);
        dist = verts_r - elec_locs(ei,:);
        dist = sum(dist.^2,2);
        [~, di] = min(dist);
        elec_locs_proj(ei,:) = verts_r(di,:);
        if nargout>2
            dist = verts - elec_locs(ei,:);
            dist = sum(dist.^2,2);
            [~, di] = min(dist);
            vt_inds(ei) = di;
        end
    end
    elec_rgb = elec_rgb/255;
end
