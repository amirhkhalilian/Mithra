function handle = PlotElecOnBrain(obj, ElecLocs, varargin)
    % function to plot electrodes on a brain surface
    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'ElecColor', 'r');
    addParameter(p, 'ElecAlpha', nan);
    addParameter(p, 'flag_AddFigure', true);
    addParameter(p, 'radius', 1.5);
    addParameter(p, 'FaceAlpha', 1);
    addParameter(p, 'BrainColor', [0.7, 0.7, 0.7]);
    addParameter(p, 'Brain_handle', []);
    addParameter(p, 'cmap', 'hot');
    addParameter(p, 'clim', []);
    % parse inputs
    parse(p, varargin{:});
    ElecColor = p.Results.ElecColor;
    ElecAlpha = p.Results.ElecAlpha;
    flag_AddFigure = p.Results.flag_AddFigure;
    radius = p.Results.radius;
    FaceAlpha = p.Results.FaceAlpha;
    BrainColor = p.Results.BrainColor;
    Brain_handle = p.Results.Brain_handle;
    cmap = p.Results.cmap;
    clim = p.Results.clim;

    if isnan(ElecAlpha); ElecAlpha = ones(size(ElecColor,1),1); end

    if isempty(Brain_handle)
        % Plot the brain surface
        handle = obj.PlotBrainSurface('flag_AddFigure', flag_AddFigure,...
                                      'FaceAlpha', FaceAlpha,...
                                      'BrainColor', BrainColor);
    else
        handle = Brain_handle;
        axes(handle.Parent);
        hold on;
    end

    % Plot electrodes
    for i = 1:size(ElecLocs,1)
        if (size(ElecColor,1) == 1) && (size(ElecColor,2)==3)
            plotSpheres(ElecLocs(i,1), ElecLocs(i,2), ElecLocs(i,3),...
                        radius, ElecColor, ElecAlpha);
        elseif (ndims(ElecColor)==2) && (size(ElecColor,1)==size(ElecLocs,1))...
                && (size(ElecColor,2)==3)
            % ElecColor has rgb values for each elec
            plotSpheres(ElecLocs(i,1), ElecLocs(i,2), ElecLocs(i,3),...
                        radius, ElecColor(i,:), ElecAlpha(i));
        elseif (ndims(ElecColor)==2) && (size(ElecColor,1)==size(ElecLocs,1))...
                && (size(ElecColor,2)==1)
            % ElecColor has the value for each elec... convert to rgb before plot
            if isempty(clim); clim=[min(ElecColor(:)), max(ElecColor(:))]; end
            [ElecColor, rgb_ind, cmap] = obj.value2color(ElecColor, cmap, clim);
            plotSpheres(ElecLocs(i,1), ElecLocs(i,2), ElecLocs(i,3),...
                        radius, ElecColor(i,:), ElecAlpha(i));
            colormap(gca, cmap); colorbar(); caxis(clim);
        else
            error('ElecColor size invalid');
        end
    end
end

function [shand]=plotSpheres(spheresX, spheresY, spheresZ,...
                             spheresRadius, col, face_alpha)
    spheresRadius = ones(length(spheresX),1).*spheresRadius;
    % set up unit sphere information
    numSphereFaces = 10;
    [unitSphereX, unitSphereY, unitSphereZ] = sphere(numSphereFaces);
    % set up basic plot
    sphereCount = length(spheresRadius);
    % for each given sphere, shift the scaled unit sphere by the
    % location of the sphere and plot
    for i=1:sphereCount
        sphereX = spheresX(i) + unitSphereX*spheresRadius(i);
        sphereY = spheresY(i) + unitSphereY*spheresRadius(i);
        sphereZ = spheresZ(i) + unitSphereZ*spheresRadius(i);
        shand = surface(sphereX, sphereY, sphereZ,...
                        'FaceColor',col,'EdgeColor','none',...
                        'FaceAlpha', face_alpha,...
                        'FaceLighting', 'gouraud',...
                        'AmbientStrength',0.55);
    end
end
 

