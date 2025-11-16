function handle = PlotElecConnectivityOnBrain(obj, ElecLocs, Adjacency, varargin)
    % function to plot electrodes on a brain surface
    p = inputParser;
    p.KeepUnmatched = true;
    addParameter(p, 'ElecColor', 'k');
    addParameter(p, 'LineColor', 'b');
    addParameter(p, 'LineWidth', 2);
    addParameter(p, 'flag_AddFigure', true);
    addParameter(p, 'radius', 1.5);
    addParameter(p, 'FaceAlpha', 0.5);
    addParameter(p, 'LineRatio', 1);
    addParameter(p, 'BrainColor', [0.8, 0.8, 0.8]);
    addParameter(p, 'cmap', 'hsv');
    addParameter(p, 'clim', []);
    % parse inputs
    parse(p, varargin{:});
    ElecColor = p.Results.ElecColor;
    LineWidth = p.Results.LineWidth;
    LineColor = p.Results.LineColor;
    flag_AddFigure = p.Results.flag_AddFigure;
    radius = p.Results.radius;
    FaceAlpha = p.Results.FaceAlpha;
    BrainColor = p.Results.BrainColor;
    LineRatio = p.Results.LineRatio;
    cmap = p.Results.cmap;
    clim = p.Results.clim;

    % Plot the brain surface
    handle = obj.PlotBrainSurface('flag_AddFigure', flag_AddFigure,...
                                  'FaceAlpha', FaceAlpha,...
                                  'BrainColor', BrainColor);

    % convert adjacency matrix to X,Y,Z line plots
    [X,Y,Z, A, elecs_2_use] = adjacency_mat_to_lines(Adjacency, ElecLocs);

    % use plot3 to plot lines in between connections
    hold on;
    [cnct_colors, rgb_inds, cmap] = obj.value2color(A, cmap, clim);
    for  k = 1:size(X,2)
        plot3(X(:,k), Y(:,k), Z(:,k),...
              'Color', cnct_colors(k,:),...
              'LineWidth', LineRatio*abs(A(k)));
    end

    ElecLocs = ElecLocs(elecs_2_use,:);

    % Plot electrodes 
    % Note: different from PlotElecOnBrain -- can only have a single color 
    % for all electrodes
    for i = 1:size(ElecLocs,1)
        if size(ElecColor,1) == 1
            plotSpheres(ElecLocs(i,1), ElecLocs(i,2), ElecLocs(i,3),...
                        radius, ElecColor);
        else
            error('ElecColor size invalid');
        end
    end

end

function [shand]=plotSpheres(spheresX, spheresY, spheresZ,...
                             spheresRadius, col)
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
                        'AmbientStrength',0.7);
    end
end
 
function [X,Y,Z,Vals,elecs_inds] = adjacency_mat_to_lines(adj, coords)
    % function to convert the adjacency matrix to XYZ values 
    % such that plot3 would plot lines in between them
    % ref: matlab gplot 

    % assuming adj matrix is symmetric and undirected

    [i,j] = find(triu(adj,1));
    elecs_inds = unique([i;j]);
    [~, p] = sort(max(i,j));
    i = i(p);
    j = j(p);

    X = [coords(i,1), coords(j,1)]';
    Y = [coords(i,2), coords(j,2)]';
    Z = [coords(i,3), coords(j,3)]';
    Vals = adj(sub2ind(size(adj),i,j))';
end
