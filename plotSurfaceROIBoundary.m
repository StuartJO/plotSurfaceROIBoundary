function [p, boundary_plot, BOUNDARY] = plotSurfaceROIBoundary(surface,varargin)
%% Plots brain surface mesh with boundaries between parcels
% This function plots a colored patch object and demarcates parcels/ROIs on that
% surface. This requires the input of the surface mesh, the ROI IDs of each
% vertex, and the data to be used to colour each vertex.
% 
% This function plots 3 separate objects: (i) the main coloured patch object,
% which adjusts to the colormap of the axes on which it is plotted; (ii) another
% patch object which is used to plot the 'medial wall', and is a fixed color;
% (iii) a patch or line object to demarcate parcel boundaries (and is also a
% fixed color). Although the data values of (i) are fixed, the colors can be
% easily changed by adjust the axes e.g. with `colormap` or `caxis`. The color
% of (ii) and (iii) is best decided at function runtime, but can be adjusted
% post hoc by accessing the data in the return arguments if needed.
% 
%% Syntax  
%  plotSurfaceROIBoundary(surface,vertex_id,data,boundary_method,cmap,linewidth)
% 
%  plotSurfaceROIBoundary(___,Name,Value)
% 
%  p = plotSurfaceROIBoundary(___)
%  [p,boundary_plot] = plotSurfaceROIBoundary(___)
%  [p,boundary_plot,BOUNDARY] = plotSurfaceROIBoundary(___)
% 
% 
%% Examples 
%   figure; plotSurfaceROIBoundary(struct('vertices', [cos([(0:6)';(0:6)']), sin([(0:6)';(0:6)']), [(1:7)';(0:6)']], 'faces', [[1,2,8]+(0:5)';[2,9,8]+(0:5)'])); view([45 45]);   
%   figure; plotSurfaceROIBoundary(struct('vertices', [cos([(0:6)';(0:6)']), sin([(0:6)';(0:6)']), [(1:7)';(0:6)']], 'faces', [[1,2,8]+(0:5)';[2,9,8]+(0:5)']), [], 1:14); view([45 45]);
%   figure; plotSurfaceROIBoundary(struct('vertices', [cos([(0:6)';(0:6)']), sin([(0:6)';(0:6)']), [(1:7)';(0:6)']], 'faces', [[1,2,8]+(0:5)';[2,9,8]+(0:5)']), [0;0;1;1;1;0;0;0;1;1;1;1;0;0], 1:14); view([45 45]);
%   figure; plotSurfaceROIBoundary(struct('vertices', [cos([(0:6)';(0:6)']), sin([(0:6)';(0:6)']), [(1:7)';(0:6)']], 'faces', [[1,2,8]+(0:5)';[2,9,8]+(0:5)']), [1;1;2;2;2;3;3;1;2;2;2;2;3;3], [1 3 2]); view([45 45]);
% 
% See the examples in the `demo_*.m` files. 
% 
% 
%% Input arguments  
%  surface - structure with fields `vertices` and `faces`
%  `vertices` should be a three column matrix (V x 3) containing xyz coordinates
%  of the surface mesh. `faces` should be a three column matrix containing the
%  vertex coordinates used to make each face.
% 
%  vertex_id - ROI allocations of each vertex (V x 1 vector)
%  Vertices indexed `0` will be treated as the 'medial wall' and will be
%  colored according to `unknown_color`.
% 
%  data - data to plot for each vertex or ROI (V x 1 or R x 1 vector)
%  Data indexed `NaN` will be treated as the 'medial wall' and will be colored
%  according to `unknown_color`.
% 
%  boundary_method - method to delineate ROIs ('faces' (default) | 'midpoint' | 'centroid' | 'edge_faces' | 'edge_vertices')
%
%  cmap - colormap to color data (three column matrix | function handle | string)
%
%  linewith - linewidth for plotting line between ROIs (positive scalar)
%
%  climits - caxis limits (two element matrix)
% 
%% Name-Value Arguments  
%  boundary_color - color to draw faces/lines to separate ROIs (RGB triplet)
% 
%  unknown_color - color to draw areas marked as unknown (RGB triplet)
% 
%  FaceColor - patch face coloring method ('flat' (default) | 'interp')
%  
%  patch_options - options to use when drawing patches (cell array)
%  The contents of this cell array will be unravelled and passed to PATCH
%  without modification.
%  
%  FacesFromRois - flag to color based on ROI data ('first' (default) | 'mean')
%  Each face has 3 vertices, but is generally one color. This flag determines
%  how to determine the color for each face. If set to 'first', the face will be
%  the color of its first vertex. If set to 'mean', the face will be the color
%  of the mean values of its vertices.
% 
% 
%% Output Arguments  
%  p - patch object containing colored data
%  
%  boundary_plot - struct containing boundary information with fields 'nanZeroPatch' and 'boundary' 
%  `nanZeroPatch` is a patch object of vertices indexed with 0 and faces indexes
%  with NaN. If boundary_method is set to 'faces', `boundary` is a patch object
%  of the vertices and faces along the boundary between ROIs. Otherwise,
%  `boundary` is an array of `plot3` objects drawing the lines that make up the
%  boundary. 
%  
%  BOUNDARY - the boundary of the ROIs
%  For 'faces' this will be a logical where a value of 1 indicates that face is
%  on the boundary between ROIs. Otherwise, BOUNDARY will be a cell where each
%  element contains the coordinates of the points making up a boundary, which
%  can be plotted as a line. Note that each boundary defines a continuous ROI,
%  if a ROI is made up of multiple non-continuous parts (i.e., the ROI is made
%  up of multiple unconnected sections), there will be a boundary for each of
%  those parts
% 
% 
%% See Also  
%  PATCH
%  https://au.mathworks.com/help/matlab/ref/matlab.graphics.primitive.patch-properties.html
%  findROIBoundaries, graphComponents
% 
% 
%% Authors  
%  Mehul Gajwani, Monash University, 2023


%% Prelims
ip = inputParser;
addRequired(ip, 'surface');
addOptional(ip, 'vertex_id', []);
addOptional(ip, 'data', []);
addOptional(ip, 'boundary_method', 'faces', ...
    @(x) any(strcmp({'faces','midpoint','centroid','edge_faces','edge_vertices'},x)));
addOptional(ip, 'cmap', parula);
addOptional(ip, 'linewidth', 2);
addOptional(ip, 'climits', []);

addParameter(ip, 'boundary_color', [0 0 0]);
addParameter(ip, 'unknown_color', [0.5 0.5 0.5]);
addParameter(ip, 'patch_options', {'EdgeColor','none','FaceLighting','gouraud'}, @(x) iscell(x));
addParameter(ip, 'FaceColor', 'flat', @(x) any(strcmp({'flat', 'interp'}, x)) );
addParameter(ip, 'FacesFromROIs', 'first', @(x) any(strcmp({'first', 'mean'}, x)) );

ip.parse(surface, varargin{:});


%% Parse inputs
verts = ip.Results.surface.vertices;
faces = ip.Results.surface.faces;
rois = ip.Results.vertex_id;
data = ip.Results.data;

% If data is in annot file format
if size(data,2) == 2
    [~,~,rois] = unique(rois, 'sorted');
    [~, temp] = sort(data(:,2));
    data = data(temp, 1);
end

[verts, faces, rois, data] = ...
    checkVertsFacesRoisData(verts, faces, rois, data, 'fillEmpty', true);

boundary_method         = ip.Results.boundary_method;
cmap                    = ip.Results.cmap;
linewidth               = ip.Results.linewidth;
climits                 = ip.Results.climits;
patch_options           = ip.Results.patch_options;
boundary_color          = ip.Results.boundary_color;
face_color_method       = ip.Results.FaceColor;
facesFromRois           = ip.Results.FacesFromROIs;


%% Format data
% put it at vertex level if it is at ROI level
% data for a ROI indexed 0 is set to NaN
if size(data, 1) == max(rois)
    if any(rois == 0)
        temp = [nan; data];
        data = temp(rois+1);
    else
        data = data(rois);
    end
end


%% Start by plotting faces indexed 0 or data indexed nan
% dataFaces are the faces that contain data and should NOT be plotted as part of
% the boundary or unknown areas

roiFaces = rois(faces);
nanFaces = any(ismember(faces, find(isnan(data))), 2);

if strcmp(boundary_method, 'faces')
    dataFaces = any(roiFaces, 2); % faces with any vertex on known region
else
    dataFaces = all(roiFaces, 2); % faces with all vertices on known regions
end

boundary_plot.nanZeroPatch = ...
    patch(struct('Vertices', verts, 'Faces', faces(nanFaces | ~dataFaces,:)), ...
    'FaceColor', ip.Results.unknown_color, patch_options{:});
hold on;


%% MAIN: Plot the data
if strcmp(boundary_method, 'faces')
    boundaryFaces = any(diff(roiFaces, [], 2), 2);
    BOUNDARY = boundaryFaces;
else
    boundaryFaces = zeros(size(faces,1),1);
end

if strcmp(face_color_method, 'flat')
    % --> manually create FaceVertexCData by using the first/mean vertex of each face
    if strcmp(facesFromRois, 'first')
        data = data(faces(dataFaces&~boundaryFaces,1));
    else
        data = mean(data(faces(dataFaces&~boundaryFaces,:)));
    end
    % --> otherwise, leave data at vertex level and allow MATLAB to interpolate
end

p = patch(struct('Vertices', verts, 'Faces', faces(dataFaces&~boundaryFaces,:)), ...
    'FaceVertexCData', data, 'FaceColor', face_color_method, patch_options{:});


%% Plot boundary
if strcmp(boundary_method, 'faces')
    boundary_plot.boundary = ...
        patch(struct('Vertices', verts, 'Faces', faces(boundaryFaces,:)), ...
        'FaceColor', boundary_color, patch_options{:});
else
    BOUNDARY = findROIboundaries(verts,faces,rois,boundary_method);
    if linewidth > 0
        for i = 1:length(BOUNDARY)
            boundary_plot.boundary(i) = ...
                plot3(BOUNDARY{i}(:,1), BOUNDARY{i}(:,2), BOUNDARY{i}(:,3), ...
                'Color', boundary_color, 'LineWidth',linewidth,'Clipping','off');
        end
    end
end


%% Beautify
material dull;
colormap(gca, cmap());
if ~isempty(climits); try caxis(climits); catch; end; end

% camlight(gca, 80, -10);
% camlight(gca, -80, -10);
% view([90 0]);
% axis off equal vis3d;


end
