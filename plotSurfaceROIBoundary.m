function [p,boundary_plot,new_cmap,BOUNDARY,new_climits,orig_data_climits] = plotSurfaceROIBoundary(surface,vertex_id,data,boundary_method,cmap,colorUnknownGrey,linewidth,climits)

% This script is a wrapper for findROIboundaries and makeFaceVertexCData so
% that they smoothly work together and you don't have to spend a lot of
% your own time making them work together. The code sets up the basics of
% the patch object

% Inputs:
%
% vertices = the vertices making up the surface
%
% faces = the faces of the surface
%
% vertex_id = the roi id of each vertex
%
% data = either data for each individual roi or data for each vertex.
% If you don't want any data to be displayed for a roi or vertex, set that
% value to NaN.
%
% boundary_method = 'faces', 'midpoint', or 'centroid'. 'faces'
% will find the faces which exist between ROIs and those will be coloured
% black to specify the boundary. 'midpoint' finds the edges that connected
% the vertices of two different ROIs and takes the midpoint of the egde and
% uses those coordinates to define the boundary. 'centroid' finds the faces
% which exist between ROIs and uses the centroid of those to draw the 
% coordinates that define the boundary
%
% cmap = an N*3 matrix specifying the RGB values making up the colormap to
% use
%
% colorUnknownGrey = set to 1 to colour any unknown rois (i.e. vertices
% with an id of 0) grey
%
% linewidth = the width of the boundary when using 'midpoint', or 
% 'centroid'.
%
% climits = the range to apply the colormap. This will work perfectly fine
% if the range specified is larger than the data itself or if the upper
% limit is larger. However if the lower limit is larger than the smallest
% value in data, then it may get strange. If colorUnknownGrey = 1, then
% faces/vertices with a value smaller than the limit will be coloured grey,
% or potentially black if 'faces' is used. If it is set to 0 and 'faces' is
% used, those regions will be set to black while if 'centroid' or midpoint' 
% are selected, the colormap will work appropriately). So if you really 
% need to enforce a lower limit I would suggest threshold the data in 
% advance and all should be good.
%
% Outputs:
%
% p = the patch surface object
%
% boundary_plot = a strcture containing a line object defining each
% boundary
%
% new_cmap = the new colormap being used for the figure
%
% BOUNDARY = the boundary of the rois. For 'faces' this will be a logical
% where a value of 1 indicates that face is on the boundary between ROIs.
% For 'midpoint'  or 'centroid', BOUNDARY will be a cell where each element
% contains the coordinates of the points making up the boundary, which can
% be plotted as a line. Note that each boundary defines a continuous ROIs,
% if a ROI is made up of multiple non-continuous parts (i.e., the ROI is
% made up of multiple unconnected sections), there will be a boundary for
% each of those parts
%
% new_climits = the new caxis limits for the data.
%
% orig_data_climits = the limits for data. If climits is used, this is just
% that, but if not this will give you the limits for the original data. If
% making a colorbar this should be used to set the range of that.
%
% If you want to include a colorbar but don't want it to display the
% black/grey values, you can do:
% c = colorbar;
% set(c, 'xlim', orig_data_climits);

if nargin < 7
    linewidth = 2;
end

if nargin < 8 
    climits = [nanmin(data) nanmax(data)];
end

% Extract the faces and vertices
vertices = surface.vertices;
faces = surface.faces;

% Find the boundaries of the ROIs
BOUNDARY = findROIboundaries(vertices,faces,vertex_id,boundary_method);

% Set up some options. 
switch boundary_method
    case 'faces'
    colorFaceBoundaries = 1;
% If boundaries are defined by faces, then the face
% color method for the patch object needs to be 'flat' otherwise it won't
% work
    face_color_method = 'flat';

    case {'midpoint','centroid'}
    colorFaceBoundaries = 0;
    
    if length(data) == length(vertex_id)
        % If there is data for each individual vertex, use 'interp' for the 
        % face color 
        face_color_method = 'interp';
        
    else
        % If there is data for each individual ROI, use 'flat' for the 
        % face color. Also by default colour any unlabelled ROIs as grey        
        face_color_method = 'flat';
        colorUnknownGrey = 1;
    end
    otherwise
        error('Unrecognised ''boundary_method'' option')
end

% Get the vertex or face data and colormap to use
[FaceVertexCData,new_cmap,new_climits,orig_data_climits] = makeFaceVertexCData(vertices,faces,vertex_id,data,cmap,colorFaceBoundaries,colorUnknownGrey,climits);

% Plot the surface using some preconfigured options
p = patch(surface);
set(p,'FaceVertexCData',FaceVertexCData,'EdgeColor','none','FaceColor',face_color_method,'Clipping','off');
p.FaceLighting = 'gouraud';
material dull

% Apply the colormap to the current axes
current_axes = gca;
colormap(current_axes,new_cmap)

% Apply the new color map limits
caxis(new_climits)

hold on

% Draw the boundary if 'midpoint' or 'centroid' was used.

switch boundary_method
    case {'midpoint','centroid'}
        for i = 1:length(BOUNDARY)
           boundary_plot.boundary(i) = plot3(BOUNDARY{i}(:,1), BOUNDARY{i}(:,2), BOUNDARY{i}(:,3), 'Color', 'k', 'LineWidth',linewidth,'Clipping','off');
        end
    otherwise
        boundary_plot = [];
end