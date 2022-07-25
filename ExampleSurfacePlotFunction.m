function [p_left,p_right,c] = ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label)

% This is just a example for how to plot the medial and lateral sides of a
% surface in the same plot

% Inputs:
%
% surface = a structure with two fields: vertices (the vertices making up 
% the surface) and faces (the faces of the surface)
%
% vertex_id = the roi id of each vertex
%
% data = either data for each individual roi or data for each vertex.
% If you don't want any data to be displayed for a roi or vertex, set that
% value to NaN.
%
% cmap = an N*3 matrix specifying the RGB values making up the colormap to
% use
%
% data_label = (optional) the name of the data being plotted
%
% Output
%
% p_left = patch object for the surface plotted on the left
%
% p_right = patch object for the surface plotted on the right
%
% c = colorbar object

figure('Position',[461   462   560   325])
ax_sub1 = axes('Position',[0.005 .33 .49 .66]);
p_left = plotSurfaceROIBoundary(surface,vertex_id,data,'midpoint',cmap,1);
camlight(80,-10);
camlight(-80,-10);
view([-90 0])

axis off
axis image

ax_sub2 = axes('Position',[.505 .33 .489 .66]);
p_right = plotSurfaceROIBoundary(surface,vertex_id,data,'midpoint',cmap,1);
camlight(80,-10);
camlight(-80,-10);
view([90 0])
axis off
axis image
colormap(cmap)
caxis([nanmin(data) nanmax(data)])
c = colorbar('Location','southoutside');
%set(c, 'xlim', [nanmin(data) nanmax(data)],'Position',[.1 .23 .8 .05],'FontSize',20);
set(c, 'Position',[.1 .23 .8 .05],'FontSize',20);
if exist('data_label','var')
    c.Label.String = data_label;
end