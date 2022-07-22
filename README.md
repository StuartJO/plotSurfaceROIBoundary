# plotSurfaceROIBoundary

<img src="./figures/SurfaceBoundaryPlotsExample.png" width="100%">

These scripts will plot the boundaries of a ROI on a surface. Additionally code is provided so that certain regions can be excluded from the colourmap.

There are five ways of defining a boundary here (I know right? Overkill!!), 'faces', 'midpoint', 'centroid', 'edge_vertices', and 'edge_faces'.

## How to run

The code can be run as follows:

```
p = plotSurfaceROIBoundary(surface,parc,data2plot,boundaryMethod,cmap,boundaryWidth);
```

surface is a structure with the fields 'vertices' and 'faces', which define the vertices and faces of the surface respectively.
parc is a vector with the parcellation ID of each vertex
data2plot is a vector of data to plot for each ROI or each vertex
boundaryMethod is one of 'faces', 'midpoint', 'centroid', 'edge_vertices' or 'edge_faces' (see below)
cmap is the colormap you wish to use
boundaryWidth is a scalar indicating the width of the boundary

## Types of boundary plotting methods

'faces' will find the faces which exist between ROIs and those will be coloured black to specify the boundary. 

'midpoint' finds the edges that connect the vertices of two different ROIs and takes the midpoint of the edge and uses those coordinates to define the boundary. 

'centroid' finds the faces which exist between ROIs and uses the centroid of those to draw the coordinates that define the boundary. 

'edge_vertices' finds the vertices which define the boundary of the ROI and draws the boundary along those (i.e., follows the edges between the vertices which make up the ROI boundary). 

'edge_faces' finds the edges which are part of two faces with different ROIs and draws the boundary along those.

(see the images belong if these descriptions are to confusing)

demo_plotSurfaceROIBoundary.m shows examples of how the code is used. See plotSurfaceROIBoundary.m, makeFaceVertexCData.m, and findROIboundaries.m for a more detailed description. I personally recommend using 'midpoint', visually it looks the cleanest and you can easily change the thickness/width of the boundary. 'faces' will look clean but you cannont change the thickness/width, 'centroid' doesn't look especially neat but you can change the thickness, while with 'edge_vertices' the boundaries of adjacent ROIs will leave a slight gap between them (i.e. they don't overlap) which can look a bit strange. 'edge_faces' can look highly irregular when zoomed in.

See makeFaceVertexCData.m for how to configure a colormap to exclude certain regions (e.g., regions which have no ROI information). This function works by assigning a colour directly to each face/vertex instead of getting MATLAB to automatically do it. This allows for somewhat easier manipulation (and it is easier to have multiple colour maps in a plot this way!).

example_surface_data.mat includes vertices and faces for fsaverage, along with parcellations for the Desikan-Killiany atlas, HCPMMP1 parcellation, and a random 200 (100 ROIs per hemisphere) parcellation.

## Boundary examples

The first five figures show the different types of boundary method projected over the surface where each face is coloured according to the vertex ROI IDs (each vertex is a dot coloured according to its ROI ID). While a face can be connected to multiple ROIs, each face is colored by the value of the first vertex that constitutes that face (in other words, when plotting each face appears to only be assigned to one ROI when it infact belongs to multiple). This is why for 'midpoint', 'centroid', and 'edge_vertices' the boundary drawn doesn't exactly follow the face colours. Anatomically based parcellations (like the Desikan-Killiany and HCPMMP1 parcellations) tend to produce jagged boundaries in my experience. All these figures can be replicated using demo_plotSurfaceROIBoundary.m

'faces'

<img src="./figures/ROIS_faces_flat.png" width="100%">

'midpoint'

<img src="./figures/ROIS_midpoint_flat.png" width="100%">

'centroid'

<img src="./figures/ROIS_centroid_flat.png" width="100%">

'edge_vertices'

<img src="./figures/ROIS_edge_vertices_flat.png" width="100%">

'edge_faces'

<img src="./figures/ROIS_edge_faces_flat.png" width="100%">

The next five figures show the different types of boundary method projected over the surface where each face is coloured according to an interpolated value of the vertex values (which in this case is the sulcal depth of each vertex).

'faces'

<img src="./figures/sulc_faces_flat.png" width="100%">

'midpoint'

<img src="./figures/sulc_midpoint_interp.png" width="100%">

'centroid'

<img src="./figures/sulc_centroid_interp.png" width="100%">

'edge_vertices'

<img src="./figures/sulc_edge_vertices_interp.png" width="100%">

'edge_faces'

<img src="./figures/sulc_edge_faces_interp.png" width="100%">

## Issues

If the surface of a particular region is very complex many of the boundary plotting approaches may not work especially well. My code makes some assumptions and these can be violated at times. However I have found that these issue are only noticable if you zoom all the way in. If in doubt, the 'faces' method should be fairly robust to weirdness
