# plotSurfaceROIBoundary
Plot the boundaries of ROIs on a surface

This script will plot the boundaries of a roi on a surface. See plotROIboundaries.m for detailed usage.

The script has three ways of defining a boundary, 'colorface', 'midpoint', and 'centroid'.

'colorface' will find the faces which exist between rois and those will be coloured black to specify the boundary. 'midpoint' finds the edges that connected the vertices of two different rois and takes the midpoint of the egde and uses those coordinates to define the boundary. 'centroid' finds the faces which exist between rois and uses the centroid of those to draw the coordinates that define the boundary. Examples are shown below

The first three figures show the different types of boundary method projected over the surface where each face is coloured according to the vertex roi ids

'colorface'

<img src="./figures/ROIS_colorface_flat.png" width="50%">

'midpoint'

<img src="./figures/ROIS_midpoint_flat.png" width="50%">

'centroid'

<img src="./figures/ROIS_centroid_flat.png" width="50%">

The next three figures show the different types of boundary method projected over the surface where each face is coloured according to an interpolated value of the vertex values (which in this case is the sulcal depth of each vertex).

'colorface'

<img src="./figures/sulc_colorface_flat.png" width="50%">

'midpoint'

<img src="./figures/sulc_midpoint_interp.png" width="50%">

'centroid'

<img src="./figures/sulc_centroid_interp.png" width="50%">
