function [BOUNDARY,BOUNDARY_ROI_ID,ROI_COMPONENTS] = findROIboundaries(vertices,faces,vertex_id,boundary_method)

% This script will plot the boundaries defined by some parcellation/ROIS on
% a surface projection. It can also alter the colormap so regions that do
% not have any information are coloured grey.
%
% Inputs:
%
% vertices = the vertices making up the surface
%
% faces = the faces of the surface
%
% vertex_id = the roi id of each vertex
%
% boundary_method = 'faces', 'midpoint', 'centroid', or 'edges'. 'faces'
% will find the faces which exist between ROIs and those will be coloured
% black to specify the boundary. 'midpoint' finds the edges that connected
% the vertices of two different ROIs and takes the midpoint of the edge and
% uses those coordinates to define the boundary. 'centroid' finds the faces
% which exist between ROIs and uses the centroid of those to draw the 
% coordinates that define the boundary. 'edges' finds the vertices which
% define the boundary of the ROI and uses them for the coordinates
%
% Outputs:
%
% BOUNDARY = the boundary of the rois. For 'faces' this will be a logical
% where a value of 1 indicates that face is on the boundary between ROIs.
% For 'midpoint', 'centroid' or 'edges', BOUNDARY will be a cell where each 
% element contains the coordinates of the points making up a boundary, 
% which can be plotted as a line. Note that each boundary defines a 
% continuous ROI, if a ROI is made up of multiple non-continuous parts 
% (i.e., the ROI is made up of multiple unconnected sections), there will 
% be a boundary for each of those parts
%
% BOUNDARY_ROI_ID = if 'midpoint' or 'centroid' is used, this will
% indicate to which roi each boundary belongs (as indicated by the 
% respective element in BOUNDARY e.g., BOUNDARY_ROI_ID(1) = 1 means
% BOUNDARY{1} defines a boundary for ROI 1, BOUNDARY_ROI_ID(2) = 1 means
% BOUNDARY{2} defines a boundary for ROI 1, BOUNDARY_ROI_ID(3) = 2 means
% BOUNDARY{3} defines a boundary for ROI 2 etc etc). BOUNDARY_ROI_ID will 
% be empty if 'faces' is used for boundary_method
%
% ROI_COMPONENTS = for each ROI, indicates how many componets make it up.
% In other words, this indicates if there are multiple nonspatially
% continuous parts to the ROI. This will only work if 'faces' is not 
% used
%
% Stuart Oldham, Monash University, 2020
% Thanks to the coronavirus for giving me the time to make this script

if nargin < 4
   boundary_method = 'midpoint'; 
end

switch boundary_method
    case 'faces'
           
            % Find the rois each face is connected to

            faces_roi_ids = vertex_id(faces);
            
            % Find the boundary faces 
            BOUNDARY = logical(diff(faces_roi_ids,2,2));
            
            BOUNDARY_ROI_ID = [];
            
            ROI_COMPONENTS = [];
            
    case {'midpoint','centroid'}
    
    % Replace vertex ids for faces with the corresponding roi ids
    faces_roi_ids = vertex_id(faces);

    % Calculate how many rois each face belongs to
    n_rois_per_face = sum(diff(sort(faces_roi_ids,2),1,2)~=0,2)+1;
    
    % Find all possible combinations of edges (each edge has 4
    % combinations)
    EDGESall = [[faces(:,1); faces(:,1); faces(:,2); faces(:,2); faces(:,3); faces(:,3)], ...
        [faces(:,3); faces(:,2); faces(:,1); faces(:,3); faces(:,2); faces(:,1)]];
    
    % Sort the edges so the same edge is grouped together   
    [EDGES_sorted,EDGES_sorted_indx] = sortrows(sort(EDGESall,2));

    % Get the number of edges
    nEdges = length(EDGES_sorted_indx);

    % Get the number of faces
    nFaces = size(faces,1);

    % Get a unique ID for each face
    FACE_ID = 1:nFaces;

    % Map the unique IDs to be the same size as EDGESall
    EDGES_FACE_ID_all = repmat(FACE_ID',6,1);

    % Sort EDGES_FACE_ID_all so the face ids are assigned to the correct
    % edges
    E_F_ID_sorted = EDGES_FACE_ID_all(EDGES_sorted_indx);

    % Each edge is represented four ways in EDGESall (each edge is assigned
    % to 2 faces, and the direction of the edge can be represented in 2
    % ways) and we don't know what order the are arranged. As we want to
    % find the two faces each edge is asigned to, and each edge is grouped
    % together in a block in E_F_ID_sorted, we define 4 indexes so we can
    % pull out each seperate representation of the edge
    
    ef1 = 1:4:nEdges-3;
    ef2 = 2:4:nEdges-2;
    ef3 = 3:4:nEdges-1;
    ef4 = 4:4:nEdges;
    
    % Extract one representation of an edge

    EDGES = EDGES_sorted(ef1,:);
    
    % Make a matrix where each row is an edge, and the faces for each of  
    % the representations of each edge are the columns. We sort each row,
    % so therefore the two unique edges are contained in columns 2 and 3

    EDGES_FACE_ID = sort([E_F_ID_sorted(ef1,:) E_F_ID_sorted(ef2,:) E_F_ID_sorted(ef3,:) E_F_ID_sorted(ef4,:)],2);

    EDGES_FACE_ID(:,[1 4]) = [];
    
    % Get the rois of the two vertices defining the edge

    edges_roi_id = vertex_id(EDGES);

    % Find edges which are contained within the one roi
    
    edge_same_roi_ind = diff(edges_roi_id,1,2)==0;
    
    % Find edges which lie between two rois (boundary edges)

    mult_edges_vert_id = EDGES(~edge_same_roi_ind,:);

    % Get the associated roi ids for the boundary edges
    
    mult_edges_roi_id = edges_roi_id(~edge_same_roi_ind,:);
    
    % Get the face ids for the boundary edges

    mult_edges_face_id = EDGES_FACE_ID(~edge_same_roi_ind,:);
    
    % For the face ids assigned to the boundary edges, get how many rois
    % they connect

    mult_edges_face_nrois = n_rois_per_face(mult_edges_face_id);

    % Get all the roi ids and the number of rois
    
    roi_ids = unique(vertex_id);
    nrois = length(roi_ids);

    % Potentially a roi is not spatially contiguous, so will have multiple 
    % boundaries. The code will iteratively add in boundaries because of this
    
    nBounds = 1;
    
    ROI_COMPONENTS = zeros(nrois,1);

    for i = 1:nrois
        
        % Find the boundary edges for a given roi
        
        index = find(any(mult_edges_roi_id == roi_ids(i), 2));
        
        % Get the face ids for the boundary edges for a given roi

        edge_roi_boundaryfaces = mult_edges_face_id(index,:);

        % Get all the unique faces along the roi boundary and the number of 
        % such faces 
        
        roi_faces = unique(edge_roi_boundaryfaces);
        nfaces = length(roi_faces);

        % Relabel the faces from 1:nfaces so the boundary of the roi can be
        % found easily using graph theory methods
        
        %roi_boundaryfaces_ordered = changem(edge_roi_boundaryfaces,1:nfaces,roi_faces);
        
        newval = 1:nfaces;
        oldval = roi_faces;
        roi_boundaryfaces_ordered = edge_roi_boundaryfaces;
        
        for k = 1:numel(newval)
            roi_boundaryfaces_ordered(edge_roi_boundaryfaces == oldval(k)) = newval(k);
        end
    
        % Make an edge list using the relabeled faces as nodes. Each row
        % defines a pair of faces which are connected. The edge list needs
        % to be symmetric

        edgeList = [roi_boundaryfaces_ordered; fliplr(roi_boundaryfaces_ordered)];
        
        % Get the id of each edge
        
        edgeID = [index; index];
        
        % Make an adjacency matrix where each connection defines a pair of
        % faces which are adjacent, and the connection weight is the
        % corresponding edge id that connects those two faces

        adj = sparse(edgeList(:, 1), edgeList(:, 2), edgeID, nfaces, nfaces);

        % The boundary of a ROI should form a circle. If for some reason the 
        % roi is not spatially continuous, each of the subrois should form a 
        % cycle. If not something is very wrong with your parcellation

        % Find the components of the adjacency matrix. Each component
        % represents a spatially continuous area for the current roi
        
        [comp,ROI_COMPONENTS(i)] = graphComponents(adj);
        
        % Loop over each roi component (i.e., the 

        for j = 1:ROI_COMPONENTS(i)
            
            % Find the "nodes" making up the component

            Nodes2Use = find(comp == j);
            
            % Define a start "node"

            N1 = Nodes2Use(1);
            
            % Find the neighbours of the "node"

            NodeNeighbors = find(adj( N1,:));
            
            % Pick one of those neighbours to be the end "node"

            N2 = NodeNeighbors(1);
            
            % Find the connection linking the start and end "nodes". Once
            % found, record it, then remove it from the adjacency matrix 

            FinalConnection = full(adj(N1,N2));

            adj(N1,N2) = 0;

            adj(N2,N1) = 0;
            
            % The boundary should be a cycle in the graph. If you remove
            % one connection between "node" N1 and N2, if you find the 
            % shortest path between N1 and N2, you get the ordering of
            % faces which define the boundary. You can also get the order
            % of the connections traversed to make this path "E". 

            G = graph(adj);

            [~,~,E] = shortestpathtree(G,N1,N2,'OutputForm','cell');
            
            % We know the connections traversed, so we use this index to 
            % get the original edge ids which make up the boundary. To make
            % the boundary complete we add in the one edge we originally
            % removed. Each edge defines a coordinate.

            EdgeINDS = [G.Edges.Weight(E{1}) ; FinalConnection];

            switch boundary_method

                case 'midpoint'
                    
    % The commented lines below define a quick and easy way to get the
    % boundary. However if done this way, where three rois meet the
    % boundaries of those rois will not intersect. This bothered me. So I
    % wrote more code to get around this you can see below. If this does
    % not bother you, you mad person, you can just uncomment the code below
    % and comment out the rest of the code up until the next case statement

    %         Boundary_verts = mult_edges_vert_id(EdgeINDS,:);
    %                 
    %         midpoints = (vertices(Boundary_verts(:, 1), :) + vertices(Boundary_verts(:, 2), :))./2;
    %                         
    %         BOUNDARY{nBounds} = [midpoints; midpoints(1,:)];
    %     
    %         nBounds = nBounds + 1;
    
            % For the rois boundary edges, get the number of rois the
            % corresponding faces are connected to

            Boundary_verts_nrois = mult_edges_face_nrois(EdgeINDS,:); 

            % Find the "corners" (i.e. where three rois intersect) of the
            % boundary

            boundary_corner_points = find(max(Boundary_verts_nrois,[],2)==3);

            % Check the there are corners for this boundary

            if ~isempty(boundary_corner_points)

            % Prevent the boundary from starting on a corner, this messes up
            % code later. The corners should come in pairs in EdgeINDS and 
            % the code expects as such. If the first edge but not second
            % edge in a corner edge, this means the last edge is a corner
            % edge. If so we just reorder to matrix to add the starting
            % edge to the end
                
                if boundary_corner_points(1) == 1 && boundary_corner_points(2) == 1          
                    EdgeINDS_1 = EdgeINDS(1,:);

                    EdgeINDS(1,:) = [];

                    EdgeINDS(length(EdgeINDS)+1,:) = EdgeINDS_1;

                    Boundary_verts_nrois = mult_edges_face_nrois(EdgeINDS,:);

                    boundary_corner_points = find(max(Boundary_verts_nrois,[],2)==3);
                end        
            
                % Find the corresponding faces of the rois boundary edge 

                Boundary_faces = mult_edges_face_id(EdgeINDS,:);

                % Find the corresponding vertices of the rois boundary edge 

                Boundary_verts = mult_edges_vert_id(EdgeINDS,:);

                % Calculate the midpoints of the boundary edges based on the
                % vertices which make up that edge

                midpoints = (vertices(Boundary_verts(:, 1), :) + vertices(Boundary_verts(:, 2), :))./2;

                % Create a binary mask of which edges are corners

                Boundary_corner_faces_mask = Boundary_verts_nrois==3;

                % Find the faces which make up the corners

                corner_faces = Boundary_faces;

                corner_faces(~Boundary_corner_faces_mask) = 0;

                boundary_corner_faces = max(corner_faces(boundary_corner_points,:),[],2);

                % Find the vertices of the corner faces

                boundary_corner_faces_verts = faces(boundary_corner_faces,:);

                % Get the centroid of the corner faces

                corner_coords = (vertices(boundary_corner_faces_verts(:, 3), :) + vertices(boundary_corner_faces_verts(:, 2), :) + vertices(boundary_corner_faces_verts(:, 1), :))./3;

                corner_coords2use_ind = 1:2:length(boundary_corner_points);

                corner_coords2use = corner_coords(corner_coords2use_ind,:);

                corner_coords_insert = boundary_corner_points(corner_coords2use_ind);

                % Create the matrix defining the final number of boundary
                % coordinates

                boundary_coords = zeros(size(midpoints,1)+size(corner_coords2use,1),size(midpoints,2));       

                % Find indices for the non-corner boundary edges

                addRows = ismember(1:size(midpoints,1), corner_coords_insert);

                oldDataInd = (1:size(midpoints,1)) + cumsum([0, addRows(1:end-1)]);

                % Add in the non-corner boundary edges

                boundary_coords(oldDataInd,:) = midpoints(:,:);

                % Find indices for the corner boundary edges

                newDataInd = (1:length(corner_coords_insert)) + corner_coords_insert';

                % Add in the corner boundary edges

                boundary_coords(newDataInd,:) = corner_coords2use(:,:);

            else

                % If there are no corners, just get the midpoints and use them

                Boundary_verts = mult_edges_vert_id(EdgeINDS,:);

                boundary_coords = (vertices(Boundary_verts(:, 1), :) + vertices(Boundary_verts(:, 2), :))./2;  

            end

            % The edges define a set of coordinates, and the order of the
            % defines how these coordinates need to be connected to form
            % the boundary. However the boundary needs to end where it
            % start so we make the first coordinate the final one as well
            
            BOUNDARY{nBounds} = [boundary_coords; boundary_coords(1,:)];
            
            % We don't know before hand how many distinct ROI boundaries
            % there will be, so we update it as we go.
            
            BOUNDARY_ROI_ID(nBounds) = i;

            nBounds = nBounds + 1;

            case 'centroid'
                
                % Find the pairs of faces making up the boundary

                Boundary_faces = mult_edges_face_id(EdgeINDS,:); 
                
                % Find the face that is not in the next pair of faces. This
                % is the face we need to calculate the coordinates at that
                % step

                Boundary_faces_next = [Boundary_faces(2:end,:); Boundary_faces(1,:)];

                Boundary_face_order_mask = any(Boundary_faces(:,1)==Boundary_faces_next,2);

                Boundary_face_order_masked = Boundary_faces;
                Boundary_face_order_masked([Boundary_face_order_mask ~Boundary_face_order_mask]) = 0;

                % Get the order of faces
                
                Boundary_face_order = max(Boundary_face_order_masked,[],2);
                
                % Get the vertex ids of those order faces
                
                Boundary_face_order_verts = faces(Boundary_face_order,:);
                
                % Get the face centroid to use as the coordinates for the
                % boundary

                Boundary_coords = (vertices(Boundary_face_order_verts(:, 3), :) + vertices(Boundary_face_order_verts(:, 2), :) + vertices(Boundary_face_order_verts(:, 1), :))./3;

                BOUNDARY{nBounds} = [Boundary_coords; Boundary_coords(1,:)];

                BOUNDARY_ROI_ID(nBounds) = i;

                nBounds = nBounds + 1;   

            end

        end

    end
    
    case 'edges'
        
            % Find the rois each face is connected to

            faces_roi_ids = vertex_id(faces);
            
            % Find the faces that are boundary faces

            bfaces = faces(logical(diff(faces_roi_ids,2,2)),:);
            
            % Find all possible combinations of edges that make up the
            % boundaries of ROIs
            
            BoundaryEdges = [[bfaces(:,1); bfaces(:,1); bfaces(:,2); bfaces(:,2); bfaces(:,3); bfaces(:,3)], ...
                [bfaces(:,3); bfaces(:,2); bfaces(:,1); bfaces(:,3); bfaces(:,2); bfaces(:,1)]];

            % Get the unique edges making up the boundary of each roi
            
            [EDGES_sorted] = sortrows(sort(BoundaryEdges,2));
            
            unique_boundary_edges = unique(EDGES_sorted,'rows');
            
            unique_boundary_edges_rois = vertex_id(unique_boundary_edges);
            
            % Get the number of unique ROIs
            
            roi_ids = unique(vertex_id);
            nrois = length(roi_ids);
            
            ROI_COMPONENTS = zeros(nrois,1);
            
            nBounds = 1;
            
            % Loop over each ROI
            
            for i = 1:nrois
                
                % Find the edges that connected two vertices of the same
                % ROI.
                
                index = sum(unique_boundary_edges_rois == roi_ids(i),2)==2;
                
                % Make an edge list of unique boundary edges for the
                % desired ROI
                
                edgeList_half = unique_boundary_edges(index,:);
                
                edgeList_orig = [edgeList_half; fliplr(edgeList_half)];
                
                unique_verts = unique(edgeList_orig(:));
                
                % Get the number of vertices that make up the boundary.
                
                nverts = length(unique_verts);
                
                % To make the next bit easier, relabel the vertex ids in 
                % edgeList_orig to be from 1:nverts
                
                new_vert_id = 1:nverts;
                old_vert_id = unique_verts;
                
                edgeList = edgeList_orig;

                for k = 1:numel(new_vert_id)
                    edgeList(edgeList_orig == old_vert_id(k)) = new_vert_id(k);
                end
                
                % Make a sparse adjacency matrix with the relabelled
                % vertices. if you graph this it should be a ring. This
                % ring indicates the order in which the vertices are
                % connected to define the boundary of the ROI. Basically
                % the nodes of this network represent vertices in the
                % original surface mesh
                                
                adj = sparse(edgeList(:, 1), edgeList(:, 2), 1, nverts, nverts);
                
            % Find the components of the adjacency matrix. Each component
            % represents a spatially continuous area for the current roi

            [comp,ROI_COMPONENTS(i)] = graphComponents(adj);

            % Loop over each component

            for j = 1:ROI_COMPONENTS(i)

                % Find the "nodes" making up the component

                Nodes2Use = find(comp == j);

                % Define a start "node"

                N1 = Nodes2Use(1);

                % Find the neighbours of the "node"

                NodeNeighbors = find(adj( N1,:));

                % Pick one of those neighbours to be the end "node"

                N2 = NodeNeighbors(1);

                % Find the connection linking the start and end "nodes". Once
                % found,then remove it from the adjacency matrix 

                adj(N1,N2) = 0;

                adj(N2,N1) = 0;

                % The boundary should be a cycle in the graph. If you remove
                % one connection between "node" N1 and N2, if you find the 
                % shortest path between N1 and N2, you get the ordering of
                % vertices which define the boundary

                G = graph(adj);

                TR = shortestpathtree(G,N1,N2,'OutputForm','cell');

                cycle = [TR{1} TR{1}(1)];

                % Insert the original vertex IDs, then extract the coordinates
                % of those IDs. This is the boundary of this component of
                % the ROI

                new_vert_id = 1:nverts;
                old_vert_id = unique_verts;

                boundary_verts = cycle;

                for k = 1:numel(new_vert_id)
                    boundary_verts(cycle == new_vert_id(k)) = old_vert_id(k);
                end

                BOUNDARY{nBounds} = vertices(boundary_verts,:);

                BOUNDARY_ROI_ID(nBounds) = i;

                nBounds = nBounds + 1;   
            end
                  
            end
    
    otherwise
        error('Unrecognised ''boundary_method'' option')
end

end
