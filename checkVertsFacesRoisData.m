function [verts,faces,rois,data] = checkVertsFacesRoisData(varargin)
%% Checks size and shape of data describing brain surface mesh
% Ensures that verts is V x 3, faces is F x 3, rois is V x 1, and data is V x 1
% or R x 1
%
%
%% Syntax
%  [verts,faces,rois,data] = checkVertsFacesRoisData(verts,faces,rois,data)
%
%
%% Input Arguments
%  verts - xyz coordinates of vertices (three column matrix)
%  faces - triangulation ie. IDs of vertices making up each face (three column matrix) 
%  rois - roi allocation of each vertex (V x 1 vector)
%  data - data allocated to each vertex or roi (V x 1 or R x 1 vector)
%
%
%% Name-Value Arguments
%  checkContents - flag to check values within matrices (true (default) | false)
%  If `checkContents` is set to false, only the shape of the input matrices will
%  be tested. If `checkContents` is true, then the contents will also be checked
%  e.g. that `faces` uses all the vertices and that `rois` cannot have negative
%  values. 
%
%  fillEmpty - flag to populate `rois` and `data`, if the inputs are empty (false (default) | true) 
%  If set to true, rois will be set to all ones, and data will be set to the
%  second column of `verts`.
%
%
%% Output Arguments
%  verts, faces, rois, data - transposed if needed
% 
% 
%% See Also
% plotSurfaceROIBoundary, read_vtk
% 
% 
%% Authors 
% Mehul Gajwani, Monash University, 2023
% 
% 


%% Prelims
ip = inputParser;
validationFcn = @(x) (isnumeric(x) && ismatrix(x));
addOptional(ip, 'verts', [], @(x) validationFcn(x));
addOptional(ip, 'faces', [], @(x) validationFcn(x));
addOptional(ip, 'rois',  [], @(x) validationFcn(x));
addOptional(ip, 'data',  [], @(x) validationFcn(x));

addParameter(ip, 'checkContents', true, @islogical);
addParameter(ip, 'fillEmpty', false, @islogical);

parse(ip, varargin{:});
verts = ip.Results.verts;
faces = ip.Results.faces;
rois  = ip.Results.rois;
data  = ip.Results.data;


%% Check shape
% Check second dimension and transpose if needed

matrices = {verts, faces, rois, data};
names = {'Vertices', 'Faces', 'ROIs', 'Data'};
target2ndDim = {3, 3, 1, 1};

for ii = 1:length(matrices)
    if isempty(matrices{ii})
        % skip - do nothing
    elseif size(matrices{ii}, 2) == target2ndDim{ii}
        % good - do nothing
    elseif size(matrices{ii}, 1) == target2ndDim{ii}
        matrices{ii} = matrices{ii}.';
    else
        error('%s must have a dimension of size %d', names{ii}, target2ndDim{ii});
    end
end

[verts, faces, rois, data] = deal(matrices{:});

if isallhere({verts, rois})
    assert(size(verts, 1) == size(rois, 1), ...
        'Vertices and ROIs must have the same number of points');
end


%% Fill empty, if desired
if ip.Results.fillEmpty
    if isempty(rois); rois = ones(size(verts, 1),1); end
    if isempty(data); data = verts(:, 2); end
end


%% Check contents
if ip.Results.checkContents
    if isallhere({verts, faces})
        if max(faces, [], "all") ~= size(verts, 1)
            warning('The vertices are not all mapped to the faces');
        end
    end

    if isallhere({faces}); mustBeNonnegative(faces); mustBeInteger(faces); end

    if isallhere({rois})
        mustBeNonnegative(rois); mustBeInteger(rois);

        goodRois = ( all(rois) && (length(unique(rois))==max(rois)) ) || ... % no non-zeros rois
            (~all(rois) && (length(unique(rois))==max(rois)+1) ); % some non-zero rois

        if ~goodRois
            warning('ROIs appear to not be sequential');
        end
    end

    if isallhere({data, verts}) || isallhere({data, rois})
        assert(size(data, 1) == size(verts, 1) || size(data, 1) == max(rois), ...
            'Data should be one per vertex or one per ROI');
    end
end

end % main

function out = isallhere(inp)
% returns true if inp is (i) a non-empty matrix OR (ii) a cell with ALL non-empty entries
% else returns false
if ~iscell(inp); out = ~isempty(inp); return; end
out = ~any(cellfun(@isempty, inp));
end



