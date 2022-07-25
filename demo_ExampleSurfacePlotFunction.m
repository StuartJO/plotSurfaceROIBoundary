% This script just gives some examples of how the code can be used

load('surface_data.mat')
load('example_data.mat')

%% Accentuate plotted data

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

vertex_id = lh_HCPMMP1;

data = zscore(lh_HCPMMP1_gene_pc1);

cmap = turbo(256);

data_label = 'Gene PC1';

ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label);

print('./figures/Example1.png','-dpng','-r300')

%% Plot parcellation over continuous data

surface.vertices = lh_verts;
surface.faces = lh_faces;

vertex_id = lh_rand200;

data = lh_sulc;

cmap = parula(256);

data_label = 'Sulcal depth';

ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label);

print('./figures/Example2.png','-dpng','-r300')

%% Plot parcellation over another parcellation

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

vertex_id = lh_rand500;

data = lh_Yeo_7net;
data(data==0) = NaN;

cmap = Yeo_7net_cmap;

data_label = 'Yeo network';

[~,~,c] = ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label);

% Make the colorbar look less weird
caxis([nanmin(data)-.5 nanmax(data)+.5])
c.Ticks = 1:7;

print('./figures/Example3.png','-dpng','-r300')

%% Plot data that has been thresholded

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

vertex_id = lh_rand500;

data = lh_func_grad1;
% Threshold to only display data 1 SD above the mean
data(data<(mean(data) +std(data))) = NaN;

% Any data which has a NaN value will not be assigned a colour from the
% colour map. Instead by default it will be displayed as grey.

cmap = turbo(256);

data_label = 'Functional gradient';

ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label);

print('./figures/Example4.png','-dpng','-r300')

%% Threshold data on a per region basis

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

vertex_id = lh_rand500;

data = zeros(length(unique(vertex_id))-1,1);
for i = 1:length(data)
    data(i) = mean(lh_func_grad1(vertex_id==i));
end
    
% Threshold to only display data 1 SD above the mean
data(data<(mean(data) +std(data))) = NaN;

% Any data which has a NaN value will not be assigned a colour from the
% colour map. Instead by default it will be displayed as grey.

cmap = turbo(256);

data_label = 'Functional gradient';

ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label);

print('./figures/Example5.png','-dpng','-r300')

%% Threshold data and borders
% This will only work if 'data' is plotted continuously on the surface
% Why you would want to do this, I have no idea but here you go anyway!

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;

vertex_id = lh_rand500;

data_temp = zeros(length(unique(vertex_id))-1,1);
for i = 1:length(data_temp)
    data_temp(i) = mean(lh_func_grad1(vertex_id==i));
end
    
% Threshold to only display data 1 SD above the mean
data_temp(data_temp<(mean(data_temp) +std(data_temp))) = NaN;

% Plot the data back to the surface
data = nan(size(vertex_id));
for i = 1:length(data_temp)
    data(vertex_id==i) = data_temp(i);
end

vertex_id(ismember(vertex_id,find(isnan(data_temp)))) = 0;

% If you set a regions vertex_id to 0, then it won't appear to be plotted 
% (technically is still is plotted but as an unknown region)

cmap = turbo(256);

data_label = 'Functional gradient';

ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label);

print('./figures/Example6.png','-dpng','-r300')