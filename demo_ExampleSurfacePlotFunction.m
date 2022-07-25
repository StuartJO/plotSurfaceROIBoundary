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

cmap = turbo(256);

data_label = 'Functional gradient';

ExampleSurfacePlotFunction(surface,vertex_id,data,cmap,data_label);

print('./figures/Example4.png','-dpng','-r300')
