

PARC = [lh_aparc Scha7_parcs.lh_scha100 lh_rand200 Scha7_parcs.lh_scha200 lh_HCPMMP1 lh_rand500 Scha7_parcs.lh_scha500];

for i = 1:size(PARC,2)
    Nrois = max(PARC(:,i));
    parc_dist_temp = zeros(length(PARC),1);
    for j = 1:Nrois
        parc_inds = PARC(:,i)==j;
        parc_verts = lh_inflated_verts(parc_inds,:);
        parc_Dists = squareform(pdist(parc_verts));
        mean_parc_Dists = mean(parc_Dists);
        [~,closest_vert] = min(mean_parc_Dists);
        dist2centre = parc_Dists(:,closest_vert);
       parc_dist_temp(parc_inds) = dist2centre;
    end
    parc_dist{i} = parc_dist_temp;
    parc_dist_norm{i} = parc_dist_temp./max(parc_dist_temp);
end

cmap = [lines(34); inferno(50); turbo(100); magma(100); viridis(180); plasma(250); parula(250)];

PARCnew(:,1) = PARC(:,1);
for i = 2:size(PARC,2)
PARCnew(:,i) = PARC(:,i)+max(PARCnew(:,i-1));
end

for i = 1:7
parc_dist2centre(:,i) = parc_dist{i};
end

PARCnew(PARCnew(:,1)==0,:) = 0;

parc = PARCnew(:,1);

surface.vertices = lh_inflated_verts;
surface.faces = lh_faces;
[p,b] = plotSurfaceROIBoundary(surface,parc,parc,'midpoint',cmap,2,[1 964]);

camlight(80,-10);
camlight(-80,-10);

view([270 0])

axis off
axis vis3d
axis tight
axis equal
%axis image

angles = [270:6:360 6:6:264];
COUNT = 0;
for i = [2:7 1]
   
    maxdist = max(parc_dist2centre(:,i));
    
    distThr = linspace(0,maxdist,60);
    
    for j = 1:60
        delete(b.boundary)
        addnewparc_ind = parc_dist2centre(:,i)<distThr(j);
        parc(addnewparc_ind) = PARCnew(addnewparc_ind,i);
       FaceVertexCData = makeFaceVertexCData(surface.vertices,surface.faces,parc,parc,cmap,[1 964]);
    p.FaceVertexCData = FaceVertexCData;
        view([angles(j) 0])
        BOUNDARY = findROIboundaries(surface.vertices,surface.faces,parc,'midpoint');
        for jj = 1:length(BOUNDARY)
           b.boundary(jj) = plot3(BOUNDARY{jj}(:,1), BOUNDARY{jj}(:,2), BOUNDARY{jj}(:,3), 'Color', 'k', 'LineWidth',2,'Clipping','off');
        end
        COUNT = COUNT + 1;
        pause(.1)
        print(['./GIF/',num2str(COUNT),'.png'],'-dpng')
    end
    
    
end
