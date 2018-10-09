function checkLandmarks(Target,Source,pTarget,pSource)

f = findobj('Tag','landmark_figure');
if(isempty(f))
    f = figure();
    set(f,'Tag','landmark_figure');
    title('Landmarks');
end
    h = f.CurrentAxes;

%Prealignment ensures both meshes will be well vissible    
[~,Source.vertices] = prealign(Target,Source);
offset = repmat([150, 0, 0], size(Source.vertices,1),1);
pS = Source.vertices + offset;
pT = Target.vertices(pTarget,:);
Q = pS(pSource,:)-pT;

trisurf(Target.faces, Target.vertices(:,1),Target.vertices(:,2),Target.vertices(:,3),'FaceAlpha',0.7,'FaceColor','r','Edgecolor','none'); 
hold
trisurf(Source.faces, pS(:,1),pS(:,2),pS(:,3),'FaceAlpha',0.7,'FaceColor','b','Edgecolor','none'); 
light
lighting phong;
set(gca, 'visible', 'off')
view(90,90)
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
%     patch(h, 'Vertices',Target.vertices, 'Faces', Target.faces, ...
%               'facecolor', 'r', 'EdgeColor',  'none', ...
%               'FaceAlpha', 0.6);
%          
%     hold on;
%     
%         patch(h, 'Vertices',pS, 'Faces', Source.faces, ...
%               'facecolor', 'b', 'EdgeColor',  'none', ...
%               'FaceAlpha', 0.6);
    
    scatter3(Target.vertices(pTarget,1),Target.vertices(pTarget,2),Target.vertices(pTarget,3),20,'filled','y');
    scatter3(pS(pSource,1),pS(pSource,2),pS(pSource,3),20,'filled','y');

    quiver3(pT(:,1),pT(:,2),pT(:,3),Q(:,1),Q(:,2),Q(:,3),'AutoScale','off','LineWidth',1);
    cameratoolbar('Show');
    %legend('Source', 'Target', 'Location', 'best')
end
