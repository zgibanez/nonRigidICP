
function disp_meshes_phong(Target,Source)
trisurf(Target.faces, Target.vertices(:,1),Target.vertices(:,2),Target.vertices(:,3),'FaceAlpha',0.7,'FaceColor','b','Edgecolor','none'); 
hold
if(nargin ==2)
trisurf(Source.faces, Source.vertices(:,1), Source.vertices(:,2), Source.vertices(:,3),'FaceAlpha',0.7,'FaceColor','b','Edgecolor','none');
end
light
lighting phong;
set(gca, 'visible', 'off')
view(90,90)
set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
end