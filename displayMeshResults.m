function displayMeshResults(Data)
figure();
cameratoolbar('Show');

    h1 = subplot(1,3,1);
    dispPatchMeshes(h1,Data.Source_Aligned,'b');
    axis vis3d
    
    h2 = subplot(1,3,2);
    dispPatchMeshes(h2,Data.SourceT{end},'b');
    axis vis3d

    h3 = subplot(1,3,3);
    dispPatchMeshes(h3,Data.Target,'b');
    axis vis3d


    Link = linkprop([h1,h2,h3],{'CameraPosition','CameraUpVector'});
    setappdata(gcf, 'StoreTheLink', Link);
end



function dispPatchMeshes(h,Mesh,color,alpha,edge_alpha)
if(nargin < 3)
    color = 'b';
end
if(nargin < 4)
    alpha = 0.5;
end
if(nargin < 5)
    edge_alpha = 0.0;
end
axes(h);
Mesh = rmfield(Mesh, 'normals');
    patch(Mesh, 'facecolor', color, 'EdgeColor',  'k', ...
                                    'EdgeAlpha', edge_alpha, ...
                                    'FaceAlpha', alpha);
material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
axis equal; axis manual;
%view([60,30]);
hold on;
end