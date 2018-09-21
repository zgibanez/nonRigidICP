function p = disp_meshes(Mesh1, Mesh2, Mesh3,h)
%Target, Source, NewSource
na = nargin;
figHandles = findall(groot, 'Type', 'figure');

flag = 0;
for i=1:numel(figHandles)
    if(strcmp(figHandles(i).Tag,'h_mesh'))
        figure(figHandles(i));
        arrayfun(@cla,findall(0,'type','axes'));
        flag = 1;
        break;
    end
end

%If the figure window of display_meshes is not found, create one
if (~flag)
    h = figure();
    screensize = get( groot, 'Screensize' );
    set(gcf,'Position',[screensize(3)/4 0 900 900]);
    set(h,'Tag','h_mesh');
end
% 


switch na
    case 1
    p = plotM(Mesh1,'b',0.8);
    case 2
    p(1) = plotM(Mesh1,'b');
    hold on
    p(2) = plotM(Mesh2,'r');
    case 3
    subplot(2,2,1);    
    title('Original');
    plotM(Mesh1,'b');
    hold on
    p(1) = plotM(Mesh2,'r');
    subplot(2,2,2);
    plotM(Mesh3,'r');
    title('Transformed');
    p(2) = plotM(Mesh1,'b');
    subplot(2,2,3);
    plotM(Mesh1,'r');
    subplot(2,2,4);
    plotM(Mesh3,'r');
%     legend('3DScan','Template');
end

    drawnow;

end


function p = plotM(Mesh,color,alpha)
if(nargin < 3)
    alpha = 0.5;
end

    p = patch('Vertices',Mesh.vertices, ...
              'Faces', Mesh.faces, ...
              'facecolor', color, 'EdgeColor',  'none', ...
              'FaceAlpha', 0.5);
    material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
    view([60,30]); axis equal; axis manual;
    %legend('Source', 'Target', 'Location', 'best')
end