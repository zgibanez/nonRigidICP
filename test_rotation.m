
LB = abs(D.LB_vec{1});
V = Source.vertices;
scatter3(V(:,1),V(:,2),V(:,3),1,LB);
% patch('Vertices',Source.vertices,'Faces',Source.faces,'FaceAlpha',0.3,'CData',LB);
axis equal;


% 
% a = 1;
% foo();
% b = 2;
% 
% 
% function foo()
% figure
% ax1 = subplot(2,1,1);
% [X1,Y1,Z1] = peaks;
% surf(X1,Y1,Z1)
% 
% ax2 = subplot(2,1,2);
% [X2,Y2,Z2] = peaks(10);
% surf(X2,Y2,Z2)
% 
% hlink = linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'}); 
% rotate3d on
% end