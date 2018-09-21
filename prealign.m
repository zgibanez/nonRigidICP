function [NTarget, NSource] = prealign(Target,Source,display)

if(nargin<3)
    display = false;
end

%Provides a coarse alignment of the meshes overlapping their centroids
%and aligning their PCA directions.
Nvertex = size(Source.vertices,1); %Number of vertex in Source mesh

%Align the two 3D vectors!
% Ref: https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0;
              norm(cross(A,B)) dot(A,B)  0;
              0              0           1];

FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
UU = @(Fi,G) Fi*G*inv(Fi);

%Align the centroid
PreSource = Source;
PreSource.vertices = PreSource.vertices - mean(PreSource.vertices);

PreTarget = Target;
PreTarget.vertices = Target.vertices - mean(Target.vertices);
coefT = pca(PreTarget.vertices);

iter = 3;

for i=1:iter
    coefS = pca(PreSource.vertices);

    U = UU(FFi(coefS(:,i),coefT(:,i)), GG(coefS(:,i),coefT(:,i)));
%     [warnMsg, warnId] = lastwarn;
%     if ~isempty(warnMsg)
%         pause();
%     end
%     disp(norm(coefT(:,i)-U*coefS(:,i)));

    %Apply rotation to each vertex
    for j=1:Nvertex
        PreSource.vertices(j,:) = (U*PreSource.vertices(j,:)')';
    end
    %quiver3(0,0,0,coefS(1,1)*100,coefS(2,1)*100,coefS(3,1)*100,'LineWidth',3);
    %disp_meshes(Source,PreSource,Target);
end
% disp_meshes(Source,PreSource,PreTarget);

if(display)
         figure();
     subplot(1,2,1);
     Tq = U*coefS(:,1);
         scatter3(PreSource.vertices(:,1),PreSource.vertices(:,2),PreSource.vertices(:,3));
         hold on;
     quiver3(0,0,0,coefS(1,1)*100,coefS(2,1)*100,coefS(3,1)*100,'LineWidth',3);
     hold on;
     quiver3(0,0,0,Tq(1)*100,Tq(2)*100,Tq(3)*100,'LineWidth',3);
     quiver3(0,0,0,coefT(1,1)*100,coefT(2,1)*100,coefT(3,1)*100,'LineWidth',3);
     axis equal;
         hold off;
         subplot(1,2,2);
         scatter3(PreTarget.vertices(:,1),PreTarget.vertices(:,2),PreTarget.vertices(:,3));
         hold on;
         quiver3(0,0,0,coefT(1,i)*100,coefT(2,i)*100,coefT(3,i)*100,'LineWidth',3);
         quiver3(0,0,0,coefS(1,i)*100,coefS(2,i)*100,coefS(3,i)*100,'LineWidth',3);
         hold off;
end

NSource = PreSource.vertices;
NTarget = PreTarget.vertices;

end

