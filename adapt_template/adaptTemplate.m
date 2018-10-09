function [vertices, faces] = adaptTemplate(Target,Source, pTarget, pSource)
%Performs smooth deformations on template to match the target
% Input:
% -Target,Source: Struct with fields 'vertices' and 'faces'
% -pTarget,pSource: Vector of indexes of points with known correspondences

nV = size(Source.vertices,1);

% checkLandmarks(Target,Source,pTarget,pSource);

%Laplace beltrami coeff
[Am,Bm,~,~] = LaplaceBeltramiCoefficients(Source);

%The landmark term uses Euclidean distance norm(v1h-v1)^2. Its derivate
%coefficients are 2d(vih) = 2(vi);
Al = sparse(nV,3);
Al(pSource,:) = 2;
Al = reshape(Al',1,[])';
Al = spdiags(Al,0,3*nV,3*nV);

b_land = sparse(nV,3);
b_land(pSource,:) = 2*Target.vertices(pTarget,:);
b_land = reshape(b_land',1,[])';

%Scheduled weights
w = linspace(1000,100,10);

newSource = Source;

for i=1:numel(w)
    A = Al+w(i)*Am;
    b = b_land + w(i)*Bm;
    
    V = A\b;
    newSource.vertices = reshape(V,3,[])';
%     disp_meshes(newSource);
%     pause();
end

vertices = newSource.vertices;
faces = newSource.faces;

end



