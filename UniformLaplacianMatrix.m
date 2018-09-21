function L = UniformLaplacianMatrix(Mesh)
%UNIFORMLAPLACIANMATRIX Computes uniform Laplacian Matrix of a mesh
% L = D - A;
% where D is the degree matrix and A is the adjacency matrix
%   Detailed explanation: https://en.wikipedia.org/wiki/Laplacian_matrix

V = Mesh.vertices;
nV = size(V,1);
F = Mesh.faces;

%Preallocation
L = sparse(nV,nV);


%Compute degree matrix
D = zeros(nV,1);
for i=1:nV
    [rows,~] = find(F==i);
    tempRows = F(rows,:);
    D(i) = numel(unique(tempRows))-1;
end

D = spdiags(D,0,nV,nV);
if(sum(D < 0) > 0)
    warning("Mesh presents isolated vertices.");
end

%Compute adjacency matrix (from toolBox)
A = triangulation2adjacency(Mesh.faces,Mesh.vertices);

L = D-A;
end


