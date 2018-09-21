function [A,b] = local_rigidity_term(Mesh)

%Extracts derivative coefficients od non rigidity term
nV = size(Mesh.vertices,1);
Adj = meshAdjacencyMatrix(Mesh.faces);
D = sum(Adj,2);
D = -spdiags(D,0,nV,nV);
A = (D+Adj);

A = -2*kron(A,eye(3)); %Repeat each row 3 times for the x,y and z coordinates
v_flat = reshape(Mesh.vertices',1,[])'; %Flatten the vertices matrix into a column vector

b = 2*sum(kron(D+Adj,eye(3))*v_flat,2);
%b = 2*(kron(D,ones(3,1))-bm);

end
