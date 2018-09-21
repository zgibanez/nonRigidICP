function V = compute_boundary_custom(faces)

%Retrieves boundary vertex indexes of the mesh described by faces


if(size(faces,1) < size(faces,2))
    faces = faces';
end    

nfaces = size(faces,1);
nvertex = max(max(faces)); %The maximum index coincides with the vertex number

A = sparse(nvertex,nvertex);
for i=1:nfaces
    f = faces(i,:);
    A(f(1),f(2)) = A(f(1),f(2)) + 1;
    A(f(1),f(3)) = A(f(1),f(3)) + 1;
    A(f(3),f(2)) = A(f(3),f(2)) + 1;
end
A = A+A';

%Find vertices with only one edge
[row,col,s] = find(A==1);
V = unique(row);

end