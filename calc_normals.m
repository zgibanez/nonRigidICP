function vertex_normal = calc_normals(mesh)
%CALC_NORMALS Calculate vertex normals of a mesh weighting its adjacent
%faces normals.

%First check if mesh has triangulated faces
faces = mesh.faces;
if(size(faces,2) ~= 3)
    error('Error: mesh does not have triangle faces.');
end

%% 1 Compute face normals:
nfaces = size(mesh.faces,1);
vertices = mesh.vertices;
face_normals = zeros(nfaces,3);


for i=1:nfaces
    %Every face is a triangle ABC
    a = vertices(faces(i,1),:);
    b = vertices(faces(i,2),:);
    c = vertices(faces(i,3),:);
    %The normal of the face will be (b-a)x(c-a). The magnitude is
    %proportional to the triangle area
    face_normals(i,:) = cross(b-a,c-a);
end


%% 2 Compute vertex normals by averaging adjacent face normals:
nvertices = size(vertices,1);
vertex_normal = zeros(nvertices,3);
for i= 1:nvertices
    %Find faces where this vertex appears
    [rows,~] = find(faces==i);
    
    %Compute the weighted average of normals
    if(numel(rows)>1)
        normal_dir = sum(face_normals(rows,:));
    else
        normal_dir = face_normals(rows,:);
    end
    
    try
        vertex_normal(i,:) = normal_dir / norm(normal_dir);
    catch
        vertex_normal(i,:) = [0 0 0];
    end
    %%TODO: Check incompatible normals
    if(sum(isnan(vertex_normal(i,:))) > 1)
        vertex_normal(i,:) = [0 0 0];
    end
    
%     normal_sum = 0;
%     normal_dir = [0,0,0];
%     for j=1:numel(rows)
%         normal_dir = normal_dir + face_normals(j);
%         normal_sum = normal_sum + norm(face_normals(rows(j)));
%     end
end

end

