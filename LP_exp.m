function [A,C] = LP_exp(Mesh)
    V = Mesh.vertices;
    F = Mesh.faces;
   
    %Calculate mass matrix

    A = massMatrix(V,F);
    C = weakLaplacianMatrix(V,F);
    
end

%Calculate mass matrix
function A = massMatrix(V,F)
nV = size(V,1);

%Vector of areas
Av = zeros(size(V,1),1);
for i=1:nV
    

    %Find connected faces to the vertex
    [connected_faces, ~] = find(F==i);
    nC = size(connected_faces,1); %Number of faces connected to the vertex
    temp_area = 0;
    
    %For every face connected to the vertex
    for n=1:nC
        
        %Take the vertices of every adjacent face
        neigh_vertex = F(connected_faces(n),:)~=i;
        neigh_vertex = F(connected_faces(n),neigh_vertex);
        triangle = [V(i,:); V(neigh_vertex,:)];
        P = triangle(1,:); Q = triangle(2,:); R = triangle(3,:);
        
        %Check if obtuse
        if(~isTObtuse(P,Q,R))
            temp_area = temp_area + vorAreaT(P,Q,R);
        else
            temp_area = temp_area + simAreaT(P,Q,R);
        end
    if(imag(temp_area)>0)
        warning('Imaginary component found');
    end
    end
    Av(i) = temp_area;
    
end
    %The result is a diagonal mass matrix with all the mixed VoronoiAreas
    %(notice that it is inverted as in Laplace-Beltrami: The Swiss Army Knife of Geometry Processing)
    A = spdiags(1./Av,0,nV,nV);
end

%Weak Laplacial matrix based on cotangent discretization
function C = weakLaplacianMatrix(V,F)
    nV = size(V,1);
    C = sparse(nV,nV);
    
    %For every vertex
    for i=1:nV
        
        %Find connected faces to the vertex
        [connected_faces, ~] = find(F==i);
        cF = size(connected_faces,1);
        faces = F(connected_faces,:);
        
        %For every adjacent face of the mesh
        for n=1:cF
            % Split the neighbouring faces into two matrices, one vector
            % containing the actual face and one matrix containing the rest
            rest_faces = faces;
            rest_faces(n,:) = [];
            actual_face = faces(n,:);

            %Of the actual face, find vertices different from the one being
            %studied (there should be two)
            idx = find(actual_face~=i);

            %Find the other faces that contain these two vertices
            v2 = i;
            [row,~] = find(rest_faces == actual_face(idx(1)));
            if (row ~= -1)
                
                v1_1 = rest_faces(row, (rest_faces(row,:) ~= i & rest_faces(row,:) ~= actual_face(idx(1))));
                v1_2 = actual_face(idx(2));
                v3 = actual_face(idx(1));
                tri1 = V([v1_1 v2 v3],:); %Triangle number 1
                tri2 = V([v1_2 v2 v3],:); %Triangle number 2
                c1 = cotg(tri1(2,:)-tri1(1,:),tri1(3,:)-tri1(1,:));
                c2 = cotg(tri2(2,:)-tri2(1,:),tri2(3,:)-tri2(1,:));
                C(i,v3) = (c1+c2)/2;
            end
%                 %If one triangle is isolated (has no two triangles asociated with the
%                 %two vertex , compute the uniform weight instead
%                 C(i,v3) = 1/cF;


            [row,~] = find(rest_faces == actual_face(idx(2)));
            if (row ~= -1)
                v1_1 = rest_faces(row, (rest_faces(row,:) ~= i & rest_faces(row,:) ~= actual_face(idx(2))));
                v1_2 = actual_face(idx(1));
                v3 = actual_face(idx(2));
                tri1 = V([v1_1 v2 v3],:); %Triangle number 1
                tri2 = V([v1_2 v2 v3],:); %Triangle number 2
                c1 = cotg(tri1(2,:)-tri1(1,:),tri1(3,:)-tri1(1,:));
                c2 = cotg(tri2(2,:)-tri2(1,:),tri2(3,:)-tri2(1,:));
                C(i,v3) = (c1+c2)/2;
%             else
%                 C(i,v3) = 1/cF;

            end
        end   
    end
    
    for i=1:nV
        C(i,i) = -sum(C(i,:));
    end
end

%Cotangent formula
function c = cotg(u,v)
    c = norm(dot(u,v))/norm(cross(u,v));
end

%Is the triangle obtuse?
function is_obtuse = isTObtuse(P,Q,R)


norms = [norm(P-Q) norm(P-R) norm(R-Q)];
[maxN,idx] = max(norms);
norms(idx) = [];

if(maxN^2 > sum(norms.^2))
    is_obtuse = true;
else
    is_obtuse = false;
end

end

%From: Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
function area = vorAreaT(P,Q,R)
    area =  1/8 * (norm(P-R)^2*cotg(P-Q,R-Q)+norm(P-Q)^2*cotg(P-R,Q-R));
end

% From https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
function area = simAreaT(P,Q,R)
v1 = R-P; v2 = Q-P;
N1 = norm(R-P); N2 = norm(Q-P);
sinT = sqrt(1-(dot(R-P,Q-P)/(N1*N2))^2);
%Check if the angle at P is obtuse
if(dot(Q-P,R-P)<0)
    %The mixed voronoi area of an obtuse triangle is the third of the standard
    %triangle area
%     area = (1/2*N1*N2*sinT)/2;
    area = norm(cross(v1,v2))/2;
else
%     area = (1/2*N1*N2*sinT)/4;
    area = norm(cross(v1,v2))/4;
end

end
