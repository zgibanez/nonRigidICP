function [LB, b,A,C] = LaplaceBeltramiCoefficients(Mesh)

    %Obtain the Laplace-Beltrami derivate expression for the mesh
    [A,C] = LP_exp(Mesh);
    L = A*C;
    
    %Coefficients of derivatives
    LB = 2*(L'*L);
    
    %Equations are repeated for x,y and z coordinates of vertices
    LB = kron(LB,eye(3));
    
    %The independent terms are just the same but multiplied by the actual
    %vertices of the mesh (knowns)
    b = reshape(Mesh.vertices',1,[]);
    b = LB *b';
end