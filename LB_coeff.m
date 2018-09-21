function [C,b] = LB_coeff(Source)
%Simplified expression of LB with fixed cotan weights

    %Output: 
    % C- cotan weights
    % b- independent coefficients
    
    [A,C] = LP_exp(Source);
    L = A*C;

    %Equations are repeated for x,y and z coordinates of vertices
    LB = kron(L,ones(3));
    
    % Coefficients of derivatives from w*(vh_i-v_i)^2
    c = LB*reshape(Source.vertices',1,[])';
    cs = 2*c.*c;
    b = -cs.*reshape(Source.vertices',1,[])';
    
    %Reshape so they fit with the six other params
    b = [b;zeros(6,1)];
    cs = [cs; zeros(6,1)];
    
    nV = size(cs,1);
    C = spdiags(cs,0,nV,nV);
    
    
end