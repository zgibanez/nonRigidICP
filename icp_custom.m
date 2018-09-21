function [s, R, t, err, newP] = icp_custom(Target,Source)
%Returns the cloud of 3D points after initial alignment
    
    %If Target and Source are given as struct take the vertices field.
    %If not, they are assumed to be a set of vertices
    if(isstruct(Target))
        M = Target.vertices;
        P = Source.vertices;
        %For display purposes
        SourceT = Source;
    else
        M = Target;
        P = Source;
        f = figure();
    end


    %Init output args
    s = 1;
    R = speye(size(M,1));
    t = sparse(size(M,1),1);
    newP = P;

    %Algorithm parameters
    min_err =  1;
    old_error = Inf;
    new_error = 0;
    iter = 0;
    
    % Number of points
    Np = size(P,1);
%     Nm = size(M,1);
%     dim = size(P,2);
    
    while (abs(old_error-new_error) >= min_err)
        %find correspondences
%         Y = zeros(Np,dim);
        U = knnsearch(M,newP);
        Y = M(U,:);
        [s, R , t, err] = find_alignment(newP, Y);
        
        %apply alignment and compute error
        total_err = 0;
        for i= 1:Np
            newP(i,:) = (s*R*newP(i,:)' + t')';
            e = Y(i,:) - newP(i,:);
            total_err = total_err + e*e';
        end
        
        old_error = new_error;
        new_error = total_err;
        fprintf('Iteration number %d: residual error is %f... \n',iter,(old_error-new_error));
        iter = iter+1;
        
        
        if(isstruct(Target) && isstruct(Source))
            SourceT.vertices = newP;
            disp_meshes(Target,Source,SourceT);
        else
            clf(f);
            subplot(1,2,1);
            scatter3(Target(:,1),Target(:,2),Target(:,3)); hold on;
            scatter3(Source(:,1),Source(:,2),Source(:,3)); hold on;
            axis equal;
            
            subplot(1,2,2);
            scatter3(Target(:,1),Target(:,2),Target(:,3)); hold on;
            scatter3(newP(:,1),newP(:,2),newP(:,3)); hold on;
            axis equal;
        end
        
        if err < min_err
            break;
        end
    end
    
end

function [s,R,t,err] = find_alignment(P,Y)
% Computes scaling factor, rotation and translation for the alignment
% of two sets of corresponding 3D points Y and P such that:
% 
%     Yi = s*R*Pi + t
%     
% Input :
%           P - 3xN matrix with 3D points
%           Y - 3xN matrix with 3D points
%           
%   Returns :
%           s - scaling factor
%           R - 3x3 rotation matrix
%           t - 3x1 translation vector
%           err - residual error after the alignment

%%Test if input args are correct
[Np, dim_p] = size(P);
[Ny, dim_y] = size(Y);

if(Np ~= Ny)
    error('Point sets need to have the same number of poins');
end

if(dim_p ~= 3 || dim_y ~= 3)
    error('Poins need to have 3 dimensions');
end

if(Np < 4)
    error('At least 4 pairs of points are needed to perform ICP');
end

%Number of points
N = Np;

%% Compute centroids
Mu_p = mean(P,1);
Mu_y = mean(Y,1);

% Express points relative to their centres
Pp = P - repmat(Mu_p,N,1);
Yp = Y - repmat(Mu_y,N,1);

%% Compute the optimal quaternion
Px = Pp(:,1); Py = Pp(:,2); Pz = Pp(:,3);
Yx = Yp(:,1); Yy = Yp(:,2); Yz = Yp(:,3);

Sxx = sum(Px.*Yx); Sxy = sum(Px.*Yy); Sxz = sum(Px.*Yz);
Syx = sum(Py.*Yx); Syy = sum(Py.*Yy); Syz = sum(Py.*Yz);
Szx = sum(Pz.*Yx); Szy = sum(Pz.*Yy); Szz = sum(Pz.*Yz);

Nmatrix = [ Sxx + Syy + Szz,    Syz-Szy,            -Sxz + Szx,         Sxy - Syx;
            -Szy + Syz,         Sxx - Szz - Syy,    Sxy + Syx,          Sxz + Szx;
            Szx - Sxz,          Syx + Sxy,          Syy - Szz - Sxx,    Syz + Szy;
            -Syx + Sxy,         Szx + Sxz,          Szy + Syz,          Szz - Syy - Sxx];

%Compute eigenvalues
[V, ~] = eig(Nmatrix);

q = V(:,4);

%%Compute rotation matrix
q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

Qbar = [ q0 -q1 -q2 -q3 ;
         q1  q0  q3 -q2 ;
         q2 -q3  q0  q1 ;
         q3  q2 -q1  q0];
     
Q =    [ q0 -q1 -q2 -q3 ;
         q1  q0 -q3  q2 ;
         q2  q3  q0 -q1 ;
         q3 -q2  q1  q0];
     
R = Qbar'*Q;
R = R(2:4,2:4);

%%Compute the scaling factor
Sp = 0;
D = 0;
for i=1:N
    D = D + Yp(i,:)* Yp(i,:)';
    Sp = Sp + Pp(i,:)* Pp(i,:)';
end

% s = sqrt(D/Sp);
s = 1; %%SCALING NOT WORKING AS INTENDED

%%Compute translational offset
t = (Mu_y' - s*R*Mu_p')';

%%Compute residual error
err = 0;
for i = 1:N
    d = Y(i,:)'-(s*R*P(i,:)'+t');
    err = err + norm(d);
end

end