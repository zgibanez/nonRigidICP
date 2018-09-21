function Data = nricp_custom(Target,Source, preal,rigidAl,pTarget,pSource)
%Non rigid registration algorithm to math P to Q
% input:
% source,target - triangle meshes with vertices and faces indexed
%
% output:
% Data - Struct with all the data fields
% F - faces of transformed source


if(preal)
% 1: Perform coarse alignment
   fprintf("Performing prealignment between Target and Source mesh... ");
   [Target.vertices,Source.vertices] = prealign(Target,Source);
   fprintf("done. \n");
   Data.Source_Original = Source;
end
if(rigidAl)
%% 2: Perform rigid ICP
   fprintf("Performing rigid allignment between Target and Source mesh... \n");
  [s, R, t, err, newP] = icp_custom(Target,Source);
%    Source_aligned = Source;
  Source.vertices = newP;
  fprintf("done. \n");
end 

%% 3: Adapt template to user geometry
    if(nargin == 5)
        [Source.vertices, ~] = adaptTemplate(Target,Source,pTarget,pSource);
        Data.Source_Warped = Source;
    end

%% 4: Perform non-rigid ICP
%    rigid_icp(Target,Source_aligned);
   Data = nonrigid_icp(Target,Source);
   Data.Source_Aligned = Source;
   Data.Target = Target;
 
end

function [V,F] = Cnon_rigid_ICP(Target,Source)

    %Calculate vertex normals of Source and target
    Tnormals = calc_normals(Target);
%     Snormals = calc_normals(Source);

    %Store vertex
    Svertex = Source.vertices;
    Tvertex = Target.vertices;
    nS = size(Source.vertices,1);
    
    %Stablish weights
    w1 = 10; w2 = 100;
    
    %And initial parameters
    alpha = 0; beta = 0; gamma = 0;
    txc = 0; tyc = 0; tzc = 0;
    Tc = [txc tyc tzc];
    T = [0 0 0];
    Rt = [1 -gamma beta; gamma 1 -alpha; -beta alpha 1];
    
    %Model matrices for easier implementation of perturbations
    Ra = [0 0 0; 0 0 -1; 0 1 0];
    Rb = [0 0 1; 0 0 0; -1 0 0];
    Rg = [0 -1 0; 1 0 0; 0 0 0];
    R = eye(3);
    eps = 1e-3;
    veps = 1; %a bigger epsilon for vertex measured in mm
    mu = 1e-3; %Step size for steepest descent
    
    %Preallocating variables
    ha = zeros(nS,1); hb = zeros(nS,1); hg = zeros(nS,1);
    htx = zeros(nS,1); hty = zeros(nS,1); htz = zeros(nS,1);

    err = zeros(1,nS);
    new_error = Inf;
    old_error = 0;
    
    %Energy function as described in the Non-rigid registration course
    %notes by Bouaziz, Pauly, Hao Li and Taylor (2017)
    energyfun = @(Tv,Tn,Sv,Rx,Tx) w1*norm(dot(Tn,Sv-Tv))^2 + w2*norm(Sv'-Rx*(R*Sv'+T')+Tx')^2;
    count = 0;
    
    %For visualizing results
    newSource = Source;
    
    %while( norm(old_error-new_error) < 1e-4 )
     while(1)   
        U = knnsearch(Target.vertices,newSource.vertices);
        
        %Optimize for rotational angles
        % Perturbations are added with a sum of the rotation matrix
        % and the model matrix: Rn =  R + Rx*eps
        for i =1 : nS
            ha(i) = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),eye(3)+eps*Ra,Tc) - energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),eye(3)-eps*Ra,Tc);
            hb(i) = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),eye(3)+eps*Rb,Tc) - energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),eye(3)-eps*Rb,Tc);
            hg(i) = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),eye(3)+eps*Rg,Tc) - energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),eye(3)-eps*Rg,Tc);
        end
        alpha = alpha - mu*sum(ha)/(2*eps);
        beta  = beta - mu*sum(hb)/(2*eps);
        gamma = gamma - mu*sum(hg)/(2*eps);
        Rt = [1 -gamma beta; gamma 1 -alpha; -beta alpha 1];

        R = Rt*R;
        
        %Optimize for translational offset
        % Perturbations are added with a sum of the rotation matrix
        % and the model matrix: Rn =  R + Rx*eps
        for i =1 : nS
            htx(i) = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),R, Tc+eps*[1 0 0]) - energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),R, Tc-eps*[1 0 0]);
            hty(i) = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),R, Tc+eps*[0 1 0]) - energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),R, Tc-eps*[0 1 0]);
            htz(i) = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),R, Tc+eps*[0 0 1]) - energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),R, Tc-eps*[0 0 1]);
        end
        T(1) = T(1) - mu*sum(htx)/(2*eps);
        T(2) = T(2) - mu*sum(hty)/(2*eps);
        T(3) = T(3) - mu*sum(htz)/(2*eps);
        
        %Optimize for vertex positions
        for i=1:nS
%             vperturbm = Svertex(i,:) - Tnormals(U(i),:)*veps;
%             vperturbp = Svertex(i,:) + Tnormals(U(i),:)*veps;
%             errm = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),vperturbm, R, Tc);
%             errp = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),vperturbp, R, Tc);
%             h = (errm - errp) / 2*eps;
            h = 0;
            Svertex(i,:) = (R*(Svertex(i,:)' - mu*h)+T')';
        end
        
        for i=1:nS
            err(i) = energyfun(Tvertex(U(i),:),Tnormals(U(i),:),Svertex(i,:),R,T);
        end
        new_error = sum(err);
        
        fprintf('Iteration number %d. Total error: %f \n',count, new_error);
        count = count+1;
        newSource.vertices = Svertex;
        disp_meshes(Target,Source,newSource);
    end
    
    % Return deformed source mesh
    V = 1;
    F = 1;
    
    newSource.vertices = Svertex;
    
    
end

function Data = nonrigid_icp(Target,Source)

%Prepare data output
Data.TSource = 0;
Data.correspondences = zeros(size(Source.vertices,1),1);
Data.LB_vectors = zeros(size(Source.vertices,1),3);

%Calculate vertex normals of Source and target
fprintf("Calculating Target mesh normals...");
Tnormals = calc_normals(Target);
fprintf(" done. \n");

%Placeholder for transformed source
newSource = Source;

%Store vertex
% Svertex = Source.vertices;
% Tvertex = Target.vertices;
nS = size(Source.vertices,1);

%Initialize translational and rotational component
R = eye(3);
T = [0 0 0];

% energyfun = @(Tv,Tn,Sv,Rx,Tx) w1*norm(dot(Tn,Sv-Tv))^2 + w2*norm(Sv'-Rx*(R*Sv'+T')+Tx')^2;

%Weights
%TODO: Schedule weight decay, e.g., w1= 100:-10:10
steps = 1;
w1 = 10;
%w2 = 5;
w2 = linspace(50,5,steps);
%w2 = zeros(steps,1);
w3 = linspace(50,5,steps);
%w3 = zeros(steps,1);
w4 = linspace(20,1,steps); %local rigidity weight
%w4 = zeros(steps,1);
eps = 1;
%%Original w1 = 10, w2 = 50-10; w3 = 50-5; w4 = 20-5; 
%Options for trimming correspondences
Options.normalTrim = true;
Options.boundaryTrim = true;
Options.distanceTrim = true;
Options.LBmatching = false; %Better if false
Options.NormalMatching = true; %Better if true
Data.Options = Options;

%Precalculate Laplace-Beltrami coefficient matrices
fprintf("Calculating Laplace-Beltrami vectors of Source mesh...");
% [lC,lb] = LB_coeff(Source);
[Al,bm,A0,C0] = LaplaceBeltramiCoefficients(Source);
fprintf(" done. \n");

% %Add extra columns to prevent dimension mismatches
Al = [Al zeros(nS*3,6)];
Al = [Al ; sparse(6,size(Al,2))];
bm = [bm; zeros(6,1)];

%LB coefficients for checking fixed cotan weights assumpition
% LB = zeros(nS,1);
fprintf("Calculating Laplace-Beltrami vectors of Target mesh...");
[AT,CT] = LP_exp(Target);
fprintf(" done. \n");
L_target = AT*CT;
L_source = A0*C0;
L_source0 = L_source;

%Parameter for Thikonov regularization
lambda = 100;

%Iteration counter
iter = 1;

%Figure handle for plotting
% h = figure();

%For every scheduled weight
for n=1:steps
    
    old_error = Inf;
    new_error = 10e6;
    in_iter = 0;
    
    %While error has not converged
    while(norm(old_error-new_error)> eps && in_iter<5)
        
        %Get correspondences and weight them
        fprintf("Calculating correspondences between Target and Source mesh...");
        [W0,U] = getCorrespondences(Target,newSource,Options,L_target,L_source);
        %Expand weights so there is a weight for every coordinate of each
        %point
%         rigid_vertices = not(W0(:,2).*W0(:,3)); %vertices with no correspondences
        W = kron(W0,ones(3,1));
        W = [W; zeros(6,3)];
        fprintf(" done. \n");
        
        vars = [Source.vertices, Target.vertices(U,:), Tnormals(U,:), repmat(reshape(R',1,[]),nS,1),repmat(T,nS,1)]; 

        %Define the coefficient matrix and pose the linear system Ax = b
        % x = [z1x z1y z1z ... znx zny znz alpha beta gamma tx ty tz]
        [A,b] = subs_nrICP(vars,w1,w2(n),false,W0);

        %Add membrane term to matrix
        b = b-w3(n)*W(:,1).*bm;
        A = A+w3(n)*W(:,1).*Al;
%         b = W(:,2).*b-w3(n)*W(:,1).*bm;
%         A = W(:,2).*A+w3(n)*W(:,1).*Al;

        %Add Thikonov regularization of alpha, beta gamma and translation
        %parameters: norm(alpha,beta,gamma,tx,ty,tz)
        % The derivatives add a term 2*parameter to the coefficient
        % correspondent to the parameter
        A(end-5:end,end-5:end) = A(end-5:end,end-5:end) + lambda*2*eye(6);
        
        %Add local rigidity term
        [A_lr,b_lr] = local_rigidity_term(newSource);
        A(1:3*nS,1:3*nS) = A(1:3*nS,1:3*nS) + ((1-kron(W0(:,1),ones(3,1)))*w4(n)).*A_lr;
        b(1:3*nS) = b(1:3*nS) + (1-kron(W0(:,1),ones(3,1)))*w4(n).*b_lr;
        
        %Finally, compute the solution
        solution = A\b;

        %Update source vertices
        newSource.vertices = reshape(solution(1:end-6),3,[])';
        
        %Calculate Laplace-Beltrami vectors
        [At,Ct] = LP_exp(newSource);
        L_source = At*Ct;
        
        aT = solution(end-2:end)';
        alpha = solution(end-5); beta = solution(end-4); gamma = solution(end-3);
        %Rt = getRotationMatrix(alpha, beta, gamma);
        Rt = [1 -gamma beta; gamma 1 -alpha; -beta alpha 1];
        
%         %Update rotation and translational vector;
%         T = aT;
%         R = Rt;
        
        %Update error
        old_error = new_error;
        [new_error, e] = energyfun(newSource.vertices,Source.vertices,vars(:,4:6),vars(:,7:9),Rt,R,aT,T,L_source0,bm,w1,w2(n),w3(n),w4(n),newSource.faces);

        %Update rotation matrix and translational vector
        %We ensure that the Rotation matrix is valid by using formula of "Linear Least-Squares Optimization for
        %Point-to-Plane ICP Surface Registration"

%         R = R*Rt;
%         T = T + aT';
        
%         %Visualize rigid transformation of source
%         TempSource = Source;
% 
%         for l=1:nS
%         TempSource.vertices(l,:) = (R*TempSource.vertices(l,:)')';
%         end
%         TempSource.vertices = TempSource.vertices + repmat(T,nS,1);
%         disp_meshes(Target,TempSource);
        
        %Display error
        fprintf('Current error: %f . Difference with previous iteration: %f \n',new_error,norm(old_error-new_error));

        %Display mesh result
        disp_meshes(Target,Source,newSource);
        
        %Gather data
        Data.SourceT{iter} = newSource;
        Data.correspondences(:,iter) = U;
        Data.LB_vec{iter} = At*Ct*newSource.vertices;
        Data.error(:,iter) = e;
        Data.weights(:,iter) = [w1; w2(n); w3(n)];
        
        %Increase iteration counter
        in_iter = in_iter +1;
        iter = iter+1;
       
    end
    
    fprintf('Weight schedule number %d . Reconstruction error: %f \n',n,new_error);
    
end


end


function [error, e] = energyfun(Z,Z0,P,N,aR,R,aT,T,Al,bm,w1,w2,w3,w4,f)

nV = size(Z0,1);
Al = kron(Al,eye(3));
e = zeros(3,1);
e(1) = w1*sum(sum((Z-P).^2,2));
for i=1:nV
    %error= error + w1*norm(dot(N(i,:),Z(i,:)-P(i,:)))^2 + w2*norm(Z(i,:)'-aR*(R*X(i,:)'+T')+aT)^2;
    e(2) = e(2) + w2*norm(Z(i,:)'-aR*(R*Z0(i,:)'+T')+aT)^2;
end
%Error of membrane term
e(3) = w3*norm(Al*(reshape(Z',[],1)-reshape(Z0',[],1)))^2;

%Error of local rigidity
d = reshape((Z-Z0)',1,[])';
e(4) = w4*sum(kron(meshAdjacencyMatrix(f),eye(3))*d);

fprintf(strcat('Correspondence error: ', num2str(e(1)), '\n', ...
               'RigidT error: ', num2str(e(2)), '\n', ...
               'Curvature error: ', num2str(e(3)), '\n'));
error = sum(e)/nV;
end

% function [A_lb,b] = LaplaceBeltramiTerm(Source)
% %     %Compute Laplace-Beltrami operator of source
% %     L = LP_exp(Source);
% %     L_d = diag(L);
% %      
% %     %Compute adjacency matrix
% %     A = triangulation2adjacency(Source.vertices,Source.faces);
% %     n = size(A,1);
% %     A = A - eye(n);
% %     b = sum(A,2);
% %     D = -diag(b);
% %     A = D + A;
%       
%       %Get Laplace Beltrami of vertices in Source
%       [A,C] = LP_exp(Source);
%       L = A*C;
%       Ls = L.*L;
%       nV = size(Source.vertices,1);
%       V = ones(nV,3);
%       
%       %Operator is applied on deformation vector di = vi'-vi
%       %The formula to minimize is mod(L*d)^2.
%       %Therefore, we have to minimize mod(L*(vi'-vi))^2
%       %Supposing matrix L = [a_11 a_12 ... a_1n; a_n1 ... a_nn]
%       %every term c_ij of the mod(L*d) matrix can be expressed as:
%       % c_ij = mod(a_ij*(v_i'-v_i))^2 = a_ij^2*(v_i'-v_i)^2
%       % c_ij = a_ij^2*(v_i'^2-2*v_i'*v_i+v_i^2)
%       
%       % And its derivative with respect to v_i' is
%       % d(c_ij)/d(v_i) = a_ij^2*(2*v_i') - a_ij^2*(2*v_i)
%       
%       % Coefficients dependent on v_i' = a_ij^2*(2*v_i')
%       A_lb = 2*Ls*V;
%       % Independent coefficients = - a_ij^2*(2*v_i)
%       % ATTENTION TO THE SIGN
%       b = -2*Ls*Source.vertices;
%       
% end

%Returns the coefficient matrix of the derivative of the Laplace Beltrami term
% function [A_lb,b] = LaplaceBeltrami2(Source)
%     
%     %Obtain the Laplace-Beltrami expression for the mesh
%     [A,C] = LP_exp(Source);
%     L = A*C;
%     L_sum = sum(L,1)';
%     
%     %Compute its derivative
%     nV = size(Source.vertices, 1);
%     A_lb = repmat(L_sum,1,nV)';
%     A_lb = kron(A_lb,k);
%     %And the independent terms
%     b = L_sum.*Source.vertices;
%     b = reshape(b',[],1);
% end

% function [A_lb,b] = LaplaceBeltrami2(Source)
%     
%     %Obtain the Laplace-Beltrami expression for the mesh
%     [A,C] = LP_exp(Source);
%     L = A*C;
%     L_sum = sum(L,1)';
%     
% end


%         %Calc-Laplace Beltrami for assumption testing
%         [Al,bm,A0,C0] = LaplaceBeltramiCoefficients(newSource);
%         Al = [Al zeros(nS*3,6)];
%         Al = [Al ; sparse(6,size(Al,2))];
%         bm = [bm; zeros(6,1)];
% %         Al = kron(Al,eye(3));
%         if(sum(sum(imag(A0)))>0)
%             disp('Found imaginary component on C');
%         end
%         if(sum(sum(imag(C0)))>0)
%             disp('Found imaginary component on A');
%         end
%         LB_vec = [LB_vec A0*C0*newSource.vertices];
        
% %         [lc,~] = LB_coeff(Source);
% %         LB = [LB diag(lc)];
% 
% %         LB = [LB, B];
% %         [M,C] = LP_exp(newSource);
% %         L = M*C;
% %         L = kron(L,eye(3));
% %         LBi = L*reshape(newSource.vertices',1,[])';
% %         LBi = reshape(LBi,3,[])';
% %         LBi = sqrt(sum(LBi.^2,2));
% %         LB = [LB, LBi];