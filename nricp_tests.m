function nricp_tests(Target,Source, w, preal, LTarget,LSource)

if(preal)
    [s, R, t, err, newP] = icp_custom(Target,Source);
    Source.vertices = newP;
end

if(nargin>4)
    [Source.vertices, ~] = adaptTemplate(Target,Source,LTarget,LSource);
end

%Nricp options
Options.normalTrim = true;
Options.boundaryTrim = true;
Options.distanceTrim = true;
Options.LBmatching = false; %Better if true
Options.NormalMatching = true; %Better if false

%Calculate vertex normals of Source and target
fprintf("Calculating Target mesh normals...");
Tnormals = calc_normals(Target);
fprintf(" done. \n");

%Store vertex
nS = size(Source.vertices,1);

%Initialize translational and rotational component
R = eye(3);
T = [0 0 0];

%Weights
steps = size(w,1);
w1 = w(:,1);
w2 = w(:,2);
w3 = w(:,3);
w4 = w(:,4); %local rigidity weight

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
% fprintf("Calculating Laplace-Beltrami vectors of Target mesh...");
% [AT,CT] = LP_exp(Target);
% fprintf(" done. \n");
L_target = 0;
L_source = A0*C0;
L_source0 = L_source;

%Parameter for Thikonov regularization
lambda = 100;

%Iteration counter
iter = 1;

%For every scheduled weight
for n=1:steps
    %Placeholder for transformed source
    newSource = Source;
    
    old_error = Inf;
    new_error = 10e6;
    in_iter = 0;
        
        %Get correspondences and weight them
        fprintf("Calculating correspondences between Target and Source mesh...");
        [W0,U] = getCorrespondences(Target,newSource,Options,L_target,L_source);
        %Expand weights so there is a weight for every coordinate of each
        %point
        W = kron(W0,ones(3,1));
        W = [W; zeros(6,3)];
        fprintf(" done. \n");
        
        vars = [Source.vertices, Target.vertices(U,:), Tnormals(U,:), repmat(reshape(R',1,[]),nS,1),repmat(T,nS,1)]; 

        %Define the coefficient matrix and pose the linear system Ax = b
        [A,b] = subs_nrICP(vars,w1(n),w2(n),false,W0);

        %Add membrane term to matrix
        b = b-w3(n)*W(:,1).*bm;
        A = A+w3(n)*W(:,1).*Al;

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
        Rt = [1 -gamma beta; gamma 1 -alpha; -beta alpha 1];
        
        %Update error
        old_error = new_error;
        [new_error, e] = energyfun(newSource.vertices,Source.vertices,vars(:,4:6),vars(:,7:9),Rt,R,aT,T,L_source0,bm,w1(n),w2(n),w3(n),w4(n),newSource.faces);

        %Display error
        fprintf('Current error: %f . Difference with previous iteration: %f \n',new_error,norm(old_error-new_error));

        %Display mesh result
        disp_meshes(Target,Source,newSource);
        
        %Increase iteration counter
        in_iter = in_iter +1;
        iter = iter+1;
        
        write_obj(strcat('test_w1_',num2str(w(n,1)),'_w2_',num2str(w(n,2)),'_w3_',num2str(w(n,3)),'_w4_',num2str(w(n,4)),'.obj'),newSource.vertices,newSource.faces);
    
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
d = abs(reshape((Z-Z0)',1,[])');
e(4) = w4*sum(kron(meshAdjacencyMatrix(f),eye(3))*d);

fprintf(strcat('Correspondence error: ', num2str(e(1)), '\n', ...
               'RigidT error: ', num2str(e(2)), '\n', ...
               'Curvature error: ', num2str(e(3)), '\n', ...
               'LocalRigid error: ', num2str(e(4)), '\n'));
error = sum(e)/nV;
end