function LB_plot(Mesh,LB)

LB = LB(1:3:end,:);
%Histograms of succesive iterations
divisions = 100;
[n,m] = size(LB);
d1 = [LB zeros(n,1)];
d2 = [zeros(n,1) LB];
D = abs(d2-d1);
% maxD = max(prctile(D,95)); minD = min(min(D));
maxD = 500; minD = 0;
step = (maxD-minD)/divisions;
ctr = minD:step:maxD;

for i=1:m-1
    subplot(1,m-1,i)
    histogram(D(:,i+1)./D(:,1),ctr);
    title(strcat('Iter',num2str(i)));
end

% 
% if(nargin == 1)
%    [M,C] = LP_exp(Mesh);
%    L = M*C;
%    L = kron(L,eye(3));
%    LB = L*reshape(Mesh.vertices',1,[])';
%    LBi = reshape(LB,3,[])';
%    LBi = sqrt(sum(LBi.^2,2));
%    idx = LBi < prctile(LBi,95);
%    LBi = LBi.*idx;
%    maxV = max(LBi);
%    colorVector = 10*LBi / maxV;
%    %Higher curvature values are treated differently because they distort
%    %the mean
%    colorVector = colorVector + not(idx)*10;
%    scatter3(Mesh.vertices(:,1),Mesh.vertices(:,2),Mesh.vertices(:,3),[],colorVector);
%     axis equal;
% else
%     LB = LB(1:3:end,:);
%     iter = size(LB,2);
% %     [row,~] = size(LB);
% %     stdColor = repmat(0.5,row,1);
% end


%Plot original LB
figure();
subplot(1,4,1);
Dt = D(:,1) < prctile(D(:,1),95);
colorVector = 10*(D(:,1).*Dt)/prctile(D(:,1),95) + 10*not(Dt);
scatter3(Mesh.vertices(:,1),Mesh.vertices(:,2),Mesh.vertices(:,3),[],colorVector);
axis equal;

for i=1:3
    subplot(1,4,i+1);
    idx = D(:,i+1) < prctile(D(:,i+1),95);
    D(:,i+1)  = D(:,i+1) .*idx;
    maxV = max(D(:,i+1) );
%     minV = min(LBi);
    colorVector = 10*D(:,i+1) ./maxV + 10*not(idx);
%     [~,idx] = sort(LBi,'ascend');
%     colorVector = 10*idx/row;
    
%     subplot(1,iter,i);
    scatter3(Mesh.vertices(:,1),Mesh.vertices(:,2),Mesh.vertices(:,3),[],colorVector);
    axis equal;
end


end
