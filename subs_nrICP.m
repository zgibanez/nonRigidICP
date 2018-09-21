function [A_total,b_total] = subs_nrICP(vars, w1, w2, linearization, W)
%SUBS_RICP Rearrange of partial derivative variables into coefficients of matrices A and b
%such that the linear system Ax = b can be solved

nVertex = size(vars,1);
%Preallocate memory for solution matrix
A_total = sparse(nVertex*3+6,nVertex*3+6);
b_total = zeros(nVertex*3+6,1);
count = 0;
A_ang_t = zeros(6,nVertex*3+6);

n=0;

%%Add non-rigid term derivation term
A_total(1:nVertex*3,1:nVertex*3) = A_total(1:nVertex*3,1:nVertex*3) + w1*2*speye(nVertex*3);
P = reshape(vars(1:end,4:6)',1,[])';
b_total(1:nVertex*3) = b_total(1:nVertex*3) + 2*w1*P;

for i=1:nVertex
    x1 = vars(i,1); x2 = vars(i,2); x3 = vars(i,3);
    px = vars(i,4); py = vars(i,5); pz = vars(i,6);
    nx = vars(i,7); ny = vars(i,8); nz = vars(i,9);
    r11 = vars(i,10); r12 = vars(i,11); r13 = vars(i,12);
    r21 = vars(i,13); r22 = vars(i,14); r23 = vars(i,15);
    r31 = vars(i,16); r32 = vars(i,17); r33 = vars(i,18);
    tx = vars(i,19); ty = vars(i,20); tz = vars(i,21);
    % subs(Am,symvars,vars(i,:));
    % A(((i-1)*9+1):9*i,:) = eval(Am);
    % b(((i-1)*9+1):9*i,:) = eval(bm);
    
    if(linearization)
%         A_temp = Am(x1,x2,x3,px,py,pz,nx,ny,nz,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,W(:,1).*w1,w2);
%         b_temp = bm(x1,x2,x3,px,py,pz,nx,ny,nz,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,W(:,1).*w1,w2);
        A_temp = Am(x1,x2,x3,px,py,pz,nx,ny,nz,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,w1,w2);
        b_temp = bm(x1,x2,x3,px,py,pz,nx,ny,nz,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,w1,w2);
    else
        A_temp = Am2(x1,x2,x3,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,W(i,1)*w1,w2);
        b_temp = bm2(x1,x2,x3,px,py,pz,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,W(i,1)*w1,w2);
%         A_temp = Am2(x1,x2,x3,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,w1,w2);
%         b_temp = bm2(x1,x2,x3,px,py,pz,r11,r12,r13,r21,r22,r23,r31,r32,r33,tx,ty,tz,w1,w2);
    end
    
    %Rearrange A coefficients
    A_total(((i-1)*3+1):3*i,((i-1)*3+1):3*i) = A_temp(1:3,1:3);
    A_total(((i-1)*3+1):3*i,(end-5):end) = A_temp(1:3,4:end);
    
    %Rearrange b coefficients
    b_total(((i-1)*3+1):3*i,:) = b_temp(1:3);
    b_total((end-5):end,:) = b_total((end-5):end,:) + b_temp(4:end);
    
    %Isolate terms of A that depend on alpha beta gamma tx ty and tz
    A_ang_t(1:end,((i-1)*3+1):3*i) = A_temp(4:9,1:3);
    A_ang_t(1:end,(end-5):end) = A_ang_t(1:end,(end-5):end)+A_temp(4:9,4:9);
    
    if(rem(count,1000)==0)
        fprintf(repmat('\b',1,n));
        msg = strcat('Calculating coefficent matrix... ', pbar(count,nVertex));
        fprintf(msg);
        n = numel(msg);
    end
    count = count+1;
end

A_total((end-5):end,:) = A_ang_t;
fprintf(repmat('\b',1,n));
fprintf('\n');
end
