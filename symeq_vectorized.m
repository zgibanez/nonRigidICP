clear all
N = 13000; %vertex of template
syms alpha beta gamma r11 r12 r13 r21 r22 r23 r31 r32 r33 tx ty tz
X1 = sym('x1_',[N,1],'real');
X2 = sym('x2_',[N,1],'real');
X3 = sym('x3_',[N,1],'real');
syms txc tyc tzc
assume(alpha,'real');
assume(beta,'real');
assume(gamma,'real');
assume(r11,'real');
assume(r12,'real');
assume(r13,'real');
assume(r21,'real');
assume(r22,'real');
assume(r23,'real');
assume(r31,'real');
assume(r32,'real');
assume(r33,'real');
% assume(X1,'real');
% assume(X2,'real');
% assume(X3,'real');
assume(tx,'real');
assume(ty,'real');
assume(tz,'real');
assume(txc,'real');
assume(tyc,'real');
assume(tzc,'real');


syms zx zy zz
syms nx ny nz
syms px py pz
syms w1 w2
Zx = sym('zx_',[N,1]);
Zy = sym('zy_',[N,1]);
Zz = sym('zz_',[N,1]);
Nx = sym('nx_',[N,1]);
Ny = sym('ny_',[N,1]);
Nz = sym('nz_',[N,1]);
Px = sym('px_',[N,1]);
Py = sym('py_',[N,1]);
Pz = sym('pz_',[N,1]);

assume(Zx,'real');
assume(Zy,'real');
assume(Zz,'real');
assume(Nx,'real');
assume(Ny,'real');
assume(Nz,'real');
assume(Px,'real');
assume(Py,'real');
assume(Pz,'real');
assume(w1,'real');
assume(w2,'real');


p = [Px Py Pz];
n = [Nx Ny Nz];
z = [Zx Zy Zz];

Rc = [1 -alpha beta; gamma 1 -alpha; -beta alpha 1];
R = [r11 r12 r13; r21 r22 r23; r31 r32 r33];
X = [X1 X2 X3]';
T = [tx ty tz]';
Tc = [txc tyc tzc]';

for i=1:N
    eq1 = simplify(w1*dot(n(i,:),(z(i,:)-p(i,:))^2));

    b = Rc*(R*X(:,i)+T);
    eq2 = simplify(w2*norm(z(i,:)'-b+Tc)^2);
end
vars = [Zx Zy Zz alpha beta gamma txc tyc tzc];
deq = gradient(eq1+eq2,vars);
deq = deq==0;
[Am,bm] = equationsToMatrix(deq,vars);

% syms zx2 zy2 zz2 px2 py2 pz2 xx2 xy2 xz2 nx2 ny2 nz2
% assume(zx2,'real');
% assume(zy2,'real');
% assume(zz2,'real');
% assume(nx2,'real');
% assume(ny2,'real');
% assume(nz2,'real');
% assume(px2,'real');
% assume(py2,'real');
% assume(pz2,'real');
% assume(xx2,'real');
% assume(xy2,'real');
% assume(xz2,'real');
% 
% Z2 = [zx2 zy2 zz2];
% N2 = [nx2 ny2 nz2];
% P2 = [px2 py2 pz2];
% X2 = [xx2 xy2 xz2]';
% 
% b = Rc*(R*X2+T);
% eq12 = simplify(w1*dot(N2,(Z2-P2))^2);
% eq22 = simplify(w2*norm(Z2'-b+Tc)^2);
% 
% eqt = eq1 + eq2 + eq12 + eq22;
% vars = [vars zx2 zy2 zz2];
% deqt = gradient(eqt,vars);
% [Am,bm] = equationsToMatrix(deqt,vars);