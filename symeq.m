syms alpha beta gamma r11 r12 r13 r21 r22 r23 r31 r32 r33 x1 x2 x3 tx ty tz
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
assume(x1,'real');
assume(x2,'real');
assume(x3,'real');
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
assume(zx,'real');
assume(zy,'real');
assume(zz,'real');
assume(nx,'real');
assume(ny,'real');
assume(nz,'real');
assume(px,'real');
assume(py,'real');
assume(pz,'real');
assume(w1,'real');
assume(w2,'real');


p = [px py pz];
n = [nx ny nz];
z = [zx zy zz];

Rc = [1 -gamma beta; gamma 1 -alpha; -beta alpha 1];
R = [r11 r12 r13; r21 r22 r23; r31 r32 r33];
X = [x1 x2 x3]';
T = [tx ty tz]';
Tc = [txc tyc tzc]';


eq1 = simplify(w1*dot(n,(z-p))^2);
% eq1 = simplify(norm(z-p)^2);
diff_eq1 = gradient(eq1,[zx zy zz]);
diff_eq1 = diff_eq1 == 0;
S1 = solve(diff_eq1,[zx zy zz]);
[A,b] = equationsToMatrix(diff_eq1,[zx zy zz]);

b = Rc*(R*X+T);
eq2 = simplify(w2*norm(z'-b+Tc)^2);
% diff_eq2 = gradient(eq2,[zx zy zz alpha beta gamma txc tyc tzc]);
% diff_eq2 = diff_eq2 == 0;
% S2 = solve(diff_eq2,[alpha beta gamma txc tyc tzc]);
% equationsToMatrix(diff_eq1,
vars = [zx zy zz alpha beta gamma txc tyc tzc];
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

