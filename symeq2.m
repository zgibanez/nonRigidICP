clear all;

syms zx zy zz
syms ax ay az
syms txc tyc tzc
syms alpha beta gamma

assume(alpha,'real');
assume(beta,'real');
assume(gamma,'real');
assume(zx,'real');
assume(zy,'real');
assume(zz,'real');
assume(txc,'real');
assume(tyc,'real');
assume(tzc,'real');
assume(ax,'real');
assume(ay,'real');
assume(az,'real');

syms zx2 zy2 zz2
assume(zx2,'real');
assume(zy2,'real');
assume(zz2,'real');

Tc = [txc tyc tzc]';
Z = [zx zy zz]';
Z2 = [zx2 zy2 zz2]';
% A = (Rx+t)
A = [ax ay az]';
Rc = [1 -alpha beta; gamma 1 -alpha; -beta alpha 1];

eq1 = simplify(norm(Z-Rc*A+Tc))^2;
eq2 = simplify(norm(Z2-Rc*A+Tc))^2;
eqT = eq1 + eq2;
deq1 = gradient(eq1,[zx zy zz alpha beta gamma txc tyc tzc]);
deq1 = deq1 == 0;
[Am,bm] = equationsToMatrix(deq1,[zx zy zz alpha beta gamma txc tyc tzc]);