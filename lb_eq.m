syms v1h v1 v2h v2 v3h v3 v4h v4
assume(v1h,'real');
assume(v2h,'real');
assume(v3h,'real');
assume(v1,'real');
assume(v2,'real');
assume(v3,'real');

syms a11 a12 a13
assume(a11,'real');
assume(a12,'real');
assume(a13,'real');

syms a21 a22 a23
assume(a21,'real');
assume(a22,'real');
assume(a23,'real');

e1 = a11*(v1h-v1)+a12*(v2h-v2)+a13*(v3h-v3)+a21*(v1h-v1) + a22*(v2h-v2) + a23*(v3h-v3);
e1_e1 = simplify(diff(norm(e1)^2,v1h));
e1_e2 = simplify(diff(norm(e1)^2,v2h));

[A,b] = equationsToMatrix(e1_e1,[v1h v2h v3h])