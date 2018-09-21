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
z = [zx zy zz]';

Rc = [1 -gamma beta; gamma 1 -alpha; -beta alpha 1];
R = [r11 r12 r13; r21 r22 r23; r31 r32 r33];
X = [x1 x2 x3]';
T = [tx ty tz]';
Tc = [txc tyc tzc]';

% eq1 = norm(z'+Rc(eye(3)*X+zeros(3,1))+T)^2;
eq1 = sum((z-Rc*(R*X+T)+Tc).^2,1);

vars = [zx zy zz alpha beta gamma txc tyc tzc];
deq = gradient(eq1,vars);

deq = deq==0;
[Am,bm] = equationsToMatrix(deq,vars);
matlabFunction(Am,'file','Am_rigid');
matlabFunction(bm,'file','bm_rigid');


