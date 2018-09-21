%Demo for alligning two meshes with arbitrary directions

load Source
load Target;
N = size(Target.vertices,1);

alpha = rand(1);
beta = rand(1);
gamma = rand(1);
[Rx, Ry, Rz, ~] = GetRotationMatrix([alpha beta gamma 0 0 0]);
R = Rx*Ry*Rz;
R = R(1:3,1:3);
t = rand(3,1)*100;
for i=1:N
    Target.vertices(i,:)= (R*Target.vertices(i,:)' + t)';
end

nricp_custom(Target,Source);
    