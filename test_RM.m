%Test for rotation matrix

R = rand(3);
V = rand(3,10);

c=sum((V*R.').*V,2);

nV = size(V,1);
for i=1:nV
    c2(i,:) = (R*V(i,:)')';
end

er = norm(c-c2);