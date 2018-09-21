function w= trimCorrespondences(Source,Target,Cor)

%Weigths wector, 1 if correct and 0 if incorrect
w = zeros(size(Source.vertices,1),1);

%Vector of correspondent vertices
corV = Target(Cor);

%Distance threshold
dist = norm(Source.vertices - corV);
s = std(dist);
m = mean(dist);
w_dist = (dist > (s+m)) || (dist < (s-d));

%Normal threshold
Snormals = calc_normals(Source);
Tnormals = calc_normals(Target);

end