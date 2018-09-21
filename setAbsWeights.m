function w_abs = setAbsWeights(distances,radius)
%Returns weights betwen 0 and 1 in function of abs(x)^p

d_min = min(distances);
d_norm = (distances-d_min)/radius;

p=0.5;

fL = @(x) (1-abs(x).^p);

w_abs = fL(d_norm);

end
