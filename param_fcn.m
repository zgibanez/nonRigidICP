function [outputArg1,outputArg2] = param_fcn(inputArg1,inputArg2)
%PARAM_FCN Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = inputArg1;
outputArg2 = inputArg2;

voronoiAreaT = @(P,Q,R) (1/8)*(norm(P-R)^2)
end

function area = VoronoiAreaT(P,Q,R)
    cotg = @(u,v) dot(u,v)/cross(u,v);

    if isTObtuse(P,Q,R)
        
    else
        area = 1/8 * (norm(P-R)^2*cotg(P-Q,R-Q)+norm(P-Q)^2*cotg(P-R,Q-R));
    end
end

%Returns True if the triangle is obtuse
function is_obtuse = isTObtuse(P,Q,R)

l(1) = norm(P-Q);
l(2) = norm(P-R);
l(3) = norm(R-Q);
l = sort(l);

if l(3)^2 > (l(2)+l(1)^2)
    is_obtuse = true;
else
    is_obtuse = false;
end

end