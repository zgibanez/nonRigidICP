function [SI,binSize] = computeSpinImage(Mesh, points, Options)
%%Based on the thesis of https://www.ri.cmu.edu/pub_files/pub2/johnson_andrew_1997_3/johnson_andrew_1997_3.pdf
%Input:
% Mesh - Mesh to compute the spin images. Must have a nx3 'vertices' field
% and nx3 'faces' field
% W - Width of the spin image (WxW).
% f - factor to multiply to the Mesh resolution to compute the bin size
% points - vector of indexes of points from which the SI will be computed.

%Output: 
% SI - a cell containing all the spin images of the points specified
% binSize - output binSize to be reused to compare two images

if(nargin<2)
    %A sample point
    error('At least one point should be specified');
end

if(nargin<3)
    %A sample point
    error('Options field needed');
end

if(isfield(Options,'binSize'))
    binSize = Options.binSize;
else
    
    if (~isfield(Options,'f'))
        f = 0.5;
    else
        f = Options.f;
    end
    
    binSize = f*medianEdgeLength(Mesh);
end

if (~isfield(Options,'W'))
    W = 100;
else
    W = Options.W;
end

if (~isfield(Options,'Aw'))
    Aw = pi/2;
else
    Aw = Options.Aw;
end

%The bin size is proportional to the median edge length

% binSize = 2;

%The function S0 Computes the alpha-beta map
%S0 = @(p,n,x) [sqrt(norm(x-p)^2-dot(n,x-p)^2),dot(n,x-p)];
S0 = @(p,n,x)[sqrt(sum((x-p).^2,1)-dot(n,x-p).^2);dot(n,x-p)];

v = Mesh.vertices';
nV = size(v,2);

if(isfield(Mesh,'normals'))
    n = Mesh.normals;
    %Transpose the matrix in case it has more columns than rows
    if(size(n,1) > size(n,2))
        n = n';
    end
else
    n = calc_normals(Mesh)';
end

nP = numel(points);
for i=1:nP

    P = repmat(v(:,points(i)),1,nV);
    N = repmat(n(:,points(i)),1,nV);

    %Compute the alpha-beta map
    ab = S0(P,N,v);

    %Identify points outside the supporting angle
    idx_angle = atan2(vecnorm(cross(N,n,1)),dot(N,n,1)) < Aw;
    
    %Identify points outside the support distance (out of the image width)
    idx_distance = not(ab(1,:)>W | ab(2,:) < -W/2 | ab(2,:) > W/2);
    
    %The collection of points that suffice the conditions become the valid candidates
    idx = idx_angle & idx_distance;
    ab_trimmed = ab(:,idx);

    %Transform ab coordinates in ij (image) coordinates
    ij  = floor([(W/2-ab_trimmed(2,:))./binSize;ab_trimmed(1,:)./binSize])+1;

    SI{i} = computeBilinearWeights(ab_trimmed,ij,binSize,W);
end

display = false;
if(display)
    figure();
    subplot(1,3,1);
    scatter(ab(1,:),ab(2,:));
    rectangle('Position',[0 -W/2 W W]);
    title('\alpha \beta map wo. trimming');
    axis equal;
    subplot(1,3,2);
    scatter(ab_trimmed(1,:),ab_trimmed(2,:));
    title('\alpha \beta map w. trimming');
    axis equal;
    subplot(1,3,3);
    image(SI{1});
    title('Spin image');
end

end

function SI = computeBilinearWeights(ab,ij,binSize,W)
    %Corrections based upon mistakes on the formula (2.3) thesis
    wa = W/2 - ab(2,:)-ij(1,:)*binSize;
    wb = ab(1,:)-ij(2,:)*binSize;
    
    nP = size(wa,2);
    nBins = ceil(W/binSize);
    SI = zeros(nBins,nBins);    
    for n=1:nP
        SI(ij(1,n),ij(2,n)) = SI(ij(1,n),ij(2,n)) + (1-wa(n)).*(1-wb(n));
        if(ij(1,n)<nBins&&ij(2,n)<nBins)
            SI(ij(1,n)+1,ij(2,n)) = SI(ij(1,n)+1,ij(2,n)) + (wa(n)).*(1-wb(n));
            SI(ij(1,n),ij(2,n)+1) = SI(ij(1,n),ij(2,n)+1) + (1-wa(n)).*wb(n);
            SI(ij(1,n)+1,ij(2,n)+1) = SI(ij(1,n)+1,ij(2,n)+1) + wa(n).*wb(n);
        end
    end
    
end

%Computes the median of the edge length as a measure of resolution
function l = medianEdgeLength(Mesh)
f = Mesh.faces;
v = Mesh.vertices;
nsq = @(A) sum(sqrt(A.^2),2);
lengths = [nsq(v(f(:,1),:)-v(f(:,2),:)); nsq(v(f(:,2),:)-v(f(:,3),:)); nsq(v(f(:,3),:)-v(f(:,1),:))];
l = median(lengths);
end