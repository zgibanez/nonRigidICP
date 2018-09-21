function [W,U] = getCorrespondences(Target,Source,Options,Lt,Ls,normalT)
%%Description:



%%Input
%Target: Target mesh
%Source: Source mesh
%Options: Struct with Options for trimming
%       distanceTrim: if true, apply distance trimming
%       lambda:       distance trimming std threshold
%       normalTrim:   if true, apply incompatible normal trimming
%       angT:         angle threshold in normal trimming
%       boundaryTrim: if true, ignore boundary vertices
%Lt,Ls: Laplace-Beltrami vectors of target and source, respectively
%normalsT: (Optional) Normals of target mesh. Helps speeding up
%performance.

%%Output
%W: n-by-3 matrix of correspondence weights, each column with a
%different criterion
%U: Vector of correspondence indexes

%If no options are defined, create one set
if(nargin<3)
    Options.normalTrim = true;
end

%Set default options if they do not exist
%Apply normal threshold?
if(~isfield(Options,'normalTrim'))
    Options.normalTrim = true;
end
%Ignore boundary vertices?
if(~isfield(Options,'boundaryTrim'))
    Options.boundaryTrim = false;
end
%Apply distance threshold?
if(~isfield(Options,'distanceTrim'))
    Options.distanceTrim = true;
end
%Apply LB matching?
if(~isfield(Options,'LBmatching') || nargin < 4)
    Options.LBmatching = false;
end
%Apply normal matching?
if(~isfield(Options,'NormalMatching') || nargin < 4)
    Options.NormalMatching = true;
end
%Apply distance reweighting?
if(~isfield(Options,'distanceWeight'))
    Options.distanceWeight = true;
end

%Distance threshold
if(~isfield(Options,'lambda'))
    Options.lambda = 5;
end
%Angle threshold
if(~isfield(Options,'angT'))
    Options.angT = 1;
end

%Initialize weight vector as 1 (all correspondenes are valid)
W = ones(size(Source.vertices,1),3);
U = knnsearch(Target.vertices,Source.vertices);

%Precalculate some variables that are common to two or more options
if(Options.LBmatching == true || Options.distanceTrim == true || Options.NormalMatching == true)
    distancesKNN = sqrt(sum((Source.vertices-Target.vertices(U,:)).^2,2));
    
    %Fourth spread of distances
    fs = median(maxk(distancesKNN,round(numel(distancesKNN)/2)));

%     h = histogram(distancesKNN,'NumBins',20);
%     hold on;
% %     line([z z],get(gca,'Ylim'),'Color',[1 0 0]);
% %     line([1.5*z 1.5*z],get(gca,'Ylim'),'Color',[1 0 0],'LineWidth',2);
%     line([1.5*fs 1.5*fs],get(gca,'Ylim'),'Color',[1 0 0],'LineWidth',2);
%     line([3*fs 3*fs],get(gca,'Ylim'),'Color',[1 0 0],'LineWidth',2);

end

if(Options.NormalMatching == true || Options.normalTrim == true)
        normalS = calc_normals(Source);
end


%Apply Normal matching if specified
%Normal matching matches each source point with the closest one to its
%vertex normal
if(Options.NormalMatching == true)
    point2lineDist = @(x0,x1,x2) sqrt(sum(cross(x0-x1,x0-x2,2).^2,2))./sqrt(sum((x2-x1).^2,2));
        
    %Function for checking if point is inside a sphere
    %From: https://stackoverflow.com/questions/26818772/python-checking-if-a-point-is-in-sphere-with-center-x-y-z
    inside_sphere = @(X,C,r) ( X(:,1)-C(1) ).^2 + (X(:,2)-C(2)).^2 + (X(:,3)-C(3)).^ 2 < r^2;
    
    %Radius of the sphere is twice the mean distance
%     radius = 5*mean(distancesKNN);
    radius = 4*fs; %Upper outliers
    %radius = 0.7*median(distancesKNN);
    
    nV = size(Source.vertices,1);
    d_normal = zeros(nV,1); %Set of normal distances
    
    for i=1:nV
        C = Source.vertices(i,:);
        idx = find(inside_sphere(Target.vertices,C,radius));
        if(~isempty(idx))
            %Calculate the distance from the candidates to the normal
            %directions
            x0 = Target.vertices(idx,:);
            x1 = repmat(Source.vertices(i,:),numel(idx),1);
            x2 = x1 + 2*normalS(i,:); %Starting from x1 
            distances = point2lineDist(x0,x1,x2);
            if(sum(distances<0)>=1)
                warning('Detected negative distances at calculating distances point2line');
            end
            [d_normal(i),u] = min(distances);
            U(i) = idx(u);
        else
            d_normal(i) = radius;
        end
    end
    
    if(Options.distanceWeight == true)
        W(:,1) = setAbsWeights(d_normal,radius);
    end
    
end

%Apply LB matching if is specified
%LB matching identifies the target points within a sphere centered in each
%Source vertex. The target vertices are candidates to be matched with the
%Target one, and the one with the most similar LB value is selected.
if(Options.LBmatching ==true)
    
    %Function for checkin if point is inside a sphere
    inside_sphere = @(X,C,r) ( X(:,1)-C(1) ).^2 + (X(:,2)-C(2)).^2 + (X(:,3)-C(3)).^ 2 < r^2;
    
    %Radius of the sphere is twice the mean distance
    %From https://stackoverflow.com/questions/26818772/python-checking-if-a-point-is-in-sphere-with-center-x-y-z
    radius = 3*fs;
%     radius = median(distancesKNN);
    Ltv = Lt*Target.vertices;
    Lsv = Ls*Source.vertices;
    
    nV = size(Source.vertices,1);
    for i=1:nV
        C = Source.vertices(i,:);
        idx = find(inside_sphere(Target.vertices,C,radius));
        if(~isempty(idx))
            R = repmat(Lsv(i,:),numel(idx),1);
            [~,u] = min(sqrt(sum((Ltv(idx,:)-R).^2,2)));
            U(i) = idx(u);
        end
    end
    
end

%%Check points on boundaries
if (Options.boundaryTrim == true)
    b = compute_boundary_custom(Source.faces);
    W(b,1) = 0;
end

%%Check incompatible distances
if (Options.distanceTrim == true)
    distances = sqrt(sum((Source.vertices-Target.vertices(U,:)).^2,2));
    s = std(distances);  
    m = mean(distances);
    Wd = distances < m+Options.lambda*s;
    W(:,2) = Wd.*W(:,2);
end

%%Check incompatible normals
if (Options.normalTrim == true)
    if(nargin <6)
        normalT = calc_normals(Target);
    end
    normalT = normalT(U,:);
    
    %Calculate angles between 2 vectors
    %From: https://nl.mathworks.com/matlabcentral/answers/101590-how-can-i-determine-the-angle-between-two-vectors-in-matlab
    ang = atan2(vecnorm(cross(normalS,normalT,2)')',dot(normalS,normalT,2));
    s = std(ang);  
    m = mean(ang);
    Wa = not(ang > m+s*Options.angT);
    W(:,3) = Wa.*W(:,3);
end

%Display results
display = false;
if(display)
    display_results(Source,W);
end

end

function display_results(Source,W)
    v = Source.vertices;
    v1 = Source.vertices(not(W(:,1)),:);
    v2 = Source.vertices(not(W(:,2)),:);
    v3 = Source.vertices(not(W(:,3)),:);

    figure();
    subplot(1,4,1);
    scatter3(v(:,1),v(:,2),v(:,3),'b');
    hold on;
    scatter3(v1(:,1),v1(:,2),v1(:,3),'r');
    axis equal;
    title('Boundary Trim');

    subplot(1,4,2);
    scatter3(v(:,1),v(:,2),v(:,3),'b');
    hold on;
    scatter3(v2(:,1),v2(:,2),v2(:,3),'r');
    axis equal;
    title('Distance Trim');
    
    subplot(1,4,3);
    scatter3(v(:,1),v(:,2),v(:,3),'b');
    hold on;
    scatter3(v3(:,1),v3(:,2),v3(:,3),'r');
    axis equal;
    title('Normal Trim');
    
    subplot(1,4,4);
    scatter3(v(:,1),v(:,2),v(:,3),1,(1-W(:,1)));
    axis equal;
    title('Distance Reweight');

end

