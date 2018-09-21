function [U, SI_source, SI_target] = getSpinImageCorrespondences(Source,Target,pointSource,pointTarget)

%if no points on Target are selected compute with all
if(nargin<4)
    pointTarget = 1:1:size(Target.vertices,1);
end

%Compute spin images of landmarks on Source
Options.W = 50;
[SourceSI, binSize] = computeSpinImage(Source,pointSource,Options);

%Compute spin images of candidates on Target
Options.binSize = binSize;
Target.normals = calc_normals(Target)';
Source.normals = calc_normals(Source)';
[TargetSI,~] = computeSpinImage(Target,pointTarget,Options);

%Set the parameter lambda (Section 3.1 of Andrew Edie Johnson thesis)
nS = numel(SourceSI);
nT = numel(TargetSI);
nO = zeros(nT,1); %Counter of non-zero bins
for i = 1: numel(TargetSI)
    nO(i) = numel(find(TargetSI{i}));
end
% for i = 1: numel(SourceSI)
%     nO(i) = numel(find(SourceSI{i}));
% end
lambda = median(nO)/2;

C = zeros(numel(TargetSI),numel(SourceSI));
for i=1:nS
    for j=1:nT
        j
        C(j,i)= spinImageCorr(SourceSI{i},TargetSI{j},lambda);
    end
end

%Find the largest correlation found
I = getBestMatch(C,Target);
% [~,I] = max(C);
U = pointTarget(I);
SI_source = SourceSI;
SI_target = TargetSI(U);

end

function best_match = getBestMatch(C,Target)

%N/2 size of correlations
nV = round(size(C,1)/2);
nP = size(C,2);

best_match = zeros(nP,1);

%Fourth spread formula
for i=1:nP
    %Get the extreme outliers with the third fourth spread of the data
    % (the 1/4 of the data with the highest values)
    b = 3*median(maxk(C(:,i),nV));
    outliers = C(:,i) > b;
    idx = find(outliers);
    
    %Get geometric center of outliers
    gC = sum(Target.vertices(outliers,:),1)/numel(outliers);
    %Get distances of individual points to the geometric center
    distances = sqrt(sum((Target.vertices(outliers,:)-gC).^2,2));
    maxD = max(distances);
    
    %Get normalized distances
    n_distances = distances/maxD;
    weighted_confidence = C(outliers,:)./n_distances;
    
    %Get index of best outlier
    [~,idxMax] = max(weighted_confidence);
    best_match(i) = idx(idxMax);
end


end