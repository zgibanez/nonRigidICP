%Read Obj
% targetMeshes = dir('D:\GMT\THESIS\Results\expressive\*_target.obj');
% sourceMeshes = dir('D:\GMT\THESIS\Results\expressive\*_source.obj');
% ambergMeshes = dir('D:\GMT\THESIS\Results\amberg\expressive\*.obj');
targetMeshes = dir('D:\GMT\THESIS\Results\neutral_2\*_target.obj');
sourceMeshes = dir('D:\GMT\THESIS\Results\neutral_2\*_source.obj');
ambergMeshes = dir('D:\GMT\THESIS\Results\amberg\*.obj');

N = numel(ambergMeshes);

meanA = zeros(N,1);
meanO = meanA;
names = cell(N,1);
%Perform KNN
for i=1:N
    [Target.vertices, Target.faces] = readObjBB([targetMeshes(i).folder '\' targetMeshes(i).name]);
    boundary = compute_boundary_custom(Target.faces);
    
    [SourceO.vertices, SourceO.faces] = readObjBB([sourceMeshes(i).folder '\' sourceMeshes(i).name]);
    [SourceA.vertices, SourceA.faces] = readObjBB([ambergMeshes(i).folder '\' ambergMeshes(i).name]);
    
%     disp_meshes(SourceA,Target);
    
    UA = knnsearch(Target.vertices,SourceA.vertices); 
    UO = knnsearch(Target.vertices,SourceO.vertices);
        
    dO = sqrt(sum((SourceO.vertices-Target.vertices(UO,:)).^2,2));
    dA = sqrt(sum((SourceA.vertices-Target.vertices(UA,:)).^2,2));
    
    %Remove KNN of borders
    [~,locA] = ismember(boundary,UA); locA(locA == 0) = [];
    dA(locA) = [];
    [~,locO] = ismember(boundary,UO); locO(locO == 0) = [];
    dO(locO) = [];

    meanA(i) = mean(dA);
    meanO(i) = mean(dO);
    
    [~,n,~] = fileparts(ambergMeshes(i).name);
    names{i} = n;
    
    fprintf('\n done \n');

end

T2 = table(names,meanA,meanO);

%Get mean distance