%Set output directory
outPath = 'D:\GMT\THESIS\Results\amberg\expressive\';

%Directory where the .OBJs are located
meshPath = 'D:\GMT\THESIS\Results\expressive\';
targetMeshes = dir([meshPath '*_target.obj']);

%Load Source mesh
load SourceDiego
Source.normals = calc_normals(Source);
newSource = Source;

% Specify that surface normals are available and can be used.
Options.useNormals = 0;
Options.normalWeighting = 0;
% Specify that the source deformations should be plotted.
Options.plot = 1;
%   Maximum iterations
Options.maxIter = 10;
%   Initialize with rigid ICP
Options.rigidInit = 1;

N = numel(targetMeshes);

for i=1:N
    %Check if file already exists
    if (exist([outPath targetMeshes(i).name], 'file') == 2)
        fprintf('File already exists. This mesh will be skipped. \n');
        continue;
    end
    %Load Target mesh
    Target = readObj(strcat(meshPath,targetMeshes(i).name));
    Target.normals = calc_normals(Target);
    fprintf(strcat('Performing nricp of Target mesh :',targetMeshes(i).name,'... \n'));
    
    %Perform icp to avoid local minima
    [Target.vertices,newSource.vertices] = prealign(Target,Source);
    [~, ~, ~, ~, newSource.vertices] = icp_custom(Target,newSource);
    
    %Perform nricp & save data
    newSource.vertices = nricp(newSource,Target,Options);
    write_obj([outPath targetMeshes(i).name],newSource.vertices,newSource.faces);
end