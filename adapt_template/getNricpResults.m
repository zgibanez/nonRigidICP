
%Get landmark .xls file
Landmarks = readtable('landmarkTest.xls');

%Set output directory
outPath = 'D:\GMT\THESIS\Results\amberg\';

%Directory where the .OBJs are located
meshPath = 'C:\Users\pey_l\Desktop\nonrigidICP\matlab_data\meshes_neutral\';
targetMeshes = dir([meshPath '*.obj']);

%Load Source mesh
load SourceDiego

%Get landmarks of the source mesh
[row,~] = find(strcmp(Landmarks.File,'template'));
LSource = table2array(Landmarks(row,2:end));

%Test with every mesh that is included on the expressive mesh folder
[rows,cols] = size(Landmarks);
for i=1:rows
    
    %Find if mesh has a landmark entry in the .xls
    [row,~] = find(strcmp(Landmarks.File,targetMeshes(i).name));
    
    %If so, proceed with the nricp
    if(~isempty(row))
        
        %Check if file already exists
        [~,n,~] = fileparts(targetMeshes(i).name);
        outFile = strcat(outPath,'GUIDATA_',n);
        if (exist(strcat(outFile,'.mat'), 'file') == 2)
            fprintf('File already exists. This mesh will be skipped. \n');
            continue;
        end
        %Load Target mesh
        Target = readObj(strcat(meshPath,targetMeshes(i).name));
        fprintf(strcat('Performing nricp of Target mesh :',targetMeshes(i).name,'...\n'));
        %Get landmarks from table
        LTarget = table2array(Landmarks(row,2:end));
        
        %Perform nricp & save data
        D = nricp_custom(Target,Source,1,LTarget,LSource);
        [~,n,~] = fileparts(targetMeshes(i).name);
        save(strcat(outFile),'D');
        
    end

end