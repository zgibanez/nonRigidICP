scans_folder = 'C:\Users\pey_l\Desktop\nonrigidICP\matlab_data\meshes_neutral\';

meshes = dir(strcat(scans_folder,'*.obj'));
pSource = [372 3870 6704 7512 2979 37 115 10629 9219 6986 11188 8787 4899 3402 3898 6281];
landmarks = readtable('landmarks.xls');
% load('SourceDiego');

for i=1:size(meshes,1)
    [idx,~] = find( strcmp(meshes(i).name, landmarks.File));
    
    if(~isempty(idx))
        SourceT = Source;
        pTarget = table2array(landmarks(idx,2:end));
        Target = readObj(strcat(scans_folder,meshes(i).name));
        [Target.vertices, Source.vertices] = prealign(Target,Source);
        [~,~,~,~,Source.vertices] = icp_custom(Target,Source);
        SourceT.vertices = adaptTemplate(Target,Source, pTarget, pSource);
        write_obj(strcat('source_',meshes(i).name),SourceT.vertices,SourceT.faces);
    end
end