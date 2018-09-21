function readObjFromFolder(folder)


if(~exist(strcat(folder,'\results'),'dir'))
    mkdir(strcat(folder,'\results'));
end

%Load Source mesh. Diegos one performs best now.
load('C:\Users\pey_l\Desktop\nonrigidICP\SourceDiego.mat');

objList = dir(strcat(folder,'\*.obj'));
outFolder = strcat(folder,'\results\');
for i=1:numel(objList.name)
    Target = readObj(strcat(folder,'/',objList(i).name));
    D = nricp_custom(Target,Source,1);
    Source = D.SourceT{end};
    %save(strcat(outFolder,objList.name{i},'Target'));
    save(strcat(outFolder,objList(i).name),'D');
end



end
