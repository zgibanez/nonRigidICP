
path = 'D:\GMT\THESIS\expressive\';

files = dir(strcat(path,'*.mat'));
[n,~] = size(files);

for i=1:n
    name = files(i).name;
    load(strcat(path,files(i).name));
    D = Data;
    save(strcat(path,files(i).name),'D');
end