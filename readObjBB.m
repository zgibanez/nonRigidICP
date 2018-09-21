% readObj: original code by BoffinBlogger, who writes ...
% ... I'm a big Matlab fan but, as I've posted before, file reading can be 
% really slow. I've just had to do some work where I need to read in some 
% quite large (100MB or so) Wavefront OBJ files. I can't use my previous 
% work around because OBJ files have a format that varies from line to 
% line. The format is actually very straightforward: you read the first 
% letter (or occasionally first 2 letters) of a line and that tells you 
% what the data represents, and then the rest of the line is the data 
% itself. The bits I'm interested in are the vertices (v followed by 3 
% floating point numbers), and faces (f followed by indices to the vertices
% that make up the face). I was using a line by line approach to this and 
% Matlab was clearly getting horribly bogged down in the parsing. 
% It doesn't help that the face definitions can optionally include 
% texture coordinates and a normal and this extra data is separated by / 
% characters. Each file was taking several minutes to read which was far 
% too slow. However I actually found a very good solution that speeds 
% things up by a factor of about 10. That was to read the data in all at 
% once and then use a series of regular expressions to separate up the 
% lines I was interested in. I could then parse them efficiently 
% using sscanf because I already know what is in each line.

% Here is the implementation. It's not particularly flexible and won't 
% read every OBJ file. That's because I know that my OBJ files are 
% triangular meshes so the face lines will always have a list of 3 
% vertexes and I've optimised the parsing code for that (and in my files 
% there is a normal but no texture index so that's the condition it 
% checks for first). A general implementation would be much slower 
% because Matlab doesn't like to handle lists of varying sizes and 
% faces would have to be a cell array of cell arrays rather than a nice 
% [num_triangles, 3] array.

% -----------------------------------------------------------------
% Updated by Nick Pears 2016 to read texture coordinates and normals
% and to provide some file read time information. Code speed improved
% by checking expected format first .. eg for headspace data we know that
% there will be vertex,normal and texture data in the file. An
% automatic format check could be implemented later

function [vertex, faces, vt, facest, vn, facesn] = readObjBoffinBlogger(filename)


% Modes 
verboseMode =1;
displayInfoMode = 0;

if(displayInfoMode)
    displayInfoReadObj.m
end

% Capture a start time
tStart = tic;


fid = fopen(filename);
if fid<0
    error(['Cannot open ' filename '.']);
end
[str, count] = fread(fid, [1,inf], 'uint8=>char'); 
if(verboseMode)
    fprintf('Read %d characters from %s\n', count, filename);
end
fclose(fid);

% get the lines that give vertex coordinates
vertex_lines = regexp(str,'v [^\n]*\n', 'match');
vertex = zeros(length(vertex_lines), 3);   % allocate space for vertex coordinates

% get the vertex coordinates line by line
for i = 1: length(vertex_lines)
    v = sscanf(vertex_lines{i}, 'v %f %f %f');
    vertex(i, :) = v';
end


% get the lines that give vt coordinates
vt_lines = regexp(str,'vt [^\n]*\n', 'match');
vt = zeros(length(vt_lines), 2);   % allocate space for texture coordinates

% get the texture coordinates line by line
for i = 1: length(vt_lines)
    t = sscanf(vt_lines{i}, 'vt %f %f');
    vt(i, :) = t';
end

% get the lines that give vn (normal) coordinates
vn_lines = regexp(str,'vn [^\n]*\n', 'match');
vn = zeros(length(vn_lines), 3);   % allocate space for texture coordinates

% get the texture coordinates line by line
for i = 1: length(vn_lines)
    n = sscanf(vn_lines{i}, 'vn %f %f %f');
    vn(i, :) = n';
end

% get the lines that describe the faces
face_lines = regexp(str,'f [^\n]*\n', 'match');
faces = zeros(length(face_lines), 3);    % allocate space for the faces (vertices)
facesn= zeros(length(face_lines), 3);    % allocate space for the faces (normals)
facest= zeros(length(face_lines), 3);    % allocate space for the faces (textures)

% get the faces, line by line
for i = 1: length(face_lines)
    
    % For speed, check the format that you are likely to get first
    % Check format for: vertices and textures and normals
    f = sscanf(face_lines{i}, 'f %d/%d/%d %d/%d/%d %d/%d/%d');
    if (length(f) == 9) % face
        faces(i, 1) = f(1);
        faces(i, 2) = f(4);
        faces(i, 3) = f(7);
        
        facest(i, 1) = f(2);
        facest(i, 2) = f(5);
        facest(i, 3) = f(8);
        
        facesn(i, 1) = f(3);
        facesn(i, 2) = f(6);
        facesn(i, 3) = f(9);
     
        continue
    end
    
    % Check format for: vertices with normals BUT no texture coordinates
    f = sscanf(face_lines{i}, 'f %d//%d %d//%d %d//%d');
    if (length(f) == 6) % face
        faces(i, 1) = f(1);
        faces(i, 2) = f(3);
        faces(i, 3) = f(5);
        
        facesn(i, 1) = f(2);
        facesn(i, 2) = f(4);
        facesn(i, 3) = f(6);
        
        continue
    end
    
    % Check format for: vertices only
    f = sscanf(face_lines{i}, 'f %d %d %d');
    if (length(f) == 3) % face
        faces(i, :) = f';
        continue
    end

    % Check format for: vertices and textures BUT no normals
    f = sscanf(face_lines{i}, 'f %d/%d %d/%d %d/%d');
    if (length(f) == 6) % face
        faces(i, 1) = f(1);
        faces(i, 2) = f(3);
        faces(i, 3) = f(5);
        
        facest(i, 1) = f(2);
        facest(i, 2) = f(4);
        facest(i, 3) = f(6);
        
        continue
    end
    
   
end

fileReadTime = toc(tStart);
% fDisp(['OBJ read time: ' num2str(fileReadTime) ' (s)'],verboseMode);


return
