function varargout = match_correspondences(varargin)
% MATCH_CORRESPONDENCES MATLAB code for match_correspondences.fig
%      MATCH_CORRESPONDENCES, by itself, creates a new MATCH_CORRESPONDENCES or raises the existing
%      singleton*.
%
%      H = MATCH_CORRESPONDENCES returns the handle to a new MATCH_CORRESPONDENCES or the handle to
%      the existing singleton*.
%
%      MATCH_CORRESPONDENCES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATCH_CORRESPONDENCES.M with the given input arguments.
%
%      MATCH_CORRESPONDENCES('Property','Value',...) creates a new MATCH_CORRESPONDENCES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before match_correspondences_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to match_correspondences_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help match_correspondences

% Last Modified by GUIDE v2.5 22-Aug-2018 12:19:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @match_correspondences_OpeningFcn, ...
                   'gui_OutputFcn',  @match_correspondences_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before match_correspondences is made visible.
function match_correspondences_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to match_correspondences (see VARARGIN)

% Choose default command line output for match_correspondences
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes match_correspondences wait for user response (see UIRESUME)
% uiwait(handles.figure1);
cameratoolbar('Show');

%Initialize buttons
set(handles.button_template,'Enable','off');
set(handles.button_target,'Enable','off');

% --- Outputs from this function are returned to the command line.
function varargout = match_correspondences_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popup_feature.
function popup_feature_Callback(hObject, eventdata, handles)
% hObject    handle to popup_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_feature contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_feature
highlightTemplatePoint(handles);

% --- Executes during object creation, after setting all properties.
function popup_feature_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_target.
function button_target_Callback(hObject, eventdata, handles)
% hObject    handle to button_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_select,'String','Loading target mesh...');
fileName = loadMesh(handles.axes_target,handles);
set(handles.text_select,'String','Target mesh loaded.');

%This UI element holds the target .obj name
set(handles.button_target,'UserData',fileName);

%Set callback function for left click
set(gcf, 'WindowButtonDownFcn', {@callbackClickA3DPoint, handles});


% --- Executes on button press in button_confirm.
function button_confirm_Callback(hObject, eventdata, handles)
% hObject    handle to button_confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in button_landmark.
function button_landmark_Callback(hObject, eventdata, handles)
% hObject    handle to button_landmark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname]=uigetfile({'*.xls';'*.csv';'*.xlsl'},'Select landmark table');
if(filename == 0)
    return;
end
Landmarks = readtable([pathname filename]);
hObject.UserData = Landmarks;
set(handles.button_template,'Enable','on');

%Set the names of colums as the popup menu options
set(handles.popup_feature,'String',Landmarks.Properties.VariableNames(2:end));
set(handles.text_select,'String','To continue load a template .obj file.');

%Add a vectors where all the idxs will be storaged
cols=zeros(size(Landmarks,2)-1,1);
set(handles.button_write,'UserData',cols);


% --- Executes on button press in button_template.
function button_template_Callback(hObject, eventdata, handles)
% hObject    handle to button_template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_select,'Visible','On');
set(handles.text_select,'String','Loading template mesh...');
loadMesh(handles.axes_template,handles);
set(handles.text_select,'String','Template mesh loaded.');
set(handles.button_target,'Enable','on');
markTemplatePoints(handles);
highlightTemplatePoint(handles);


function changeMesh(h, showEdges)
axes(h);
cla(gca);
Mesh = h.UserData;
if(showEdges)
          patch(h, ...
            'Faces', Mesh.faces, ...
            'Vertices', Mesh.vertices, ...
            'Facecolor', 'b', ...
            'EdgeColor',  'k', ...
            'EdgeAlpha', 0.1, ...
            'FaceAlpha', 0.5);  
else

            patch(h, ...
            'Faces', Mesh.faces, ...
            'Vertices', Mesh.vertices, ...
            'facecolor', 'b', ...
            'EdgeColor',  'none', ...
            'EdgeAlpha', 0.3, ...
            'FaceAlpha', 0.8);
end

material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
axis equal; axis manual;
view([60,30]); hold on;

function filename = loadMesh(h,handles)
[filename, pathname] = uigetfile({'*.obj';'*.mat'},'Select data');
if(filename == 0)
    return;
end

Mesh = readObj([pathname filename]);
% [vertex, faces, ~, ~, ~, ~] = readObjBB([pathname filename]);
% Mesh.vertices = vertex;
% Mesh.faces = faces;
set(h,'UserData',Mesh);
changeMesh(h,0);


function pointCloudIndex = callbackClickA3DPoint(src, eventData, handles)
% CALLBACKCLICK3DPOINT mouse click callback function for CLICKA3DPOINT
%
%   The transformation between the viewing frame and the point cloud frame
%   is calculated using the camera viewing direction and the 'up' vector.
%   Then, the point cloud is transformed into the viewing frame. Finally,
%   the z coordinate in this frame is ignored and the x and y coordinates
%   of all the points are compared with the mouse click location and the 
%   closest point is selected.
%
%   Babak Taati - May 4, 2005
%   revised Oct 31, 2007
%   revised Jun 3, 2008
%   revised May 19, 2009

if(~strcmp(get(gca,'Tag'),'axes_target'))
    return;
end
Mesh = get(gca,'UserData');
pointCloud = Mesh.vertices';

point = get(gca, 'CurrentPoint'); % mouse click position
camPos = get(gca, 'CameraPosition'); % camera position
camTgt = get(gca, 'CameraTarget'); % where the camera is pointing to

camDir = camPos - camTgt; % camera direction
camUpVect = get(gca, 'CameraUpVector'); % camera 'up' vector

% build an orthonormal frame based on the viewing direction and the 
% up vector (the "view frame")
zAxis = camDir/norm(camDir);    
upAxis = camUpVect/norm(camUpVect); 
xAxis = cross(upAxis, zAxis);
yAxis = cross(zAxis, xAxis);

rot = [xAxis; yAxis; zAxis]; % view rotation 

% the point cloud represented in the view frame
rotatedPointCloud = rot * pointCloud; 

% the clicked point represented in the view frame
rotatedPointFront = rot * point' ;

% find the nearest neighbour to the clicked point 
pointCloudIndex = dsearchn(rotatedPointCloud(1:2,:)', ... 
    rotatedPointFront(1:2));

h = findobj(gca,'Tag','pt'); % try to find the old point
selectedPoint = pointCloud(:, pointCloudIndex); 

if isempty(h) % if it's the first click (i.e. no previous point to delete)
    
    % highlight the selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'y.', 'MarkerSize', 20); 
    set(h,'Tag','pt'); % set its Tag property for later use   

else % if it is not the first click

    delete(h); % delete the previously selected point
    
    % highlight the newly selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'y.', 'MarkerSize', 20);  
    set(h,'Tag','pt');  % set its Tag property for later use
end

%Add selected vertex to landmarks
addTargetVertex(handles,pointCloudIndex);
v = get(handles.popup_feature,'Value');
nFeatures = numel(get(handles.button_write,'UserData'));
if(v < nFeatures)
    set(handles.popup_feature,'Value',v+1);
    highlightTemplatePoint(handles);
end

%Highlights feature point on template axes
function highlightTemplatePoint(handles)
landmarks = get(handles.button_landmark,'UserData');
[row, ~] = find(strcmp(landmarks.File, 'template'));
axes(handles.axes_template);

featureIdx = get(handles.popup_feature,'Value')+1;
pointCloudIndex = table2array(landmarks(row,featureIdx));
Mesh = get(gca,'UserData');
pointCloud = Mesh.vertices';
h = findobj(gca,'Tag','ptT'); % try to find the old point
selectedPoint = pointCloud(:, pointCloudIndex); 

if isempty(h) % if it's the first click (i.e. no previous point to delete)
    
    % highlight the selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'y.', 'MarkerSize', 20); 
    set(h,'Tag','ptT'); % set its Tag property for later use   

else % if it is not the first click

    delete(h); % delete the previously selected point
    
    % highlight the newly selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'y.', 'MarkerSize', 20);  
    set(h,'Tag','ptT');  % set its Tag property for later use
end

function markTemplatePoints(handles)
axes(handles.axes_template);

%Get template mesh
Mesh = get(handles.axes_template,'UserData');

%Get idx of landmarks
landmarks = get(handles.button_landmark,'UserData');
[row,~] = find(strcmp(landmarks.File,'template'));
pointsIdx = table2array(landmarks(row,2:end));

%Point them on the axes
    plot3(Mesh.vertices(pointsIdx,1), Mesh.vertices(pointsIdx,2), ...
        Mesh.vertices(pointsIdx,3), 'r.', 'MarkerSize', 15);  


% --- Executes on button press in button_write.
function button_write_Callback(hObject, eventdata, handles)
% hObject    handle to button_write (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Get landmark table
landmarkTable = get(handles.button_landmark,'UserData');
%Get target landmarks
landmarkTarget = get(handles.button_write,'UserData')';

%Construct new row with target name and landmarks
newRow{1} = get(handles.button_target,'UserData');
newRow = [newRow, num2cell(landmarkTarget)];

%Check if row with the target mesh name already exists on original table
[row,~] = find(strcmp(landmarkTable.File,newRow{1}));
if(~isempty(row))
    landmarkTable(row,:) = newRow;
else
    landmarkTable = [landmarkTable; newRow];
end
%Write the full table on a xls file
writetable(landmarkTable,'landmarkTest.xls');

%Save the table in the GUI information for future annotations
set(handles.button_landmark,'UserData',landmarkTable);
set(handles.text_select,'String', strcat('File', {' '}, newRow{1}, 'written successfully.'));

function addTargetVertex(handles,pointIdx)
targetLandmarks = get(handles.button_write,'UserData');
featureIdx = get(handles.popup_feature,'Value');
targetLandmarks(featureIdx) = pointIdx;

msg = strcat(num2str(nnz(targetLandmarks)), {' '},' of ', {' '}, num2str(numel(targetLandmarks)),  {' '}, 'landmarks annotated');
set(handles.text_landmarks,'String',msg);

%If all landmarks are written, we can save the file
if(nnz(targetLandmarks) == numel(targetLandmarks))
    set(handles.button_write,'Enable','on');
end

set(handles.button_write,'UserData',targetLandmarks);
markTargetPoints(handles);

function markTargetPoints(handles)
axes(handles.axes_target);
Mesh = get(handles.axes_target,'UserData');
V = Mesh.vertices;
targetL = get(handles.button_write,'UserData');
targetL(targetL ==0) = [];

h = findobj(gca,'Tag','Landmarks'); % try to find the old point

if(~isempty(h))
    delete(h)
end

h = plot3(V(targetL,1),V(targetL,2),V(targetL,3), ...
    'r.', 'MarkerSize', 15);
set(h,'Tag','Landmarks'); % set its Tag property for later use  
