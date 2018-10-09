function varargout = disp_correspondence(varargin)
% DISP_CORRESPONDENCE MATLAB code for disp_correspondence.fig
%      DISP_CORRESPONDENCE, by itself, creates a new DISP_CORRESPONDENCE or raises the existing
%      singleton*.
%
%      H = DISP_CORRESPONDENCE returns the handle to a new DISP_CORRESPONDENCE or the handle to
%      the existing singleton*.
%
%      DISP_CORRESPONDENCE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISP_CORRESPONDENCE.M with the given input arguments.
%
%      DISP_CORRESPONDENCE('Property','Value',...) creates a new DISP_CORRESPONDENCE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before disp_correspondence_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to disp_correspondence_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help disp_correspondence

% Last Modified by GUIDE v2.5 09-Sep-2018 11:26:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @disp_correspondence_OpeningFcn, ...
                   'gui_OutputFcn',  @disp_correspondence_OutputFcn, ...
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


% --- Executes just before disp_correspondence is made visible.
function disp_correspondence_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to disp_correspondence (see VARARGIN)

% Choose default command line output for disp_correspondence
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Disable checkboxes
cameratoolbar('Show');
set(handles.checkbox_left,'Enable','off');
set(handles.checkbox_middle,'Enable','off');
set(handles.checkbox_right,'Enable','off');

%Set callback function for left click
set(gcf, 'WindowButtonDownFcn', {@callbackClickA3DPoint, handles});


% UIWAIT makes disp_correspondence wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = disp_correspondence_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_left.
function pushbutton_left_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadMesh(handles.axes_left,handles);

% --- Executes on button press in checkbox_left.
function checkbox_left_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_left (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_left
changeMesh(handles.axes_left,hObject.Value);


% --- Executes on button press in pushbutton_middle.
function pushbutton_middle_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadMesh(handles.axes_middle,handles);

% --- Executes on button press in checkbox_middle.
function checkbox_middle_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_middle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_middle
changeMesh(handles.axes_middle,hObject.Value);


% --- Executes on button press in pushbutton_right.
function pushbutton_right_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadMesh(handles.axes_right,handles);

% --- Executes on button press in checkbox_right.
function checkbox_right_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_right (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
changeMesh(handles.axes_right,hObject.Value);
% Hint: get(hObject,'Value') returns toggle state of checkbox_right


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

% if(~strcmp(get(gca,'Tag'),'display'))
%     return;
% end
D = get(gca,'UserData');
pointCloud = D.vertices';

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
value = get(handles.checkbox_landmark,'Value');

if isempty(h) % if it's the first click (i.e. no previous point to delete)
    
    % highlight the selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'r.', 'MarkerSize', 20); 
    set(h,'Tag','pt'); % set its Tag property for later use   

else % if it is not the first click
    if(~value)
        delete(h); % delete the previously selected point
    end
    % highlight the newly selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'r.', 'MarkerSize', 20);  
    set(h,'Tag','pt');  % set its Tag property for later use

end

hS = {handles.axes_left, handles.axes_middle, handles.axes_right};

for i=1:numel(hS)
    if(hS{i}==gca) 
        continue;
    end
    M = get(hS{i},'UserData');
    if(isempty(M))
        continue;
    else
        try
            Pt = M.vertices(pointCloudIndex,:);
        catch
            writeError(handles,'Error: cannot find point across other axes');
        end
        hp = findobj(hS{i},'Tag','Pt');
        if(~isempty(hp))
            if(~value)
                delete(hp);
            end
        end 
        hp = plot3(hS{i},Pt(1), Pt(2), ...
            Pt(3), 'r.', 'MarkerSize', 20);
        set(hp,'Tag','Pt');
    end
end

set(handles.text_point,'String',strcat('Point Selected: ',num2str(pointCloudIndex)));


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



function loadMesh(h,handles)
[filename, pathname] = uigetfile({'*.obj';'*.mat'},'Select data');
if(filename == 0)
    return;
end
% [~,~,ext] = fileparts(strcat(pathname,filename));
% if(strcmp(ext,'.mat'))
%     Mesh = load(strcat(pathname,filename));
% end

Mesh = readObj([pathname filename]);
set(h,'UserData',Mesh);
changeMesh(h,0);


function writeError(handles,msg)
set(handles.error_msg,'String',msg);


% --- Executes on button press in checkbox_landmark.
function checkbox_landmark_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_landmark (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_landmark
