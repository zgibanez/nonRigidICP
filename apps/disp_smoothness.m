function varargout = disp_smoothness(varargin)
% DISP_SMOOTHNESS MATLAB code for disp_smoothness.fig
%      DISP_SMOOTHNESS, by itself, creates a new DISP_SMOOTHNESS or raises the existing
%      singleton*.
%
%      H = DISP_SMOOTHNESS returns the handle to a new DISP_SMOOTHNESS or the handle to
%      the existing singleton*.
%
%      DISP_SMOOTHNESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISP_SMOOTHNESS.M with the given input arguments.
%
%      DISP_SMOOTHNESS('Property','Value',...) creates a new DISP_SMOOTHNESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before disp_smoothness_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to disp_smoothness_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help disp_smoothness

% Last Modified by GUIDE v2.5 26-Jul-2018 14:53:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @disp_smoothness_OpeningFcn, ...
                   'gui_OutputFcn',  @disp_smoothness_OutputFcn, ...
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


% --- Executes just before disp_smoothness is made visible.
function disp_smoothness_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to disp_smoothness (see VARARGIN)

cameratoolbar('Show');
set(handles.showSmoothnessLeft, 'Enable','Off');
set(handles.showSmoothnessRight, 'Enable','Off');

%Link rotation between plots
% Link = linkprop([handles.axesLeft,handles.axesRight],{'CameraPosition','CameraUpVector'});
% setappdata(gcf, 'StoreTheLink', Link);

% Choose default command line output for disp_smoothness
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes disp_smoothness wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = disp_smoothness_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadButton1.
function LoadButton1_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadMesh(handles.axesLeft);
set(handles.showSmoothnessLeft, 'Enable','On');

% --- Executes on button press in LoadButton2.
function LoadButton2_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
loadMesh(handles.axesRight);
set(handles.showSmoothnessRight, 'Enable','On');

function loadMesh(h)
[filename, pathname] = uigetfile({'.obj'},'Select data');
if(filename == 0)
    return;
end
Mesh = readObj([pathname filename]);
Mesh.S = calcSmoothness(Mesh);
set(h,'UserData',Mesh);
changeMesh(h,0);


function changeMesh(h, showSmooth)
cla(h);
Mesh = h.UserData;
if(~showSmooth)
    color = 'b';
          patch(h, ...
            'Faces', Mesh.faces, ...
            'Vertices', Mesh.vertices, ...
            'Facecolor', color, ...
            'EdgeColor',  'k', ...
            'EdgeAlpha', 0.1, ...
            'FaceAlpha', 0.5);  
else

    c = log(Mesh.S);
    color = mat2gray(c);
              patch(h, ...
            'Faces', Mesh.faces, ...
            'Vertices', Mesh.vertices, ...
            'FaceVertexCData', color, ...
            'Facecolor', 'interp', ...
            'EdgeColor',  'k', ...
            'EdgeAlpha', 0.1, ...
            'FaceAlpha', 0.5);  
end

axis equal; axis manual;

%Calculates smoothness
function S =calcSmoothness(Mesh)
[A,C] = LP_exp(Mesh);
S1 = sqrt(sum((A*C*Mesh.vertices).^2,2));

MeshS = smoothpatch(Mesh);
[A,C] = LP_exp(MeshS);
S2 = sqrt(sum((A*C*MeshS.vertices).^2,2));

S = abs(S1-S2);

% --- Executes on button press in showSmoothnessLeft.
function showSmoothnessLeft_Callback(hObject, eventdata, handles)
% hObject    handle to showSmoothnessLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showSmoothnessLeft
changeMesh(handles.axesLeft,hObject.Value);

% --- Executes on button press in showSmoothnessRight.
function showSmoothnessRight_Callback(hObject, eventdata, handles)
% hObject    handle to showSmoothnessRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
changeMesh(handles.axesRight,hObject.Value);
% Hint: get(hObject,'Value') returns toggle state of showSmoothnessRight
