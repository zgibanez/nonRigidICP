function varargout = disp_normals(varargin)
% DISP_NORMALS MATLAB code for disp_normals.fig
%      DISP_NORMALS, by itself, creates a new DISP_NORMALS or raises the existing
%      singleton*.
%
%      H = DISP_NORMALS returns the handle to a new DISP_NORMALS or the handle to
%      the existing singleton*.
%
%      DISP_NORMALS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISP_NORMALS.M with the given input arguments.
%
%      DISP_NORMALS('Property','Value',...) creates a new DISP_NORMALS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before disp_normals_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to disp_normals_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help disp_normals

% Last Modified by GUIDE v2.5 09-Jul-2018 09:32:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @disp_normals_OpeningFcn, ...
                   'gui_OutputFcn',  @disp_normals_OutputFcn, ...
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


% --- Executes just before disp_normals is made visible.
function disp_normals_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to disp_normals (see VARARGIN)

%Meshes
handles.SourceMesh = [];
handles.TargetMesh = [];

%Current data to display is a generic Mesh
% axeshandle = findobj('Tag','patchPlot');
% m = readObj('sample.obj');
% m = rmfield(m,'normals');
% handles.mesh = m;
% patch(axeshandle, handles.mesh, 'facecolor', 'blue', 'EdgeColor',  'none', ...
%               'FaceAlpha', 0.5);
%     material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
%     view([60,30]); axis equal; axis manual;

%Grey out unusable options
set(handles.CheckSource,'Enable','off');
set(handles.CheckTarget,'Enable','off');
set(handles.checkbox_module,'Value',0);

%Show camera toolbar
cameratoolbar('Show');


% Choose default command line output for disp_normals
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes disp_normals wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = disp_normals_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in CheckSource.
function CheckSource_Callback(hObject, eventdata, handles)
% hObject    handle to CheckSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckSource


% --- Executes on button press in CheckTarget.
function CheckTarget_Callback(hObject, eventdata, handles)
% hObject    handle to CheckTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CheckTarget


% --- Executes on button press in LoadBttnSource.
function LoadBttnSource_Callback(hObject, eventdata, handles)
% hObject    handle to LoadBttnSource (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path] = uigetfile('*.obj','C:\Users\pey_l\Desktop\nonrigidICP\test_meshes');

if(sum(file ~= 0) && sum(path~= 0))
    
    msginfo(handles,'Reading mesh file...');
    m = readObj(strcat(path,file));
    hObject.UserData{1} = rmfield(m,'normals');
    
      
    %Draw difference
%     axes(handles.diffPlot);
%     h = gca;
%     h = figure();
%     v = m.vertices;
%     f = 50;
%     N = (normals-normalsLB);
%     N_matlab = vertexNormal(m.vertices,m.faces);
% %     quiver3(v(:,1),v(:,2),v(:,3),v(:,1)+f*normals(:,1),v(:,2)+f*normals(:,2),v(:,3)+f*normals(:,3),'Color','r','LineWidth',2);
%     quiver3(v(:,1),v(:,2),v(:,3),v(:,1)+f*N(:,1),v(:,2)+f*N(:,2),v(:,3)+f*N(:,3),'Color','r','LineWidth',2,'AutoScale','off');
%     hold on;
% %     quiver3(v(:,1),v(:,2),v(:,3),v(:,1)+f*normalsLB(:,1),v(:,2)+f*normalsLB(:,2),v(:,3)+f*normalsLB(:,3),'Color','b','LineWidth',2);
% %      quiver3(v(:,1),v(:,2),v(:,3),v(:,1)+f*normalsLB(:,1),v(:,2)+f*normalsLB(:,2),v(:,3)+f*normalsLB(:,3),'Color','b','LineWidth',2);
% 
% %     ang = atan2(vecnorm(cross(normals,normalsLB,2)')',dot(normals,normalsLB,2));
%     dv = sqrt(sum(N.*N,2));
% %     dv = normr(normals-normalsLB);
%     scatter3(v(:,1),v(:,2),v(:,3),dv*30);
%     axis equal;
%     cameratoolbar('show');
    
    %Link rotation of plots    
    set(handles.CheckSource,'Enable','on')
    handles.LoadBttnSource.String = 'Load';
    set(handles.LoadBttnSource,'Enable','on')
    
    plotMeshes(handles);
end

function plotMeshes(handles)

m = get(handles.LoadBttnSource,'UserData');
m = m{1};
%Draw normal directions
axes(handles.patchPlot);
hn = gca;
cla;
%     normals = drawNormals(hn,hObject.UserData{1});
normals = drawMatlabNormals(hn,m);
    
c = get(handles.checkbox_module,'Value');
axes(handles.patchPlotLB);
hl = gca;
cla; %Clear axes

    %Draw LB vector
%    normals = drawNormals(hl,hObject.UserData{1});
    normalsLB = drawLB(hl,m,c);

Link = linkprop([hn,hl],{'CameraPosition','CameraUpVector'});
setappdata(gcf, 'StoreTheLink', Link);


%Displays a message in the info box
function msginfo(handles,message)
set(handles.text_info,'String',message);

function Nmat = drawMatlabNormals(h,Mesh)
v = Mesh.vertices;
Nmat = vertexNormal(Mesh.vertices,Mesh.faces);
f = 10; %multiplication factor for extending arrows
% quiver3(h,v(:,1),v(:,2),v(:,3),v(:,1)+f*Hvn(:,1),v(:,2)+f*Hvn(:,2),v(:,3)+f*Hvn(:,3),'Color','b','LineWidth',2);
quiver3(h,v(:,1),v(:,2),v(:,3),v(:,1)+f*Nmat(:,1),v(:,2)+f*Nmat(:,2),v(:,3)+f*Nmat(:,3),'Color','b','LineWidth',2);

hold on;

%And patch plot
patch(h,'vertices',Mesh.vertices,'faces', Mesh.faces, ... 
                                'facecolor', 'b', ... 
                                'EdgeColor',  'k', ...
                                'EdgeAlpha', 0.8, ...
                                'FaceAlpha', 0.4);
                            
axis equal;

function Hvn = drawLB(h,Mesh,plotModule)
[A,C] = LP_exp(Mesh);
L = A*C;
v = Mesh.vertices;
Hvn = L*v;
nHvn = normr(Hvn); %Normalization for visualization;
f = 50; %multiplication factor for extending arrows

if(plotModule)
    markSizes = sqrt(sum(Hvn.*Hvn,2))+1;
    scatter3(h,v(:,1),v(:,2),v(:,3),markSizes);
else
    % quiver3(h,v(:,1),v(:,2),v(:,3),v(:,1)+f*Hvn(:,1),v(:,2)+f*Hvn(:,2),v(:,3)+f*Hvn(:,3),'Color','b','LineWidth',2);
    quiver3(h,v(:,1),v(:,2),v(:,3),v(:,1)+f*nHvn(:,1),v(:,2)+f*nHvn(:,2),v(:,3)+f*nHvn(:,3),'Color','b','LineWidth',2);

    hold on;

    %And patch plot
    patch(h,'vertices',Mesh.vertices,'faces', Mesh.faces, ... 
                                    'facecolor', 'b', ... 
                                    'EdgeColor',  'k', ...
                                    'EdgeAlpha', 0.8, ...
                                    'FaceAlpha', 0.4);
end                            
axis equal;

function nvn = drawNormals(h,Mesh)
    if(~isfield(Mesh,'vertex_normals'))
        disp('Mesh does not have a "vertex normal" field. Calculating normals...');
        Mesh.vertex_normals = calc_normals(Mesh);
    end

    %Display quiver plot
    v = Mesh.vertices;
    vn = Mesh.vertex_normals;
    nvn = normr(vn);
    f = 11; %multiplication factor
    quiver3(h,v(:,1),v(:,2),v(:,3),v(:,1)+f*vn(:,1),v(:,2)+f*vn(:,2),v(:,3)+f*vn(:,3),'Color','b','LineWidth',2,'AutoScale','on');
    hold on;

    %And patch plot
    patch(h,'vertices',Mesh.vertices,'faces', Mesh.faces, ... 
                                    'facecolor', 'b', ... 
                                    'EdgeColor',  'k', ...
                                    'EdgeAlpha', 0.8, ...
                                    'FaceAlpha', 0.4);

%     vn2 = vertexNormal(Mesh.vertices,Mesh.faces);
%     quiver3(v(:,1),v(:,2),v(:,3),v(:,1)+4*vn2(:,1),v(:,2)+4*vn2(:,2),v(:,3)+4*vn2(:,3),'LineWidth',1);
%     for i = 1:n
%         quiver3(v(i,1),v(i,2),v(i,3),v(i,1)+vn(i,1),v(i,2)+vn(i,2),v(i,3)+vn(i,3));
%     end
    axis equal;

% --- Executes on button press in LoadBttnTarget.
function LoadBttnTarget_Callback(hObject, eventdata, handles)
% hObject    handle to LoadBttnTarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_module.
function checkbox_module_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_module
plotMeshes(handles);