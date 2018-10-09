function varargout = disp_data(varargin)
% DISP_DATA MATLAB code for disp_data.fig
%      DISP_DATA, by itself, creates a new DISP_DATA or raises the existing
%      singleton*.
%
%      H = DISP_DATA returns the handle to a new DISP_DATA or the handle to
%      the existing singleton*.
%
%      DISP_DATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISP_DATA.M with the given input arguments.
%
%      DISP_DATA('Property','Value',...) creates a new DISP_DATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before disp_data_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to disp_data_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help disp_data

% Last Modified by GUIDE v2.5 15-Aug-2018 10:55:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @disp_data_OpeningFcn, ...
                   'gui_OutputFcn',  @disp_data_OutputFcn, ...
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


% --- Executes just before disp_data is made visible.
function disp_data_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to disp_data (see VARARGIN)


%Reset axes
cla(handles.display,'reset');
cla(handles.target_display,'reset');
cla(handles.source_display,'reset');

% Choose default command line output for disp_data
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Set default strings for pop-up menu
set(handles.popup_display,'String',{'Plot type...','Plot3','Patch','Curvature patch'});
set(handles.popup_point,'String',{'Select Point...'});

%Disable Source and Target radio buttons
set(handles.target_toggle,'Enable','off');
set(handles.source_toggle,'Enable','off');
set(handles.source_toggle,'Value',1);
%Hide unusable options
set(handles.display_panel,'Visible','Off');
%Set the Plot panel invisible
set(handles.plot_panel,'Visible','Off');
%Disable the plot panel buttons
set(handles.checkbox_LBm,'Enable','Off');

%Add the dummy option to the popup_point UI control
set(handles.popup_point,'UserData',0);

% UIWAIT makes disp_data wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = disp_data_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in source_toggle.
function source_toggle_Callback(hObject, eventdata, handles)
% hObject    handle to source_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of source_toggle
changeMesh(handles);

% --- Executes on button press in target_toggle.
function target_toggle_Callback(hObject, eventdata, handles)
% hObject    handle to target_toggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of target_toggle
changeMesh(handles);

% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'.mat'},'Select data');
if(filename == 0)
    return;
end

%Reset previous figures
cla(handles.display,'reset');
cla(handles.target_display,'reset');
cla(handles.source_display,'reset');

data = load([pathname filename]);
set(handles.load_button,'UserData',data);

%Set the Options panel visible
set(handles.display_panel,'Visible','On');

%Set the Plot panel visible
set(handles.plot_panel,'Visible','On');

%Set the slider properties and idx
setSlider(handles,data);
set(handles.iter_number,'String',1);

%Set the Target plot
dispPatchMeshes(handles.source_display,data.D.Source_Aligned,'b',0.8);
hs = gca;
set(hs,'Tag','source_display');
% set(ha,'ButtonDownFcn',{@copyRotationCallback,handles});
% vertices = data.D.Source_Aligned.vertices;
% plot3(vertices(:,1), vertices(:,2), vertices(:,3), 'r.'); 

%Set the Source plot
dispPatchMeshes(handles.target_display,data.D.Target,'b',0.8);
ht = gca;
set(ht,'Tag','target_display');

%Link rotation between plots
Link = linkprop([hs,ht],{'CameraPosition','CameraUpVector'});
setappdata(gcf, 'StoreTheLink', Link);
% set(ha,'ButtonDownFcn',{@copyRotationCallback,handles});
% vertices = data.D.Target.vertices;
% plot3(vertices(:,1), vertices(:,2), vertices(:,3), 'b.'); 

changeMesh(handles);
% %Add the click callback routine
% axes(handles.mesh_display)
% axes(handles.mesh_display);
% h = gcf;
% vertices = data.D.SourceT{1}.vertices;
% plot3(vertices(:,1), vertices(:,2), vertices(:,3), 'b.');
% ha = gca;
% set(ha,'Tag','display');
% cameratoolbar('Show');
% hold on;
% set(h, 'WindowButtonDownFcn', {@callbackClickA3DPoint, handles, vertices'});


function changeMesh(handles)
idx = handles.iter_slider.Value;
axes(handles.display);
contents = get(handles.popup_display,'String');
selection = contents{get(handles.popup_display,'Value')};
h = gcf;
cla; %Clear previous figure
data = handles.load_button.UserData;
vertices = data.D.SourceT{idx}.vertices;
ha = gca;
switch selection
    case 'Plot3'
        plot3(vertices(:,1), vertices(:,2), vertices(:,3), 'b.');
        %ha = gca;
%         cameratoolbar('Show');
        hold on;
    case 'Patch'
        %ha = gca;
%         cameratoolbar('Show');
        iter = get(handles.iter_slider,'Value');
        if(handles.source_toggle.Value)
            if(handles.checkbox_showedges.Value == false)
                dispPatchMeshes(ha,data.D.SourceT{iter},'b',0.3,0.5);
            else
                dispPatchMeshes(ha,data.D.SourceT{iter},'b',0.8,0);
            end
        end
        if(handles.target_toggle.Value)
           dispPatchMeshes(ha,data.D.Target,'r',0.3);
        end
    case 'Curvature patch'
        dispPatchMeshes(ha,data.D.Target,'r',0.3);
    otherwise
        return;

end

set(ha,'Tag','display');
cameratoolbar('Show');
set(h, 'WindowButtonDownFcn', {@callbackClickA3DPoint, handles, vertices'});

        



% function changeIDX(h_fig, h_idx, data)
%     set(handles.idx, 'string', callbackClickA3DPoint([], [], data.vertices'));
% end


%--- Set sliders properties after loading the data ---%
function setSlider(handles,data)

iterations = numel(data.D.SourceT);

h =handles.iter_slider;
set(h, 'Min', 1);
set(h, 'Max', iterations);
set(h, 'Value', 1);
set(h, 'SliderStep', [1/(iterations-1) , 1/(iterations-1) ]);

% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pointCloudIndex = callbackClickA3DPoint(src, eventData, handles, pointCloud)
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

if(~strcmp(get(gca,'Tag'),'display'))
    return;
end

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
        selectedPoint(3,:), 'r.', 'MarkerSize', 20); 
    set(h,'Tag','pt'); % set its Tag property for later use   

else % if it is not the first click

    delete(h); % delete the previously selected point
    
    % highlight the newly selected point
    h = plot3(selectedPoint(1,:), selectedPoint(2,:), ...
        selectedPoint(3,:), 'r.', 'MarkerSize', 20);  
    set(h,'Tag','pt');  % set its Tag property for later use

end

%Highlight point on the Source and Template plot
iter = get(handles.iter_slider,'Value');
U = handles.load_button.UserData.D.correspondences(:,iter);
Pt = handles.load_button.UserData.D.Target.vertices(U(pointCloudIndex),:);

% set(handles.idx, 'string', pointCloudIndex);
changePointIdx(handles,pointCloudIndex,U(pointCloudIndex));

axes(handles.source_display); hold on;
hp = findobj(gca,'Tag','PtS');
if(~isempty(hp))
    delete(hp)
end
hp = plot3(Pt(1), Pt(2), ...
        Pt(3), 'y.', 'MarkerSize', 20);
set(hp,'Tag','PtS'); hold off;

axes(handles.target_display); hold on;
hp = findobj(gca,'Tag','PtH');
if(~isempty(hp))
    delete(hp)
end
hp = plot3(Pt(1), Pt(2), ...
        Pt(3), 'y.', 'MarkerSize', 20);
set(hp,'Tag','PtH');

% fprintf('you clicked on point number %d\n', pointCloudIndex);

%Changes the selected point
function changePointIdx(handles,idx,U)
s = strcat('Source point selected:',num2str(idx));
s2 = strcat('Corr. target point:',num2str(U));
set(handles.idx, 'String',s);
set(handles.idxT, 'String',s2);
set(handles.idx, 'UserData',idx);

%Display patch mesh at graph handle h
function dispPatchMeshes(h,Mesh,color,alpha,edge_alpha)
if(nargin < 4)
    alpha = 0.5;
end
if(nargin < 5)
    edge_alpha = 0.0;
end
axes(h);
    patch('Vertices', Mesh.vertices, 'Faces', Mesh.faces, ...
            'facecolor', color, 'EdgeColor',  'k', ...
            'EdgeAlpha', edge_alpha, ...
            'FaceAlpha', alpha);
material dull; light; grid on; xlabel('x'); ylabel('y'); zlabel('z');
axis equal; axis manual;
view([60,30]);
hold on;


% --- Executes during object creation, after setting all properties.
function idx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to idx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on slider movement.
function iter_slider_Callback(hObject, eventdata, handles)
% hObject    handle to iter_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val=round(hObject.Value);
hObject.Value=val;

%Update Point plot
UpdatePointPlot(handles);

%Set the iteration number
set(handles.iter_number,'String',num2str(val));
changeMesh(handles);



% --- Executes during object creation, after setting all properties.
function iter_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iter_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popup_display.
function popup_display_Callback(hObject, eventdata, handles)
% hObject    handle to popup_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_display contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_display
contents = cellstr(get(hObject,'String'));
selection = contents{get(hObject,'Value')};

switch selection
    case 'Plot3'
        set(handles.target_toggle,'Enable','off');
        set(handles.source_toggle,'Enable','off');
    case 'Patch'
        set(handles.target_toggle,'Enable','on');
        set(handles.source_toggle,'Enable','on');
    case 'Curvature patch'
        set(handles.target_toggle,'Enable','on');
        set(handles.source_toggle,'Enable','on');
    otherwise
        return;
       
end
changeMesh(handles);

% --- Executes during object creation, after setting all properties.
function popup_display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on mouse press over axes background.
function copyRotationCallback(src, eventData, handles)
% hObject    handle to target_display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[az,el] = view(gca);
t = src.Tag;
if(t == 'source_display')
    axes(handles.target_display);
else
    axes(handles.source_display);
end
view([az, el]);


% --- Executes on button press in addpoint_button.
function addpoint_button_Callback(hObject, eventdata, handles)
% hObject    handle to addpoint_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = get(handles.idx,'UserData');
if(isempty(idx))
    set(handles.idx,'String',strcat('No point selected'));
    return;
else
    hObject.UserData = [hObject.UserData idx];
end

%Add point to the list

point_list = get(handles.popup_point,'String');
if(~ismember(point_list,num2str(idx)))
    point_list_data = get(handles.popup_point,'UserData');
    %Add checkbox information
    point_list_data = [point_list_data;0];
    set(handles.popup_point,'UserData',point_list_data);
    point_list{end+1} = num2str(idx);
else
    set(handles.idx,'String',strcat({'Point'},{' '},num2str(idx),' is already added'));
end
set(handles.popup_point,'String',point_list);


% --- Executes on selection change in popup_point.
function popup_point_Callback(hObject, eventdata, handles)
% hObject    handle to popup_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_point contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_point
% contents = cellstr(get(hObject,'String'));
selection = get(hObject,'Value');

if(selection == 1)
    set(handles.checkbox_LBm,'Enable','Off');
else
    %Enable all checkboxes in the subpanel
    set(handles.checkbox_LBm,'Enable','On');
    
    %Load switches configuration
    switches = get(hObject,'UserData');
    switches = switches(selection,:);
    set(handles.checkbox_LBm,'Value',switches(1));
end

% --- Executes during object creation, after setting all properties.
function popup_point_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_point (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in delete_button.
function delete_button_Callback(hObject, eventdata, handles)
% hObject    handle to delete_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switches = get(handles.popup_point,'UserData');
contents = cellstr(get(handles.popup_point,'String'));
selection = get(handles.popup_point,'Value');

%Delete the selected point from the popup_menu of the list
if(selection~=1)
    set(handles.popup_point,'Value',1);
    switches(selection,:) = [];
    contents(selection) = [];
    set(handles.popup_point,'UserData',switches);
    set(handles.popup_point,'String',contents);
else
    return;
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in checkbox_w1.
function checkbox_w1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_w1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_w1
UpdatePointPlot(handles);

% --- Executes on button press in checkbox_w2.
function checkbox_w2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_w2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_w2
UpdatePointPlot(handles);

% --- Executes on button press in checkbox_w3.
function checkbox_w3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_w3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_w3
UpdatePointPlot(handles);

% --- Executes on button press in checkbox_e1.
function checkbox_e1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_e1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_e1
UpdatePointPlot(handles);

% --- Executes on button press in checkbox_e2.
function checkbox_e2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_e2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_e2
UpdatePointPlot(handles);

% --- Executes on button press in checkbox_e3.
function checkbox_e3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_e3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_e3
UpdatePointPlot(handles);


function UpdatePointPlot(handles)
axes(handles.point_plot);
cla;
data = get(handles.load_button,'UserData');
data = data.D;
% h = findobj('Tag','h_point');
h = gca;
hold on;
if(get(handles.checkbox_w1,'Value'))
    plot(h,data.weights(1,:));
end

if(get(handles.checkbox_w2,'Value'))
    plot(h,data.weights(2,:));
end

if(get(handles.checkbox_w3,'Value'))
    plot(h,data.weights(3,:));
end

if(get(handles.checkbox_w4,'Value'))
    try
        plot(h,data.weights(4,:));
        hold on;
    catch
    end
end

if(get(handles.checkbox_e1,'Value'))
    plot(h,data.error(1,:));
end

if(get(handles.checkbox_e2,'Value'))
    plot(h,data.error(2,:));
end

if(get(handles.checkbox_e3,'Value'))
    plot(h,data.error(3,:));
end

if(get(handles.checkbox_e4,'Value'))
    try
        plot(h,data.error(4,:));
        hold on;
    catch
    end
end

if(get(handles.checkbox_E,'Value'))
    plot(h,sum(data.error,1));
end

%%Update plot according to switch configuration in popup_point UI
switches = get(handles.popup_point,'UserData');
contents = cellstr(get(handles.popup_point,'String'));
%Get module difference of the selected vertices
if(numel(contents)>1)
    for i=2:numel(contents)
        if(switches(i,1) == 0)
            continue;
        end
        idx = str2double(contents{i});
        %Gather the LB vector of the indexed point for all iterations
        iterations = numel(data.LB_vec);
        LB_v = zeros(iterations,3);
        for j=1:iterations
            LB_v(j,:) = data.LB_vec{j}(idx,:);
        end
        LB_mod = getLBMod(LB_v);
        plot(h,LB_mod);
        hold on;
    end
end

%Plot vertical line on the iteration
iter = get(handles.iter_slider,'Value');
% line(h,[iter iter],get(gca,'ylim'));
line([iter iter],get(gca,'ylim'),'LineWidth',2,'LineStyle',':');
drawnow;

function LB_mod = getLBMod(LB_v)
LB1 = LB_v(1,:);
LB_v = LB_v-LB1;
LB_mod = sqrt(sum(LB_v.^2, 2));



% set(h,'Tag','h_point');


% --- Executes on button press in checkbox_LBm.
function checkbox_LBm_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_LBm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_LBm
% contents = cellstr(get(handles.popup_point,'String'));

%Register the change in switch configuration to the popup_point UserData
selection = get(handles.popup_point,'Value');
switches = get(handles.popup_point,'UserData');
switches(selection,1) = hObject.Value;
set(handles.popup_point,'UserData',switches);

%Update the plot graph
UpdatePointPlot(handles);


% --- Executes during object creation, after setting all properties.
function display_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate display


% --- Executes on button press in checkbox_showedges.
function checkbox_showedges_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showedges (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showedges
changeMesh(handles);


% --- Executes on button press in checkbox_w4.
function checkbox_w4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_w4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_w4
UpdatePointPlot(handles);


% --- Executes on button press in checkbox_e4.
function checkbox_e4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_e4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_e4
UpdatePointPlot(handles);


% --- Executes on button press in checkbox_E.
function checkbox_E_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_E
UpdatePointPlot(handles);
