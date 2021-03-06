function varargout = system2DLayoutDiagram(varargin)
% SYSTEM2DLAYOUTDIAGRAM MATLAB code for system2DLayoutDiagram.fig
%      SYSTEM2DLAYOUTDIAGRAM, by itself, creates a new SYSTEM2DLAYOUTDIAGRAM or raises the existing
%      singleton*.
%
%      H = SYSTEM2DFLAYOUTDIAGRAM returns the handle to a new SYSTEM2DLAYOUTDIAGRAM or the handle to
%      the existing singleton*.
%
%      SYSTEM2DLAYOUTDIAGRAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SYSTEM2DLAYOUTDIAGRAM.M with the given input arguments.
%
%      SYSTEM2DLAYOUTDIAGRAM('Property','Value',...) creates a new SYSTEM2DLAYOUTDIAGRAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before system2DLayoutDiagram_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to system2DLayoutDiagram_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help system2DLayoutDiagram

% Last Modified by GUIDE v2.5 27-Oct-2013 09:04:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @system2DLayoutDiagram_OpeningFcn, ...
                   'gui_OutputFcn',  @system2DLayoutDiagram_OutputFcn, ...
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


% --- Executes just before system2DLayoutDiagram is made visible.
function system2DLayoutDiagram_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to system2DLayoutDiagram (see VARARGIN)

% Choose default command line output for system2DLayoutDiagram
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes system2DLayoutDiagram wait for user response (see UIRESUME)
% uiwait(handles.system2DLayoutDiagram);



% --- Outputs from this function are returned to the command line.
function varargout = system2DLayoutDiagram_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when panelSetting2DLayout is resized.
function panelSetting2DLayout_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to panelSetting2DLayout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function popSurfaceIndex_Callback(hObject, eventdata, handles)
% hObject    handle to popSurfaceIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popSurfaceIndex as text
%        str2double(get(hObject,'String')) returns contents of popSurfaceIndex as a double
currentOpticalSystem  = saveTheOpticalSystem ();
surfIndex = get(hObject,'Value');

if surfIndex == 1 % enter new coating
    % enable coating name text box
    set(handles.txtCoatingName,'String','None','Enable','On'); 
elseif surfIndex > 1
    % write the coating name on the text box and disable it
    coating = currentOpticalSystem.SurfaceArray(surfIndex).Coating;
    set(handles.txtCoatingName,'String',coating.Name,'Enable','Off');
else
    
end

% --- Executes during object creation, after setting all properties.
function popSurfaceIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popSurfaceIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtNumberOfRay_Callback(hObject, eventdata, handles)
% hObject    handle to txtNumberOfRay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtNumberOfRay as text
%        str2double(get(hObject,'String')) returns contents of txtNumberOfRay as a double


% --- Executes during object creation, after setting all properties.
function txtNumberOfRay_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtNumberOfRay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function popField_Callback(hObject, eventdata, handles)
% hObject    handle to popField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popField as text
%        str2double(get(hObject,'String')) returns contents of popField as a double
currentOpticalSystem  = saveTheOpticalSystem ();
maxFldIndexDefined =  currentOpticalSystem.NumberOfFieldPoints;
fldIndex = get(hObject,'Value')-1;
if fldIndex > maxFldIndexDefined
   set(hObject,'Value',maxFldIndexDefined+1);
   fldIndex = maxFldIndexDefined;
end

% --- Executes during object creation, after setting all properties.
function popField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txtMaxWavLen_Callback(hObject, eventdata, handles)
% hObject    handle to txtMaxWavLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtMaxWavLen as text
%        str2double(get(hObject,'String')) returns contents of txtMaxWavLen as a double


% --- Executes during object creation, after setting all properties.
function txtMaxWavLen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtMaxWavLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function btnSetting_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to btnSetting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panel2DLayout,'Visible','Off')
set(handles.panelSetting2DLayout,'Visible','On')



% --- Executes on button press in btnCancel.
function btnCancel_Callback(hObject, eventdata, handles)
% hObject    handle to btnCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function txtCoatingName_Callback(hObject, eventdata, handles)
% hObject    handle to txtCoatingName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txtCoatingName as text
%        str2double(get(hObject,'String')) returns contents of txtCoatingName as a double


% --- Executes during object creation, after setting all properties.
function txtCoatingName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txtCoatingName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function popWavelength_Callback(hObject, eventdata, handles)
% hObject    handle to popWavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of popWavelength as text
%        str2double(get(hObject,'String')) returns contents of popWavelength as a double
currentOpticalSystem  = saveTheOpticalSystem ();
maxWavelengthDefined =  currentOpticalSystem.NumberOfWavelengths;
wavIndex = get(hObject,'Value')-1;
if wavIndex > maxWavelengthDefined
   set(hObject,'Value',maxWavelengthDefined+1);
   wavIndex = maxWavelengthDefined;
end

% --- Executes during object creation, after setting all properties.
function popWavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popWavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnOk.
function btnOk_Callback(hObject, eventdata, handles)
% hObject    handle to btnOk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

update2DLayoutGraph(handles);



function rayPathMatrix = computeRayPathMatrix(currentOpticalSystem,wavInd,fldInd,pupilSamplingPoints,pupilZLocation,objThick)
rayPathMatrix = [];
Wavelen =  currentOpticalSystem.WavelengthMatrix(wavInd,1); 
FieldPoint =  [currentOpticalSystem.FieldPointMatrix(fldInd,1:2),-objThick];  

if abs(currentOpticalSystem.SurfaceArray(1).Thickness) > 10^10 %object at infinity 
         
    angX = FieldPoint(1)*pi/180;
    angY = FieldPoint(2)*pi/180;                

    %convert field angle to ray direction as in Zemax
    dz = sqrt(1/((tan (angX))^2+(tan (angY))^2+1));
    dx = dz*tan (angX);
    dy = dz*tan (angY);    
    
    commonRayDirectionCosine=[dx,dy,dz];
    [ initialRayBundlePositions ] = computeInitialRayBundlePositions(...
        commonRayDirectionCosine,pupilSamplingPoints,pupilZLocation,objThick);
    nRay=size(initialRayBundlePositions);
    nRay = nRay(1);
    initialRayBundleDirections =  repmat(commonRayDirectionCosine,nRay,1);    
else
    [ initialRayBundleDirections ] = computeInitialRayBundleDirections...
    (FieldPoint,pupilSamplingPoints);
    nRay=size(initialRayBundleDirections);
    nRay = nRay(1);
    initialRayBundlePositions =  repmat(FieldPoint,nRay,1);
end

for ri = 1:1:nRay
initialRayBundle(ri) = Ray;
initialRayBundle(ri).Position = initialRayBundlePositions(ri,:);
initialRayBundle(ri).Direction = initialRayBundleDirections(ri,:);

initialRayBundle(ri).JonesVector = [0,0;0,0];
initialRayBundle(ri).Polarized = 0;
initialRayBundle(ri).Wavelength = Wavelen;
end

%===============RAYTRACE For Bundle of Ray========================
tracedRayCounter = 0;
for rayIndex = 1:1:nRay    
    rayTraceResult = currentOpticalSystem.tracePolarizedRay(initialRayBundle(rayIndex));
%     if isempty(rayTraceResult)
%     initialRayBundle(rayIndex).Vignated = 1;
%     continue;
%     end
%     initialRayBundle(rayIndex).Vignated = 0;
    if rayTraceResult.NoIntersectionPoint || rayTraceResult.TotalInternalReflection
        continue;
    end
    tracedRayCounter = tracedRayCounter+1;
    polarizedRayTracerResult(tracedRayCounter)=rayTraceResult;
end

%display the layout
try
    nRay = length(polarizedRayTracerResult);
catch
    
end
for rayIndex = 1:1:nRay     
    try
%         if ~isempty(polarizedRayTracerResult(rayIndex).RayIntersectionPoint)
%            rayPathMatrix(:,:,rayIndex) =  polarizedRayTracerResult(rayIndex).RayIntersectionPoint;       
%         end
       rayPathMatrix(:,:,rayIndex) =  polarizedRayTracerResult(rayIndex).RayIntersectionPoint;
    catch
    end
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function btnUpdate_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to btnUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update2DLayoutGraph(handles);

function update2DLayoutGraph(handles)
global INF_OBJ_Z;
cla(handles.axes2DLayout,'reset')
currentOpticalSystem = saveTheOpticalSystem();

fldIndex = get(handles.popField);
fldIndex = fldIndex.Value - 1;

wavIndex = get(handles.popWavelength);
wavIndex = wavIndex.Value - 1;

% fldIndex = get(popField,'Tag','Value')-1;
% wavIndex = get(popWavelength,'Tag','Value')-1; 

if abs(currentOpticalSystem.SurfaceArray(1).Thickness) > 10^10 %object at infinity
    objThick = INF_OBJ_Z;
else
    objThick  = currentOpticalSystem.SurfaceArray(1).Thickness;
end

Polarized = 0;
% NumberOfRays = str2num(get(txtNumberOfRay,'Tag','String'));
NumberOfRays = get(handles.txtNumberOfRay);
NumberOfRays = str2num(NumberOfRays.String) - 1;

PupSampling = 5; % Tangential plane sampling
pupilRadius = (currentOpticalSystem.getEntrancePupilDiameter)/2;
pupilZLocation = (currentOpticalSystem.getEntrancePupilLocation);
[ pupilSamplingPoints ] = computePupilSamplingPoints(NumberOfRays,pupilZLocation,pupilRadius,PupSampling,objThick);

totalRayPathMatrix = [];
if fldIndex == 0
    for fldInd=1:1:currentOpticalSystem.NumberOfFieldPoints
        if wavIndex == 0
            for wavInd=1:1:currentOpticalSystem.NumberOfWavelengths
                rayPathMatrix = computeRayPathMatrix(currentOpticalSystem,wavInd,fldInd,pupilSamplingPoints,pupilZLocation,objThick);
                totalRayPathMatrix = cat(3,totalRayPathMatrix,rayPathMatrix);
            end
        else
            wavInd = wavIndex;
%              for fldInd=1:1:currentOpticalSystem.NumberOfFieldPoints
                rayPathMatrix = computeRayPathMatrix(currentOpticalSystem,wavInd,fldInd,pupilSamplingPoints,pupilZLocation,objThick);
                totalRayPathMatrix = cat(3,totalRayPathMatrix,rayPathMatrix);            
%              end
        end
    end
else
    fldInd = fldIndex;
    if wavIndex == 0
        for wavInd=1:1:currentOpticalSystem.NumberOfWavelengths
            rayPathMatrix = computeRayPathMatrix(currentOpticalSystem,wavInd,fldInd,pupilSamplingPoints,pupilZLocation,objThick);
            totalRayPathMatrix = cat(3,totalRayPathMatrix,rayPathMatrix);
        end
    else
        wavInd = wavIndex;
        rayPathMatrix = computeRayPathMatrix(currentOpticalSystem,wavInd,fldInd,pupilSamplingPoints,pupilZLocation,objThick);
        totalRayPathMatrix = cat(3,totalRayPathMatrix,rayPathMatrix);
    end
end

set(handles.panel2DLayout,'Visible','On')
set(handles.panelSetting2DLayout,'Visible','Off')

if PupSampling == 5 % the ray are meridional rays and so show 2D layout
     % plot2DLayout(currentOpticalSystem,totalRayPathMatrix,handles.axes2DLayout);
     currentOpticalSystem.plot2DLayout(handles.axes2DLayout,totalRayPathMatrix);
end

% --------------------------------------------------------------------
function btnGraph_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to btnGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panel2DLayout,'Visible','On')
set(handles.panelSetting2DLayout,'Visible','Off')
