function varargout = princeGUI(varargin)
% PRINCEGUI MATLAB code for princeGUI.fig
%      PRINCEGUI, by itself, creates a new PRINCEGUI or raises the existing
%      singleton*.
%
%      H = PRINCEGUI returns the handle to a new PRINCEGUI or the handle to
%      the existing singleton*.
%
%      PRINCEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PRINCEGUI.M with the given input arguments.
%
%      PRINCEGUI('Property','Value',...) creates a new PRINCEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before princeGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to princeGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help princeGUI

% Last Modified by GUIDE v2.5 01-Nov-2017 16:12:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @princeGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @princeGUI_OutputFcn, ...
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


% --- Executes just before princeGUI is made visible.
function princeGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to princeGUI (see VARARGIN)

% Choose default command line output for princeGUI
%handles.output = hObject;
% handles.output.skipcomparison = 0;
% handles.output.skipalignment = 0;
% handles.output.skipinteractions = 0;
% handles.output.skipcomplexes = 0;
% handles.output.skipcomplexes = 0;
% handles.output.fastgaussbuild = 0;
% handles.output.fastcomparison = 0;
% handles.output.skipgaussbuild = 0;

% Get handle for GUI figure
handles.GUIhandle = hObject;

% Change appearance
set(handles.figure1,'Name', 'Enter experiment details');
% set(handles.text20,'fontsize',10)
% set(handles.text21,'fontsize',10)
% set(handles.text22,'fontsize',10)
% set(handles.text23,'fontsize',10)
% set(handles.text24,'fontsize',10)
% set(handles.radiobutton1,'fontsize',9)
% set(handles.radiobutton2,'fontsize',9)
% set(handles.radiobutton3,'fontsize',9)
% set(handles.radiobutton4,'fontsize',9)
% set(handles.radiobutton5,'fontsize',9)
% set(handles.radiobutton6,'fontsize',9)
% set(handles.radiobutton7,'fontsize',9)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes princeGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = princeGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = '';

if not(isempty(handles))
    handles.output.majorproteingroupsfile = handles.edit9.String;
    handles.output.corumfile = handles.edit10.String;
    handles.output.treatmentcondition = ['condition' handles.edit12.String];
    if not(isempty(handles.edit11.String))
        handles.output.notreatmentcondition = ['condition' handles.edit11.String];
    else
        handles.output.notreatmentcondition = handles.edit11.String;
    end
    handles.output.desiredPrecision = handles.edit13.Value;
    handles.output.skipgaussbuild = handles.radiobutton1.Value;
    handles.output.skipalignment = handles.radiobutton2.Value;
    handles.output.skipcomparison = handles.radiobutton3.Value;
    handles.output.skipinteractions = handles.radiobutton4.Value;
    handles.output.skipcomplexes = handles.radiobutton5.Value;
    handles.output.fastgaussbuild = handles.radiobutton6.Value;
    handles.output.fastcomparison = handles.radiobutton7.Value;
    
    % Get default command line output from handles structure
    varargout{1} = handles.output;
    
    delete(hObject);
end




function handles = edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = get(hObject,'string');
I = ismember(tmp,['1' '2' '3' '4' '5' '6' '7' '8' '9' '0']);
if sum(I)==0
    warndlg('Non-treatment condition must be an integer.');
end
tmp = str2double(tmp(I));
set(hObject,'string',tmp,'value',tmp);
%handles.output.notreatmentcondition = ['condition' num2str(tmp) '.csv'];
% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = get(hObject,'string');
I = ismember(tmp,['1' '2' '3' '4' '5' '6' '7' '8' '9' '0']);
if sum(I)==0
    warndlg('Treatment condition must be an integer.');
end
tmp = str2double(tmp(I));
set(hObject,'string',tmp,'value',tmp);
%handles.output.treatmentcondition = ['condition' num2str(tmp) '.csv'];
% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = get(hObject,'string');
I = ismember(tmp,['1' '2' '3' '4' '5' '6' '7' '8' '9' '0' '.']);
tmp = str2double(tmp(I));
if isempty(tmp)
    warndlg('Desired precision must be a number.');
end
if tmp>1
    tmp = tmp/100;
end
if tmp>1
    warndlg('Desired precision too big. Enter a number between 0 and 100.');
    tmp = 0.5;
end
set(hObject,'string',[num2str(tmp*100) '%'],'value',tmp);
%handles.output.desiredPrecision = tmp;
% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = get(hObject,'string');

% trim leading and trailing white space
tmp = strtrim(tmp);

set(hObject,'string',tmp);
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function handles = edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = get(hObject,'string');

% trim leading and trailing white space
tmp = strtrim(tmp);

set(hObject,'string',tmp);
% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function handles = pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[mpgfile, mpgpath] = uigetfile('*.csv');

% set to 0 if user presses cancel. delete this.
if mpgfile==0
    mpgfile = '';
end
if mpgpath==0
    mpgpath = '';
end

handles.majorproteingroupsfile = [mpgpath mpgfile];
set(handles.edit9,'string',mpgfile);


% --- Executes on button press in pushbutton6.
function handles = pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[corfile, corpath] = uigetfile('*');

% set to 0 if user presses cancel. delete this.
if corfile==0
    corfile = '';
end
if corpath==0
    corpath = '';
end

handles.corumfile = [corpath corfile];
set(handles.edit10,'string',corfile);


% --- Executes on button press in radiobutton6.
function handles = radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output.fastgaussbuild = get(hObject,'value');
% Hint: get(hObject,'Value') returns toggle state of radiobutton6


% --- Executes on button press in radiobutton7.
function handles = radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output.fastcomparison = get(hObject,'value');
% Hint: get(hObject,'Value') returns toggle state of radiobutton7


% --- Executes on button press in radiobutton1.
function handles = radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output.skipgaussbuild = get(hObject,'value');
% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function handles = radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output.skipcomparison = get(hObject,'value');
% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function handles = radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output.skipalignment = get(hObject,'value');
% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in radiobutton4.
function handles = radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.output.skipinteractions = get(hObject,'value');
% Hint: get(hObject,'Value') returns toggle state of radiobutton4


% --- Executes on button press in radiobutton5.
function handles = radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get(hObject,'value')
handles.output.skipcomplexes = get(hObject,'value');
% Hint: get(hObject,'Value') returns toggle state of radiobutton5


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% is major protein groups field empty?
Iempty = zeros(1,5);
warnString = cell(1,5);
%defaultString = cell(1,5);
tmp = handles.edit9.String;
tmp = strtrim(tmp);
if isempty(tmp)
    Iempty(1) = 1;
    warnString{1} = [handles.text20.String ': no file\n'];
    %defaultString{1} = 'No file';
    handles.edit9.String = '';
end

% is corumfile field empty?
tmp = handles.edit10.String;
tmp = strtrim(tmp);
if isempty(tmp)
    Iempty(2) = 1;
    warnString{2} = [handles.text21.String ': allComplexes.txt\n'];
    %defaultString{2} = 'allComplexes.txt';
    handles.edit10.String = 'allComplexes.txt';
end

% is treatment condition field empty?
tmp = handles.edit12.String;
tmp = strtrim(tmp);
if isempty(tmp)
    Iempty(3) = 1;
    warnString{3} = 'Treatment condition: condition1\n';
    %defaultString{3} = 'condition1';
    handles.edit12.String = '1';
    handles.edit12.Value = 1;
end

% is no-treatment condition field empty?
tmp = handles.edit11.String;
tmp = strtrim(tmp);
if isempty(tmp)
    Iempty(4) = 1;
    warnString{4} = 'Non-treatment condition: condition2 (if condition2.csv exists)\n';
    %defaultString{4} = 'condition2';
    handles.edit11.String = '2';
    handles.edit11.Value = 2;
end

% is desired precision field empty?
tmp = handles.edit13.String;
tmp = strtrim(tmp);
if isempty(tmp)
    Iempty(5) = 1;
    warnString{5} = 'Interaction precision: 50%%\n';
    %defaultString{5} = '50%';
    handles.edit13.String = '50%';
    handles.edit13.Value = 0.5;
end

if sum(Iempty)>0
    warnString(cellfun('isempty', warnString)) = [];
    %defaultString(cellfun('isempty', defaultString)) = [];
    hh = warndlg(sprintf(strjoin(['Some files and/or parameters were not set!\n\n' ...
        'Attempting to use the following default values (this might cause issues!):\n\n' ...
        warnString])));
    uiwait(hh);
end

uiresume;
