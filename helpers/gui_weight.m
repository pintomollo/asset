function varargout = gui_weight(varargin)
% GUI_WEIGHT_EXPORT M-file for gui_weight_export.fig
%      GUI_WEIGHT_EXPORT, by itself, creates a new GUI_WEIGHT_EXPORT or raises the existing
%      singleton*.
%
%      H = GUI_WEIGHT_EXPORT returns the handle to a new GUI_WEIGHT_EXPORT or the handle to
%      the existing singleton*.
%
%      GUI_WEIGHT_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_WEIGHT_EXPORT.M with the given input arguments.
%
%      GUI_WEIGHT_EXPORT('Property','Value',...) creates a new GUI_WEIGHT_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_weight_export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_weight_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

if (nargin == 0)
  fname = {'1056-14-inc_med-060711_1_'};
else
  fname = varargin(1);

  if (any(fname{1} == '*'))
    files = dir(fname{1});
    fname = cell(length(files), 1);

    for i=1:length(fname)
      fname{i} = files(i).name;
    end
  end
end
  
  kymo = load(fname{1});
  [domain, ruffles, theta] = gather_quantification(kymo.mymovie, kymo.opts);
  domain = imnorm(domain);
  kymo.opts = load_parameters(kymo.opts, 'domain_center.txt');
  kymo.opts.quantification.weights.filt = ruffles;
  weights = kymo.opts.quantification.weights;

  curr_indx = 1;

  haxes = gui_weight_export_LayoutFcn;
  update_display;

function update_display

  img = weight_symmetry(domain, weights);
  imagesc(img, 'Parent', haxes);
  display(num2str([weights.alpha weights.beta weights.gamma weights.delta]));

end

function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

weights.alpha = get(hObject,'Value');
update_display;

end

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
weights.beta = get(hObject,'Value');
update_display;

end

% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

weights.gamma = get(hObject,'Value');
update_display;

end

% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

weights.delta = get(hObject,'Value');
update_display;

end
% --- Executes on slider movement.
function choice_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

indx = get(hObject,'Value');

if (indx ~= curr_indx)
  strings = get(hObject,'String');
  kymo = load(strings{indx});
  [domain, ruffles, theta] = gather_quantification(kymo.mymovie, kymo.opts);
  domain = imnorm(domain);
  weights.filt = ruffles;

  curr_indx = indx;
  update_display;
end

end

% --- Creates and returns a handle to the GUI figure. 
function h2 = gui_weight_export_LayoutFcn
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

h1 = figure(...
'MenuBar','none',...
'Name','gui_weight',...
'NumberTitle','off',...
'Units','characters',...
'Position',[103.714285714286 29.1428571428571 112.142857142857 32.2857142857143],...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',[],...
'Visible','on');

h2 = axes(...
'Parent',h1,...
'Units','characters',...
'Position',[5.14285714285714 3 67.8571428571428 27.0714285714286],...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode',get(0,'defaultaxesCameraPositionMode'),...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'Tag','axes1');

%'Callback',@(hObject,eventdata)gui_weight_export('slider1_Callback',hObject,eventdata,guidata(hObject)),...

h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@slider1_Callback,...
'Position',[78.4285714285714 27.6428571428571 30.4285714285714 1.78571428571429],...
'String',{  'Slider' },...
'Style','slider',...
'Min', 0, ...
'Max', 1, ...
'Value', weights.alpha, ...
'Tag','slider1');

h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Position',[78.4285714285714 22.5 30.4285714285714 1.78571428571429],...
'Callback',@slider2_Callback,...
'String',{  'Slider' },...
'Min', 0, ...
'Max', 1, ...
'Style','slider',...
'Value', weights.beta, ...
'Tag','slider2');

h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@slider3_Callback,...
'Position',[78.4285714285714 17.2857142857143 30.4285714285714 1.78571428571429],...
'String',{  'Slider' },...
'Min', 0, ...
'Max', 1, ...
'Style','slider',...
'Value', weights.gamma, ...
'Tag','slider3');

h10 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback',@slider4_Callback,...
'Position',[78.4285714285714 12.5714285714286 30.4285714285714 1.78571428571429],...
'String',{  'Slider' },...
'Min', 0, ...
'Max', 1, ...
'Style','slider',...
'Value', weights.delta, ...
'Tag','slider4');

h11 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Position',[78.4285714285714 26.7142857142857 7.42857142857143 0.928571428571428],...
'String','alpha',...
'Style','text',...
'Tag','text1');

h12 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Position',[78.4285714285714 21.7857142857143 7.42857142857143 0.928571428571428],...
'String','beta',...
'Style','text',...
'Tag','text2');

h13 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Position',[78.4285714285714 16.7142857142857 7.42857142857143 0.928571428571428],...
'String','gamma',...
'Style','text',...
'Tag','text3');

h14 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Position',[78.4285714285714 11.4285714285714 7.42857142857143 0.928571428571428],...
'String','delta',...
'Style','text',...
'Tag','text4');

h15 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Position',[78.4285714285714 6.4285714285714 35 0.928571428571428],...
'Style','popupmenu',...
'Callback',@choice_Callback,...
'String', fname, ...
'Tag','choice');

end

end