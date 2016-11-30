function f=LiveFeed
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

%  Create and then hide the UI as it is being constructed.
f = figure('Visible','off','Position',[360,500,1000,680]);

hax1 = axes('Parent',f,'Units','pixels','Position',[100,100,200,500]);
hax2 = axes('Parent',f,'Units','pixels','Position',[400,100,200,500]);
hax3 = axes('Parent',f,'Units','pixels','Position',[700,100,200,500]);
hax1.Box='on';
hax2.Box='on';
hax3.Box='on';
set(hax1,'tag','hax1');
% Initialize the UI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
hax1.Units = 'normalized';
hax2.Units = 'normalized';
hax3.Units = 'normalized';
hrecon.Units = 'normalized';
hstop.Units = 'normalized';
hbrowse.Units = 'normalized';


% Construct the components.
% Control Buttons
hbrowse    = uicontrol('Parent',f,'Style','pushbutton',...
    'String','Browse','Position',[5,580,70,25],...
    'Callback',@browsebutton_Callback);
hrecon    = uicontrol('Parent',f,'Style','pushbutton',...
    'String','Reconstruct','Position',[5,550,70,25],...
    'Callback',@reconbutton_Callback,...
    'Tag','StartButton');
hPAF = uicontrol('Parent',f,'Style','pushbutton',...
    'String','PAF','Position',[5,520,70,25],...
    'Callback',@PAFbutton_Callback,...
    'Tag','StopButton');

%Flow Settings
hspeedtext = uicontrol('Parent',f,'Style','text',...
    'String','Flow Speed (ml/min)','Position',[405,650,70,25]);
hspeed = uicontrol('Parent',f,'Style','edit',...
    'String','0.05','Position',[405,620,70,25],...
    'Callback',@speed_Callback);
hangletext = uicontrol('Parent',f,'Style','text',...
    'String','Angle (deg)','Position',[505,650,70,25]);
hangle = uicontrol('Parent',f,'Style','edit',...
    'String','60','Position',[505,620,70,25],...
    'Callback',@angle_Callback);
hPRFtext = uicontrol('Parent',f,'Style','text',...
    'String','Pulse Separation (ms)','Position',[595,650,80,25]);
hPRF = uicontrol('Parent',f,'Style','edit',...
    'String','20','Position',[605,620,70,25],...
    'Callback',@PRF_Callback);

%Other
hcontrasttext = uicontrol('Parent',f,'Style','text',...
    'String','Contrast','Position',[100,625,200,25]);
hcontrast1 = uicontrol('Parent',f,'Style','slide',...
    'min',1,'value',1,'max',100,'Position',[100,610,200,25],...
    'Callback',{@contrast_Callback,hax1,hax2});
hOrdertext = uicontrol('Parent',f,'Style','text',...
    'String','Shift Initial Frame','Position',[90,650,100,25]);
hOrderToggle = uicontrol('Parent',f,'Style','togglebutton',...
    'min',0,'max',1,'Position',[185,659,15,15],...
    'Callback',@shiftOrder_Callback);



% Assign the a name to appear in the window title.
f.Name = 'PAF GUI';
% Move the window to the center of the screen.
movegui(f,'center')
% Make the window visible.
f.Visible = 'on';
speed_Callback(hspeed);
angle_Callback(hangle);
PRF_Callback(hPRF);
shiftOrder_Callback


%%
% Push button callbacks.

    function reconbutton_Callback(hObject,eventdata)
        rfObj=evalin('base','rfObj');
        reconObj=ReconFcn(rfObj);
                reconObj.wallfilter;
%         reconObj.highpass(5,2);
        reconObj.plot(f,hax2);
        assignin('base','reconObj',reconObj);

    end

    function PAFbutton_Callback(hObject,eventdata)
        handles = guidata(hObject);
        
        reconObj=evalin('base','reconObj');
        reconObj.wallfilter;
        flw.theta=handles.theta;
        flw.tube_diameter=580E-3;
        flw.PRF=handles.PRF;%evalin('base','PRF');        
        flw.rate=handles.rate;
%         xcObj=evalin('base','xcObj');
        shift_order=handles.shift_order;
        xcObj=EnsembleCorrelateFcn(reconObj, flw,...
            'ShiftOrder', shift_order,'Truncate', true);
        assignin('base','xcObj',xcObj);
        xcObj.plot(f,hax3)
        %         reconObj.highpass(10,2);

    end

    function browsebutton_Callback(hObject,eventdata)
        % Display contour plot of the currently selected data.
        hObject.UserData = uigetdir();
        rfObj = UltrasonixImport('Path',hObject.UserData);
        rfObj.detrend;
        rfObj.del=0;
        rfObj.plot(f,hax1);
        assignin('base','rfObj',rfObj); 
    end

    function speed_Callback(hObject,eventdata)
        handles = guidata(hObject);
        
        handles.rate=str2double(get(hObject,'string'));
        
        guidata(hObject, handles)
    end
    function angle_Callback(hObject,eventdata)
        handles = guidata(hObject);
        
        handles.theta=str2double(get(hObject,'string'));
        
        guidata(hObject, handles)
    end
    function PRF_Callback(hObject,eventdata,handles)
        handles = guidata(hObject);
        
        PRF=1000/(str2double(get(hObject,'string')));
        handles.PRF=PRF;
        
        guidata(hObject, handles)
    end
    function contrast_Callback(varargin)
        % Callback for the slider.
        [h ax1 ax2] = deal(varargin{[1;3;4]});  % Get the calling handle and structure.
        caxis(ax1,[-100+get(h,'val') 100-get(h,'val')]);
        caxis(ax2,[-20+get(h,'val')/5.0 20-get(h,'val')/5.0]);
    end
    function shiftOrder_Callback(hObject,eventdata)
        handles = guidata(hObject);
        
        handles.shift_order=get(hObject,'val');
        
        guidata(hObject, handles)
    end
end