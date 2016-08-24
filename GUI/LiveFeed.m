function f=LiveFeed
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

%  Create and then hide the UI as it is being constructed.
f = figure('Visible','off','Position',[360,500,1000,900]);

% Construct the components.
hbrowse    = uicontrol('Parent',f,'Style','pushbutton',...
             'String','Browse','Position',[300,800,70,25],...
             'Callback',@browsebutton_Callback);
hstart    = uicontrol('Parent',f,'Style','pushbutton',...
             'String','Start','Position',[400,800,70,25],...
             'Callback',@startbutton_Callback,...
             'Tag','StartButton');
hstop = uicontrol('Parent',f,'Style','pushbutton',...
             'String','Stop','Position',[500,800,70,25],...
             'Callback',@stopbutton_Callback,...
             'Tag','StopButton');
         
hax1 = axes('Parent',f,'Units','pixels','Position',[100,100,200,700]);
hax2 = axes('Parent',f,'Units','pixels','Position',[400,100,200,700]);
hax3 = axes('Parent',f,'Units','pixels','Position',[700,100,200,700]);
hax1.Box='on';
hax2.Box='on';
hax3.Box='on';

% Initialize the UI.
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
hax1.Units = 'normalized';
hax2.Units = 'normalized';
hax3.Units = 'normalized';
hstart.Units = 'normalized';
hstop.Units = 'normalized';
hbrowse.Units = 'normalized';

% Assign the a name to appear in the window title.
f.Name = 'PAF GUI';
% Move the window to the center of the screen.
movegui(f,'center')
% Make the window visible.
f.Visible = 'on';

%%
% Push button callbacks.

    function startbutton_Callback(hObject,eventdata)
        [PAobj, ~, ~] = LoadFiles(pathname);
        while get(hObject,'Value')
            sinc_data=abs(sinc_data.^1.1);
            surf(hax1,sinc_data);
            drawnow;
            pause(0.1);
        end
    end

    function stopbutton_Callback(hObject,eventdata)
        set(hstart,'Value',0)        
    end

    function browsebutton_Callback(hObject,eventdata)
        % Display contour plot of the currently selected data.
        pathname = uigetdir();
    end
end