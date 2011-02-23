function progressbar(x,whichbar, title, position)
%WAITBAR Display wait bar.
%   H = WAITBAR(X,'message', property, value, property, value, ...)
%   creates and displays a waitbar of fractional length X.  The
%   handle to the waitbar figure is returned in H.
%   X should be between 0 and 1.  Optional arguments property and
%   value allow to set corresponding waitbar figure properties.
%   Property can also be an action keyword 'CreateCancelBtn', in
%   which case a cancel button will be added to the figure, and
%   the passed value string will be executed upon clicking on the
%   cancel button or the close figure button.
%
%   WAITBAR(X) will set the length of the bar in the most recently
%   created waitbar window to the fractional length X.
%
%   WAITBAR(X,H) will set the length of the bar in waitbar H
%   to the fractional length X.
%
%   WAITBAR(X,H,'message') will update the message text in
%   the waitbar figure, in addition to setting the fractional
%   length to X.
%
%   WAITBAR is typically used inside a FOR loop that performs a
%   lengthy computation.  
%
%   Example:
%       h = waitbar(0,'Please wait...');
%       for i=1:1000,
%           % computation here %
%           waitbar(i/1000,h)
%       end
%
%   See also DIALOG, MSGBOX.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1.23.4.13 $  $Date: 2007/09/27 22:48:16 $

if nargin>=2
    if isnumeric(whichbar)
        hfig = whichbar;

        p = findobj(hfig,'Type','patch','-and','Tag','ProgressBar');
        l = findobj(hfig,'Type','line','-and','Tag','ProgressBar');
        if(isempty(p) || isempty(l))
          type=2;
        else
  
          if(x<0)
            delete(findobj(hfig,'Tag','ProgressBar'));
            return;
          end
          
          type=1;
        end
    else
      error('MATLAB:progressbar:InvalidArguments', 'Input arguments not valid.');
    end
else
    error('MATLAB:progressbar:InvalidArguments', 'Input arguments not valid.');
end

x = max(0,min(100*x,100));
try 
switch type
    case 1,  % waitbar(x)    update

        xpatch = [0 x x 0];
        set(p,'XData',xpatch)
        xline = get(l,'XData');
        set(l,'XData',xline);

        if nargin>2,
            % Update waitbar title:
            hAxes = findobj(hfig,'type','axes','-and','Tag','ProgressBar');
            hTitle = get(hAxes,'title');
            set(hTitle,'string',title);
        end
    case 2,  % waitbar(x,name)  initialize
       
        if nargin<=3
          position=[.025 .025 .95 .03];
        end

        h = axes('XLim',[0 100],...
            'YLim',[0 1],...
            'Position',position,...
            'XTickMode','manual',...
            'YTickMode','manual',...
            'XTick',[],...
            'YTick',[],...
            'XTickLabelMode','manual',...
            'XTickLabel',[],...
            'YTickLabelMode','manual',...
            'YTickLabel',[], ...
            'Tag','ProgressBar', ...
            'Parent', hfig);

        tHandle=get(h,'title');
        set(tHandle,'String', title);

        xpatch = [0 x x 0];
        ypatch = [0 0 1 1];
        xline = [100 0 0 100 100];
        yline = [0 0 1 1 0];
        
        patch(xpatch,ypatch,'r','EdgeColor','r','EraseMode','normal','Tag','ProgressBar','Parent',h);
        l = line(xline,yline,'EraseMode','none','Tag','ProgressBar','Parent',h);
        set(l,'Color',get(gca,'XColor'));
    
end  % case

catch
       close(findobj(hfig,'Tag','ProgressBar'));
       error('MATLAB:waitbar:InvalidArguments','Improper arguments for progressbar');
end
drawnow;
