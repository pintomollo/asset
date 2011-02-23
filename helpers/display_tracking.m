function handles = display_tracking(trackings, varargin)

  [h, indx, color, name, width, new_size, center, orientation] = parse_input(varargin{:}); 

  if (isempty(h) || ~ishandle(h))
    h = figure;
    h = axes;
  end

  handles = display_really(trackings, h, indx, color, name, width, new_size, center, orientation);

  return;
end

function handles = display_really(trackings, h, indx, color, name, width, new_size, center, orientation)

  handles = [];

  if (isstruct(trackings))
    if (isfield(trackings, 'splines'))

      if (isfield(trackings, 'fname'))
        [tokens,junk]=regexp(trackings.fname,'(.+[-_])?([^-_\.]+)(\..+)','tokens');
        name = ['...-' tokens{1}{2} ];
      end

      if (~isempty(center))
        spline = realign(trackings.splines(indx(1), indx(2)), new_size, center, orientation);
      else
        spline = trackings.splines(indx(1), indx(2));
      end
      %pts = fnplt(trackings.splines(indx(1), indx(2)));
      %handles = display_tracking(h, pts, indx, color, name, width);
      handles = myplot(spline, 'Color', color, 'LineWidth', width, 'Parent', h, 'DisplayName', name);

    elseif (isfield(trackings, 'breaks'))

      if (~isempty(center))
        spline = realign(trackings, new_size, center, orientation);
      else
        spline = trackings;
      end
      %pts = fnplt(trackings);
      %handles = display_tracking(h, pts, indx, color, name, width);
      handles = myplot(spline, 'Color', color, 'LineWidth', width, 'Parent', h, 'DisplayName', name);

    else

      if (isfield(trackings, 'child'))
        nchild = length(trackings.child);
      else
        nchild = 0;
      end
      if (isfield(trackings, 'files'))
        nfiles = length(trackings.files);
      else
        nfiles = 0;
      end
      if (isfield(trackings, 'mean'))
        nmean = 1;
      else
        nmean = 0;
      end
      if (isfield(trackings, 'name'))
        name = trackings.name;

        if (isempty(name))
          name = 'Mean';
        end
      end

      nlines = nchild + nfiles + nmean;

      if (size(color, 1) < nlines)
        if (~isempty(color))
          tmp = color(end,:);
          tmp = rgb2hsv(tmp);
          incr = tmp(2) / (nlines+1);

          tmp = repmat(tmp, nlines, 1);
          tmp(:,2) = tmp(:,2) - [nlines:-1:1]' * incr;
          color = hsv2rgb(tmp);
        else
          color = hsv(6* ceil(nlines / 6));
        end
      end

      for i=1:nlines
        if (i <= nchild)
          htmp = display_really(trackings.child(i), h, indx, color(i,:), name, width, new_size, center, orientation);
        elseif (i <= nfiles + nchild)
          htmp = display_really(trackings.files(i-nchild), h, indx, color(i,:), name, width, new_size, center, orientation);
        else
          htmp = display_really(trackings.mean(indx(1), indx(2)), h, indx, color(i,:), name, 2, new_size, center, orientation);
          if (isfield(trackings, 'errors'))
            set(htmp, 'UserData', squeeze(trackings.errors(indx(1),indx(2),:)));
          end
        end

        handles = [handles htmp];
      end

      if (nfiles > 0)
        linkprop(handles(nchild+1:end), 'Visible');
      end
    end
  else
    trackings

    if (isempty(color))
      color = hsv2rgb([rand(1) 1 1]);
    end

    if (size(trackings, 2) > 2)
      trackings = trackings.';
    end

    handles = line(trackings(:,1), trackings(:,2), 'Color', color, 'LineWidth', width, 'Parent', h, 'DisplayName', name);
  end

  return;
end

function [h, indx, color, name, width, new_size, center, orientation] = parse_input(varargin)

  h = [];
  indx = [];
  color = zeros(0,3);
  name = ' ';
  width = 1;
  new_size = [];
  center = [];
  orientation = 0;

  if (nargin > 0)
    for i=1:length(varargin)
      var = varargin{i};
      type = get_type(var);
      switch type
        case 'num'
          if (numel(var) == 1)
            if (ishandle(var))
              h = var;
            elseif (isfloat(var) | (var == 0))
              orientation = var;
            elseif (isempty(indx))
              indx = var;
            else
              width = var;
            end
          elseif (numel(var) == 2)
            if (isempty(new_size))
              new_size = var;
            elseif (isempty(center))
              center = var;
            end
          elseif (size(var, 2) == 3)
            color = var;
          end
        case 'char'
          name = var;
      end
    end
  end

  if (isempty(indx))
    indx = [1 1];
  end

  return;
end
