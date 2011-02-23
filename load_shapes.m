function [shapes, groups] = load_shapes(fname)

  if (isempty(fname))
    [fname,dirpath] = uigetfile('*.shapes', 'Load a .SHAPES file');
    fname = [dirpath fname];
  elseif (isstruct(fname) & isfield(fname, 'fname'))
    fname = fname.fname;
  end

  groups = [];
  frame_shift = 0;

  fid = fopen(fname,'rt');
  if (fid<0)
    fname = [pwd fname];
    fid = fopen(fname,'rt');
    if (fid<0)
      shapes = {};

      return;
    end
  end

  nshapes = 0;

  line = fgetl(fid);
  status = 'header';
  while ischar(line)
    
    switch status
      case 'header'
        shift = regexp(line,'FrameShift=(\d+)','tokens');
        if (~isempty(shift))
          frame_shift = str2num(shift{1}{1});
        end
        find_groups = regexp(line,'Groups=([\w,]+)','tokens');
        if (~isempty(find_groups))
          groups = find_groups{1};
          groups = regexp(groups,',','split');
          groups = groups{1};
          groups = groups(1,1:end-1);
        else
          if(length(line)==0)
            status = 'find';
          end
        end
      case 'find'
        shift = regexp(line,'FrameShift=(\d+)','tokens');
        if (~isempty(shift))
          frame_shift = str2num(shift{1}{1});
        end
        type = regexp(line,'type=(\w+)ROI','tokens');
        if (~isempty(type))
          group_indx = 0;
          slice_indx = 0;
          status = 'locate';
        end
      case 'locate'
        group = regexp(line,'^group=(\w+)','tokens');
        if(~isempty(group))
          group_indx = identify_group(group{1},groups);
        else
          slice = regexp(line,'^in slice=(\d+)','tokens');
          if(~isempty(slice))
            slice_indx = str2num(char(slice{1})) + 1;
          end
        end
        
        if (group_indx>0 & slice_indx>0)
          status = 'points';
          x = [];
          y = [];
        end
      case 'points'
        if (length(line)==0)
          status = 'find';
          good_indx = ~(any(isnan([x y]),2));

          x = x(good_indx);
          y = y(good_indx);
          [x, y] = poly2cw(x, y);

          shapes{group_indx,slice_indx + frame_shift} = [x y];
        else
          if (strncmp(line,'x=',2))
            x = [x; str2num(line(3:end))];
          elseif (strncmp(line,'y=',2))
            y = [y; str2num(line(3:end))];
          end
        end
    end

    line = fgetl(fid);
  end

  if (strncmp(status, 'points', 6))
    good_indx = ~(any(isnan([x y]),2));

    x = x(good_indx);
    y = y(good_indx);
    [x, y] = poly2cw(x, y);

    shapes{group_indx,slice_indx + frame_shift} = [x y];
  end

  return;
end

function indx = identify_group(group,groups)

  n = length(group);
  indx = 0;

  for i=1:length(groups)
    if (strncmp(char(groups(i)),group,n))
      indx = i;
      return;
    end
  end

  return;
end
