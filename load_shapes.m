function [shapes, groups] = load_shapes(fname)
% LOAD_SHAPES loads .shapes files and outputs the coordinates of the tracking in cells.
%
%   SHAPES = LOAD_SHAPES(FNAME) loads the content of the .shapes FNAME file into the GxN
%   cell SHAPES. G being the number of "groups" (i.e. the different types of trackigns),
%   while N is the number of frames.
%
%   [SHAPES, GROUPS] = LOAD_SHAPES(...) also returns the strings identifying the GROUPS.
%
%   [...] = LOAD_SHAPES(MYSTRUCT) loads the content of the .shapes file specified in
%   the "fname" field.
%
% SHAPES FILES are text file used to store segmentations of images created by Albert 
%   Cardona for his ImageJ plugin "A 3D editing plugin" (http://www.mcdb.ucla.edu/
%   research/hartenstein/software/imagej/). They have a very simple structure, namely 
%   a header that defines the content of the file, followed by the x and y coordinates
%   of each segmentation, separated by empty lines. More formally (pseudo-context free
%   grammar, parenthesized lines are optional):
%
%   SHAPES: HEADER
%           //
%           BODY
%
%   HEADER: Groups=NAMES
%           (FrameShift=INDEX)
%           MORE
%
%   BODY:   SEGMENT
%           //
%           BODY
%
%   SEGMENT:(FrameShift=INDEX)
%           type=TYPE
%           group=NAME
%           in slice=INDEX
%           POINTS
%
%   POINTS: x=POSITION
%           y=POSITION
%           POINTS
%
% where NAMES is a coma-separated list of the NAME of the different groups of segmentations,
% INDEX is the integer value referring to the segmented frame, TYPE is usually "bezier"
% but it can take different values depending on the segmentation strategy, POSITION is
% the real-value coordinate of the point and MORE represent any non-empty line of text
% which is not considered by this loading function.
%
% Gonczy & Naef labs, EPFL
% Simon Blanchoud
% 12.04.2011

  % Handle different input values
  if (nargin < 1 || isempty(fname))
    % Load a file if none input
    [fname,dirpath] = uigetfile('*.shapes', 'Load a .SHAPES file');
    fname = fullfile(dirpath, fname);

  % Or if it's a structure, try to extract a file from it
  elseif (isstruct(fname) & isfield(fname, 'fname'))
    fname = fname.fname;
  end

  % Initiallize some parameters
  groups = [];
  shapes = {};
  frame_shift = 0;

  % Open the file, line by line
  fid = fopen(fname,'rt');

  % If it does not work, try the absolute path instead
  if (fid<0)
    fname = absolutepath(fname);
    fid = fopen(fname,'rt');

    % If it still does not work, we give up
    if (fid<0)
      shapes = {};
      warning(['File ' fname ' does not exist, skipping it.'])

      return;
    end
  end

  % Read the first line, initiallize the state of the parser
  line = fgetl(fid);
  status = 'header';

  % Loop until we have read the whole file
  while ischar(line)
    % Switch between the different states of the parser
    switch status

      % We start by the header
      case 'header'

        % If there is a shift, we store it
        shift = regexp(line,'FrameShift=(\d+)','tokens');
        if (~isempty(shift))
          frame_shift = str2double(shift{1}{1});
        end

        % Look for the definition of the groups
        find_groups = regexp(line,'Groups=([\w,]+)','tokens');
        if (~isempty(find_groups))

          % If we found it, we separate the group names using the commas 
          groups = find_groups{1};
          groups = regexp(groups,',','split');
          groups = groups{1};

          % Remove the last one if it is empty due to an additional comma
          if (isempty(groups{1, end}))
            groups = groups(1,1:end-1);
          end
        else

          % When we find an empty line, the header is finished, we now look for segmetnations
          if(length(line)==0)
            status = 'find';
          end
        end

      % We now look for the "header" of the segmentation
      case 'find'
        % Maybe there is another shift
        shift = regexp(line,'FrameShift=(\d+)','tokens');
        if (~isempty(shift))
          frame_shift = str2double(shift{1}{1});
        end

        % We look for the definition of the segmentation type
        type = regexp(line,'type=(\w+)','tokens');
        if (~isempty(type))
          % When we found it, we can look for the points
          group_indx = 0;
          slice_indx = 0;
          status = 'locate';
        end

      % Here we look for the appartenance of the group of points defining a segmentation
      case 'locate'
        % We need to knoe the group they belong to
        group = regexp(line,'^group=(\w+)','tokens');
        if(~isempty(group))
          [junk, group_indx] = ismember(group{1}, groups);
        else

          % And in which slice they are
          slice = regexp(line,'^in slice=(\d+)','tokens');
          if(~isempty(slice))
            slice_indx = str2double(char(slice{1})) + 1;
          end
        end
        
        % As soon as we have both information, we can start looking for the points
        if (group_indx>0 & slice_indx>0)
          status = 'points';
          x = [];
          y = [];
        end

      % Finally we load the points individually
      case 'points'

        % An empty line finishes the segmentation
        if (length(line)==0)
          % We look for a new segmentation
          status = 'find';
          good_indx = ~(any(isnan([x y]),2));

          % We remove the bad points
          x = x(good_indx);
          y = y(good_indx);
          [x, y] = poly2cw(x, y);

          % And store them
          shapes{group_indx,slice_indx + frame_shift} = [x y];
        else

          % We look for x and y points
          if (strncmp(line,'x=',2))
            x = [x; str2double(line(3:end))];
          elseif (strncmp(line,'y=',2))
            y = [y; str2double(line(3:end))];
          end
        end
    end

    % Read a new line
    line = fgetl(fid);
  end

  % If the file finished during points acquisition, we need to store them
  if (strncmp(status, 'points', 6))
    good_indx = ~(any(isnan([x y]),2));

    x = x(good_indx);
    y = y(good_indx);
    [x, y] = poly2cw(x, y);

    % And store them
    shapes{group_indx,slice_indx + frame_shift} = [x y];
  end

  % Close the .shapes file
  fclose(fid);

  return;
end
