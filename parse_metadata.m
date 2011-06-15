function metadata = parse_metadata(fname, opts)

  metadata = [];
  mymovie = [];
  z_pos = [];
  if (isstruct(fname))
    mymovie = fname;
    fname = mymovie.experiment;
  end

  if (~exist(fname, 'file'))
    if (exist([fname '.txt'], 'file'))
      fname = [fname '.txt'];
    elseif (exist(['Movies' filesep fname], 'file'))
      fname = ['Movies' filesep fname];
    elseif (exist(['Movies' filesep fname '.txt'], 'file'))
      fname = ['Movies' filesep fname '.txt'];
    elseif (exist(['Movies' filesep 'metadata' filesep fname], 'file'))
      fname = ['Movies' filesep 'metadata' filesep fname];
    elseif (exist(['Movies' filesep 'metadata' filesep fname '.txt'], 'file'))
      fname = ['Movies' filesep 'metadata' filesep fname '.txt'];
    else
      fname = [];
    end
  end
  
  if (~isempty(fname))
    % Open the file, line by line
    fid = fopen(fname,'rt');

    % If it does not work, try the absolute path instead
    if (fid<0)
      fname = absolutepath(fname);
      fid = fopen(fname,'rt');
    end

    if (~isempty(mymovie) & ~isempty(mymovie.dic))
      z_pos = NaN(1, length(mymovie.dic.eggshell));
    end

    % Read the first line
    line = fgetl(fid);
    index = -1;
    pos = NaN;
    % Loop until we have read the whole file
    while ischar(line)

      if (any(line == '{'))
        if (index > 0)
          z_pos(index) = pos;
        end
        index = -1;
        pos = NaN;
      else
        tokens = regexp(line, '^\s*"(.+)": "?(.*?)"?,?$','tokens');
        if (~isempty(tokens))
          tokens = tokens{1};
        end
        
        if (numel(tokens) == 2)
          switch tokens{1}
            case 'Z-um'
              pos = str2double(tokens{2});
            case 'Frame'
              index = str2double(tokens{2}) + 1;
            otherwise
              disp(['Parsing incomplete ...']);
          end
        end
      end
    
      % Read a new line
      line = fgetl(fid);
    end

    fclose(fid);

    metadata.z_pos = z_pos;
  end

  if (~isempty(mymovie))
    mymovie.metadata = metadata;
    metadata = mymovie;
  end

  return;
end
