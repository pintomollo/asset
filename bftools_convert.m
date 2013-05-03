function [newfile, policy] = bftools_convert(fname, policy, opts)

  % Create a null policy if it is not provided
  if (nargin == 2)
    opts = policy;
    policy = 0;
  end

  % If there is no movie file, stop here
  if(isempty(fname) | ~ischar(fname))
    return;
  end

  % We need the absolute path for Java to work properly
  fname = absolutepath(fname);

  if (exist(fname, 'file') ~= 2)
    error('ASSET:InvalidID', ['File ' fname ' does not exist']);
  end

  % Split the filename using the provided pattern to change its extension to OME-TIFF
  [tokens,junk]=regexp(fname, opts.file_regexpr, 'tokens');
  name = tokens{1}{1};
  suffix = tokens{1}{2};
  ext = tokens{1}{3};

  % Remove the extension
  name = [name suffix];

  % Identify the filename VS the path
  [slash] = findstr(name, filesep);
  if(length(slash)>0)
    slash = slash(end) + 1;
  else
    slash = 1;
  end
  
  % Creat the fancy name for display (otherwise it thinks they are LaTeX commands)
  printname = strrep(name(slash:end),'_','\_');

  curdir = pwd;
  %cmd_path = which('bftools')
  cmd_path = which('bfconvert.bat');
  if (isempty(cmd_path))
    infos = imfinfo(fname);

    if (strncmpi(infos.Format, 'tif', 3))
      warning('ASSET:lociMissing','The LOCI command line tools are not present but the file seems to be a TIF already, trying without conversion...\nYou can download them from http://www.loci.wisc.edu/bio-formats/downloads');

      newfile = fname;
      newfile = relativepath(newfile);

      return;
    else
      error('ASSET:lociMissing', 'The LOCI command line tools are not present !\nPlease download them from http://www.loci.wisc.edu/bio-formats/downloads');
    end
  end
  [mypath, junk] = fileparts(cmd_path);

  cd(mypath);

  if (ispc)
    cmd_name = ['"' fname '"'];
    [res, infos] = system(['showinf.bat -nopix -nometa ' cmd_name]);
  else
    cmd_name = strrep(fname,' ','\ ');
    [res, infos] = system(['./showinf -nopix -nometa ' cmd_name]);
  end
  if (res ~= 0)
    cd(curdir);
    error(infos);
  end
  
  format = regexp(infos, 'file format\s*\[([ -\w]+)\]', 'tokens');
  is_rgb = regexp(infos, 'RGB\s*=\s*(\w+)', 'tokens');

  if (isempty(format) | isempty(is_rgb))
    cd(curdir);
    error('ASSET:lociFormat', ['The metadata does not present the expected information: ''file format'' and ''RGB'' :\n\n' infos]);
  end

  format = format{1}{1};
  is_rgb = strncmp(is_rgb{1}{1}, 'true', 4);
  
  if (strncmpi(format,'OME-TIFF',8) | (length(format) >= 19 & strncmpi(format([1 8 14 19]), 'TIFF', 4)) | strncmpi(format, 'TIFF', 4))

    % If it's already an OME-TIFF file, we do not need to convert it, so stop here
    newfile = fname;
    cd(curdir);
    newfile = relativepath(newfile);

    return;
  end

%  if (is_rgb)
%    warning('ASSET:RGB', ['ASSET does not support RGB recordings, thus only the first channel of ' printname ' will be kept.\nIf it is not a grayscale recording, data might be lost !']);
%  end
  if (is_rgb)
    % We create an OME-TIFF file
    %%newname = [name '_%c.ome.tiff'];
    disp('Should handle RGB here')
    newname = [name '.ome.tiff'];
  else
    % We create an OME-TIFF file
    newname = [name '.ome.tiff'];
  end

  % If the file already exists, we have several possibilities (c.f. help section)
  if(exist(newname,'file'))

    % We initially do not know what to do
    answer = 0;

    % If there is no provided policy, we need to ask the user
    if (policy == 0)

      % We do not accept "empty" answers
      while (answer == 0)
        answer = menu(['The OME-TIFF version of ' printname ' already exists, overwrite it ?'],'Yes', 'Yes to All','No','No to All');
      end

    % Otherwise, we use the provided policy
    else
      answer = policy;
    end

    % Store the answer in case we need to output it
    policy = answer;

    % If it's a "to All" answer, act as it's a standard yes/no
    if (mod(answer, 2) == 0)
      answer = answer - 1;

    % Otherwise, reset the policy for the potential next file
    else
      policy = 0;
    end

    % Act accorindly to the policy
    switch answer

      % Delete the current files (did not dare overwriting it directly)
      case 1
        delete(newname);

      % Otherwise we can stop here
      case 3
        % Store the new name
        newfile = newname;
        cd(curdir);
        newfile = relativepath(newfile);

        return;
    end
  end
  
  if (ispc)
    cmd_newname = ['"' newname '"'];
    [res, infos] = system(['bfconvert.bat -separate ' cmd_name ' ' cmd_newname]);
  else
    cmd_newname = strrep(newname,' ','\ ');
    [res, infos] = system(['./bfconvert -separate ' cmd_name ' ' cmd_newname]);
  end
  if (res ~= 0)
    cd(curdir);
    error(infos);
  end

  newfile = newname;
  cd(curdir);
  newfile = relativepath(newfile);

  return;
end
