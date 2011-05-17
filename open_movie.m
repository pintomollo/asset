function [mymovie] = open_movie(mymovie, expr_name, opts)

  if (nargin < 1)
    mymovie = get_struct('mymovie',1);
    expr_name = '';
  elseif (nargin < 2)
    expr_name = '';
  end

  if(~isstruct(mymovie) | length(mymovie) == 0 | (isempty(mymovie.dic) & isempty(mymovie.eggshell & isempty(mymovie.cortex))))
    if(ischar(mymovie))
      indx = strfind(mymovie, filesep);

      if (~isempty(indx))
        curdir = cd;
        cd(mymovie(1:indx(end)-1));
        dirpath = [pwd filesep];
        cd(curdir);
        fname = mymovie(indx(end)+1:end);
      else
        fname = mymovie;
        dirpath = pwd;
      end
    else

      % In case a subfolder name Movies exists, move into it
      curdir = '';
      if(exist('Movies', 'dir'))
        curdir = cd;
        cd('Movies');
      elseif(exist(['..' filesep 'Movies'], 'dir'))
        curdir = cd;
        cd(['..' filesep 'Movies']);
      end

      disp('[Select a movie file]');
      %[fnames,path] = uigetfile({'*.*'}, 'Load a file','MultiSelect','on');
      [fname,dirpath] = uigetfile({'*.*'}, ['Load a movie (' expr_name ')']);

      if(~isempty(curdir))
        cd(curdir);
      end

      if (length(fname)==0  ||  isequal(dirpath,0))
        disp(['No movie selected']);
        return;
      end  
    end

    [tokens,junk]=regexp(fname, opts.file_regexpr,'tokens');
    name = tokens{1}{1};
    suffix = tokens{1}{2};
    ext = tokens{1}{3};

    if (strncmp(ext, '.mat', 4))
      load(fname);

      return;
    end

    if (~isempty(name))
      files = dir([dirpath name '*' ext]);

      for i = length(files):-1:1
        [tokens,junk]=regexp(files(i).name,opts.file_regexpr,'tokens');

        if (~strcmp(name, tokens{1}{1}))
          files(i) = [];
        end
      end
    else
      name = suffix;
      suffix = '';

      files = struct('name', [name ext]);
    end

    nchannels = length(files);
    channels = get_struct('channel',[1 nchannels]);
    %channels = repmat(struct('file','','color',zeros(1,3),'type','data','detrend',false,'hot_pixels',false,'fname','','min', Inf, 'max', -Inf, 'metadata', ''),1,nchannels);

    policy = 0;
    for i=1:nchannels
      channels(i).file = [dirpath files(i).name];
      channels(i).file = relativepath(channels(i).file);
      %[channels(i).fname, policy] = convert_movie(channels(i).file, 'Uncompressed', policy);
      [channels(i).fname, policy] = convert_movie(channels(i).file, 'LZW', policy, opts);
      %channels(i) = bfopen(channels(i));
    end

    %keyboard

    mymovie = identify_channels(channels, opts);

    mymovie = rescale_movie(mymovie, true);
    mymovie.experiment = name;
  end

  return;
end
