function [mymovie,trackings] = RECOS(varargin)
  %RECOS  REference COordinate System automated mapping of C. elegans embryos
  %
  %  mymovie = RECOS  returns
  %
  %  mymovie = RECOS(mymovie)  returns
  %
  %  mymovie = RECOS(mymovie, CONDS)  returns
  %
  %  CONDS may be 

  if (~exist('USE_MM', 'var'))
    global USE_MM;
    if (exist('mmversion') == 3)
      USE_MM = true;
    else
      USE_MM = false;
    end
  end

  [mymovie, trackings, opts] = parse_input(varargin{:});
  expr_name = '';

  if (isempty(opts))
    return;
  end

  if (opts.debug)
    keyboard;
  end

  display('Parameters loaded !');

  if (opts.measure_performances || ~strncmp(opts.do_ml, 'none', 4))
    [trackings, expr_name, updated] = import_trackings(trackings, opts);

    if (updated && opts.auto_save)
      save(mymovie.experiment, 'mymovie', 'trackings');
    end
  end

  display('Trackings loaded !');

  mymovie = open_movie(mymovie, expr_name);

  display('Movie loaded !');

  if (isempty(mymovie.experiment))
    return;
  else
    disp(['--Movie: ' mymovie.experiment]);
  end

  switch (opts.segmentation_type)
    case {'dic', 'all'}
      [imgsize nframes] = size_data(mymovie.dic);
    case 'markers'
      switch (opts.ml_type)
        case {'eggshell', 'all'}
          [imgsize nframes] = size_data(mymovie.eggshell);
        case 'cortex'
          [imgsize nframes] = size_data(mymovie.cortex);
      end
    otherwise
      error 'None of the expected field are present in ''mymovie''';
  end

  opts.segmentation_parameters = set_image_size(opts.segmentation_parameters, imgsize);

  if (~strncmp(opts.do_ml, 'none', 4))
    display('ML starts !');

    opts = find_parameters(mymovie, trackings, opts);
  end

  if (opts.auto_save)
    save(mymovie.experiment, 'mymovie', 'trackings');
  end

  if (opts.segment)
    [mymovie, updated] = segment_movie(mymovie, opts);

    if (opts.auto_save && updated)
      save(mymovie.experiment, 'mymovie', 'trackings');
    end
  end

  if (opts.export_movie)
    export_movie(mymovie, opts);
  end
  
  if (opts.measure_performances)

    [mymovie, trackings] = analyse_segmentation(mymovie, trackings, opts);

    if (opts.auto_save)
      save(mymovie.experiment, 'mymovie', 'trackings');
    end
  end

  if (length(opts.application) > 0)
    mymovie = perform_application(mymovie, opts);

    if (opts.auto_save)
      save(mymovie.experiment, 'mymovie', 'trackings');
    end
  end

  if (opts.normalize)
    mymovie = carth2RECOS(mymovie, opts);

    if (opts.auto_save)
      save(mymovie.experiment, 'mymovie', 'trackings');
    end
  end

  % Really save at least once
  if (~opts.auto_save)
    save(mymovie.experiment, 'mymovie', 'trackings');
  end

  beep;beep;beep;

  return;
end

function [mymovie, trackings, opts] = parse_input(varargin)

  mymovie = [];
  trackings = [];
  opts = [];
  
  if (nargin > 0 && isstruct(varargin{1}) && ~isfield(varargin{1}, 'segmentation_type'))
    mymovie = varargin{1};
    varargin(1) = [];
  elseif (nargin > 0 && ischar(varargin{1}))

    if (~isempty(findstr(varargin{1}, '*')))
      datas = dir(varargin{1});
      for i=1:length(datas)
        RECOS(datas(i).name, varargin{2:end});
      end

      return;
    end

    done = false;
    real_mat = true;
    waiting_time = 0;
    while (~done)
      try
        if (strncmp(varargin{1}(end-3:end), '.mat', 4) && exist(varargin{1}, 'file'))
          load(varargin{1});
          varargin(1) = [];
        elseif (exist([varargin{1} '.mat'], 'file'))
          load([varargin{1} '.mat']);
          varargin(1) = [];
        else
          real_mat = false;
          mymovie = get_struct('mymovie',1);
        end

        done = true;
      catch ME
        if (real_mat && waiting_time < 600)
          nsecs = rand(1) * 30;
          waiting_time = waiting_time + nsecs;

          disp(['MAT file is currently used, trying again in ' str2double(nsecs) 's']);
          pause(nsecs);
        else
          error(ME);
        end
      end
    end
  else
    mymovie = get_struct('mymovie',1);
  end

  if (~exist('trackings', 'var'))
    trackings = get_struct('tracking');
  end

  if (nargin > 0 && ~isempty(varargin) && isstruct(varargin{1}))
    opts = varargin{1};
    varargin(1) = [];
  else
    opts = get_struct('RECOS',1);
  end

  npairs = length(varargin) / 2;
  if (npairs ~= floor(npairs))
    error 'Properties pairs must come in PAIRS.';
  end

  for i = 1:npairs
    if (isfield(opts, varargin{2*i - 1}))
      opts.(varargin{2*i - 1}) = varargin{2*i};
    else
      error(['Property ''' varargin{2*i -1} ''' does not exist.']);
    end
  end

  opts = load_parameters(opts);

  return;
end
