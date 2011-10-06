function params = learn_separation(mymovie, nneighs, opts)

  centers = mymovie.dic.centers;
  axes_length = mymovie.dic.axes_length;
  orientations = mymovie.dic.orientations;
  eggshell = mymovie.dic.eggshell;
  cortex = mymovie.dic.cortex;
  neighbors = mymovie.dic.neighbors;
  [nframes, imgsize] = size_data(mymovie.dic);

  parameters = opts.segmentation_parameters.dic;
  mymovie.dic.parameters = parameters;

  img = [];

  frames = randperm(nframes);
  frames = frames(1:10);

      opt = cmaes('defaults');
      opt.MaxFunEvals = 100000;
      opt.TolFun = 1e-25;
      opt.LBounds = zeros(4,1);
      opt.SaveFilename = '';
      opt.EvalParallel = 'yes';
      opt.LogPlot = 'false';

      opt.LogFilenamePrefix = ['ML-' mymovie.experiment];

  params = [1/3, pi/20, 15, 0.05].';
  nparams = length(params);

  for nimg = frames
    img = imnorm(double(load_data(mymovie.dic,nimg)));

    img = imadm_mex(img);
    thresh = graythresh(img);
    img = (img > thresh*0.5*(max(img(:))) );

    [params, fval, ncoutns, stopflag, out] = cmaes(@error_function, params, 0.15, opt);
  end

  return;

  function err_all = error_function(p_all)
    
    [curr_nparams, nevals] = size(p_all);
    flip = false;
    if (nevals == nparams)
      flip = true;
      p_all = p_all.';
      nevals = curr_nparams; 
    end
    err_all = NaN(1, nevals);

    for i = 1:nevals
      new_p = p_all(:, i);


      [ellipse, estim] =  split_cells(img, false, opts, new_p);
      err_all(i) = abs(size(ellipse, 1) - nneighs);
    end

    if (flip)
      err_all = err_all.';
    end

    return;
  end
end
