function [mymovie, mymovies] = combine_movies(mymovie, method)

  if (nargin == 1)
    method = 'mean';
  end
  repl_str = 'MERGE';

  if (ischar(mymovie))
    fnames = dir(mymovie);
    mymovie = {};
    indx = 1;
    for i=1:length(fnames)
      if (isempty(strfind(fnames(i).name, repl_str)))
        mymovie(indx) = {fnames(i).name};
        indx = indx + 1;
      end
    end
  end

  if (iscell(mymovie))
    for i=1:length(mymovie)
      tmp = load(mymovie{i});
      mymovies(i) = tmp.mymovie;
    end
  elseif (isstruct(mymovie))
    mymovies = mymovie;
  end

  if (strncmp(method, 'combine', 7))
    mymovie = mymovies;

    return;
  end

  nmovies = length(mymovies);
  timings = zeros(nmovies, 3);
  fields = false(nmovies, 3);

  new_name = '';

  for i=1:nmovies
    if (i==1)
      new_name = mymovies(i).experiment;
    else
      new_name = common_substring(new_name, mymovies(i).experiment, repl_str);
    end

    if (isfield(mymovies(i), 'dic') & isfield(mymovies(i).dic, 'cytokinesis'))
      fields(i,1) = true;
    end
    if (isfield(mymovies(i), 'markers') & isfield(mymovies(i).markers, 'cytokinesis'))
      fields(i,2) = true;
    end
    if (isfield(mymovies(i), 'data') & isfield(mymovies(i).data, 'centrosomes') & any(fields(i,1:2)))
      fields(i,3) = true;
    end

    nframes = 0;
    if (fields(i,1) & isfield(mymovies(i).dic, 'fname'))
      [junk, nframes] = size_data(mymovies(i).dic);
      zero = mymovies(i).dic.cytokinesis;
    elseif (fields(i,2))
      if (isfield(mymovies(i), 'eggshell') & isfield(mymovies(i).eggshell, 'fname'))
        [junk, nframes] = size_data(mymovies(i).eggshell);
      elseif (isfield(mymovies(i), 'cortex') & isfield(mymovies(i).cortex, 'fname'))
        [junk, nframes] = size_data(mymovies(i).cortex);
      elseif (isfield(mymovies(i), 'markers') & isfield(mymovies(i).markers, 'fname'))
        [junk, nframes] = size_data(mymovies(i).markers);
      end
      zero = mymovies(i).markers.cytokinesis;
    end
    if (nframes == 0)
      disp('Error, cannot find any suitable field to measure time');
      beep;
      keyboard;
    else
      timings(i,1:3) = [zero 1-zero nframes-zero];
    end

    mymovies(i) = carth2RECOS(mymovies(i));
  end
  range = [min(timings(:,2)):max(timings(:,3))];
  timings(:,4) = timings(:,1) - find(range == 0);

  mymovie = mymovies(1); 
  mymovie.experiment = new_name;
  myfields = any(fields, 1);

  if (myfields(1))
    mymovie.dic.cytokinesis = find(range == 0);

    mymovie.dic.fname = '';
    mymovie.dic.centers = [];
    mymovie.dic.axes_length = [];
    mymovie.dic.orientations = [];
    mymovie.dic.eggshell = [];
    mymovie.dic.cortex = [];
    mymovie.dic.ruffles = [];
    mymovie.dic.warpers = get_struct('warper', [1 length(range)]);
  end
  if (myfields(2))
    mymovie.markers.cytokinesis = find(range == 0);

    mymovie.markers.fname = '';
    mymovie.markers.centers = [];
    mymovie.markers.axes_length = [];
    mymovie.markers.orientations = [];
    mymovie.markers.eggshell = [];
    mymovie.markers.cortex = [];
    mymovie.markers.ruffles = [];
    mymovie.markers.warpers = get_struct('warper', [1 length(range)]);
  end
  if (myfields(3))
    mymovie.data.centrosomes = [];
  end

  for indx = 1:length(range)
    splines = get_struct('spline', [nmovies 2]);
    centers = NaN(2,2,nmovies);

    for i = 1:nmovies
      tmp_indx = indx + timings(i,4);

      if (fields(i,1) & range(indx) >= timings(i,2) & range(indx) <= timings(i,3))
        splines(i,1) = create_spline(mymovies(i).dic.cortex(tmp_indx).warped);
      end
      if (fields(i,2) & range(indx) >= timings(i,2) & range(indx) <= timings(i,3))
        splines(i,2) = create_spline(mymovies(i).markers.cortex(tmp_indx).warped);
      end
      if (fields(i,3) & range(indx) >= timings(i,2) & range(indx) <= timings(i,3))
        centers(:,:,i) = mymovies(i).data.centrosomes(tmp_indx).warped;
      end
    end

    if (myfields(1))
      tmp_spline = mean_spline(splines(:,1));
      
      if (~isempty(tmp_spline.breaks))
        mymovie.dic.cortex(indx).warped = fnval(tmp_spline, tmp_spline.breaks).';
        mymovie.dic.cortex(indx).all = splines(:,1);
      end
    end
    if (myfields(2))
      tmp_spline = mean_spline(splines(:,2));
      if (~isempty(tmp_spline.breaks))
        mymovie.markers.cortex(indx).warped = fnval(tmp_spline, tmp_spline.breaks).';
        mymovie.markers.cortex(indx).all = splines(:,1);
      end
    end
    if (myfields(3))
      [means stds] = mymean(centers, 3);
      mymovie.data.centrosomes(indx).warped = [means stds];
      mymovie.data.centrosomes(indx).all = centers;
    end
  end

  return;
end
