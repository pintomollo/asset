function [times, names] = get_manual_timing(mymovie, opts)

  times = NaN(1, 3);

  if (nargin == 0)
    fname = 'all';
    mymovie = [];
  else
    if (ischar(mymovie))
      fname = mymovie;
    elseif (isfield(mymovie, 'experiment'))
      fname = mymovie.experiment;
    else
      return;
    end
  end

  if (exist('timings.txt', 'file'))
    fid = fopen('timings.txt');
    values = textscan(fid, '%s %u %f %u');
    fclose(fid);
  elseif (exist(['Config' filesep 'timings.txt'], 'file'))
    fid = fopen(['Config' filesep 'timings.txt']);
    values = textscan(fid, '%s %u %f %u');
    fclose(fid);
  elseif (exist(['celegans-analysis' filesep 'Config' filesep 'timings.txt'], 'file'))
    fid = fopen(['celegans-analysis' filesep 'Config' filesep 'timings.txt']);
    values = textscan(fid, '%s %u %f %u');
    fclose(fid);
  else
    warning('Cannot find the timings text file');

    return;
  end

  names = values{1};
  times = cat(2, values{2:end});

  if (~isempty(mymovie))
    index = find(strncmp(fname, names, length(fname)), 1);

    if (isempty(index))
      times = NaN(1, 3);
      names = fname;
    else
      times = times(index, :);
      names = names{index};
    end
  end

  times = double(times);
  times(times == 0) = NaN;

  return;
end
