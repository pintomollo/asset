%function display_mcmc_results(fname, merge, varying)
function display_mcmc_results(data, func)

  %{
  if (nargin == 1)
    merge = false;
  end

  if (ischar(fname))
    if (~isempty(findstr(fname, '*')))
      files = dir(fname);
      [path, name, ext] = fileparts(fname);

      if (merge)
        obj_val = -1;
        if (nargin == 3 & isnumeric(varying))
          obj_val = varying;
        end
        for i=1:length(files)
          if (~isempty(files(i).name))
            [header, data, varying] = parse_mcmc_results(fullfile(path, files(i).name));
            nparams = sum(varying);

            if (obj_val < 0 | nparams == obj_val)
              display(['Loaded ' files(i).name]);

              for j=i+1:length(files)
                [tmp_header, tmp_data, tmp_varying] = parse_mcmc_results(fullfile(path, files(j).name));
                if (sum(tmp_varying) == nparams)
                  display(['Loaded ' files(j).name]);

                  tmp_header(1,:) = NaN;
                  tmp_data(1,:) = NaN;

                  header = [header; tmp_header];
                  data = [data; tmp_data];
                  files(j).name = '';
                end
              end

              display_mcmc_results(header, data, varying);
            end
          end
        end
      else
        for i=1:length(files)
          display_mcmc_results(fullfile(path, files(i).name));
        end
      end

      return;
    end

    [header, data, varying] = parse_mcmc_results(fname);

  elseif (isnumeric(fname))
    if (nargin ==3)
      header = fname;
      data = merge;
    elseif (nargin == 2)
      varying = true(1, size(data, 2));
    else
      return;
    end
  else
    return;
  end
  %}

  if (nargin == 1)
    func = @mean;
  end

  header = data(:,1);
  %data = data(:,2:end).^2;
  data = data(:,2:end);
  %disp('Squared !');

  %varying = true(1, size(data,2));
  merge = false;
  %nparams = sum(varying);
  nparams = size(data,2);

  if (size(header, 2) == 4)
    valids = (header(:, 4) ~= 1);

    errs = data(~valids, :);

    header = header(valids, :);
    data = data(valids, :);
  end

  %best = data(1,varying);
  %data = data(2:end, varying);
  %data = data(2:end, varying);
  %%header = header(2:end, :);
  ndims = size(data, 2);

  nbins = 32;

  gaps = any(isnan(data), 2);
  gaps_index = find(gaps).';

  tmp_data = data(~gaps, :);
  %edges = sort(tmp_data);
  %edges = edges([1:ceil(end/nbins):end end],:);

  %discards = (diff(header(:,3)) == 1);
  %path = data(~discards, :);
  path_index = ceil(([1:nbins]/nbins)*size(data, 1));
  path_index = unique([path_index 1 gaps_index gaps_index+1 gaps_index-1]);
  path = data(path_index,:);

  for n=1:ndims
    for m=n+1:ndims

      x = unique(data(:, n));
      y = unique(data(:, m));

      last_x = x(end);
      last_y = y(end);

      x = x([1:ceil(end/nbins):end end],:);
      y = y([1:ceil(end/nbins):end end],:);

      if (x(end)==x(end-1))
        x(end) = x(end) + (x(end-1)-x(end-2));
      else
        x(end+1) = x(end) + (x(end)-x(end-1));
      end

      if (y(end)==y(end-1))
        y(end) = y(end) + (y(end-1)-y(end-2));
      else
        y(end+1) = y(end) + (y(end)-y(end-1));
      end

      %dx = min(diff(x));
      %dy = min(diff(y));

      %full_x = [x(1):dx:x(end)];
      %full_x = [full_x, full_x(end)+dx];
      %full_y = [y(1):dy:y(end)];
      %full_y = [full_y, full_y(end)+dy];

      %x(end) = full_x(end);
      %y(end) = full_y(end);

      %[X,Y] = meshgrid(full_x, full_y);
     % [full_map, junk, junk, areas] = histcn([X(:) Y(:)], x, y);

      figure;hold on;

      [map, pos_edges, xs, pos] = histcn(data(:, [n m]), x, y);
      [coords, junk, index] = unique(pos, 'rows');

      score_map = NaN(size(map)+1);
      params_map = cell(size(map));
      centers = NaN(size(coords));

      for i=1:size(coords, 1)
        if (any(coords(i, :) == 0))
          continue;
        end

        currents = (index==i);

        tmp_score = header(currents, 1);
        %[min_val, min_indx] = min(tmp_score);
        %%[min_val, min_indx] = func(tmp_score);
        %%score_map(coords(i, 1), coords(i, 2)) = min_val;
        score_map(coords(i, 1), coords(i, 2)) = func(tmp_score);
        %score_map(coords(i, 1)==areas(:,1) & coords(i, 2) == areas(:,2)) = min_val;

        %%tmp_params = data(currents, :);
        %%params_map{coords(i, 1), coords(i, 2)} = tmp_params(min_indx, :);
        %%centers(i,:) = tmp_params(min_indx, [n m]);
      end

      %pcolor(x,y,[[exp(score_map.') zeros(length(y)-1, 1)]; zeros(1,length(x))]);
      %pcolor(x,y,exp(score_map.'));
      pcolor(x,y,score_map.');
      set(gca, 'XScale', 'log', 'YScale', 'log');

      whitebg('k');
%      scatter(errs(:, n), errs(:, m), '+b');
      %pcolor(xs{:}, exp(score_map));
      %imagesc(exp(score_map));
%      plot(path(:, n), path(:, m), 'w');
      %%scatter(centers(:, 1), centers(:, 2), 'xk');
      %%scatter(data([1 gaps_index+1], n), data([1 gaps_index+1], m), '^g');
      %%scatter(best(1, n), best(1, m), 'oy', 'filled');
      %axis equal;
      axis([x([1 end]); y([1 end])].');
      mtit(num2str(ndims));
      xlabel(num2str(n))
      ylabel(num2str(m))
    end
  end

  return;
end
