function display_domain_results(res_num)

  close all;

  [timing, names] = get_manual_timing();
  low_temp = strncmp('1056-14', names, 7);
  low_data = [];
  high_data = [];

  %% TODO
  %% Example segmentation
  %% AVG domain profile vwith std -- stretch bewtween PC & Cytok
  %% RE-align full quantif & avg
  %% ---avg profile of front
  %% Z-reconstruct
  %% DIC segmentation for timing
  %% ---Simulations
  
  switch(res_num)
    case 1
      plot_timings(timing, low_temp);
    case 2
      plot_example();
    case 3
      low_data = plot_all_domains(timing(low_temp, :), names(low_temp));
      high_data = plot_all_domains(timing(~low_temp, :), names(~low_temp));
    case 4
      if (isempty(low_data) | isempty(high_data))
        %realign_domains(timing(low_temp, :), names(low_temp));
        realign_domains(timing(~low_temp, :), names(~low_temp));
      else
        realign_domains(low_data);
        realign_domains(high_data);
      end
  end

  return;
end

function realign_domains(timing, names)

  if (nargin == 2)
    data = plot_all_domains(timing, names);
  else
    data = timing;
  end

  all_domains = [];
  bounds = zeros(1, 2);
  sides = zeros(1, 2);

  all_centers = [];
  all_firsts = [];

  for i=1:size(data, 1)

    tmp = load(data{i, 1});
    [img, pos, center]=align_domain(tmp.mymovie, data{i, 3}, tmp.opts); 
    img = imnorm(img);

    %center = find(pos == 0);
    resolution = median(diff(pos));

    centers = data{i, 3}(:, 1) - center;
    %plot(centers);hold on;
    %plot(data{i,3}(:,1));hold off;
    %keyboard

    times = data{i, 2};

    if (isnan(times(3)))
      continue;
    end

    all_firsts(end+1) = centers(find(~isnan(centers), 1, 'first'));

    min_pos = (1 - center) * resolution;
    max_pos = (length(pos) - center) * resolution;
    
    if (min_pos < sides(1))
      nmin = round((sides(1) - min_pos) / resolution);
      sides(1) = min_pos;
    else
      nmin = 0;
    end
    if (max_pos > sides(2))
      nmax = round((max_pos - sides(2)) / resolution);
      sides(2) = max_pos;
    else
      nmax = 0;
    end     

    first = 1-times(3);
    last = size(img, 1) - times(3);

    if (first < bounds(1))
      nfirst = bounds(1) - first;
      bounds(1) = first;
    else
      nfirst = 0;
    end
    if (last > bounds(2))
      nlast = last - bounds(2);
      bounds(2) = last;
    else
      nlast = 0;
    end     

    if (i==1)
      all_domains = img;
      all_centers = centers;
    else

      [nrows, ncols, ndomains] = size(all_domains);

        all_centers = [NaN(nfirst, ndomains);  all_centers; NaN(nlast, ndomains)];
        first_indx = first - bounds(1) + 1;
        last_indx = bounds(2) - last;

        all_centers(:, end+1) = NaN;
        all_centers(first_indx:first_indx+size(img,1)-1, end) = centers;



      all_domains = cat(2, NaN(nrows, nmin, ndomains), all_domains, NaN(nrows, nmax, ndomains));
      ncols = size(all_domains, 2);
      all_domains = cat(1, NaN(nfirst, ncols, ndomains), all_domains, NaN(nlast, ncols, ndomains));

      shift_indx = round((min_pos - sides(1))/resolution) + 1;

      all_domains(first_indx:first_indx+size(img,1)-1, shift_indx:shift_indx+size(img,2)-1, end+1) = img;



    end

  end
  all_domains = mymean(all_domains, 3);
  figure;imagesc(all_domains);

    all_centers = all_centers * resolution;
    [m, s] = mymean(all_centers, 2);

    figure;
    plot([bounds(1):bounds(2)], all_centers.', 'b');
    hold on;
    plot([bounds(1):bounds(2)], m, 'k', 'LineWidth', 2);
    plot([bounds(1):bounds(2)], m+s, 'k', 'LineWidth', 2);
    plot([bounds(1):bounds(2)], m-s, 'k', 'LineWidth', 2);

    figure;
    hist(all_firsts * resolution)

  return;
end

function data = plot_all_domains(timings, names)

  manuals = dir('*-Manual-DP.mat');
  data = cell(0, 4);
  for i=1:length(manuals)
    tmp_name = manuals(i).name(1:end-18);
    index = find(strncmp(tmp_name, names, length(tmp_name)), 1);

    if (~isempty(index))
      tmp = load(manuals(i).name);
      path = tmp.path;
      tmp = load(names{index});
      [junk, pos] = gather_quantification(tmp.mymovie, tmp.opts);
      data(end+1,:) = {names{index}, timings(index, :), path, pos};
    end
  end

  for i=1:size(timings, 2)
    aligned_pos = [];
    bounds = zeros(1, 2);
    for j = 1:size(data, 1)
      path = data{j, 3} * tmp.opts.pixel_size;
      times = data{j, 2};

      if (isnan(times(i)))
        continue;
      end
      
      first = 1-times(i);
      last = size(path, 1) - times(i);

      if (first < bounds(1))
        nfirst = bounds(1) - first;
        bounds(1) = first;
      else
        nfirst = 0;
      end
      if (last > bounds(2))
        nlast = last - bounds(2);
        bounds(2) = last;
      else
        nlast = 0;
      end     
      if (j==1)
        aligned_pos = path(:,2);
      else
        ncols = size(aligned_pos, 2);
        aligned_pos = [NaN(nfirst, ncols);  aligned_pos; NaN(nlast, ncols)];
        first_indx = first - bounds(1) + 1;
        last_indx = bounds(2) - last;

        aligned_pos(:, end+1) = NaN;
        aligned_pos(first_indx:end-last_indx, end) = path(:,2);
      end
    end
    for i=1:size(aligned_pos, 2)

      first = find(~isnan(aligned_pos(:,i)), 1, 'first');
      last = find(~isnan(aligned_pos(:,i)), 1, 'last')-2;

      aligned_pos(1:first-1, i) = 0;
      aligned_pos(last+1:end, i) = aligned_pos(last,i);
    end
    [m, s] = mymean(aligned_pos, 2);

    figure;
    plot([bounds(1):bounds(2)], aligned_pos.', 'b');
    hold on;
    plot([bounds(1):bounds(2)], m, 'k', 'LineWidth', 2);
    plot([bounds(1):bounds(2)], m+s, 'k', 'LineWidth', 2);
    plot([bounds(1):bounds(2)], m-s, 'k', 'LineWidth', 2);
  end

  return;
end

function plot_timings(timing, groups)

  values = 10*abs(diff(timing(:, [1:end 1]), [], 2));
  all_groups = ones(size(values,1), 1) * [1:size(values, 2)];

  indexes = unique(groups(:)).';
  for i = indexes
    tmp_val = values(groups == i, :);
    tmp_group = all_groups(groups == i, :);

    figure;boxplot(tmp_val(:), tmp_group(:));
    title(num2str(i));
  end

  for i=1:size(values, 2)
    [h,p] = myttest(values(:, i), groups);
    fprintf(1, 'Column %d: %f\n', i, p(1));
  end

  return;
end

function plot_example

  %names = dir('1056-*_.mat');

  %for i=1:length(names)
  %try
  %fname = names(i).name(1:end-4)

  fname = '1056-14-inc_med-080811_2_';
  tmp = load(fname);
  mymovie = tmp.mymovie;
  path = mymovie.data.domain;

  img = gather_quantification(mymovie, tmp.opts);
  [h,w] = size(img);
  
  figure;imagesc(img);hold on;
  plot(path(:, 1)*w,[1:h], 'k');
  plot((path(:, 1) + path(:,2))*w, [1:h], 'k');
  plot((path(:, 1) - path(:,2))*w, [1:h], 'k');

  tmp = load([fname(1:end-1) '-DP.-Manual-DP.mat']);
  path = tmp.path;
  plot(path(:, 1),[1:h], 'g');
  plot((path(:, 1) + path(:,2)), [1:h], 'g');
  plot((path(:, 1) - path(:,2)), [1:h], 'g');

  %  pause
  %  catch ME
  %    continue;
  %  end
  %end

  return;
end
