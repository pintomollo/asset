function [res, vars, pos, time] = combine_domains(mymovies, min_num, thresh)
  %close all;
  %fname = regexprep(mymovies, '\*', '');

  if (nargin == 1)
    min_num = 3;
    thresh = 0.775;
  elseif (nargin < 3)
    thresh = 0.775;
  end

  if (ischar(mymovies))
    if (~isempty(findstr(mymovies, '*')))
      tmp = dir(mymovies); 
      mymovies = cell(1, length(tmp));

      for i=1:length(tmp)
        mymovies{i} = tmp(i).name;
      end
    else
      mymovies = {mymovies}; 
    end
  elseif (~iscell(mymovies))
    tmp = mymovies;
    mymovies = {};
    mymovies{1} = tmp;
  end

  nmovies = length(mymovies);
  %tmp_quant = cell(1, nmovies);
  %cytos = NaN(1, nmovies);
  %opts = get_struct('ASSET');
  %paths = cell(1, nmovies);;

  %for i=1:nmovies
  %  load(mymovies{i});

    %if (isfield(mymovie.data, 'domain') & ~isempty(mymovie.data.domain))
    %  tmp_path = mymovie.data.domain;
    %else
    %  continue;
    %end

  %  try
  %  img = gather_quantification(mymovie, opts);
  %  catch ME
  %    continue;
  %  end
%%[img,~,cytos(i)] = gather_quantification(mymovie, opts);

    %tmp_path = mymovie.data.domain;
  %  [h,w] = size(img);

    %imagesc(imnorm(img));
    %hold on
    %plot(tmp_path(:, 1)*w,[1:h], 'k');
    %plot(((tmp_path(:,1) + tmp_path(:, 2))*w)-1,[1:h], 'k');
    %plot(((tmp_path(:,1) - tmp_path(:, 2))*w)+1,[1:h], 'k');

%%plot([1 size(img,2)],[cytos(i) cytos(i)], 'k');

   % imwrite(imnorm(img), [mymovies{i}(1:end-5) '-DP.png']);
    %hold off

    %plot(tmp_path(:,1)*w * opts.pixel_size, 'k');
    %hold on;
    %plot(((tmp_path(:,2) *w) - 1) * opts.pixel_size, 'r');
    %saveas(gca, [mymovies{i}(1:end-5) '-domain.png']);
    %hold off

    %if (isfield(mymovie, 'metadata') & isfield(mymovie.metadata, 'timing'))
    %  mymovie.metadata.timing
    %else
    %  'none'
    %end

    %paths{i} = tmp_path;
  %end

  %res{1} = mymovies;
  %res{2} = paths;
  %res{3} = cytos;

  %return;

  %nframes = last - first + 1;
  %cyto = 1-first;
  %res = NaN(nframes, size(tmp_quant{i}, 2), nmovies);
  %orig_res = res;
  %maxint = double(intmax('uint16'));

  %y_label = [1:nframes]-cyto;
  all_indx = 0;
  all_cytok = 0;
  tmp_res = [];
  res = [];
  vars = [];
  count = [];

  for i=1:nmovies
    load(mymovies{i});
    try
      [domain, ruffles, junk] = gather_quantification(mymovie, opts);
      domain = imnorm(domain);
    catch ME
      continue;
    end
    opts = load_parameters(opts, 'domain_center.txt');
    opts.quantification.weights.filt = ruffles;
    mymovie.data.domain = dynamic_programming(domain, opts.quantification.params, @weight_symmetry, opts.quantification.weights, opts);
    [domain, ruffles, pos, indx] = align_domain(mymovie, opts);
    %[domain, pos, indx] = align_domain(mymovie, opts);
    boundary = min(indx - 1, size(domain, 2) - indx);
    time = get_manual_timing(mymovie, opts);
    if (all(isnan(time)))
      warning(['No valid timing available for ' mymovie.experiment ', skipping.'])

      continue;
    end

    opts = load_parameters(opts, 'domain_expansion.txt');
    fraction = domain_expansion(domain, indx, time(end), opts);
    align_time = find(fraction >  thresh, 1, 'first');
    if (isempty(align_time))
      align_time = time(end);
    end
    maintenance = time(end) - align_time;
    time = align_time;
    %time = time(end);

    %figure;imagesc(domain);
    %hold on;plot([1 size(domain, 2)], [align_time align_time], 'k');
    %title(mymovie.experiment);

    center = mymovie.data.domain(1:time+maintenance);
    center = center - center(end);
    %pos = pos(1:time+maintenance) - pos(end);
    fraction = fraction(1:time+maintenance);

    domain = domain(1:time+maintenance,[-boundary:boundary]+indx);
    %domain = ruffles(:,[-boundary:boundary]+indx);
    indx = boundary + 1;
    valids = ~isnan(domain);
    %valids = valids(1:time+maintenance, :);

    if (isempty(tmp_res))
      all_indx = indx;
      all_align = time;
      all_maintenance = maintenance;
      tmp_res = domain(1:time+maintenance, :);
      tmp_res(~valids) = 0;
      vars = tmp_res.^2;
      res = tmp_res;
      count = double(valids);
      all_valids = valids;

%      figure;
%      subplot(1,2,1)
%      plot([-time:maintenance-1], center.');
%      hold on;
%
 %     subplot(1,2,2)
 %     plot([-time:maintenance-1], fraction.');
 %     hold on;
    else
      if (boundary >= all_indx)
        tmp_res = padarray(tmp_res, [0 boundary - all_indx + 1], 0, 'both');
        count = padarray(count, [0 boundary - all_indx + 1], 0, 'both');
        vars = padarray(vars, [0 boundary - all_indx + 1], 0, 'both');
        all_valids = padarray(all_valids, [0 boundary - all_indx + 1], 0, 'both');
        all_indx = boundary + 1;
      end
      if (time > all_align)
        tmp_res = padarray(tmp_res, [time - all_align 0], 0, 'pre');
        count = padarray(count, [time - all_align 0], 0, 'pre');
        vars = padarray(vars, [time - all_align 0], 0, 'pre');
        all_valids = padarray(all_valids, [time - all_align 0], 0, 'pre');
        all_align = time;
      end
      %if (maintenance < all_maintenance)
      %  all_maintenance = maintenance;
      %end
      %ends = size(domain, 1) - time;
      %all_ends = size(tmp_res, 1) - all_align;

      if (maintenance > all_maintenance)
        tmp_res = padarray(tmp_res, [maintenance - all_maintenance 0], 0, 'post');
        count = padarray(count, [maintenance - all_maintenance 0], 0, 'post');
        vars = padarray(vars, [maintenance - all_maintenance 0], 0, 'post');
        all_valids = padarray(all_valids, [maintenance - all_maintenance 0], 0, 'post');
        all_maintenance = maintenance;
      end

      res = tmp_res ./ count;

      width = all_indx - boundary;
      start = all_align - time + 1;
      ends = all_maintenance - maintenance;

      valids = all_valids([start:end-ends], [width:end-width+1]) & valids(1:time+maintenance, :);

      window = res([start:end-ends], [width:end-width+1]);
      c = robustfit(domain(valids), window(valids));
      domain = c(1) + c(2)*domain;
      domain(~valids) = 0;

      tmp_res([start:end-ends], [width:end-width+1]) = tmp_res([start:end-ends], [width:end-width+1]) + domain(1:time+maintenance, :);
      vars([start:end-ends], [width:end-width+1]) = vars([start:end-ends], [width:end-width+1]) + domain(1:time+maintenance, :).^2;
      count([start:end-ends], [width:end-width+1]) = count([start:end-ends], [width:end-width+1]) + double(valids);
      all_valids([start:end-ends], [width:end-width+1]) = valids(1:time+maintenance, :);


%      subplot(2,1,1)
%      hold on;
%      plot([-time:maintenance-1], center.');
%      subplot(2,1,2)
%      hold on;
%      plot([-time:maintenance-1], fraction.');
    %imagesc(res);
    %drawnow;
    end

    %orig_res((cyto-cytos(i))+[1:tmp_frame],:,i) = tmp_quant{i}*(minmax(i,2) - minmax(i,1))/maxint + minmax(i,1);;

    %imagesc(orig_res(:,:,i));
    %hold on
    %set(gca,'YTickLabel',y_label)
    %saveas(gca, [mymovies{i}(1:end-5) '.png']);
    %hold off
    %imagesc(imnorm(orig_res(:,:,i), [], [], 'row'));
    %hold on
    %set(gca,'YTickLabel',y_label)
    %saveas(gca, [mymovies{i}(1:end-5) '-scaled.png']);
    %hold off
    %imagesc(res(:,:,i));
    %saveas(gca, [mymovies{i}(1:end-5) '-norm.png']);
    %imagesc(res(:,:,i));
    %saveas(gca, [mymovies{i}(1:end-5) '-norm.png']);
  end
  vars = vars ./ count;
  res = tmp_res ./ count;
  
  all_valids = all_valids & count >= min_num;

  valids = any(all_valids, 2);
  first_frame = find(valids, 1, 'first');
  last_frame = find(valids, 1, 'last');

  valids = any(all_valids, 1);
  
  %vars = sqrt(vars - res.^2) ./ sqrt(count);
  vars = sqrt(vars - res.^2);

  first = all_indx - find(valids(1:all_indx), 1, 'first');
  last = find(valids(all_indx:end), 1, 'last') - 1;
  boundary = min(first, last);

  %first = find(max(count, [], 2) >= min_num, 1, 'first');

  res(~all_valids) = NaN;
  vars(~all_valids) = NaN;

  res = (res([first_frame:last_frame],[-boundary:boundary]+all_indx));
  vars = (vars([first_frame:last_frame],[-boundary:boundary]+all_indx));
  %res = (res([first:all_align],[-boundary:boundary]+all_indx));
  %vars = (vars([first:all_align],[-boundary:boundary]+all_indx));


  pos = [-boundary:boundary] * median(diff(pos));

  %[val1, errs] = mymean(res, 3);
  %imagesc(val1);
  %saveas(gca, [fname '-norm.png']);

  %[val1, errs] = mymean(orig_res, 3);
  %imagesc(val1);
  %saveas(gca, [fname '-scaled.png']);

  return;

  npts = size(res, 2);
  theta = [0:2*pi/npts:2*pi];
  theta = theta(1:end-1);
  all_thetas = repmat(theta,[1 1 nmovies]);
  all_thetas = all_thetas(:);

  hfig = figure;

  hax = axes;
    set(hax, ...
    'FontSize',18, ...
    'NextPlot','add', ...
    'Tag', 'axes',...
    'Parent',hfig, ...
    'Ylim', [0 1], ...
    'Xlim', [0 2*pi]);

  for i=1:nframes
    tmp_res = res(i,:,:);
    errorbar(theta, val1(i,:), errs(i,:), 'b');
    hold on;
    scatter(all_thetas, tmp_res(:), 'k');
    plot(theta, val1(i,:), 'r', 'LineWidth', 2);
    text(0.1,min(res(:))+0.05,[num2str(10*(i-cyto)) 's (' num2str(sum(~isnan(res(i,1,:)))) ')']);
    hold off;
    ylim([0 1]);
    xlim([0 2*pi])
    drawnow
    movie(i) = getframe(hfig);
  end

  movie2avi(movie,[fname '-norm.avi'],'FPS',15);
  res = orig_res;

  for i=1:nmovies
    res(:,:,i) = res(:,:,i)*(minmax(i,2) - minmax(i,1)) + minmax(i,1);
  end
  [val2, errs] = mymean(res, 3);
  figure;imagesc(val2);
  saveas(gca, [fname '-scale.png']);

  hfig = figure;

      hax = axes;
    set(hax, ...
    'FontSize',18, ...
    'NextPlot','add', ...
    'Tag', 'axes',...
    'Parent',hfig, ...
    'Ylim', [min(res(:)) max(res(:))], ...
    'Xlim', [0 2*pi]);

  for i=1:nframes
    tmp_res = res(i,:,:);
    errorbar(theta, val2(i,:), errs(i,:), 'b');
    hold on;
    scatter(all_thetas, tmp_res(:), 'k');
    plot(theta, val2(i,:), 'r', 'LineWidth', 2);
    text(0.1,min(res(:))+3,[num2str(10*(i-cyto)) 's (' num2str(sum(~isnan(res(i,1,:)))) ')']);
    hold off;
    ylim([min(res(:)) max(res(:))]);
    xlim([0 2*pi])
    drawnow
    movie(i) = getframe(hfig);
  end
  movie2avi(movie,[fname '-scale.avi'],'FPS',15);

  return;
end
