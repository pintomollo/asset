function [res, pos] = combine_domains(mymovies, min_num)
  %close all;
  %fname = regexprep(mymovies, '\*', '');

  if (nargin == 1)
    min_num = 2;
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
  count = [];

  for i=1:nmovies
    load(mymovies{i});
    try
      domain = imnorm(gather_quantification(mymovie, opts));
    catch ME
      continue;
    end
    opts = load_parameters(opts, 'domain_center.txt');
    mymovie.data.domain = dynamic_programming(domain, opts.quantification.params, @weight_symmetry, opts.quantification.weights, opts);
    [domain, pos, indx] = align_domain(mymovie, opts);
    boundary = min(indx - 1, size(domain, 2) - indx);
    time = get_manual_timing(mymovie, opts);
    if (all(isnan(time)))
      warning(['No valid timing available for ' mymovie.experiment ', skipping.'])

      continue;
    end
    time = time(end);
    domain = domain(:,[-boundary:boundary]+indx);
    indx = boundary + 1;
    valids = ~isnan(domain);

    if (i == 1)
      all_indx = indx;
      all_cytok = time;
      tmp_res = domain;
      tmp_res(~valids) = 0;
      res = tmp_res;
      count = double(valids);
    else
      if (boundary >= all_indx)
        tmp_res = padarray(tmp_res, [0 boundary - all_indx + 1], 0, 'both');
        count = padarray(count, [0 boundary - all_indx + 1], 0, 'both');
        all_indx = boundary + 1;
      end
      if (time > all_cytok)
        tmp_res = padarray(tmp_res, [time - all_cytok 0], 0, 'pre');
        count = padarray(count, [time - all_cytok 0], 0, 'pre');
        all_cytok = time;
      end
      ends = size(domain, 1) - time;
      all_ends = size(tmp_res, 1) - all_cytok;

      if (ends > all_ends)
        tmp_res = padarray(tmp_res, [ends - all_ends 0], 0, 'post');
        count = padarray(count, [ends - all_ends 0], 0, 'post');
        all_ends = ends;
      end

      res = tmp_res ./ count;

      width = all_indx - boundary;
      start = all_cytok - time + 1;
      ends = all_ends - ends;

      window = res([start:end-ends], [width:end-width+1]);
      c = robustfit(domain(valids), window(valids));
      domain = c(1) + c(2)*domain;
      domain(~valids) = 0;

      tmp_res([start:end-ends], [width:end-width+1]) = tmp_res([start:end-ends], [width:end-width+1]) + domain;
      count([start:end-ends], [width:end-width+1]) = count([start:end-ends], [width:end-width+1]) + double(valids);


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
  res = tmp_res ./ count;
  valids = any(isnan(res), 1);

  first = all_indx - find(valids(1:all_indx), 1, 'last');
  last = find(valids(all_indx:end), 1, 'first') - 1;
  boundary = min(first, last);

  first = find(max(count, [], 2) >= min_num, 1, 'first');

  res = imnorm(res([first:all_cytok],[-boundary:boundary]+all_indx));
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
