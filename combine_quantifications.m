function [res] = combine_quantifications(mymovies)
  close all;
  fname = regexprep(mymovies, '\*', '');

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
  end

  nmovies = length(mymovies);
  %tmp_quant = cell(1, nmovies);
  cytos = NaN(1, nmovies);
  opts = get_struct('ASSET');

  res = [];
  first = 0;
  last = 0;
  minmax = NaN(nmovies, 2);

  for i=1:nmovies
    load(mymovies{i});
    [tmp_quant, ~, cytos(i)] = gather_quantification(mymovie, opts);

    if (1-cytos(i) < first)
      first = 1-cytos(i);
    end
    %if (size(tmp_quant{i}, 1) - cytos(i) > last)
    %  last = size(tmp_quant{i}, 1) - cytos(i);
    %end
    minmax(i,:) = [mymovie.data.min mymovie.data.max]; 


    y_label = [1:5:size(tmp_quant, 1)];

    imagesc(tmp_quant);
    hold on
    set(gca, 'YTick', y_label,'YTickLabel',y_label - cytos(i))
    saveas(gca, [mymovies{i}(1:end-5) '.png']);
    hold off
    imagesc(imnorm(tmp_quant, [], [], 'row'));
    hold on
    set(gca,'YTick',y_label,'YTickLabel',y_label - cytos(i))
    saveas(gca, [mymovies{i}(1:end-5) '-scaled.png']);
    hold off
  end

  return;

  nframes = last - first + 1;
  cyto = 1-first;
  res = NaN(nframes, size(tmp_quant{i}, 2), nmovies);
  orig_res = res;
  maxint = double(intmax('uint16'));

  y_label = [1:nframes]-cyto;

  for i=1:nmovies
    tmp_frame = size(tmp_quant{i}, 1);
    orig_res((cyto-cytos(i))+[1:tmp_frame],:,i) = tmp_quant{i}*(minmax(i,2) - minmax(i,1))/maxint + minmax(i,1);;
    %res(:,:,i) = imnorm(orig_res(:,:,i));

    imagesc(orig_res(:,:,i));
    hold on
    set(gca,'YTickLabel',y_label)
    saveas(gca, [mymovies{i}(1:end-5) '.png']);
    hold off
    imagesc(imnorm(orig_res(:,:,i), [], [], 'row'));
    hold on
    set(gca,'YTickLabel',y_label)
    saveas(gca, [mymovies{i}(1:end-5) '-scaled.png']);
    hold off
    %imagesc(res(:,:,i));
    %saveas(gca, [mymovies{i}(1:end-5) '-norm.png']);
    %imagesc(res(:,:,i));
    %saveas(gca, [mymovies{i}(1:end-5) '-norm.png']);
  end

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
