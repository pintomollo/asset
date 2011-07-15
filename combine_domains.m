function [res] = combine_domains(mymovies)
  close all;
  %fname = regexprep(mymovies, '\*', '');

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
  cytos = NaN(1, nmovies);
  opts = get_struct('ASSET');
  paths = cell(1, nmovies);;

  for i=1:nmovies
    load(mymovies{i});
    [img,~,cytos(i)] = gather_quantification(mymovie, opts);
    %tmp_path = mymovie.data.domain;

    [h,w] = size(img);

    imagesc(imnorm(img));
    hold on
    %plot(tmp_path(:, 1)*w,[1:h], 'k');
    %plot(tmp_path(:, 2)*w,[1:h], 'k');
    plot([1 size(img,2)],[cytos(i) cytos(i)], 'k');

    saveas(gca, [mymovies{i}(1:end-5) '-DP.png']);
    hold off

    %paths{i} = tmp_path;
  end

  %res{1} = mymovies;
  %res{2} = paths;
  %res{3} = cytos;

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
