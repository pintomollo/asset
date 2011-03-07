function display_domain(mymovie, opts, removed)
  
  if (ischar(mymovie))
    if (~isempty(findstr(mymovie, '*')))
      tmp = dir(mymovie); 
      mymovies = cell(1, length(tmp));

      for i=1:length(tmp)
        load(tmp(i).name);
        display_domain(mymovie, opts);
      end

      return;
    else
      load(mymovie);
    end
  end

  [data, theta, cyto] = gather_quantification(mymovie, opts);

  if (nargin > 2 & ~isempty(removed))
    data(removed, :) = NaN;
  end

  ticks = [1:50:length(theta) length(theta)];

  figure;imagesc(imnorm(data, [], [], 'r'));
  colormap(hot);
  hold on;
  plot([1 length(theta)], cyto+[0 0],'b', 'LineWidth', 2)
  set(gca, 'XTick', ticks, 'XTickLabel', theta(ticks));

  title(mymovie.experiment);

  saveas(gca, [mymovie.experiment 'domain.png']);

  return;
end
