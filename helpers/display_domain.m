function display_domain(mymovie, opts, removed)
  
  if (ischar(mymovie))
    if (~isempty(findstr(mymovie, '*')))
      tmp = dir(mymovie); 
      mymovies = cell(1, length(tmp));

      close all;

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

  ticks = [1:100:length(theta) length(theta)];

  %imagesc(imnorm(data, [], [], 'r'));
  %colormap(jet);
  %hold on;
  %plot([1 length(theta)], cyto+[0 0],'b', 'LineWidth', 2)
  %set(gca, 'XTick', ticks, 'XTickLabel', theta(ticks));

  %title(mymovie.experiment);

  %saveas(gca, [mymovie.experiment '-norm.png']);

  imagesc(imnorm(data));
  colormap(jet);
  hold on;
  plot([1 length(theta)], cyto+[0 0],'k', 'LineWidth', 2)
  set(gca, 'XTick', ticks, 'XTickLabel', theta(ticks));

  title(mymovie.experiment);

  saveas(gca, [mymovie.experiment '.png']);
  hold off;

  return;
end
