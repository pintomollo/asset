function analyse_domains(results)

  if (nargin == 0 | isempty(results))
    results = combine_domains('1056-24-*.mat');
  end

  nmovies = length(results{1});
  nframes = NaN(nmovies, 1);
  centered = cell(nmovies, 1);

  first = 0;
  last = 0;
  cytos = results{3};

  for i=1:nmovies
    nframes(i) = size(results{2}{i}, 1);
    if (1-cytos(i) < first)
      first = 1-cytos(i);
    end
    if (cytos(i)-nframes(i) > last)
      last = cytos(i)-nframes(i);
    end

    inv = (results{2}{i}(:,1) > results{2}{i}(:,2));
    results{2}{i}(inv,1) = results{2}{i}(inv,1) - 0.5;
    results{2}{i}(inv,1) = results{2}{i}(inv,1) + 0.5;
    results{2}{i} = results{2}{i} - 0.5;

    center = mean(results{2}{i}, 2);
    centered{i} = [center results{2}{i}(:,2) - center];

    plot([1-cytos(i) nframes(i)-cytos(i)], [0 0], 'g');
    hold on;
    plot([1-cytos(i):nframes(i)-cytos(i)], results{2}{i}(:,1), 'k');
    plot([1-cytos(i):nframes(i)-cytos(i)], results{2}{i}(:,2), 'k');
    plot([1-cytos(i):nframes(i)-cytos(i)], center, 'b');
    plot([1-cytos(i):nframes(i)-cytos(i)], centered{i}(:,2), 'r');

    saveas(gca, [results{1}{i}(1:end-5) '-domain.png']);
    hold off
  end

  datas = NaN(nmovies, last - first + 1, 2);
  cyto = 1-first;

  for i=1:nmovies
    datas(i, [1-cytos(i)-first+1:1-cytos(i)-first+nframes(i)],1) = centered{i}(:,1);
    datas(i, [1-cytos(i)-first+1:1-cytos(i)-first+nframes(i)],2) = centered{i}(:,2);
  end

  avg = squeeze(mymean(datas, 1));
  last = size(datas, 2) - cyto;
  
  plot([first last], [0 0], 'g');
  hold on;
  plot([first:last], squeeze(datas(:,:,1)).', 'b');
  plot([first:last], avg(:,1), 'c');
  hold off;
  saveas(gca, 'merge-centers.png');

  plot([first last], [0 0], 'g');
  hold on;
  plot([first:last], squeeze(datas(:,:,2)).', 'r');
  plot([first:last], avg(:,2), 'm');
  saveas(gca, 'merge-domains.png');

  return;
end
