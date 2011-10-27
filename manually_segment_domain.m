function manually_segment_domains(fname)

  if (ischar(fname))
    if (~isempty(findstr(fname, '*')))
      tmp = dir(fname); 
      mymovies = cell(1, length(tmp));

      for i=1:length(tmp)
        mymovies{i} = tmp(i).name;
      end
    else
      mymovies = {fname}; 
    end
  elseif (~iscell(fname))
    mymovies = {};
    mymovies{1} = fname;
  end

  nmovies = length(mymovies);

  for i=1:nmovies
    mymovie = mymovies{i};

    if (ischar(mymovie))
      domain = imread(mymovie);
    elseif (isnumeric(mymovie))
      domain = mymovie;
      mymovie = 'unkown_name.png'
    end
    [h, w] = size(domain);

    imagesc(domain);
    [bw, xi, yi] = roipoly();

    if (numel(xi) < 10)
      continue;
    end

    start = find(yi == min(yi), 1);
    if (start ~= 1)
      xi = [xi(start:end); xi(1:start-1)];
      yi = [yi(start:end); yi(1:start-1)];
    end
    start = ceil(yi(1));

    dy = [0; sign(diff(yi))];
    ends = find(dy == -dy(2), 1, 'first');
    dist = hypot(diff(xi(ends-1:ends+1)), diff(yi(ends-1:ends+1)));
  
    if (dist(1) > dist(2))
      ends = ends-1;
    end

    half1 = [xi(1:ends) yi(1:ends)];
    half2 = [xi([ends+1:end 1]) yi([ends+1:end 1])];

    if (dy(2) < 0)
      half1 = half1([end:-1:1], :);
    else
      half2 = half2([end:-1:1], :);
    end

    half1(end, 2) = h;
    half2(end, 2) = h;

    indexes = [start:h].';
    x1 = interp1q(half1(:, 2), half1(:, 1), indexes);

    x2 = interp1q(half2(:, 2), half2(:, 1), indexes);

    mean_pos = mean([x1,x2], 2);
    width = abs((x1 - x2) / 2);

    path = NaN(h, 2);
    path(start:end,:) = [mean_pos, width];

    %imagesc(domain)
    %hold on;
    %plot(path(:, 1), [1:h], 'k');
    %plot(path(:, 1) + path(:,2), [1:h], 'k');
    %plot(path(:, 1) - path(:,2), [1:h], 'k');
    %hold off;

    %keyboard

    save([mymovie(1:end-3) '-Manual-DP.mat'], 'domain', 'path');
  end

  return;
end
