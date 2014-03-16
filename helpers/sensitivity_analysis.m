%function test7(files, aligns)
%function test7(files, aligns)
function [conf_int] = sensitivity_analysis(pts, name)

  % To use conf_int:
  %tmp = conf_int([4 5 12 13], :)
  %tmp(1,:) = tmp(1,:) * tmp(4,1) 
  %tmp(3,:) = tmp(3,:) * tmp(2,1)
  %max(abs(tmp(:,[2 3]) - tmp(:,[1 1])), [], 2)

  if (nargin < 2)
    name = '';
  end

  legend = {'D_A', 'k_{A+}', 'k_{A-}', 'k_{AP}', '\alpha','\rho_A','\psi', 'L', ...
            'D_P', 'k_{P+}', 'k_{P-}', 'k_{PA}', '\beta', '\rho_P', '\nu', '\gamma'};

  if (isstruct(pts))

    nrates = size(pts.rate, 2);
    noffset = size(pts.offset, 2);
    nenergy = size(pts.energy, 2);
    nvisc = size(pts.viscosity, 2);
    nflow = size(pts.flow, 2);
    nscale = size(pts.flow_scaling, 2);

    pts = [pts.score pts.rate pts.offset pts.energy pts.viscosity pts.flow pts.flow_scaling];

    legend = [legend(1:nrates), ... 
              cellstr([repmat('\delta_', noffset, 1), num2str([1:noffset].')]).', ...
              cellstr([repmat('E_', nenergy, 1), num2str([1:nenergy].')]).', ...
              cellstr([repmat('\eta_', nvisc, 1), num2str([1:nvisc].')]).', ...
              legend(end-(nflow+nscale-1):end)];

    legend = legend(~cellfun('isempty', legend));

  else
    nshifts = size(pts,2) - length(legend) - 1;
    shifts = cellstr([repmat('\delta_', nshifts, 1), num2str([1:nshifts].')]);
    legend = [legend(1:end-1) shifts.' legend(end)];
  end

  is_fixed = any(isnan(pts), 1) | all(bsxfun(@eq, pts, pts(1,:)), 1);
  pts = pts(:, ~is_fixed);
  legend = legend(~is_fixed(2:end));

  hfig = figure;
  haxes = axes('Parent', hfig, 'NextPlot', 'add', 'box', 'on');
  hinset = axes('Parent', hfig, 'Position', [0.2 0.65 0.25 0.25], 'box', 'on', 'NextPlot', 'add');

  npts = length(unique(pts(:,2)));
  sep_sampling = NaN(npts, size(pts,2)-1);
  pos_sampling = NaN(npts, size(pts,2)-1);
  mscore = NaN;
  counter = zeros(1, size(pts,2));

  %pts = sortrows(pts);
  indx = 1;

  std_values = median(pts(:,2:end));
  for i=1:size(pts,1)
    var_indx = find(pts(i,2:end)~=std_values, 1);

    if (isempty(var_indx))
      mscore = pts(i,1);
    else
      counter(var_indx) = counter(var_indx)+1;
      sep_sampling(counter(var_indx), var_indx) = pts(i,1);
      pos_sampling(counter(var_indx), var_indx) = pts(i,var_indx+1);
    end
  end
  sep_sampling(end,:) = mscore;
  pos_sampling(end,:) = std_values;

  conf_int = NaN(length(std_values), 3);
  conf_int(:,1) = std_values(:);

  pos_sampling = (pos_sampling);
  tmp_std = std_values;
  tmp_std(tmp_std == 0) = 1;
  rel_sampling = bsxfun(@rdivide, pos_sampling, tmp_std);

  for i=1:size(sep_sampling,2)
    [rel_sampling(:,i), indxs] = sort(rel_sampling(:,i));
    pos_sampling(:,i) = pos_sampling(indxs, i);
    sep_sampling(:,i) = sep_sampling(indxs, i);

    pos = rel_sampling(:,i);
    likeli = sep_sampling(:,i);

    is_neg = any(pos <= 0);

    if (is_neg || std_values(i) == 0)
      center = find(pos == 0);
    else
      center = find(pos == 1);
    end

    [mval, mindx] = min(likeli);
    ci = (likeli <= likeli(center) + 1.92);
    left_bound = find(ci, 1, 'first')-1;
    right_bound = find(ci, 1, 'last')+1;

    if (isempty(left_bound) || left_bound < 1)
      left_bound = 1;
      conf_int(i, 2) = -Inf;
    else
      conf_int(i, 2) = pos_sampling(left_bound, i);
    end
    if (isempty(right_bound) || right_bound > length(ci))
      right_bound = length(ci);
      conf_int(i, 3) = Inf;
    else
      conf_int(i, 3) = pos_sampling(right_bound, i);
    end
    right_bound = min(right_bound, length(ci));
    ci = pos([left_bound right_bound]);

    rindx = -2;
    vals = roundn(pos_sampling([left_bound right_bound], i), rindx);
    while(vals(1) ==  vals(2))
      rindx = rindx-1;
      vals = roundn(pos_sampling([left_bound right_bound], i), rindx);
    end

    if (is_neg)
      zoom = ((pos >= ci(1)*5) & pos <= ci(2)*5);
    else
      zoom = ((pos >= ci(1)*0.75) & pos <= ci(2)*1.2);
    end
    cla(haxes)
    xlim(haxes, [min(pos(zoom)) max(pos(zoom))]);
    if (is_neg)
      set(haxes, 'XScale', 'linear');
    else
      set(haxes, 'XScale', 'log');
    end

    plot(haxes, pos(zoom), likeli(zoom), 'LineWidth', 2, 'Color', [83 83 83] / 255);

    ylims = ylim(haxes);
    xlims = xlim(haxes);
    plot(haxes, pos(center)*[1 1],ylims, 'k');

    %if (mindx~=center)
      scatter(haxes, pos(mindx), mval, 72, [215 25 28]/255, 'filled');
      text(pos(mindx), mval, num2str(pos(mindx)), 'Parent', haxes);
    %end

    m = mean(ci);
    s = ci(1) - m;

    errorbarxy(haxes, m, mean(ylims), s, [], {'k', 'k', 'k'});
    text(ci, [1 1]*(mean(ylims)+0.05*range(ylims)), num2str(vals), 'HorizontalAlignment', 'center', 'Parent', haxes);
    %plot(pos([left_bound, right_bound]),likeli([left_bound, right_bound]), 'g');

    set(haxes,'XTickLabel',roundn(get(haxes,'xtick'), -2));
    xlabel(haxes, legend{i})

    %axes(hinset);

    cla(hinset);
    xlim(hinset, [pos(1) pos(end)]);
    if (is_neg)
      set(hinset, 'XScale', 'linear');
    else
      set(hinset, 'XScale', 'log');
    end

    rectangle('Parent', hinset, 'Position', [xlims(1) ylims(1) diff(xlims) diff(ylims)], 'EdgeColor', 'none', 'FaceColor', [217 217 217]/255);
    plot(pos, likeli, 'LineWidth', 2, 'Color', [83 83 83] / 255);
    ylims = ylim(hinset);
    plot(hinset, pos(center)*[1 1], ylims, 'k');

    set(hinset,'XTickLabel',roundn(get(hinset,'xtick'), -2));

    %{
    subplot(1,2,1)
    if (is_neg)
      plot(pos, likeli);
    else
      semilogx(pos, likeli);
    end
    hold on;
    if (is_neg)
      plot([0 0],[mval max(likeli)], 'k');
    else
      plot([1 1],[mval max(likeli)], 'k');
    end
    if (mindx~=center)
      scatter(pos(mindx), mval, 'r');
    end
    plot(pos([left_bound, right_bound]),likeli([left_bound, right_bound]), 'g');
    title(num2str(pos_sampling([left_bound, right_bound], i)))

    zoom = ((pos >= pos(left_bound)*0.75) & pos <= pos(right_bound)*1.2);
    hold off;

    subplot(1,2,2)
    if (is_neg)
      plot(pos(zoom), likeli(zoom));
    else
      semilogx(pos(zoom), likeli(zoom));
    end
    hold on;
    if (is_neg)
      plot([0 0],[mval max(likeli(zoom))], 'k');
    else
      plot([1 1],[mval max(likeli(zoom))], 'k');
    end
    if (mindx~=center)
      scatter(pos(mindx), mval, 'r');
    end
    plot(pos([left_bound, right_bound]),likeli([left_bound, right_bound]), 'g');
    hold off;
    %}

    plot2svg(['PNG/' name '_sensitivity_' num2str(i) '.svg'], hfig);
  end

  return;

  for i=1:length(files)
    load(files{i});
    figure;imagesc(gather_quantification(mymovie, opts));
    title(files{i});
  end

  return;

  nfiles = length(files);
  [timing, names] = get_manual_timing;

  all_times = NaN(nfiles, 4);

  for i=1:nfiles
    goods = ismember(names, files{i});
    all_times(i, 1:3) = timing(goods, :);
    all_times(i, end) = aligns(i);
  end

  keyboard

  return;
end
